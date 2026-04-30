"""
Profile likelihood analysis for practical identifiability.

For each parameter p_i, we fix p_i = p_i* + delta (varying delta over a
grid), re-optimise the remaining 11 parameters, and record the achieved
minimum cost. A non-identifiable parameter yields a FLAT profile: the cost
does not increase when the parameter is moved far from its true value,
because the other parameters can compensate.

This is the gold standard for practical identifiability:
  - Parabolic profile -> well identified
  - Flat profile      -> structurally or practically non-identifiable
  - Asymmetric / multi-modal profile -> optimization issues

Reference: Raue et al. 2009 "Structural and practical identifiability
analysis of partially observed dynamical models..." (Bioinformatics).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import contextlib
import io
import time


PARAM_NAMES = [
    "lom_1", "lom_2", "lom_3",
    "phio_1", "phio_2", "phio_3",
    "Fom_1", "Fom_2", "Fom_3",
    "lst_1", "lst_2", "lst_3",
]
PARAM_UNITS = ["m"] * 3 + ["rad"] * 3 + ["N"] * 3 + ["m"] * 3


def profile_one_parameter(
    param_index,
    true_params,
    grid,
    data,
    lower_band,
    upper_band,
    skeleton_num,
    unknown_parameters,
    casadi_function,
    optimization_nlp,
    silent=True,
):
    """For ONE parameter p_i, compute the profile cost on the given grid.

    At each grid point:
      1. Build a modified NLP where p_i is clamped (lb = ub = grid value).
      2. Solve the NLP on the remaining 11 parameters.
      3. Record the achieved cost.

    Returns
    -------
    grid : np.ndarray
    costs : np.ndarray of the same length as grid, NaN for failed runs
    """
    costs = np.full(len(grid), np.nan)

    # Start each profile point from the true params as warm start
    x0_base = np.asarray(true_params, dtype=float).copy()

    for k, value in enumerate(grid):
        lb_k = np.array(lower_band, dtype=float).copy()
        ub_k = np.array(upper_band, dtype=float).copy()

        # Clamp the profiled parameter
        lb_k[param_index] = value
        ub_k[param_index] = value

        # Initial guess: true params but with the profiled component set
        x0 = x0_base.copy()
        x0[param_index] = value
        # Also make sure x0 respects the other bounds
        x0 = np.clip(x0, lb_k, ub_k)

        try:
            if silent:
                with contextlib.redirect_stdout(io.StringIO()):
                    params_est = optimization_nlp(
                        data, x0, lb_k, ub_k,
                        skeleton_num, true_params,
                        unknown_parameters, casadi_function,
                    )
            else:
                params_est = optimization_nlp(
                    data, x0, lb_k, ub_k,
                    skeleton_num, true_params,
                    unknown_parameters, casadi_function,
                )
            if params_est is None:
                continue
            # The cost is not returned by optimization_nlp directly;
            # we capture it by running a lightweight cost evaluation.
            # Simpler: monkey-patch optimization_nlp to return the cost.
            # For now, we recompute an approximate cost by re-running in a
            # way that just reports. Users can adapt this.
            # Here, we re-solve with tight bounds and the cost will be
            # printed; we parse it from stdout.
            costs[k] = _evaluate_cost_after_fit(
                data, params_est, skeleton_num,
                unknown_parameters, casadi_function,
            )
        except Exception as e:
            print(f"  [profile {PARAM_NAMES[param_index]}] "
                  f"grid point {k+1}/{len(grid)} failed: {e}")

    return costs


def _evaluate_cost_after_fit(
    data, params_fit, skeleton_num,
    unknown_parameters, casadi_function,
):
    """Evaluate the NLP cost given a fitted parameter vector.

    This re-constructs the cost expression with the states solved by the
    equilibrium rootfinder and evaluates it numerically.
    """
    from casadi import SX, vertcat, sum1
    import casadi as ca

    w_torque = 1.0
    w_length = 0.005
    w_angle = (1.0 / 180.0) * np.pi

    # Rootfinder template
    x_rf = SX.sym("x_rf", 9)
    a_rf = SX.sym("a_rf", 3)
    mtu_rf = SX.sym("mtu_rf", 3)
    p_rf = SX.sym("p_rf", 12)
    eq = casadi_function["equilibrium_error_all_muscle"](
        x_rf, vertcat(a_rf, mtu_rf, p_rf)
    )
    rf = ca.rootfinder(
        "rf", "newton",
        {"x": x_rf, "p": vertcat(a_rf, mtu_rf, p_rf), "g": eq},
        {"error_on_fail": False},
    )

    total_cost = 0.0
    n_trials = data.shape[1]
    for trial in range(n_trials):
        d = data[:, trial]
        q_trial = [0, 0, 0, 0] + list(d[1:3])
        a_trial = d[3:6]
        fl_meas = np.abs(d[6:9])
        pa_meas = d[9:12]
        tl_meas = np.abs(d[12:15])
        mtu = casadi_function["get_mtu_length"](
            q_trial + list(skeleton_num)
        )
        x0 = np.concatenate([fl_meas, pa_meas, tl_meas])
        try:
            states = np.array(
                rf(x0, vertcat(a_trial, mtu, params_fit))
            ).flatten()
        except Exception:
            return np.nan

        neuromusculo = np.concatenate([a_trial, q_trial, skeleton_num])
        all_states_num = np.concatenate([neuromusculo, states])
        torque_sim = float(
            casadi_function["get_joint_moment"](all_states_num, params_fit)
        )

        e_t = d[0] - torque_sim
        e_f = fl_meas - states[0:3]
        e_p = pa_meas - states[3:6]

        total_cost += w_torque * e_t ** 2
        total_cost += float(np.sum(w_length * e_f ** 2))
        total_cost += float(np.sum(w_angle * e_p ** 2))
    return total_cost


def run_profile_likelihood(
    data,
    true_params,
    lower_band,
    upper_band,
    skeleton_num,
    unknown_parameters,
    casadi_function,
    optimization_nlp,
    out_dir,
    n_grid=11,
    rel_range=0.30,
    param_indices=None,
):
    """Profile every parameter (or the selected subset) and save plots.

    Parameters
    ----------
    n_grid : int
        Number of grid points per parameter.
    rel_range : float
        Relative range around the true value (e.g. 0.3 -> +/-30%).
    param_indices : list[int] or None
        If None, all 12 parameters are profiled.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    true_params = np.asarray(true_params, dtype=float).flatten()
    if param_indices is None:
        param_indices = list(range(12))

    print("=" * 72)
    print("  Profile likelihood analysis")
    print("=" * 72)
    print(f"  n_grid          = {n_grid} points per parameter")
    print(f"  rel_range       = +/- {rel_range*100:.0f}% of true value")
    print(f"  parameters      = {[PARAM_NAMES[i] for i in param_indices]}")
    print()

    all_profiles = {}
    t_start = time.time()

    for pi in param_indices:
        name = PARAM_NAMES[pi]
        true_val = true_params[pi]

        # Build grid
        if abs(true_val) > 1e-10:
            delta = rel_range * abs(true_val)
        else:
            delta = 0.1  # fallback for true_val ~ 0
        lo = max(lower_band[pi], true_val - delta)
        hi = min(upper_band[pi], true_val + delta)
        grid = np.linspace(lo, hi, n_grid)

        print(f"  Profiling {name}  ({true_val:.4g}, range [{lo:.4g},{hi:.4g}])...")
        t0 = time.time()
        costs = profile_one_parameter(
            pi, true_params, grid,
            data, lower_band, upper_band,
            skeleton_num, unknown_parameters, casadi_function,
            optimization_nlp,
        )
        dt = time.time() - t0
        print(f"    done in {dt:.1f} s (costs range: "
              f"[{np.nanmin(costs):.2e}, {np.nanmax(costs):.2e}])")

        all_profiles[name] = {"grid": grid, "costs": costs,
                              "true_value": true_val, "index": pi}

    t_total = time.time() - t_start
    print(f"\n  Total profiling time: {t_total:.1f} s")

    # Save and plot
    _save_profiles_csv(all_profiles, out / "profile_likelihood.csv")
    _plot_profiles(all_profiles, out / "profile_likelihood.png")

    print(f"\n  Outputs saved to: {out.resolve()}")
    return all_profiles


def _save_profiles_csv(profiles, out_path):
    rows = []
    for name, d in profiles.items():
        for g, c in zip(d["grid"], d["costs"]):
            rows.append({
                "parameter": name,
                "grid_value": g,
                "cost": c,
                "true_value": d["true_value"],
            })
    pd.DataFrame(rows).to_csv(out_path, index=False)


def _plot_profiles(profiles, out_path):
    n = len(profiles)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes = axes.flatten() if n > 1 else [axes]

    for ax, (name, d) in zip(axes, profiles.items()):
        grid = d["grid"]
        costs = d["costs"]
        true_val = d["true_value"]
        ax.plot(grid, costs, "o-", color="C0")
        ax.axvline(true_val, color="C3", linestyle="--", alpha=0.7,
                   label="true value")
        ax.set_xlabel(f"{name}")
        ax.set_ylabel("min cost (other params free)")
        ax.grid(alpha=0.3)

        # Flatness score: ratio max/min cost
        valid = ~np.isnan(costs)
        if valid.sum() > 2:
            cmin = np.nanmin(costs[valid])
            cmax = np.nanmax(costs[valid])
            if cmin > 0:
                ratio = cmax / cmin
                ax.set_title(f"{name}  (max/min = {ratio:.1f})")
            else:
                ax.set_title(name)
        else:
            ax.set_title(name)

    for ax in axes[len(profiles):]:
        ax.set_visible(False)

    fig.suptitle("Profile likelihood: cost vs clamped parameter value",
                 fontsize=14)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

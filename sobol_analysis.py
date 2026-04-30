"""
Sobol global sensitivity analysis of the identified parameters.

Whereas local tools (Hessian, profile likelihood) probe the cost LANDSCAPE
at/near the optimum, Sobol indices quantify how much of the variance of
a model OUTPUT (here: the simulated torque / fiber length / pennation)
is attributable to each parameter, across the full admissible parameter
space.

  - First-order index S1_i: fraction of variance explained by parameter i
    alone.
  - Total-order index ST_i: fraction of variance to which parameter i
    contributes, INCLUDING its interactions with other parameters.

Low ST means that the output is insensitive to this parameter -> hard to
identify. Large gap (ST - S1) means interactions dominate.

Dependencies
------------
- SALib (pip install SALib)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import time

try:
    from SALib.sample import sobol as sobol_sample
    from SALib.analyze import sobol as sobol_analyze
    SALIB_AVAILABLE = True
except ImportError:
    SALIB_AVAILABLE = False


PARAM_NAMES = [
    "lom_1", "lom_2", "lom_3",
    "phio_1", "phio_2", "phio_3",
    "Fom_1", "Fom_2", "Fom_3",
    "lst_1", "lst_2", "lst_3",
]


def simulate_model(
    params,
    data,
    skeleton_num,
    casadi_function,
):
    """Run the forward model for a given parameter vector.

    Returns a scalar summary of the prediction (here: the sum of squared
    torques across trials). Alternative summaries can be implemented:
    mean torque, mean fiber length, etc.
    """
    from casadi import SX, vertcat
    import casadi as ca

    n_trials = data.shape[1]

    # Rootfinder for equilibrium
    x_rf = SX.sym("x_rf", 9)
    a_rf = SX.sym("a_rf", 3)
    mtu_rf = SX.sym("mtu_rf", 3)
    p_rf = SX.sym("p_rf", 12)
    eq = casadi_function["equilibrium_error_all_muscle"](
        x_rf, vertcat(a_rf, mtu_rf, p_rf)
    )
    rf = ca.rootfinder(
        "rf_sobol", "newton",
        {"x": x_rf, "p": vertcat(a_rf, mtu_rf, p_rf), "g": eq},
        {"error_on_fail": False},
    )

    torque_pred = np.full(n_trials, np.nan)

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
                rf(x0, vertcat(a_trial, mtu, params))
            ).flatten()
            if np.any(np.isnan(states)):
                continue
        except Exception:
            continue

        neuromusculo = np.concatenate([a_trial, q_trial, skeleton_num])
        all_states_num = np.concatenate([neuromusculo, states])
        try:
            torque_pred[trial] = float(
                casadi_function["get_joint_moment"](all_states_num, params)
            )
        except Exception:
            continue

    # Return RMS torque as the scalar output
    valid = ~np.isnan(torque_pred)
    if valid.sum() == 0:
        return np.nan
    return float(np.sqrt(np.mean(torque_pred[valid] ** 2)))


def run_sobol_analysis(
    data,
    lower_band,
    upper_band,
    skeleton_num,
    casadi_function,
    out_dir,
    n_base=256,
    calc_second_order=False,
    seed=42,
):
    """Sobol analysis of the model output.

    Parameters
    ----------
    n_base : int
        Saltelli base sample size. Total evaluations ~ n_base * (2D + 2)
        for first+total indices, ~ n_base * (2D + 2) + n_base * D*(D-1)/2
        if calc_second_order=True. With D=12, n_base=256 -> ~6700 model
        evaluations.

    Returns
    -------
    Si : dict returned by SALib.analyze.sobol.
    """
    if not SALIB_AVAILABLE:
        raise ImportError("SALib is required: pip install SALib")

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)

    print("=" * 72)
    print("  Sobol global sensitivity analysis")
    print("=" * 72)
    print(f"  n_base           = {n_base}")
    print(f"  parameters       = 12")
    mult = (2 * 12 + 2) if not calc_second_order else (12 + 2)
    total = n_base * mult
    print(f"  estimated evals  = {total}")
    print(f"  output quantity  = RMS torque across trials")
    print()

    problem = {
        "num_vars": 12,
        "names": PARAM_NAMES,
        "bounds": [[lb, ub] for lb, ub in zip(lower_band, upper_band)],
    }

    X = sobol_sample.sample(
        problem, n_base,
        calc_second_order=calc_second_order,
        seed=seed,
    )
    print(f"  Sampling done: {X.shape[0]} parameter sets to evaluate.")
    print(f"  Evaluating model on all samples... (this will take a while)")

    Y = np.zeros(X.shape[0])
    t0 = time.time()
    n_fail = 0
    for i in range(X.shape[0]):
        Y[i] = simulate_model(X[i], data, skeleton_num, casadi_function)
        if np.isnan(Y[i]):
            n_fail += 1
        if (i + 1) % max(1, X.shape[0] // 20) == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (X.shape[0] - (i + 1))
            print(f"    {i+1}/{X.shape[0]}  "
                  f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining, "
                  f"{n_fail} failures)")

    print(f"\n  Model evaluation done ({n_fail} failed).")

    # Handle NaNs by replacing with median
    nan_mask = np.isnan(Y)
    if nan_mask.any():
        Y[nan_mask] = np.nanmedian(Y)
        print(f"  WARNING: {nan_mask.sum()} NaN values replaced with median.")

    Si = sobol_analyze.analyze(
        problem, Y,
        calc_second_order=calc_second_order,
        print_to_console=False,
        seed=seed,
    )

    _print_sobol_table(Si, calc_second_order)
    _save_sobol_csv(Si, out / "sobol_indices.csv", calc_second_order)
    _plot_sobol(Si, out / "sobol_indices.png")

    print(f"\n  Outputs saved to: {out.resolve()}")
    return Si


def _print_sobol_table(Si, calc_second_order):
    print("\n  Sobol indices (RMS torque output):")
    print("  " + "-" * 65)
    print(f"  {'Parameter':<10} {'S1':>10} {'S1_conf':>10} "
          f"{'ST':>10} {'ST_conf':>10}")
    print("  " + "-" * 65)
    order = np.argsort(Si["ST"])[::-1]
    for i in order:
        print(f"  {PARAM_NAMES[i]:<10} "
              f"{Si['S1'][i]:>10.4f} {Si['S1_conf'][i]:>10.4f} "
              f"{Si['ST'][i]:>10.4f} {Si['ST_conf'][i]:>10.4f}")
    print("  " + "-" * 65)


def _save_sobol_csv(Si, out_path, calc_second_order):
    df = pd.DataFrame({
        "parameter": PARAM_NAMES,
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })
    df.to_csv(out_path, index=False)


def _plot_sobol(Si, out_path):
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(12)
    width = 0.4

    s1 = np.asarray(Si["S1"])
    st = np.asarray(Si["ST"])
    s1_err = np.asarray(Si["S1_conf"])
    st_err = np.asarray(Si["ST_conf"])

    ax.bar(x - width / 2, s1, width, yerr=s1_err, label="S1 (first order)",
           capsize=3, color="C0", alpha=0.8)
    ax.bar(x + width / 2, st, width, yerr=st_err, label="ST (total order)",
           capsize=3, color="C3", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(PARAM_NAMES, rotation=45, ha="right")
    ax.set_ylabel("Sobol index")
    ax.set_title("Global sensitivity of RMS torque to muscle-tendon parameters")
    ax.axhline(0.05, color="gray", linestyle=":", alpha=0.5,
               label="5% threshold")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

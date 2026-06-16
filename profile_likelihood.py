"""
Profile likelihood analysis for practical identifiability.

[PATCH param_config]
--------------------
Le vecteur de paramètres n'est plus figé à 12 : il correspond aux seuls
paramètres 'sym' déclarés dans param_config, dans l'ordre de param_index
(sortie de useful.get_initial_guess). On profile donc n_opt paramètres,
pas 12. Les paramètres 'fixed' sont codés en dur dans les casadi_function
et ne sont jamais profilés.

Pour chaque paramètre p_i du vecteur optimisé, on fixe p_i = p_i* + delta
(grille), on ré-optimise les autres (n_opt - 1) paramètres, et on enregistre
le coût minimal atteint. Profil plat -> non identifiable.

Reference: Raue et al. 2009 (Bioinformatics).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import contextlib
import io
import time

from param_layout import Layout


def profile_one_parameter(
    param_pos,
    true_params,
    grid,
    data,
    lower_band,
    upper_band,
    skeleton_num,
    unknown_parameters,
    casadi_function,
    optimization_nlp,
    param_index,
    layout,
    silent=True,
):
    """Profil de coût d'UN paramètre (à la position param_pos) sur la grille.

    À chaque point de grille :
      1. NLP modifié où le paramètre profilé est verrouillé (lb = ub = valeur).
      2. Résolution sur les (n_opt - 1) autres paramètres.
      3. Enregistrement du coût atteint.
    """
    name = layout.names[param_pos]
    costs = np.full(len(grid), np.nan)
    x0_base = np.asarray(true_params, dtype=float).copy()

    for k, value in enumerate(grid):
        lb_k = np.array(lower_band, dtype=float).copy()
        ub_k = np.array(upper_band, dtype=float).copy()

        # Verrouille le paramètre profilé
        lb_k[param_pos] = value
        ub_k[param_pos] = value

        x0 = x0_base.copy()
        x0[param_pos] = value
        x0 = np.clip(x0, lb_k, ub_k)

        try:
            if silent:
                with contextlib.redirect_stdout(io.StringIO()):
                    params_est = optimization_nlp(
                        data, x0, lb_k, ub_k,
                        skeleton_num, true_params,
                        unknown_parameters, casadi_function,
                        param_index,
                    )
            else:
                params_est = optimization_nlp(
                    data, x0, lb_k, ub_k,
                    skeleton_num, true_params,
                    unknown_parameters, casadi_function,
                    param_index,
                )
            if params_est is None:
                continue
            costs[k] = _evaluate_cost_after_fit(
                data, params_est, skeleton_num,
                unknown_parameters, casadi_function,
            )
        except Exception as e:
            print(f"  [profile {name}] "
                  f"grid point {k+1}/{len(grid)} failed: {e}")

    return costs


def _evaluate_cost_after_fit(
    data, params_fit, skeleton_num,
    unknown_parameters, casadi_function,
):
    """Évalue le coût NLP pour un vecteur de paramètres ajusté.

    NB: params_fit contient uniquement les 'sym' (taille n_opt) ; les 'fixed'
    sont déjà internes aux casadi_function, donc rien ne change ici hormis le
    fait que params_fit n'est plus supposé de taille 12.
    """
    from casadi import SX, vertcat
    import casadi as ca

    n_opt = unknown_parameters.shape[0]

    w_torque = 1.0
    w_length = 0.005
    w_angle = (1.0 / 180.0) * np.pi

    x_rf = SX.sym("x_rf", 9)
    a_rf = SX.sym("a_rf", 3)
    mtu_rf = SX.sym("mtu_rf", 3)
    p_rf = SX.sym("p_rf", n_opt)
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
    param_index,
    out_dir,
    n_grid=11,
    rel_range=0.30,
    param_names=None,
):
    """Profile chaque paramètre 'sym' (ou un sous-ensemble) et sauve les plots.

    Parameters
    ----------
    param_index : dict
        Sortie de get_initial_guess ; pilote taille, ordre et NLP.
    param_names : list[str] or None
        Sous-ensemble à profiler par NOM (ex: ['lom_1', 'Fom_2']).
        Si None, tous les paramètres 'sym' sont profilés.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    layout = Layout.from_param_index(param_index)
    true_params = np.asarray(true_params, dtype=float).flatten()
    assert true_params.shape[0] == layout.n_opt, (
        f"true_params ({true_params.shape[0]}) incohérent avec "
        f"n_opt={layout.n_opt} déduit de param_index"
    )

    if param_names is None:
        positions = list(range(layout.n_opt))
    else:
        positions = [layout.position(nm) for nm in param_names]

    print("=" * 72)
    print("  Profile likelihood analysis")
    print("=" * 72)
    print(f"  {layout.summary()}")
    print(f"  n_grid          = {n_grid} points per parameter")
    print(f"  rel_range       = +/- {rel_range*100:.0f}% of true value")
    print(f"  parameters      = {[layout.names[p] for p in positions]}")
    print()

    all_profiles = {}
    t_start = time.time()

    for pos in positions:
        name = layout.names[pos]
        true_val = true_params[pos]

        if abs(true_val) > 1e-10:
            delta = rel_range * abs(true_val)
        else:
            delta = 0.1
        lo = max(lower_band[pos], true_val - delta)
        hi = min(upper_band[pos], true_val + delta)
        grid = np.linspace(lo, hi, n_grid)

        print(f"  Profiling {name}  ({true_val:.4g}, range [{lo:.4g},{hi:.4g}])...")
        t0 = time.time()
        costs = profile_one_parameter(
            pos, true_params, grid,
            data, lower_band, upper_band,
            skeleton_num, unknown_parameters, casadi_function,
            optimization_nlp, param_index, layout,
        )
        dt = time.time() - t0
        print(f"    done in {dt:.1f} s (costs range: "
              f"[{np.nanmin(costs):.2e}, {np.nanmax(costs):.2e}])")

        all_profiles[name] = {"grid": grid, "costs": costs,
                              "true_value": true_val, "index": pos}

    t_total = time.time() - t_start
    print(f"\n  Total profiling time: {t_total:.1f} s")

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
    ncols = min(4, max(1, n))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, (name, d) in zip(axes, profiles.items()):
        grid = d["grid"]
        costs = d["costs"]
        costs = np.log(costs)
        true_val = d["true_value"]
        ax.plot(grid, costs, "o-", color="C0")
        ax.axvline(true_val, color="C3", linestyle="--", alpha=0.7,
                   label="true value")
        ax.set_xlabel(f"{name}")
        ax.set_ylabel("log min cost (other params free)")
        ax.grid(alpha=0.3)

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

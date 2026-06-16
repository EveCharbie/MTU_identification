"""
Sobol global sensitivity analysis of the identified parameters.

[PATCH param_config]
--------------------
Le nombre de variables Sobol n'est plus 12 mais n_opt = nombre de paramètres
'sym' (déduit de param_index / get_initial_guess). Les bornes (lower_band,
upper_band) doivent être de taille n_opt et alignées sur l'ordre de
param_index. Les paramètres 'fixed' sont dans les casadi_function et ne
participent pas à l'analyse.

  - S1_i : fraction de variance expliquée par le paramètre i seul.
  - ST_i : fraction de variance à laquelle i contribue, interactions incluses.

Dependencies: SALib (pip install SALib)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import time

from param_layout import Layout

try:
    from SALib.sample import sobol as sobol_sample
    from SALib.analyze import sobol as sobol_analyze
    SALIB_AVAILABLE = True
except ImportError:
    SALIB_AVAILABLE = False


def simulate_model(
    params,
    data,
    skeleton_num,
    casadi_function,
):
    """Forward model pour un vecteur de paramètres (taille n_opt).

    Renvoie un résumé scalaire (RMS du couple sur les essais).
    """
    from casadi import SX, vertcat
    import casadi as ca

    n_trials = data.shape[1]
    n_opt = len(params)

    x_rf = SX.sym("x_rf", 9)
    a_rf = SX.sym("a_rf", 3)
    mtu_rf = SX.sym("mtu_rf", 3)
    p_rf = SX.sym("p_rf", n_opt)
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
    param_index,
    out_dir,
    n_base=256,
    calc_second_order=False,
    seed=42,
):
    """Analyse de Sobol de la sortie du modèle.

    Parameters
    ----------
    param_index : dict
        Sortie de get_initial_guess ; fixe le nombre/ordre des variables.
    n_base : int
        Taille de base Saltelli. Total ~ n_base * (2*D + 2), D = n_opt.
    """
    if not SALIB_AVAILABLE:
        raise ImportError("SALib is required: pip install SALib")

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    layout = Layout.from_param_index(param_index)
    D = layout.n_opt

    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)
    assert lower_band.shape[0] == D and upper_band.shape[0] == D, (
        f"bornes (taille {lower_band.shape[0]}) incohérentes avec "
        f"n_opt={D} déduit de param_index"
    )

    mult = (2 * D + 2) if not calc_second_order else (D + 2)
    total = n_base * mult

    print("=" * 72)
    print("  Sobol global sensitivity analysis")
    print("=" * 72)
    print(f"  {layout.summary()}")
    print(f"  n_base           = {n_base}")
    print(f"  parameters       = {D}")
    print(f"  estimated evals  = {total}")
    print(f"  output quantity  = RMS torque across trials")
    print()

    problem = {
        "num_vars": D,
        "names": layout.names,
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

    _print_sobol_table(Si, layout)
    _save_sobol_csv(Si, out / "sobol_indices.csv", layout)
    _plot_sobol(Si, out / "sobol_indices.png", layout)

    print(f"\n  Outputs saved to: {out.resolve()}")
    return Si


def _print_sobol_table(Si, layout):
    print("\n  Sobol indices (RMS torque output):")
    print("  " + "-" * 65)
    print(f"  {'Parameter':<10} {'S1':>10} {'S1_conf':>10} "
          f"{'ST':>10} {'ST_conf':>10}")
    print("  " + "-" * 65)
    order = np.argsort(Si["ST"])[::-1]
    for i in order:
        print(f"  {layout.names[i]:<10} "
              f"{Si['S1'][i]:>10.4f} {Si['S1_conf'][i]:>10.4f} "
              f"{Si['ST'][i]:>10.4f} {Si['ST_conf'][i]:>10.4f}")
    print("  " + "-" * 65)


def _save_sobol_csv(Si, out_path, layout):
    df = pd.DataFrame({
        "parameter": layout.names,
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })
    df.to_csv(out_path, index=False)


def _plot_sobol(Si, out_path, layout):
    D = layout.n_opt
    fig, ax = plt.subplots(figsize=(max(8, D * 0.8), 6))
    x = np.arange(D)
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
    ax.set_xticklabels(layout.names, rotation=45, ha="right")
    ax.set_ylabel("Sobol index")
    ax.set_title("Global sensitivity of RMS torque to muscle-tendon parameters")
    ax.axhline(0.05, color="gray", linestyle=":", alpha=0.5,
               label="5% threshold")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

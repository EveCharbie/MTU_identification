"""
Conditioning analysis of the identification NLP.

Computes the Hessian of the cost at the true parameters (local curvature
of the cost landscape along each parameter direction). Eigen-analysis of
this Hessian reveals:
  - well-identified directions (large eigenvalues -> steep curvature)
  - sloppy directions (small eigenvalues -> flat valleys, poor identifiability)
  - the parameter combinations that form those directions (eigenvectors)

Outputs
-------
- Console: eigenvalues sorted, condition number, dominant parameter
  combinations for the N sloppiest directions.
- `hessian_eigenvalues.png`: scree plot of eigenvalues (log scale).
- `hessian_sloppy_directions.png`: bar plots of the loadings of the
  worst eigenvectors.
- `conditioning_summary.csv`: raw eigenvalues and eigenvectors.

Key references
--------------
- Gutenkunst et al. 2007, "Universally Sloppy Parameter Sensitivities..."
- Chis et al. 2011, "Structural Identifiability of Systems Biology Models"
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import casadi as ca
from pathlib import Path


PARAM_NAMES = [
    "lom_1", "lom_2", "lom_3",
    "phio_1", "phio_2", "phio_3",
    "Fom_1", "Fom_2", "Fom_3",
    "lst_1", "lst_2", "lst_3",
]


def build_cost_function_params_only(
    data,
    skeleton_num,
    unknown_parameters,
    casadi_function,
):
    """Build a CasADi Function cost(p) that depends ONLY on the 12 parameters.

    For each trial, the muscle states (fiber, pennation, tendon) are solved
    via a rootfinder on the equilibrium equations. The cost is then the sum
    over trials of the residuals on torque, fiber, pennation (same weights
    as the NLP).

    This flattens the augmented NLP into a reduced NLP on the parameters
    alone, which is what we need for a clean Hessian on the physical
    parameters.

    Returns
    -------
    cost_func : ca.Function([params12] -> [scalar cost])
    """
    from casadi import SX, vertcat, sum1

    n_trials = data.shape[1]
    n_muscle = 3

    # Weights — must match the NLP's cost definition
    w_torque = 1.0
    w_length = 0.005
    w_angle = (1.0 / 180.0) * np.pi

    # Parameters as symbolic
    p_sym = SX.sym("p", 12)

    # Build a rootfinder template: given (activations, mtu_length, params),
    # solve for (fiber, pennation, tendon) such that equilibrium = 0.
    x_rf = SX.sym("x_rf", 9)
    a_rf = SX.sym("a_rf", 3)
    mtu_rf = SX.sym("mtu_rf", 3)
    params_rf = SX.sym("params_rf", 12)
    p_for_equilibrium = vertcat(a_rf, mtu_rf, params_rf)

    equilibrium_expr = casadi_function["equilibrium_error_all_muscle"](
        x_rf, p_for_equilibrium
    )
    rf_problem = {
        "x": x_rf,
        "p": vertcat(a_rf, mtu_rf, params_rf),
        "g": equilibrium_expr,
    }
    rf = ca.rootfinder("rf_states", "newton", rf_problem,
                       {"error_on_fail": False})

    cost = SX(0.0)

    for trial in range(n_trials):
        d = data[:, trial]
        mesured_torque = d[0]
        q_trial = [0, 0, 0, 0] + list(d[1:3])
        a_trial = d[3:6]
        fl_meas = np.abs(d[6:9])
        pa_meas = d[9:12]
        tl_meas = np.abs(d[12:15])

        musculo = q_trial + list(skeleton_num)
        mtu_length = casadi_function["get_mtu_length"](musculo)

        # Initial guess: use the measured values
        x0_guess = np.concatenate([fl_meas, pa_meas, tl_meas])

        # Solve equilibrium for these params (symbolic in p_sym)
        states = rf(x0_guess, vertcat(a_trial, mtu_length, p_sym))
        fiber_k = states[0:3]
        pennat_k = states[3:6]

        # Simulated torque
        neuromusculo = np.concatenate([a_trial, q_trial, skeleton_num])
        all_states = vertcat(SX(neuromusculo.tolist()), states)
        torque_sim = casadi_function["get_joint_moment"](all_states, p_sym)

        # Residuals
        e_t = mesured_torque - torque_sim
        e_f = fl_meas - fiber_k
        e_p = pa_meas - pennat_k

        cost = cost + w_torque * e_t ** 2
        cost = cost + sum1(w_length * e_f ** 2)
        cost = cost + sum1(w_angle * e_p ** 2)

    return ca.Function("cost_params_only", [p_sym], [cost])


def build_cost_from_augmented(
    data,
    skeleton_num,
    unknown_parameters,
    casadi_function,
):
    """Alternative: build the FULL augmented cost (params + per-trial states)
    and return the Function along with the layout of the decision vector.

    Use this when the rootfinder approach in `build_cost_function_params_only`
    fails (e.g. equilibrium has multiple solutions, rootfinder divergence).
    The Hessian will then be taken over the full 12 + 9*n_trials vector and
    we marginalise on the 12x12 parameter block.
    """
    from casadi import SX, vertcat, sum1

    n_trials = data.shape[1]
    n_muscle = 3

    w_torque = 1.0
    w_length = 0.005
    w_angle = (1.0 / 180.0) * np.pi

    p_sym = SX.sym("p", 12)

    w_vars = [p_sym]
    cost = SX(0.0)
    g_list = []

    for trial in range(n_trials):
        d = data[:, trial]
        mesured_torque = d[0]
        q_trial = [0, 0, 0, 0] + list(d[1:3])
        a_trial = d[3:6]
        fl_meas = np.abs(d[6:9])
        pa_meas = d[9:12]
        tl_meas = np.abs(d[12:15])

        musculo = q_trial + list(skeleton_num)
        mtu_length = casadi_function["get_mtu_length"](musculo)

        fl_k = SX.sym(f"fl_{trial}", 3)
        pa_k = SX.sym(f"pa_{trial}", 3)
        tl_k = SX.sym(f"tl_{trial}", 3)
        w_k = vertcat(fl_k, pa_k, tl_k)
        w_vars.append(w_k)

        k_eq = vertcat(a_trial, mtu_length, p_sym)
        g_list.append(
            casadi_function["equilibrium_error_all_muscle"](w_k, k_eq)
        )

        neuromusculo = np.concatenate([a_trial, q_trial, skeleton_num])
        all_states = vertcat(SX(neuromusculo.tolist()), w_k)
        torque_sim = casadi_function["get_joint_moment"](all_states, p_sym)

        e_t = mesured_torque - torque_sim
        e_f = fl_meas - fl_k
        e_p = pa_meas - pa_k

        cost = cost + w_torque * e_t ** 2
        cost = cost + sum1(w_length * e_f ** 2)
        cost = cost + sum1(w_angle * e_p ** 2)

    w_full = vertcat(*w_vars)
    g_full = vertcat(*g_list)
    return {
        "w": w_full,
        "cost": cost,
        "g": g_full,
        "p_sym": p_sym,
        "n_params": 12,
        "n_states": 9 * n_trials,
        "n_constraints": 9 * n_trials,
    }


def compute_reduced_hessian_via_constraint_projection(augmented):
    """Compute the 12x12 Hessian of the cost on the manifold defined by g=0,
    evaluated analytically at the true parameters and the solved states.

    This uses the constraint Jacobian to project out the state directions.
    Given the augmented cost J(p, x) with equality constraint g(p, x) = 0
    (x being per-trial states), the reduced Hessian on p is obtained by:

        dx/dp = -(dg/dx)^{-1} (dg/dp)
        H_reduced = d^2 J / dp^2  + (dx/dp)^T @ d^2 J / dx^2 @ (dx/dp)
                    + cross terms ...

    This is more cleanly computed via Schur complement of the KKT system,
    but here we use CasADi automatic differentiation directly on the
    Lagrangian evaluated at x*(p).
    """
    raise NotImplementedError(
        "Reduced Hessian via KKT projection is deferred. "
        "The rootfinder approach in build_cost_function_params_only is used "
        "as the primary path."
    )


def hessian_analysis(
    data,
    true_params,
    skeleton_num,
    unknown_parameters,
    casadi_function,
    out_dir,
    method="rootfinder",
):
    """Main entry point.

    Parameters
    ----------
    data : (15, n_trials)
        Noise-free synthetic dataset.
    true_params : (12,)
        Ground-truth parameters used to generate `data`.
    method : str
        'rootfinder' -> reduced cost via equilibrium rootfinder (cleanest).
        'augmented'  -> fallback using KKT projection (not yet implemented).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    true_params = np.asarray(true_params, dtype=float).flatten()

    print("=" * 72)
    print("  Conditioning analysis: Hessian of the reduced cost")
    print("=" * 72)
    print(f"  method          = {method}")
    print(f"  true_params     = {true_params}")
    print()

    if method == "rootfinder":
        cost_func = build_cost_function_params_only(
            data, skeleton_num, unknown_parameters, casadi_function
        )

        # Sanity: cost at true params should be ~0
        f_true = float(cost_func(true_params))
        print(f"  cost at true_params : {f_true:.4e} (expected ~0)")

        # Build Hessian function symbolically
        p_sym = cost_func.sx_in(0)
        H_expr = ca.hessian(cost_func(p_sym), p_sym)[0]
        H_func = ca.Function("H", [p_sym], [H_expr])
        H = np.array(H_func(true_params))
    else:
        raise NotImplementedError(method)

    # Symmetrise (numerical)
    H = 0.5 * (H + H.T)

    # Eigen-analysis
    eigvals, eigvecs = np.linalg.eigh(H)  # ascending order
    print(f"\n  Eigenvalue range :")
    print(f"    lambda_min = {eigvals[0]:.4e}")
    print(f"    lambda_max = {eigvals[-1]:.4e}")
    if eigvals[0] > 0:
        kappa = eigvals[-1] / eigvals[0]
        print(f"    condition number kappa = {kappa:.4e}")
        if kappa > 1e8:
            print("    -> SEVERELY ILL-CONDITIONED (kappa > 1e8)")
        elif kappa > 1e4:
            print("    -> poorly conditioned (kappa > 1e4)")
    else:
        print(f"    negative eigenvalue detected ({eigvals[0]:.2e})")
        print("    -> not a local minimum, or numerical issue")

    # Report the sloppiest and stiffest directions
    _report_directions(eigvals, eigvecs, "sloppiest (worst-identified)",
                       ascending=True, n=4)
    _report_directions(eigvals, eigvecs, "stiffest (best-identified)",
                       ascending=False, n=3)

    # Save raw data
    df_eigvals = pd.DataFrame({
        "index": range(len(eigvals)),
        "eigenvalue": eigvals,
        "log10_eigenvalue": np.log10(np.abs(eigvals) + 1e-300),
    })
    df_eigvals.to_csv(out / "hessian_eigenvalues.csv", index=False)

    df_eigvecs = pd.DataFrame(
        eigvecs, index=PARAM_NAMES,
        columns=[f"v{i+1}" for i in range(len(eigvals))]
    )
    df_eigvecs.to_csv(out / "hessian_eigenvectors.csv")

    # Plots
    _plot_scree(eigvals, out / "hessian_eigenvalues.png")
    _plot_sloppy_directions(eigvals, eigvecs,
                            out / "hessian_sloppy_directions.png",
                            n_directions=3)
    _plot_eigenvector_heatmap(eigvecs, eigvals,
                              out / "hessian_eigenvectors_heatmap.png")

    print(f"\n  Outputs saved to: {out.resolve()}")

    return {
        "H": H,
        "eigvals": eigvals,
        "eigvecs": eigvecs,
        "cost_func": cost_func,
    }


def _report_directions(eigvals, eigvecs, label, ascending, n):
    idx = np.argsort(eigvals)
    if not ascending:
        idx = idx[::-1]
    idx = idx[:n]

    print(f"\n  {n} {label} directions :")
    for rank, i in enumerate(idx):
        lam = eigvals[i]
        v = eigvecs[:, i]
        # Dominant components
        loads = np.abs(v)
        order = np.argsort(loads)[::-1]
        top = order[:4]
        desc = ", ".join(
            f"{PARAM_NAMES[k]}({v[k]:+.2f})" for k in top
        )
        print(f"    [{rank+1}] lambda = {lam:+.3e}  |  {desc}")


def _plot_scree(eigvals, out_path):
    fig, ax = plt.subplots(figsize=(8, 5))
    pos = eigvals[eigvals > 0]
    neg = -eigvals[eigvals <= 0]

    x_pos = np.where(eigvals > 0)[0] + 1
    x_neg = np.where(eigvals <= 0)[0] + 1

    if len(pos) > 0:
        ax.semilogy(x_pos, pos, "o-", color="C0",
                    label="positive eigenvalues")
    if len(neg) > 0:
        ax.semilogy(x_neg, neg, "s", color="C3",
                    label="|negative eigenvalues|")

    ax.set_xlabel("eigenvalue index (ascending)")
    ax.set_ylabel("|eigenvalue|")
    ax.set_title("Hessian scree plot (cost function, 12 parameters)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_sloppy_directions(eigvals, eigvecs, out_path, n_directions=3):
    idx = np.argsort(eigvals)[:n_directions]
    fig, axes = plt.subplots(1, n_directions, figsize=(4 * n_directions, 5),
                             sharey=True)
    if n_directions == 1:
        axes = [axes]
    for ax, i in zip(axes, idx):
        v = eigvecs[:, i]
        colors = ["C0" if x >= 0 else "C3" for x in v]
        ax.barh(range(12), v, color=colors, alpha=0.7, edgecolor="k")
        ax.set_yticks(range(12))
        ax.set_yticklabels(PARAM_NAMES)
        ax.axvline(0, color="k", linewidth=0.5)
        ax.set_title(f"lambda = {eigvals[i]:.2e}")
        ax.grid(alpha=0.3, axis="x")
    fig.suptitle(f"{n_directions} sloppiest eigen-directions (parameter loadings)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_eigenvector_heatmap(eigvecs, eigvals, out_path):
    # Sort columns by eigenvalue ascending (sloppiest first)
    order = np.argsort(eigvals)
    V = eigvecs[:, order]
    lam_sorted = eigvals[order]

    fig, ax = plt.subplots(figsize=(11, 7))
    im = ax.imshow(V, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_yticks(range(12))
    ax.set_yticklabels(PARAM_NAMES)
    ax.set_xticks(range(len(lam_sorted)))
    ax.set_xticklabels(
        [f"{l:.1e}" for l in lam_sorted],
        rotation=45, ha="right", fontsize=8
    )
    ax.set_xlabel("eigenvalues (ascending -> sloppiest on the left)")
    ax.set_title("Hessian eigenvectors (columns, rows = parameters)")
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            if abs(V[i, j]) > 0.25:
                ax.text(j, i, f"{V[i, j]:+.2f}",
                        ha="center", va="center", fontsize=7,
                        color="white" if abs(V[i, j]) > 0.6 else "black")
    fig.colorbar(im, ax=ax, label="loading")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

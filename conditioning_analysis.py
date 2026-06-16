"""
Conditioning analysis (local identifiability via the Hessian / FIM).

[PATCH param_config]
--------------------
On construit la matrice de sensibilité (jacobienne du résidu par rapport aux
paramètres 'sym') et la Hessienne de Gauss-Newton H ≈ Jᵀ J sur les n_opt
paramètres optimisés (déduits de param_index). L'analyse spectrale de H donne :

  - les valeurs propres : une valeur propre quasi nulle = direction
    (combinaison de paramètres) non identifiable localement ;
  - le conditionnement κ = λ_max / λ_min : κ grand (> 1e8 typiquement) =
    problème mal posé ;
  - le vecteur propre associé à la plus petite valeur propre : indique
    QUELS paramètres se compensent.

La jacobienne est calculée par différences finies autour de true_params,
en réutilisant le même forward model (rootfinder d'équilibre) que les autres
analyses. Cohérent avec sobol_analysis.simulate_model mais renvoie ici le
vecteur complet des résidus (et non un scalaire).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

from param_layout import Layout


def _forward_residual(params, data, skeleton_num, casadi_function):
    """Vecteur des résidus (torque, fibre, pennation) empilés sur les essais.

    Renvoie un vecteur 1D ; NaN si le rootfinder échoue sur un essai.
    """
    from casadi import SX, vertcat
    import casadi as ca

    n_opt = len(params)
    n_trials = data.shape[1]

    x_rf = SX.sym("x_rf", 9)
    a_rf = SX.sym("a_rf", 3)
    mtu_rf = SX.sym("mtu_rf", 3)
    p_rf = SX.sym("p_rf", n_opt)
    eq = casadi_function["equilibrium_error_all_muscle"](
        x_rf, vertcat(a_rf, mtu_rf, p_rf)
    )
    rf = ca.rootfinder(
        "rf_cond", "newton",
        {"x": x_rf, "p": vertcat(a_rf, mtu_rf, p_rf), "g": eq},
        {"error_on_fail": False},
    )

    # Pondérations cohérentes avec _evaluate_cost_after_fit (poids = 1/sigma)
    w_torque = 1.0
    w_length = np.sqrt(0.005)
    w_angle = np.sqrt((1.0 / 180.0) * np.pi)

    res = []
    for trial in range(n_trials):
        d = data[:, trial]
        q_trial = [0, 0, 0, 0] + list(d[1:3])
        a_trial = d[3:6]
        fl_meas = np.abs(d[6:9])
        pa_meas = d[9:12]
        tl_meas = np.abs(d[12:15])
        mtu = casadi_function["get_mtu_length"](q_trial + list(skeleton_num))

        x0 = np.concatenate([fl_meas, pa_meas, tl_meas])
        try:
            states = np.array(rf(x0, vertcat(a_trial, mtu, params))).flatten()
            if np.any(np.isnan(states)):
                return None
        except Exception:
            return None

        neuromusculo = np.concatenate([a_trial, q_trial, skeleton_num])
        all_states_num = np.concatenate([neuromusculo, states])
        torque_sim = float(
            casadi_function["get_joint_moment"](all_states_num, params)
        )

        res.append(w_torque * (d[0] - torque_sim))
        res.extend(list(w_length * (fl_meas - states[0:3])))
        res.extend(list(w_angle * (pa_meas - states[3:6])))

    return np.asarray(res, dtype=float)


def _jacobian_fd(params, data, skeleton_num, casadi_function, rel_step=1e-6):
    """Jacobienne du résidu par rapport aux n_opt paramètres (diff. finies)."""
    params = np.asarray(params, dtype=float)
    r0 = _forward_residual(params, data, skeleton_num, casadi_function)
    if r0 is None:
        raise RuntimeError("Forward model failed at true_params.")
    m, n = r0.size, params.size
    J = np.zeros((m, n))
    for j in range(n):
        h = rel_step * max(abs(params[j]), 1.0)
        pj = params.copy(); pj[j] += h
        rj = _forward_residual(pj, data, skeleton_num, casadi_function)
        if rj is None:
            J[:, j] = np.nan
        else:
            J[:, j] = (rj - r0) / h
    return J


def hessian_analysis(
    data,
    true_params,
    skeleton_num,
    unknown_parameters,
    casadi_function,
    param_index,
    out_dir,
):
    """Analyse de conditionnement local au point true_params."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    layout = Layout.from_param_index(param_index)
    true_params = np.asarray(true_params, dtype=float).flatten()
    assert true_params.shape[0] == layout.n_opt, (
        f"true_params ({true_params.shape[0]}) != n_opt={layout.n_opt}"
    )

    print("=" * 72)
    print("  Conditioning analysis (Gauss-Newton Hessian / FIM)")
    print("=" * 72)
    print(f"  {layout.summary()}")

    J = _jacobian_fd(true_params, data, skeleton_num, casadi_function)
    H = J.T @ J                       # Hessien de Gauss-Newton ≈ FIM
    eigvals, eigvecs = np.linalg.eigh(H)   # H symétrique
    eigvals = np.clip(eigvals, 0, None)    # nettoyer le bruit numérique

    lam_max = eigvals[-1]
    lam_min = eigvals[0]
    kappa = lam_max / lam_min if lam_min > 0 else np.inf

    print(f"  eig(H) min/max   = {lam_min:.3e} / {lam_max:.3e}")
    print(f"  cond(H) = κ      = {kappa:.3e}  (log10 = "
          f"{np.log10(kappa) if np.isfinite(kappa) else np.inf:.2f})")

    # Direction la moins identifiable : vecteur propre de lam_min
    weak_dir = eigvecs[:, 0]
    print("\n  Direction la moins contrainte (|composantes| triées) :")
    order = np.argsort(np.abs(weak_dir))[::-1]
    for i in order:
        print(f"    {layout.names[i]:<10} {weak_dir[i]:+.3f}")

    # Sauvegardes
    pd.DataFrame({
        "eigenvalue": eigvals[::-1],
        "rank": np.arange(1, layout.n_opt + 1),
    }).to_csv(out / "conditioning_eigenvalues.csv", index=False)

    pd.DataFrame(eigvecs, index=layout.names,
                 columns=[f"v{i+1}" for i in range(layout.n_opt)]
                 ).to_csv(out / "conditioning_eigenvectors.csv")

    _plot_spectrum(eigvals, kappa, out / "conditioning_spectrum.png")
    _plot_weak_direction(weak_dir, layout, out / "conditioning_weak_dir.png")

    print(f"\n  Outputs saved to: {out.resolve()}")
    return {
        "H": H, "eigenvalues": eigvals, "eigenvectors": eigvecs,
        "cond": kappa, "weak_direction": weak_dir, "names": layout.names,
    }


def _plot_spectrum(eigvals, kappa, out_path):
    fig, ax = plt.subplots(figsize=(8, 5))
    ev = np.sort(eigvals)[::-1]
    ax.semilogy(range(1, len(ev) + 1), np.maximum(ev, 1e-30), "o-")
    ax.set_xlabel("eigenvalue rank")
    ax.set_ylabel("eigenvalue (log)")
    title = "Hessian spectrum"
    if np.isfinite(kappa):
        title += f"  (κ = {kappa:.1e})"
    ax.set_title(title)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_weak_direction(weak_dir, layout, out_path):
    fig, ax = plt.subplots(figsize=(max(8, layout.n_opt * 0.7), 5))
    ax.bar(range(layout.n_opt), weak_dir, color="C3", alpha=0.8)
    ax.set_xticks(range(layout.n_opt))
    ax.set_xticklabels(layout.names, rotation=45, ha="right")
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_ylabel("component")
    ax.set_title("Least-identifiable direction (smallest eigenvalue)")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

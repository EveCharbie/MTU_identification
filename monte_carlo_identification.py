"""
Monte-Carlo analysis of the noise robustness of the muscle-tendon
parameter identification NLP.

[PATCH param_config]
--------------------
Le vecteur estimé n'est plus de taille 12 : il correspond aux paramètres
'sym' (déduits de param_index / get_initial_guess). Toutes les statistiques
et figures s'adaptent à n_opt. Le scatter φo–Fom n'est tracé que pour les
muscles dont φo ET Fom sont tous deux 'sym'.

Le bruit est injecté sur les DONNÉES (15 lignes), indépendamment de
param_config : add_gaussian_noise est inchangé.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import time
import contextlib
import io

from param_layout import Layout, PROPERTY_LABEL


CONFIG = {
    "n_mc": 60,
    "sigma_torque": 0.1,       # N.m
    "sigma_fiber": 0.000,      # m
    "sigma_pennation": np.deg2rad(0),  # rad
    "sigma_tendon": 0.00,     # m
    "seed": 42,
    "output_dir": "mc_results",
    "silent_ipopt": True,
    "random_x0": True,
    "x0_noise_relative": 0.2,
}


# ============================================================
# Monte-Carlo core
# ============================================================

def add_gaussian_noise(data_clean, cfg, rng):
    """Inject Gaussian noise on measured channels of `data` (inchangé).

    data layout (15 x n_trials): voir docstring d'origine.
    Indépendant de param_config (porte sur les mesures, pas les paramètres).
    """
    n_trials = data_clean.shape[1]
    d = data_clean.copy()

    d[0, :]     += rng.normal(0.0, cfg["sigma_torque"],    size=n_trials)
    d[6:9, :]   += rng.normal(0.0, cfg["sigma_fiber"],     size=(3, n_trials))
    d[9:12, :]  += rng.normal(0.0, cfg["sigma_pennation"], size=(3, n_trials))
    d[12:15, :] += rng.normal(0.0, cfg["sigma_tendon"],    size=(3, n_trials))

    return d


def perturb_initial_guess(x0, rel_noise, rng):
    x0 = np.asarray(x0, dtype=float)
    factor = 1.0 + rng.normal(0.0, rel_noise, size=x0.shape)
    return x0 * factor


def run_monte_carlo(
    data_clean,
    initial_guess,
    lower_band,
    upper_band,
    skeleton_num,
    muscle_tendon_parameters_num,
    unknown_parameters,
    casadi_function,
    optimization_nlp,
    param_index,
    cfg=None,
):
    """Campagne Monte-Carlo complète. Renvoie (estimates, meta)."""
    cfg = cfg or CONFIG
    rng = np.random.default_rng(cfg["seed"])
    n_mc = cfg["n_mc"]

    estimates = []
    meta = []
    t0_global = time.time()

    for i in range(n_mc):
        data_noisy = add_gaussian_noise(data_clean, cfg, rng)

        if cfg["random_x0"]:
            x0_i = perturb_initial_guess(
                initial_guess, cfg["x0_noise_relative"], rng
            )
            x0_i = np.clip(x0_i, lower_band, upper_band)
        else:
            x0_i = np.asarray(initial_guess, dtype=float)

        t0 = time.time()
        success = True
        params_est = None

        try:
            if cfg["silent_ipopt"]:
                with contextlib.redirect_stdout(io.StringIO()):
                    params_est = optimization_nlp(
                        data_noisy, x0_i, lower_band, upper_band,
                        skeleton_num, muscle_tendon_parameters_num,
                        unknown_parameters, casadi_function,
                        param_index,
                    )
            else:
                params_est = optimization_nlp(
                    data_noisy, x0_i, lower_band, upper_band,
                    skeleton_num, muscle_tendon_parameters_num,
                    unknown_parameters, casadi_function,
                    param_index,
                )
        except Exception as e:
            success = False
            print(f"  MC {i+1:>3}/{n_mc}: FAILED — {type(e).__name__}: {e}")

        wall = time.time() - t0

        if success and params_est is not None and not np.any(np.isnan(params_est)):
            estimates.append(np.asarray(params_est).flatten())
            meta.append({"run": i, "success": True, "wall_s": wall})
            print(f"  MC {i+1:>3}/{n_mc}: OK  ({wall:5.2f} s)")
        else:
            meta.append({"run": i, "success": False, "wall_s": wall})

    total_wall = time.time() - t0_global
    n_success = len(estimates)
    print(
        f"\nMonte-Carlo complete: {n_success}/{n_mc} successful "
        f"in {total_wall:.1f} s (mean {total_wall/max(n_mc,1):.2f} s/run)"
    )
    return np.asarray(estimates), pd.DataFrame(meta)


# ============================================================
# Statistical analysis
# ============================================================

def summarise_results(estimates, true_params, layout):
    true_params = np.asarray(true_params).flatten()
    mean = estimates.mean(axis=0)
    std = estimates.std(axis=0, ddof=1)
    bias = mean - true_params
    denom = np.where(np.abs(true_params) > 0, np.abs(true_params), 1.0)
    rel_bias = bias / denom * 100
    rel_std = std / denom * 100

    n = estimates.shape[0]
    ci_half = 1.96 * std / np.sqrt(n)

    return pd.DataFrame({
        "parameter": layout.names,
        "unit": layout.units,
        "true": true_params,
        "mean": mean,
        "bias": bias,
        "rel_bias_%": rel_bias,
        "std": std,
        "rel_std_%": rel_std,
        "ci95_low": mean - ci_half,
        "ci95_high": mean + ci_half,
    })


def print_summary_table(df):
    print("\n" + "=" * 100)
    print(f"{'Parameter':<8} {'Unit':<5} {'True':>12} {'Mean':>12} "
          f"{'Bias':>12} {'Std':>12} {'RelStd%':>9}")
    print("-" * 100)
    for _, row in df.iterrows():
        print(f"{row['parameter']:<8} {row['unit']:<5} "
              f"{row['true']:>12.4g} {row['mean']:>12.4g} "
              f"{row['bias']:>12.3g} "
              f"{row['std']:>12.3g} {row['rel_std_%']:>9.2f}")
    print("=" * 100)


# ============================================================
# Plots
# ============================================================

def plot_parameter_distributions(estimates, true_params, layout, out_path):
    """Boxplots de l'erreur relative, groupés par propriété 'sym' présente."""
    true_params = np.asarray(true_params).flatten()
    denom = np.where(np.abs(true_params) > 0, np.abs(true_params), 1.0)
    rel_err = (estimates - true_params) / denom * 100

    groups = layout.group_slices()  # [(label, i0, i1, unit), ...]
    n_groups = len(groups)
    ncols = min(2, n_groups)
    nrows = int(np.ceil(n_groups / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, (label, i0, i1, unit) in zip(axes, groups):
        sub = rel_err[:, i0:i1]
        bp = ax.boxplot(sub, labels=[f"M{k+1}" for k in range(i1 - i0)],
                        showmeans=True, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("#b3cde0")
        ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
        ax.set_title(f"Relative error on {label}")
        ax.set_ylabel("Error (%)")
        ax.grid(alpha=0.3)

    for ax in axes[n_groups:]:
        ax.set_visible(False)

    fig.suptitle("Monte-Carlo: parameter identification relative error")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_error_correlation(estimates, true_params, layout, out_path):
    """Heatmap des corrélations entre erreurs (taille n_opt x n_opt)."""
    true_params = np.asarray(true_params).flatten()
    err = estimates - true_params
    corr = np.corrcoef(err.T)
    D = layout.n_opt

    fig, ax = plt.subplots(figsize=(max(8, D * 0.7), max(7, D * 0.65)))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(D))
    ax.set_yticks(range(D))
    ax.set_xticklabels(layout.names, rotation=45, ha="right")
    ax.set_yticklabels(layout.names)

    for i in range(D):
        for j in range(D):
            txt_color = "white" if abs(corr[i, j]) > 0.5 else "black"
            ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center",
                    color=txt_color, fontsize=8)

    fig.colorbar(im, ax=ax, label="Pearson correlation")
    ax.set_title("Correlation of parameter estimation errors")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_scatter_phio_Fom(estimates, true_params, layout, out_path):
    """Scatter erreur(φo) vs erreur(Fom) par muscle.

    Tracé uniquement pour les muscles où φo ET Fom sont 'sym'. Si aucun
    muscle ne satisfait la condition, la figure n'est pas produite.
    """
    true_params = np.asarray(true_params).flatten()
    err = estimates - true_params

    # Muscles ayant les deux propriétés dans le vecteur optimisé
    pairs = []  # (muscle_id, pos_phio, pos_Fom)
    for m in range(1, 4):
        n_phi, n_fom = f"phio_{m}", f"Fom_{m}"
        if n_phi in layout.index_map and n_fom in layout.index_map:
            pairs.append((m, layout.position(n_phi), layout.position(n_fom)))

    if not pairs:
        print("  [scatter φo–Fom] aucun muscle avec φo ET Fom 'sym' — figure ignorée.")
        return False

    fig, axes = plt.subplots(1, len(pairs), figsize=(5 * len(pairs), 4.5))
    axes = np.atleast_1d(axes)
    for ax, (m, pp, pf) in zip(axes, pairs):
        ax.scatter(err[:, pp], err[:, pf], alpha=0.6, s=30,
                   edgecolor="k", linewidth=0.5)
        ax.axhline(0, color="k", linewidth=0.5, alpha=0.5)
        ax.axvline(0, color="k", linewidth=0.5, alpha=0.5)
        r = np.corrcoef(err[:, pp], err[:, pf])[0, 1]
        ax.set_title(f"Muscle {m}  (r = {r:+.2f})")
        ax.set_xlabel("φo error (rad)")
        ax.set_ylabel("Fom error (N)")
        ax.grid(alpha=0.3)

    fig.suptitle("φo — Fom identifiability coupling (per muscle)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ============================================================
# Main driver
# ============================================================

def main(
    data_clean,
    initial_guess,
    lower_band,
    upper_band,
    skeleton_num,
    muscle_tendon_parameters_num,
    unknown_parameters,
    casadi_function,
    optimization_nlp,
    param_index,
    cfg=None,
):
    """Orchestre l'analyse Monte-Carlo.

    param_index : dict (sortie de get_initial_guess) — fixe taille/ordre/NLP.
    Toutes les array (12,) d'origine sont désormais de taille n_opt.
    """
    cfg_ = dict(CONFIG)
    if cfg:
        cfg_.update(cfg)

    out = Path(cfg_["output_dir"])
    out.mkdir(parents=True, exist_ok=True)

    layout = Layout.from_param_index(param_index)

    print("=" * 72)
    print("  Monte-Carlo noise-robustness analysis")
    print("=" * 72)
    print(f"  {layout.summary()}")
    print(f"  n_mc           = {cfg_['n_mc']}")
    print(f"  sigma_torque   = {cfg_['sigma_torque']} N.m")
    print(f"  sigma_fiber    = {cfg_['sigma_fiber']*1000:.2f} mm")
    print(f"  sigma_pennat.  = {np.rad2deg(cfg_['sigma_pennation']):.2f} deg")
    print(f"  sigma_tendon   = {cfg_['sigma_tendon']*1000:.2f} mm")
    print(f"  random_x0      = {cfg_['random_x0']}")
    print(f"  output_dir     = {out.resolve()}")
    print("=" * 72)

    estimates, meta = run_monte_carlo(
        data_clean, initial_guess, lower_band, upper_band,
        skeleton_num, muscle_tendon_parameters_num,
        unknown_parameters, casadi_function,
        optimization_nlp, param_index, cfg_,
    )

    if estimates.size == 0:
        print("No successful run — aborting analysis.")
        return None

    summary = summarise_results(estimates, muscle_tendon_parameters_num, layout)
    print_summary_table(summary)

    df_raw = pd.DataFrame(estimates, columns=layout.names)
    df_raw.to_csv(out / "mc_estimates.csv", index=False)
    summary.to_csv(out / "mc_summary.csv", index=False)
    meta.to_csv(out / "mc_meta.csv", index=False)

    plot_parameter_distributions(
        estimates, muscle_tendon_parameters_num, layout,
        out / "mc_distributions.png"
    )
    plot_error_correlation(
        estimates, muscle_tendon_parameters_num, layout,
        out / "mc_correlations.png"
    )
    plot_scatter_phio_Fom(
        estimates, muscle_tendon_parameters_num, layout,
        out / "mc_phio_Fom_coupling.png"
    )

    print(f"\nAll outputs saved to: {out.resolve()}")
    return {"estimates": estimates, "summary": summary, "meta": meta}

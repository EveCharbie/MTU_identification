"""
Monte-Carlo analysis of the noise robustness of the muscle-tendon
parameter identification NLP.

For each Monte-Carlo run, Gaussian noise is injected into the synthetic
measurements (torque, fiber length, pennation angle, tendon length), the
NLP is re-solved, and the estimated parameters are stored.

The script outputs:
- A summary table (true value, mean, bias, std, relative std, 95% CI)
- Distribution plots of the estimated parameters
- A correlation heatmap between parameter errors (reveals φo/Fom coupling)
- A CSV with all Monte-Carlo realisations for further analysis

Dependencies: numpy, pandas, matplotlib, scipy (only used for chi2 test)
Requires `optimization_nlp` and all its arguments to be available in scope.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import time
import contextlib
import io


# ============================================================
# CONFIGURATION — adapt to your experimental setup
# ============================================================

CONFIG = {
    # Number of Monte-Carlo realisations
    # 30 is a minimum; 100+ is recommended for stable statistics
    "n_mc": 50,

    # Noise standard deviations (adapt to your instruments)
    # Default values are typical for ultrasound + isokinetic dynamometer
    "sigma_torque": 0.5,       # N.m
    "sigma_fiber": 0.002,      # m  (2 mm — ultrasound fiber length)
    "sigma_pennation": np.deg2rad(2.0),  # rad (2°)
    "sigma_tendon": 0.003,     # m  (3 mm — ultrasound tendon length)

    # Reproducibility
    "seed": 42,

    # Output directory (will be created if missing)
    "output_dir": "mc_results",

    # Silence the verbose IPOPT output inside the MC loop
    "silent_ipopt": True,

    # Randomise the initial guess at each MC run (robust check)
    # If False, use the same `initial_guess` each time.
    "random_x0": True,
    "x0_noise_relative": 0.10,  # 10% relative noise on x0 if random_x0=True
}

PARAM_NAMES = [
    "lom_1", "lom_2", "lom_3",
    "phio_1", "phio_2", "phio_3",
    "Fom_1", "Fom_2", "Fom_3",
    "lst_1", "lst_2", "lst_3",
]
PARAM_UNITS = ["m"] * 3 + ["rad"] * 3 + ["N"] * 3 + ["m"] * 3
PARAM_GROUPS = ["ℓom"] * 3 + ["φo"] * 3 + ["Fom"] * 3 + ["ℓst"] * 3


# ============================================================
# Monte-Carlo core
# ============================================================

def add_gaussian_noise(data_clean, cfg, rng):
    """Inject Gaussian noise on measured channels of `data`.

    data layout (rows, 15 x n_trials):
        [0]     torque
        [1:3]   joint angles       (NOT noised — these are inputs, not measures)
        [3:6]   activations        (NOT noised — assumed control inputs)
        [6:9]   fiber lengths
        [9:12]  pennation angles
        [12:15] tendon lengths
    """
    n_trials = data_clean.shape[1]
    d = data_clean.copy()

    d[0, :]     += rng.normal(0.0, cfg["sigma_torque"],    size=n_trials)
    d[6:9, :]   += rng.normal(0.0, cfg["sigma_fiber"],     size=(3, n_trials))
    d[9:12, :]  += rng.normal(0.0, cfg["sigma_pennation"], size=(3, n_trials))
    d[12:15, :] += rng.normal(0.0, cfg["sigma_tendon"],    size=(3, n_trials))

    return d


def perturb_initial_guess(x0, rel_noise, rng):
    """Multiplicative Gaussian perturbation on the initial guess."""
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
    cfg=None,
):
    """Run the full Monte-Carlo campaign.

    Returns
    -------
    results : np.ndarray, shape (n_success, 12)
        Estimated parameters for each successful MC run.
    mc_meta : pd.DataFrame
        Per-run metadata: success flag, wall time, seed.
    """
    cfg = cfg or CONFIG
    rng = np.random.default_rng(cfg["seed"])
    n_mc = cfg["n_mc"]

    estimates = []
    meta = []

    t0_global = time.time()

    for i in range(n_mc):
        # Fresh noise realisation
        data_noisy = add_gaussian_noise(data_clean, cfg, rng)

        # Fresh initial guess (optional)
        if cfg["random_x0"]:
            x0_i = perturb_initial_guess(
                initial_guess, cfg["x0_noise_relative"], rng
            )
            # Clip to bounds to avoid starting infeasible
            x0_i = np.clip(x0_i, lower_band, upper_band)
        else:
            x0_i = np.asarray(initial_guess, dtype=float)

        t0 = time.time()
        success = True
        params_est = None

        try:
            if cfg["silent_ipopt"]:
                # Redirect stdout to silence IPOPT's verbose banner
                with contextlib.redirect_stdout(io.StringIO()):
                    params_est = optimization_nlp(
                        data_noisy, x0_i, lower_band, upper_band,
                        skeleton_num, muscle_tendon_parameters_num,
                        unknown_parameters, casadi_function,
                    )
            else:
                params_est = optimization_nlp(
                    data_noisy, x0_i, lower_band, upper_band,
                    skeleton_num, muscle_tendon_parameters_num,
                    unknown_parameters, casadi_function,
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

def summarise_results(estimates, true_params):
    """Return a DataFrame with per-parameter statistics."""
    true_params = np.asarray(true_params).flatten()
    mean = estimates.mean(axis=0)
    std = estimates.std(axis=0, ddof=1)
    bias = mean - true_params
    rel_bias = bias / np.where(np.abs(true_params) > 0, np.abs(true_params), 1.0) * 100
    rel_std = std / np.where(np.abs(true_params) > 0, np.abs(true_params), 1.0) * 100

    # 95% confidence interval of the mean estimator (normal approx)
    n = estimates.shape[0]
    ci_half = 1.96 * std / np.sqrt(n)

    df = pd.DataFrame({
        "parameter": PARAM_NAMES,
        "unit": PARAM_UNITS,
        "true": true_params,
        "mean": mean,
        "bias": bias,
        "rel_bias_%": rel_bias,
        "std": std,
        "rel_std_%": rel_std,
        "ci95_low": mean - ci_half,
        "ci95_high": mean + ci_half,
    })
    return df


def print_summary_table(df):
    """Pretty-print the summary to the terminal."""
    print("\n" + "=" * 100)
    print(f"{'Parameter':<8} {'Unit':<5} {'True':>12} {'Mean':>12} "
          f"{'Bias':>12} {'RelBias%':>10} {'Std':>12} {'RelStd%':>9}")
    print("=" * 100)
    for _, row in df.iterrows():
        print(f"{row['parameter']:<8} {row['unit']:<5} "
              f"{row['true']:>12.4g} {row['mean']:>12.4g} "
              f"{row['bias']:>+12.3g} {row['rel_bias_%']:>+10.2f} "
              f"{row['std']:>12.3g} {row['rel_std_%']:>9.2f}")
    print("=" * 100)


# ============================================================
# Plots
# ============================================================

def plot_parameter_distributions(estimates, true_params, out_path):
    """Boxplots of relative estimation error per parameter, grouped by type."""
    true_params = np.asarray(true_params).flatten()
    rel_err = (estimates - true_params) / np.abs(true_params) * 100  # %

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    groups = [("ℓom", 0, 3, "%"), ("φo", 3, 6, "%"),
              ("Fom", 6, 9, "%"), ("ℓst", 9, 12, "%")]

    for ax, (name, i0, i1, unit) in zip(axes.flat, groups):
        sub = rel_err[:, i0:i1]
        bp = ax.boxplot(sub, labels=[f"M{k+1}" for k in range(i1 - i0)],
                        showmeans=True, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("#b3cde0")
        ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
        ax.set_title(f"Relative error on {name}")
        ax.set_ylabel(f"Error ({unit})")
        ax.grid(alpha=0.3)

    fig.suptitle("Monte-Carlo: parameter identification relative error")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_error_correlation(estimates, true_params, out_path):
    """Heatmap of Pearson correlations between parameter errors.

    A strong correlation between φo and Fom for the same muscle is the
    classical identifiability signature of the Hill model.
    """
    true_params = np.asarray(true_params).flatten()
    err = estimates - true_params
    corr = np.corrcoef(err.T)

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(12))
    ax.set_yticks(range(12))
    ax.set_xticklabels(PARAM_NAMES, rotation=45, ha="right")
    ax.set_yticklabels(PARAM_NAMES)

    # Annotate
    for i in range(12):
        for j in range(12):
            txt_color = "white" if abs(corr[i, j]) > 0.5 else "black"
            ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center",
                    color=txt_color, fontsize=8)

    fig.colorbar(im, ax=ax, label="Pearson correlation")
    ax.set_title("Correlation of parameter estimation errors")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_scatter_phio_Fom(estimates, true_params, out_path):
    """Scatter plot of φo error vs Fom error for each muscle (per-muscle
    panel). Reveals the compensation direction.
    """
    true_params = np.asarray(true_params).flatten()
    err = estimates - true_params

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for m, ax in enumerate(axes):
        ax.scatter(err[:, 3 + m], err[:, 6 + m], alpha=0.6, s=30,
                   edgecolor="k", linewidth=0.5)
        ax.axhline(0, color="k", linewidth=0.5, alpha=0.5)
        ax.axvline(0, color="k", linewidth=0.5, alpha=0.5)

        # Correlation annotation
        r = np.corrcoef(err[:, 3 + m], err[:, 6 + m])[0, 1]
        ax.set_title(f"Muscle {m+1}  (r = {r:+.2f})")
        ax.set_xlabel("φo error (rad)")
        ax.set_ylabel("Fom error (N)")
        ax.grid(alpha=0.3)

    fig.suptitle("φo — Fom identifiability coupling (per muscle)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


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
    cfg=None,
):
    """Orchestrate the Monte-Carlo analysis.

    Parameters
    ----------
    data_clean : np.ndarray (15, n_trials)
        Noise-free reference dataset (e.g. synthetic ground truth).
    initial_guess, lower_band, upper_band : array-like (12,)
        As in optimization_nlp.
    skeleton_num, unknown_parameters, casadi_function : as in optimization_nlp.
    muscle_tendon_parameters_num : array-like (12,)
        Ground-truth parameters used to generate `data_clean`.
    optimization_nlp : callable
        The NLP wrapper function (expected signature matches yours).
    cfg : dict, optional
        Configuration dict (see CONFIG above). Merged with defaults.
    """
    cfg_ = dict(CONFIG)
    if cfg:
        cfg_.update(cfg)

    out = Path(cfg_["output_dir"])
    out.mkdir(parents=True, exist_ok=True)

    # --- Banner
    print("=" * 72)
    print("  Monte-Carlo noise-robustness analysis")
    print("=" * 72)
    print(f"  n_mc           = {cfg_['n_mc']}")
    print(f"  sigma_torque   = {cfg_['sigma_torque']} N.m")
    print(f"  sigma_fiber    = {cfg_['sigma_fiber']*1000:.2f} mm")
    print(f"  sigma_pennat.  = {np.rad2deg(cfg_['sigma_pennation']):.2f} deg")
    print(f"  sigma_tendon   = {cfg_['sigma_tendon']*1000:.2f} mm")
    print(f"  random_x0      = {cfg_['random_x0']}")
    print(f"  output_dir     = {out.resolve()}")
    print("=" * 72)

    # --- Run
    estimates, meta = run_monte_carlo(
        data_clean, initial_guess, lower_band, upper_band,
        skeleton_num, muscle_tendon_parameters_num,
        unknown_parameters, casadi_function,
        optimization_nlp, cfg_,
    )

    if estimates.size == 0:
        print("No successful run — aborting analysis.")
        return None

    # --- Statistics
    summary = summarise_results(estimates, muscle_tendon_parameters_num)
    print_summary_table(summary)

    # --- Save raw data
    df_raw = pd.DataFrame(estimates, columns=PARAM_NAMES)
    df_raw.to_csv(out / "mc_estimates.csv", index=False)
    summary.to_csv(out / "mc_summary.csv", index=False)
    meta.to_csv(out / "mc_meta.csv", index=False)

    # --- Plots
    plot_parameter_distributions(
        estimates, muscle_tendon_parameters_num,
        out / "mc_distributions.png"
    )
    plot_error_correlation(
        estimates, muscle_tendon_parameters_num,
        out / "mc_correlations.png"
    )
    plot_scatter_phio_Fom(
        estimates, muscle_tendon_parameters_num,
        out / "mc_phio_Fom_coupling.png"
    )

    print(f"\nAll outputs saved to: {out.resolve()}")
    print(f"  - mc_estimates.csv  (raw {estimates.shape[0]} x 12 realisations)")
    print(f"  - mc_summary.csv    (per-parameter statistics)")
    print(f"  - mc_meta.csv       (per-run metadata)")
    print(f"  - mc_distributions.png")
    print(f"  - mc_correlations.png")
    print(f"  - mc_phio_Fom_coupling.png")

    return {"estimates": estimates, "summary": summary, "meta": meta}


# ============================================================
# Usage example — adapt to how you load your setup
# ============================================================

if __name__ == "__main__":
    # This block is illustrative. In practice you likely have a main
    # script that already builds `data`, `initial_guess`, the CasADi
    # functions, etc. Import `main` from this module and pass them in.
    #
    # from your_setup_module import (
    #     data, initial_guess, lower_band, upper_band, skeleton_num,
    #     muscle_tendon_parameters_num, unknown_parameters, casadi_function
    # )
    # from your_nlp_module import optimization_nlp
    #
    # results = main(
    #     data_clean=data,
    #     initial_guess=initial_guess,
    #     lower_band=lower_band,
    #     upper_band=upper_band,
    #     skeleton_num=skeleton_num,
    #     muscle_tendon_parameters_num=muscle_tendon_parameters_num,
    #     unknown_parameters=unknown_parameters,
    #     casadi_function=casadi_function,
    #     optimization_nlp=optimization_nlp,
    #     cfg={"n_mc": 50},   # override defaults if needed
    # )
    raise SystemExit(
        "Import `main` from this module and call it from your setup script.\n"
        "See the commented example at the bottom of the file."
    )
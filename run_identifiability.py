"""
Orchestrator: run all identifiability analyses in sequence.

[PATCH param_config]
--------------------
setup doit maintenant fournir `param_index` (sortie de get_initial_guess) et
`param_config`. Toutes les analyses en dérivent le layout (n_opt, ordre, noms)
et ne travaillent que sur les paramètres 'sym'. Les bornes / true_params /
unknown_parameters sont de taille n_opt.

Analyses : conditioning (Hessien/FIM local), profile likelihood (identifiabilité
pratique), Sobol (sensibilité globale).

Usage:
    python run_identifiability.py
    # ou
    from run_identifiability import run_all
    run_all(setup, analyses=["conditioning", "profile"])
"""

from __future__ import annotations

from pathlib import Path
import argparse

from conditioning_analysis import hessian_analysis
from profile_likelihood import run_profile_likelihood
from sobol_analysis import run_sobol_analysis
from param_layout import Layout


OUTPUT_ROOT = Path("identifiability_results")


def load_setup():
    """Build the arguments needed by each analysis.

    Doit renvoyer un dict avec EXACTEMENT ces clés (tailles = n_opt) :
      - data               : np.ndarray (15, n_trials), NOISE-FREE
      - initial_guess      : np.ndarray (n_opt,)
      - lower_band         : np.ndarray (n_opt,)
      - upper_band         : np.ndarray (n_opt,)
      - skeleton_num       : iterable, géométrie squelette
      - true_params        : np.ndarray (n_opt,) — valeurs vraies des 'sym'
      - unknown_parameters : CasADi SX (n_opt,) symbolique
      - casadi_function    : dict de CasADi Functions ('fixed' déjà en dur)
      - optimization_nlp   : callable, wrapper NLP (signature avec param_index)
      - param_index        : dict, sortie de useful.get_initial_guess
      - param_config       : dict {prop: 'sym'|'fixed'} (pour traçabilité)

    Exemple d'obtention de initial_guess/bornes/param_index :
        initial_guess, upper_band, lower_band, param_index = \
            useful.get_initial_guess(
                muscle_tendon_parameters_num, data_train, param_config,
                'measured', verbose=True)
    """
    raise NotImplementedError(
        "Adapt load_setup() to return data, bounds, CasADi functions, "
        "param_index and param_config."
    )


def run_all(
    setup,
    analyses=("conditioning", "profile", "sobol"),
    profile_kwargs=None,
    sobol_kwargs=None,
):
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    profile_kwargs = profile_kwargs or {}
    sobol_kwargs = sobol_kwargs or {}
    results = {}

    # Vérification de cohérence du setup une seule fois, en amont
    layout = Layout.from_param_index(setup["param_index"])
    print(f"\n[setup] {layout.summary()}")
    for key in ("initial_guess", "lower_band", "upper_band", "true_params"):
        v = setup[key]
        assert len(v) == layout.n_opt, (
            f"setup['{key}'] de taille {len(v)} != n_opt={layout.n_opt}"
        )

    if "conditioning" in analyses:
        print("\n" + "#" * 72)
        print("# 1/3 Conditioning analysis")
        print("#" * 72)
        out = OUTPUT_ROOT / "conditioning"
        results["conditioning"] = hessian_analysis(
            data=setup["data"],
            true_params=setup["true_params"],
            skeleton_num=setup["skeleton_num"],
            unknown_parameters=setup["unknown_parameters"],
            casadi_function=setup["casadi_function"],
            param_index=setup["param_index"],
            out_dir=out,
        )

    if "profile" in analyses:
        print("\n" + "#" * 72)
        print("# 2/3 Profile likelihood analysis")
        print("#" * 72)
        out = OUTPUT_ROOT / "profile_likelihood"
        results["profile"] = run_profile_likelihood(
            data=setup["data"],
            true_params=setup["true_params"],
            lower_band=setup["lower_band"],
            upper_band=setup["upper_band"],
            skeleton_num=setup["skeleton_num"],
            unknown_parameters=setup["unknown_parameters"],
            casadi_function=setup["casadi_function"],
            optimization_nlp=setup["optimization_nlp"],
            param_index=setup["param_index"],
            out_dir=out,
            **profile_kwargs,
        )

    if "sobol" in analyses:
        print("\n" + "#" * 72)
        print("# 3/3 Sobol global sensitivity analysis")
        print("#" * 72)
        out = OUTPUT_ROOT / "sobol"
        results["sobol"] = run_sobol_analysis(
            data=setup["data"],
            lower_band=setup["lower_band"],
            upper_band=setup["upper_band"],
            skeleton_num=setup["skeleton_num"],
            casadi_function=setup["casadi_function"],
            param_index=setup["param_index"],
            out_dir=out,
            **sobol_kwargs,
        )

    print("\n" + "=" * 72)
    print(f"  All analyses complete. Results in {OUTPUT_ROOT.resolve()}")
    print("=" * 72)
    return results


def main(setup=None):
    parser = argparse.ArgumentParser(
        description="Run identifiability analyses on the MTU NLP."
    )
    parser.add_argument("--analyses", nargs="+",
                        default=["conditioning", "profile", "sobol"],
                        choices=["conditioning", "profile", "sobol"])
    parser.add_argument("--n-grid", type=int, default=11)
    parser.add_argument("--rel-range", type=float, default=0.30)
    parser.add_argument("--sobol-n-base", type=int, default=256)
    args = parser.parse_args()

    if setup is None:
        setup = load_setup()

    run_all(
        setup,
        analyses=tuple(args.analyses),
        profile_kwargs={"n_grid": args.n_grid, "rel_range": args.rel_range},
        sobol_kwargs={"n_base": args.sobol_n_base},
    )


if __name__ == "__main__":
    main()

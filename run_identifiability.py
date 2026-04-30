"""
Orchestrator: run all identifiability analyses in sequence.

Analyses performed
------------------
1. Conditioning analysis (Hessian eigendecomposition)
   -> local identifiability at the true parameters.

2. Profile likelihood
   -> practical identifiability: for each parameter, how much the cost
      increases when that parameter is forced away from its true value,
      with the others re-optimised.

3. Sobol global sensitivity analysis
   -> which parameters drive the variance of the model prediction (torque)
      over the whole admissible space.

Usage
-----
Adapt the `load_setup()` function to your own setup (data loading, CasADi
function creation, bounds, etc.). Then run:

    python run_identifiability.py

Or import selectively:

    from run_identifiability import run_all
    run_all(..., analyses=["conditioning", "profile"])

Outputs are saved under `identifiability_results/{analysis_name}/`.
"""

from __future__ import annotations

from pathlib import Path
import argparse

# Local modules (must be in the same directory)
from conditioning_analysis import hessian_analysis
from profile_likelihood import run_profile_likelihood
from sobol_analysis import run_sobol_analysis


OUTPUT_ROOT = Path("identifiability_results")


def load_setup():
    """Build the arguments needed by each analysis.

    Replace this with your own loading logic: it must produce a dict with
    exactly these keys:
      - data             : np.ndarray (15, n_trials), NOISE-FREE
      - initial_guess    : np.ndarray (12,)
      - lower_band       : np.ndarray (12,)
      - upper_band       : np.ndarray (12,)
      - skeleton_num     : iterable with the skeleton geometry
      - true_params      : np.ndarray (12,) ground-truth parameters
      - unknown_parameters : CasADi SX (12,) symbolic
      - casadi_function  : dict of CasADi Functions
      - optimization_nlp : callable, the NLP wrapper
    """

    raise NotImplementedError(
        "Adapt load_setup() to return your data, bounds, and CasADi functions."
    )


def run_all(
    setup,
    analyses=("conditioning", "profile", "sobol"),
    profile_kwargs=None,
    sobol_kwargs=None,
):
    """Run the requested analyses.

    Parameters
    ----------
    setup : dict
        The dict returned by load_setup().
    analyses : tuple of str
        Subset of {"conditioning", "profile", "sobol"}.
    profile_kwargs : dict
        Extra keyword arguments for run_profile_likelihood
        (e.g. n_grid=15, rel_range=0.5, param_indices=[3,4,5]).
    sobol_kwargs : dict
        Extra keyword arguments for run_sobol_analysis
        (e.g. n_base=128 for a quicker run).
    """
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    profile_kwargs = profile_kwargs or {}
    sobol_kwargs = sobol_kwargs or {}
    results = {}

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
            out_dir=out,
            **sobol_kwargs,
        )

    print("\n" + "=" * 72)
    print(f"  All analyses complete. Results in {OUTPUT_ROOT.resolve()}")
    print("=" * 72)
    return results


def main(setup):
    parser = argparse.ArgumentParser(
        description="Run identifiability analyses on the MTU NLP."
    )
    parser.add_argument(
        "--analyses",
        nargs="+",
        default=["conditioning", "profile", "sobol"],
        choices=["conditioning", "profile", "sobol"],
        help="Which analyses to run.",
    )
    parser.add_argument(
        "--n-grid", type=int, default=11,
        help="Profile likelihood: grid points per parameter.",
    )
    parser.add_argument(
        "--rel-range", type=float, default=0.30,
        help="Profile likelihood: relative range around the true value.",
    )
    parser.add_argument(
        "--sobol-n-base", type=int, default=256,
        help="Sobol: base sample size (total ~ n_base * 26 for 12 params).",
    )
    args = parser.parse_args()

    # setup = load_setup()

    run_all(
        setup,
        analyses=tuple(args.analyses),
        profile_kwargs={
            "n_grid": args.n_grid,
            "rel_range": args.rel_range,
        },
        sobol_kwargs={
            "n_base": args.sobol_n_base,
        },
    )


if __name__ == "__main__":
    main()

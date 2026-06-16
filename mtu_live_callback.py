"""
mtu_live_callback.py
====================

Live plot pour l'identification des parametres muscle-tendon (NLP), compatible
Windows / PyCharm. Architecture THREAD (pas de multiprocessing) :
  - matplotlib tourne dans le thread PRINCIPAL (obligatoire sous Windows),
  - le solveur IPOPT tourne dans un thread de travail,
  - le callback pousse les donnees via une queue thread-safe.

Aucune dependance a `ocp`, `g_names`, ni a la machinerie states/controls.

Affiche deux fenetres mises a jour a chaque iteration IPOPT :
  1. "IPOPT output"   : f, max|g|, inf_pr, inf_du (echelle log)
  2. "MTU Parameters" : un sous-plot par parametre, une courbe par muscle,
                        en VALEURS PHYSIQUES (de-scalees par `scale`).

Scaling
-------
Le NLP optimise des variables sans dimension `up ~ O(1)` avec
`up_phys = up * scale` et `scale = initial_guess`. Le vecteur `x` recu par le
callback contient donc, sur ses `n_opt` premiers elements, les valeurs SCALEES.
Pour l'affichage physique on multiplie par `scale`.

Utilisation (recommandee) — un seul appel qui lance l'optim + le plot
--------------------------------------------------------------------
    from mtu_live_callback import run_with_live_plot

    sol = run_with_live_plot(
        nlp={'x': w, 'f': j, 'g': g},
        x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg,
        opts_ipopt=opts_ipopt,              # SANS iteration_callback (ajoute auto)
        grad_f_func=grad_f_func,
        grad_g_func=grad_g_func,
        param_index=param_index,
        scale=scale,                        # = initial_guess
        lower_band=lower_band,
        upper_band=upper_band,
        initial_guess=initial_guess,
        muscles=('ta', 'sol', 'gast'),
        n_trials=n_trials,                  # active la decomposition de inf_du
    )
    w_opt = sol['x'].full().flatten()

IMPORTANT (Windows) : appelle ceci depuis ton bloc principal protege :
    if __name__ == '__main__':
        ...
        sol = run_with_live_plot(...)
"""

import threading
import queue as _queue

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import casadi as cas


PARAM_LABELS = {
    "l0m": r"$\ell_0^m$ (optimal fiber length) [m]",
    "phi0": r"$\phi_0$ (pennation angle) [rad]",
    "f0m": r"$F_0^m$ (max isometric force) [N]",
    "km": r"$k_m$ (passive muscle stiffness)",
    "lst": r"$\ell_s^t$ (tendon slack length) [m]",
    "kt": r"$k_t$ (tendon stiffness)",
}


def _to_1d(v):
    return np.asarray(v).flatten()


# ===========================================================================
#  Figures
# ===========================================================================

DUAL_BLOCK_COLORS = {
    "params": "tab:red",
    "fiber": "tab:blue",
    "pennation": "tab:green",
    "tendon": "tab:orange",
}
DUAL_BLOCK_LABELS = {
    "params": "params MTU",
    "fiber": "fiber length",
    "pennation": "pennation",
    "tendon": "tendon length",
}


def _create_ipopt_plot():
    fig, axs = plt.subplots(4, 1, num="IPOPT output", figsize=(6, 8))
    for ax, lab in zip(axs, ["f", "max|g|", "inf_pr", "inf_du"]):
        ax.set_ylabel(lab, fontweight="bold")
        ax.grid(True, alpha=0.3)
        ax.set_yscale("log")
    axs[3].set_xlabel("iteration")

    # f, max|g|, inf_pr : une courbe noire chacun
    plots = [axs[i].plot([0], [1], "-", marker=".", color="k")[0] for i in range(3)]

    # inf_du : total (noir pointille) + une courbe par bloc de variables
    du_lines = {}
    du_lines["total"] = axs[3].plot([], [], "--", color="k", lw=1.2,
                                    label="total")[0]
    for name in ("params", "fiber", "pennation", "tendon"):
        du_lines[name] = axs[3].plot(
            [], [], "-", marker=".", ms=3, color=DUAL_BLOCK_COLORS[name],
            label=DUAL_BLOCK_LABELS[name])[0]
    axs[3].legend(fontsize=7, ncol=2, loc="upper right")

    fig.tight_layout()
    return fig, plots, du_lines, axs


def _update_ipopt_plot(plots, du_lines, axes, hist):
    # f, max|g|, inf_pr
    for i, k in enumerate(["f", "max_g", "inf_pr"]):
        y = hist[k]
        plots[i].set_data(range(len(y)), y)
        axes[i].set_xlim(0, max(len(y), 1))
        ypos = [v for v in y if v is not None and v > 0]
        if ypos:
            axes[i].set_ylim(min(ypos) * 0.5, max(ypos) * 2)

    # inf_du : total + blocs
    n = len(hist["inf_du"])
    xs = range(n)
    du_lines["total"].set_data(xs, hist["inf_du"])
    all_pos = [v for v in hist["inf_du"] if v and v > 0]
    has_blocks = len(hist["du_params"]) == n and n > 0
    if has_blocks:
        for name in ("params", "fiber", "pennation", "tendon"):
            y = hist["du_" + name]
            du_lines[name].set_data(xs, y)
            all_pos += [v for v in y if v and v > 0]
    else:
        for name in ("params", "fiber", "pennation", "tendon"):
            du_lines[name].set_data([], [])
    axes[3].set_xlim(0, max(n, 1))
    if all_pos:
        axes[3].set_ylim(min(all_pos) * 0.5, max(all_pos) * 2)


def _create_param_plot(param_index, lower_band, upper_band, initial_guess, muscles):
    param_names = list(param_index.keys())
    n = len(param_names)
    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))
    fig, axs = plt.subplots(nrows, ncols, num="MTU Parameters",
                            figsize=(4.8 * ncols, 3.2 * nrows))
    axs = np.atleast_1d(axs).reshape(-1)

    lb, ub, w0 = _to_1d(lower_band), _to_1d(upper_band), _to_1d(initial_guess)
    cmap = get_cmap("viridis")
    n_m = len(muscles)
    colors = {m: cmap(k / max(n_m - 1, 1)) for k, m in enumerate(muscles)}

    plots, axes = {}, {}
    for i, pname in enumerate(param_names):
        ax = axs[i]
        axes[pname] = ax
        lbs, ubs = [], []
        for m in muscles:
            idx = param_index[pname][m]
            lbs.append(float(lb[idx])); ubs.append(float(ub[idx]))
            plots[(pname, m)] = ax.plot([], [], "-", marker=".",
                                        color=colors[m], label=m)[0]
            ax.axhline(float(w0[idx]), color=colors[m], ls=":", lw=0.8, alpha=0.6)
        lo, hi = min(lbs), max(ubs)
        ax.axhline(lo, color="grey", ls="--", lw=0.7)
        ax.axhline(hi, color="grey", ls="--", lw=0.7)
        span = hi - lo
        pad = 0.1 * span if span > 0 else max(abs(hi), 1.0) * 0.1
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel("iteration")
        ax.set_title(PARAM_LABELS.get(pname, pname), fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=n_m)
    for j in range(n, len(axs)):
        axs[j].axis("off")
    fig.tight_layout()
    return fig, plots, axes, param_names


def _update_param_plot(plots, axes, param_names, muscles, hist,
                       param_index, x, scale, n_opt):
    x = _to_1d(x); scale = _to_1d(scale)
    phys = x[:n_opt] * scale[:n_opt]
    for pname in param_names:
        ax = axes[pname]
        n_iter, bits = 0, []
        for m in muscles:
            idx = param_index[pname][m]
            val = float(phys[idx])
            hist[(pname, m)].append(val)
            y = hist[(pname, m)]; n_iter = len(y)
            plots[(pname, m)].set_data(range(n_iter), y)
            bits.append(f"{m}={val:.4g}")
        ax.set_xlim(0, max(n_iter, 1))
        ax.set_title(f"{PARAM_LABELS.get(pname, pname)}\n" + "  ".join(bits),
                     fontsize=8)


# ===========================================================================
#  Callback CasADi (thread-safe : pousse dans une queue)
# ===========================================================================

class MtuOnlineCallback(cas.Callback):
    """
    Callback IPOPT dedie au NLP MTU. Pousse les donnees dans `out_queue`
    (queue.Queue thread-safe) consommee par le thread principal (matplotlib).
    """

    def __init__(self, nx, ng, grad_f_func, grad_g_func, n_opt, out_queue,
                 n_trials=None, name="mtu_cb"):
        cas.Callback.__init__(self)
        self.nx = int(nx)
        self.ng = int(ng)
        self.grad_f_func = grad_f_func
        self.grad_g_func = grad_g_func
        self.n_opt = int(n_opt)
        self.out_queue = out_queue

        # ---- Masques d'indices des blocs de w, pour decomposer inf_du ----
        # Structure de w :
        #   [0 : n_opt]                      -> parametres MTU (up)
        #   puis par essai t, un bloc de 9 :
        #     [n_opt + 9t      : n_opt+9t+3] -> fiber_length
        #     [n_opt + 9t + 3  : n_opt+9t+6] -> pennation_angle
        #     [n_opt + 9t + 6  : n_opt+9t+9] -> tendon_length
        self.dual_blocks = None
        if n_trials is not None:
            n_trials = int(n_trials)
            base = self.n_opt
            fib = np.zeros(self.nx, dtype=bool)
            pen = np.zeros(self.nx, dtype=bool)
            ten = np.zeros(self.nx, dtype=bool)
            for t in range(n_trials):
                s = base + 9 * t
                fib[s + 0:s + 3] = True
                pen[s + 3:s + 6] = True
                ten[s + 6:s + 9] = True
            par = np.zeros(self.nx, dtype=bool)
            par[0:self.n_opt] = True
            self.dual_blocks = {
                "params": par,
                "fiber": fib,
                "pennation": pen,
                "tendon": ten,
            }

        self.out_queue = out_queue
        self.construct(name, {})

    @staticmethod
    def get_n_in():
        return cas.nlpsol_n_out()

    @staticmethod
    def get_n_out():
        return 1

    @staticmethod
    def get_name_in(i):
        return cas.nlpsol_out(i)

    @staticmethod
    def get_name_out(_):
        return "ret"

    def get_sparsity_in(self, i):
        n = cas.nlpsol_out(i)
        if n == "f":
            return cas.Sparsity.scalar()
        elif n in ("x", "lam_x"):
            return cas.Sparsity.dense(self.nx)
        elif n in ("g", "lam_g"):
            return cas.Sparsity.dense(self.ng)
        return cas.Sparsity(0, 0)

    def eval(self, arg):
        a = {s: arg[i] for i, s in enumerate(cas.nlpsol_out())}
        x = np.array(a["x"]).flatten()
        f = float(np.array(a["f"]).flatten()[0])
        g = np.array(a["g"]).flatten()
        lam_x = np.array(a["lam_x"]).flatten()
        lam_g = np.array(a["lam_g"]).flatten()

        inf_pr = float(np.max(np.abs(g))) if g.size else 0.0
        grad_f = np.array(self.grad_f_func(x)).flatten()
        grad_g_lam = np.array(self.grad_g_func(x) @ lam_g).flatten() if g.size else 0.0 * grad_f
        r = grad_f + grad_g_lam - lam_x          # residu dual complet (par variable)
        inf_du = float(np.max(np.abs(r)))

        # Decomposition de inf_du par bloc de variables (max|r| sur chaque bloc)
        du_blocks = {}
        if self.dual_blocks is not None:
            absr = np.abs(r)
            for name, mask in self.dual_blocks.items():
                du_blocks[name] = float(np.max(absr[mask])) if mask.any() else 1e-12
                du_blocks[name] = max(du_blocks[name], 1e-12)

        self.out_queue.put({
            "x": x, "f": f,
            "max_g": max(inf_pr, 1e-12),
            "inf_pr": max(inf_pr, 1e-12),
            "inf_du": max(inf_du, 1e-12),
            "du_blocks": du_blocks,
        })
        return [0]


# ===========================================================================
#  Runner : lance l'optim dans un thread, le plot dans le thread principal
# ===========================================================================

def run_with_live_plot(nlp, x0, lbx, ubx, lbg, ubg, opts_ipopt,
                       grad_f_func, grad_g_func,
                       param_index, scale, lower_band, upper_band, initial_guess,
                       muscles=("ta", "sol", "gast"),
                       n_trials=None,
                       save_on_finish=True):
    """
    Lance le solveur IPOPT (thread de travail) et affiche le live plot
    (thread principal). Retourne le dict `sol` de CasADi a la fin.

    opts_ipopt : NE PAS y mettre 'iteration_callback' (ajoute automatiquement).
    n_trials   : nombre d'essais. Si fourni, le graphe inf_du est decompose
                 en 4 courbes colorees (params MTU / fiber / pennation / tendon).
                 Si None, inf_du reste une seule courbe (total).
    """
    n_opt = int(_to_1d(initial_guess).shape[0])
    w = nlp["x"]
    g = nlp.get("g", cas.SX.zeros(0))
    nx = w.shape[0]
    ng = g.shape[0]

    data_q = _queue.Queue()
    callback = MtuOnlineCallback(nx, ng, grad_f_func, grad_g_func, n_opt, data_q,
                                 n_trials=n_trials)

    opts = dict(opts_ipopt)
    opts["iteration_callback"] = callback
    solver = cas.nlpsol("solver", "ipopt", nlp, opts)

    result = {}
    done = threading.Event()

    def _solve():
        try:
            result["sol"] = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
        except Exception as e:  # noqa
            result["error"] = e
        finally:
            done.set()

    worker = threading.Thread(target=_solve, daemon=True)

    # --- Figures dans le thread principal ---
    plt.ion()
    ip_fig, ip_plots, ip_du_lines, ip_axes = _create_ipopt_plot()
    p_fig, p_plots, p_axes, p_names = _create_param_plot(
        param_index, lower_band, upper_band, initial_guess, muscles)

    ip_hist = {k: [] for k in ("f", "max_g", "inf_pr", "inf_du")}
    for name in ("params", "fiber", "pennation", "tendon"):
        ip_hist["du_" + name] = []
    p_hist = {(p, m): [] for p in p_names for m in muscles}

    worker.start()

    nb = 0
    while not done.is_set() or not data_q.empty():
        drained = False
        while not data_q.empty():
            a = data_q.get()
            for k in ("f", "max_g", "inf_pr", "inf_du"):
                ip_hist[k].append(a[k])
            du = a.get("du_blocks") or {}
            for name in ("params", "fiber", "pennation", "tendon"):
                if name in du:
                    ip_hist["du_" + name].append(du[name])
            _update_ipopt_plot(ip_plots, ip_du_lines, ip_axes, ip_hist)
            _update_param_plot(p_plots, p_axes, p_names, muscles, p_hist,
                               param_index, a["x"], scale, n_opt)
            nb += 1
            drained = True
        if drained:
            ip_fig.canvas.draw_idle()
            p_fig.canvas.draw_idle()
        plt.pause(0.1)

    # Dernier rafraichissement
    ip_fig.canvas.draw_idle()
    p_fig.canvas.draw_idle()
    if save_on_finish:
        ip_fig.savefig("ipopt_output_final.png", dpi=100)
        p_fig.savefig("mtu_params_final.png", dpi=100)

    if "error" in result:
        plt.ioff()
        raise result["error"]

    print(f"[live plot] {nb} iterations tracees. Fenetres laissees ouvertes.")
    plt.ioff()
    plt.show()  # bloque a la fin pour laisser les fenetres ouvertes
    return result["sol"]
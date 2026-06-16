"""
param_layout.py
===============

Couche unique de "layout" partagée par toutes les analyses d'identifiabilité
(conditioning, profile likelihood, Sobol, Monte-Carlo).

Avant : chaque script supposait 12 paramètres dans un ordre fixe
        (l0m_1..3, phi0_1..3, f0m_1..3, lst_1..3).

Maintenant : le vecteur optimisé a une taille variable n_opt = 3 * n_sym et
        un ordre dicté par `param_config` / `param_index` (sortie de
        useful.get_initial_guess). Seuls les paramètres déclarés 'sym' sont
        dans le vecteur ; les 'fixed' sont codés en dur dans les
        casadi_function et ne sont JAMAIS analysés.

Ce module transforme `param_index` (ou directement `param_config`) en :
  - names    : noms ordonnés du vecteur optimisé (ex: ['l0m_1', 'l0m_2', ...])
  - units    : unité par composante
  - groups   : groupe (propriété MTU) par composante
  - n_opt    : taille du vecteur
  - sym_props: liste ordonnée des propriétés 'sym'

Aucune analyse ne doit plus écrire "12" ou indexer par position absolue :
elle passe par `Layout`.
"""

from __future__ import annotations

import numpy as np


# ------------------------------------------------------------------
# Métadonnées statiques des propriétés MTU (indépendantes de sym/fixed)
# ------------------------------------------------------------------


PROPERTY_ORDER = ["l0m", "phi0", "f0m", "lst","kt"]


PROPERTY_UNIT = {
    "l0m": "m",
    "phi0": "rad",
    "f0m": "N",
    "lst": "m",
    "kt": "[-]",
}

# Étiquette d'affichage (pour titres de figures / manuscrit)
PROPERTY_LABEL = {
    "l0m": "ℓom",
    "phi0": "φo",
    "f0m": "f0m",
    "lst": "ℓst",
    "kt": "kt",
}

N_MUSCLE = 3


class Layout:
    """Décrit le vecteur de paramètres optimisés à partir de param_index.

    Le format attendu de param_index est :
        { 'l0m':  {1: 0, 2: 1, 3: 2},
          'phi0': {1: 3, 2: 4, 3: 5},
          ... }
    c.-à-d. {param_name: {muscle_id: position_dans_le_vecteur_optimisé}}.

    Si ton get_initial_guess renvoie un format légèrement différent, seul
    `from_param_index` est à ajuster ; tout le reste du module est générique.
    """

    def __init__(self, names, units, groups, sym_props, index_map):
        self.names = list(names)
        self.units = list(units)
        self.groups = list(groups)
        self.sym_props = list(sym_props)
        self.index_map = dict(index_map)   # {name: position}
        self.n_opt = len(self.names)

    # -- Constructeurs -------------------------------------------------

    @classmethod
    def from_param_index(cls, param_index):
        """Construit le layout depuis le dict param_index de get_initial_guess.

        On reconstruit l'ordre du vecteur en triant par la position (idx)
        renvoyée par get_initial_guess, ce qui rend le module robuste à
        l'ordre d'itération du dict.
        """
        flat = []  # (position, name, prop)
        for prop, muscle_map in param_index.items():
            if prop not in PROPERTY_UNIT:
                raise KeyError(
                    f"Propriété inconnue dans param_index: '{prop}'. "
                    f"Attendu parmi {list(PROPERTY_UNIT)}."
                )
            for muscle_id, pos in muscle_map.items():
                name = f"{prop}_{muscle_id}"
                flat.append((int(pos), name, prop))

        flat.sort(key=lambda t: t[0])  # ordre du vecteur optimisé

        positions = [p for p, _, _ in flat]
        if positions != list(range(len(positions))):
            raise ValueError(
                f"Les positions de param_index ne sont pas 0..n-1 contiguës: "
                f"{positions}. Vérifie get_initial_guess."
            )

        names = [n for _, n, _ in flat]
        props = [pr for _, _, pr in flat]
        units = [PROPERTY_UNIT[pr] for pr in props]
        groups = [PROPERTY_LABEL[pr] for pr in props]
        sym_props = [p for p in PROPERTY_ORDER if p in param_index]
        index_map = {n: i for i, n in enumerate(names)}
        return cls(names, units, groups, sym_props, index_map)

    @classmethod
    def from_param_config(cls, param_config):
        """Construit le layout directement depuis param_config.

        Pratique si tu veux dériver le layout sans rappeler get_initial_guess
        (ex: pré-visualiser ce qui sera analysé). Reproduit la convention
        "3 muscles par propriété 'sym'", dans l'ordre PROPERTY_ORDER.
        """
        sym_props = [p for p in PROPERTY_ORDER if param_config.get(p) == "sym"]
        names, units, groups = [], [], []
        for prop in sym_props:
            for m in range(1, N_MUSCLE + 1):
                names.append(f"{prop}_{m}")
                units.append(PROPERTY_UNIT[prop])
                groups.append(PROPERTY_LABEL[prop])
        index_map = {n: i for i, n in enumerate(names)}
        return cls(names, units, groups, sym_props, index_map)

    # -- Accès pratiques ----------------------------------------------

    def __len__(self):
        return self.n_opt

    def position(self, name):
        """Position d'un paramètre nommé dans le vecteur optimisé."""
        return self.index_map[name]

    def group_slices(self):
        """Renvoie [(label, i0, i1, unit), ...] pour les plots groupés.

        Itère sur les propriétés 'sym' présentes, dans l'ordre du vecteur.
        Gère le cas où une propriété aurait moins de 3 muscles (sécurité).
        """
        slices = []
        i = 0
        for prop in self.sym_props:
            count = sum(1 for g in self.groups[i:] if g == PROPERTY_LABEL[prop])
            # nombre de composantes consécutives de cette propriété
            j = i
            while j < self.n_opt and self.groups[j] == PROPERTY_LABEL[prop]:
                j += 1
            slices.append((PROPERTY_LABEL[prop], i, j, PROPERTY_UNIT[prop]))
            i = j
        return slices

    def summary(self):
        return (
            f"Layout: n_opt={self.n_opt}, "
            f"sym={self.sym_props}, "
            f"names={self.names}"
        )

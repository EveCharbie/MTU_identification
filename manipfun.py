# function to work on the data from neuro-musculo-skeletal (emg and echo driven)
#################################################################################
# lib
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path
import matplotlib.pyplot as plt
import os
from datetime import datetime

# function
def select_folder_and_get_files():
    root = tk.Tk()
    root.withdraw()

    folder = filedialog.askdirectory(
        title="select the interest folder"
    )

    if not folder:
        return None, None, None

    folder_path = Path(folder)
    results = {
        'osim': None,
        'train': None,
        'test': None
    }

    # find .osim
    osim_files = list(folder_path.glob('*.osim'))
    if osim_files:
        results['osim'] = osim_files[0]

    # Chercher Excel (train et test)
    excel_files = list(folder_path.glob('*.xlsx')) + \
                  list(folder_path.glob('*.xls')) + \
                  list(folder_path.glob('*.csv'))

    for excel_file in excel_files:
        filename_lower = excel_file.name.lower()
        if 'train' in filename_lower and not results['train']:
            results['train'] = excel_file
        elif 'test' in filename_lower and not results['test']:
            results['test'] = excel_file

    # return paht as string
    return (str(results['osim']) if results['osim'] else None,
            str(results['train']) if results['train'] else None,
            str(results['test']) if results['test'] else None)

def split_path_name(full_path):
    """
    Avoir le nom et le path
        Args:
        full_path (str):
            [full_path]
    """
    full_path = Path(full_path)
    folder_name = str(full_path.parent)
    file_name = full_path.name
    return folder_name, file_name

def import_data_from_excel(full_path):
    """
    Importe les données depuis un fichier Excel

    Args:
        full_path (str): Chemin complet du fichier Excel

    Returns:
        data (np.ndarray): Shape (15, ntrials) contenant :
                          [ankle_torque, q_knee, q_ankle,
                           a_tibialis, a_soleus, a_gastrocnemius,
                           fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
                           pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius,
                           tendon_length_tibialis, tendon_length_soleus, tendon_length_gastrocnemius]
    """

    # Vérifier que le fichier existe
    if not Path(full_path).exists():
        raise FileNotFoundError(f"Le fichier {full_path} n'existe pas")

    # Charger le fichier
    print(f"Chargement : {full_path}")
    try:
        if full_path.endswith('.csv'):
            df = pd.read_csv(full_path)
        else:
            df = pd.read_excel(full_path)
    except Exception as e:
        raise ValueError(f"Erreur lecture fichier : {e}")

    # Mapping des colonnes Excel vers les noms attendus
    column_mapping = {
        'τ_ankle': 'ankle_torque',
        'q_1': 'q_knee',
        'q_2': 'q_ankle',
        'a_Ta': 'a_tibialis',
        'lf_Ta^': 'fiber_length_tibialis',
        'φ_Ta': 'pennation_angle_tibialis',
        'a_sol': 'a_soleus',
        'lf_sol^': 'fiber_length_soleus',
        'φ_sol': 'pennation_angle_soleus',
        'a_Gast': 'a_gastrocnemius',
        'lf_Gast^': 'fiber_length_gastrocnemius',
        'φ_Gast': 'pennation_angle_gastrocnemius'
    }

    # Vérifier colonnes manquantes
    missing_columns = [col for col in column_mapping.keys() if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Colonnes manquantes : {missing_columns}")

    # Renommer les colonnes
    df = df.rename(columns=column_mapping)

    # Ajouter colonnes tendon_length (si pas dans Excel, initialiser à NaN)
    if 'tendon_length_tibialis' not in df.columns:
        df['tendon_length_tibialis'] = np.nan
    if 'tendon_length_soleus' not in df.columns:
        df['tendon_length_soleus'] = np.nan
    if 'tendon_length_gastrocnemius' not in df.columns:
        df['tendon_length_gastrocnemius'] = np.nan

    # Ordre des colonnes pour la sortie
    output_columns = [
        'ankle_torque',
        'q_knee', 'q_ankle',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius',
        'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
        'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
        'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
    ]

    # Extraire et transposer
    data = df[output_columns].values.T  # Shape (15, ntrials)

    # Inverser le signe du torque (première ligne) - changement de sens du repère
    data[0, :] = data[0, :]
    data[1, :] = data[1, :]
    data[2, :] = data[2, :]

    # ANGLES : deg to rad
    # q_knee (ligne 1) et q_ankle (ligne 2)
    data[1, :] = np.deg2rad(data[1, :])
    data[2, :] = np.deg2rad(data[2, :])

    # pennation angles (lignes 9, 10, 11)
    data[9, :] = np.deg2rad(data[9, :])  # pennation_angle_tibialis
    data[10, :] = np.deg2rad(data[10, :])  # pennation_angle_soleus
    data[11, :] = np.deg2rad(data[11, :])  # pennation_angle_gastrocnemius

    # LONGUEURS : cm to m (diviser par 100)
    # fiber_length (lignes 6, 7, 8)
    data[6, :] = data[6, :] / 100  # fiber_length_tibialis
    data[7, :] = data[7, :] / 100  # fiber_length_soleus
    data[8, :] = data[8, :] / 100  # fiber_length_gastrocnemius
    # tendon_length (lignes 12, 13, 14)
    data[12, :] = data[12, :] / 100  # tendon_length_tibialis
    data[13, :] = data[13, :] / 100  # tendon_length_soleus
    data[14, :] = data[14, :] / 100  # tendon_length_gastrocnemius

    print(f"✓ Shape : {data.shape}")
    print(f"  Conversions appliquées :")
    print(f"    - Torque : inversé")
    print(f"    - Angles (q, φ) : deg → rad")
    print(f"    - Longueurs (lf, lst) : cm → m")

    return data

def add_tendon_length_to_data(data, skeleton_num, casadi_function):
    """
    Ajoute les longueurs de tendon aux données en utilisant les fonctions CasADi

    Args:
        data (np.ndarray): Shape (15, ntrials) contenant les données
        skeleton_num (np.ndarray): Paramètres du squelette OpenSim
        casadi_function (dict): Dictionnaire contenant les fonctions CasADi

    Returns:
        data (np.ndarray): Shape (15, ntrials) avec tendon_length calculées
    """
    from math import cos

    # Indices dans data
    idx_q_knee = 1
    idx_q_ankle = 2
    idx_fiber_length_ta = 6
    idx_fiber_length_sol = 7
    idx_fiber_length_gast = 8
    idx_pennation_ta = 9
    idx_pennation_sol = 10
    idx_pennation_gast = 11
    idx_tendon_ta = 12
    idx_tendon_sol = 13
    idx_tendon_gast = 14

    # Extraire q_knee et q_ankle
    q_knee = data[idx_q_knee, :]
    q_ankle = data[idx_q_ankle, :]

    # Extraire longueurs de fibres et angles de pennation
    fiber_length_ta = data[idx_fiber_length_ta, :]
    fiber_length_sol = data[idx_fiber_length_sol, :]
    fiber_length_gast = data[idx_fiber_length_gast, :]

    pennation_ta = data[idx_pennation_ta, :]
    pennation_sol = data[idx_pennation_sol, :]
    pennation_gast = data[idx_pennation_gast, :]

    ntrials = data.shape[1]

    # Initialiser les arrays pour les longueurs de tendon
    tendon_length_ta = np.zeros(ntrials)
    tendon_length_sol = np.zeros(ntrials)
    tendon_length_gast = np.zeros(ntrials)

    # Boucle sur tous les essais
    for trial in range(ntrials):
        # Construire le vecteur d'état musculo-squelettique
        q = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        q[4] = np.deg2rad(q_knee[trial])  # q5
        q[5] = np.deg2rad(q_ankle[trial])  # q6

        pennation_ta_rad = np.deg2rad(pennation_ta[trial])
        pennation_sol_rad = np.deg2rad(pennation_sol[trial])
        pennation_gast_rad = np.deg2rad(pennation_gast[trial])

        # Concaténer avec les paramètres du squelette
        musculoskeletal_states_num = np.concatenate((q, skeleton_num))

        # Calculer les longueurs MTU (muscle-tendon unit)
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_num)
        mtu_length_ta = float(mtu_length[0])
        mtu_length_sol = float(mtu_length[1])
        mtu_length_gast = float(mtu_length[2])

        # Calculer les longueurs de tendon
        # tendon_length = mtu_length - cos(pennation_angle) * fiber_length
        tendon_length_ta[trial] = mtu_length_ta - cos(pennation_ta_rad) * fiber_length_ta[trial]
        tendon_length_sol[trial] = mtu_length_sol - cos(pennation_sol_rad) * fiber_length_sol[trial]
        tendon_length_gast[trial] = mtu_length_gast - cos(pennation_gast_rad) * fiber_length_gast[trial]

    # Ajouter les longueurs de tendon calculées aux données
    data[idx_tendon_ta, :] = tendon_length_ta
    data[idx_tendon_sol, :] = tendon_length_sol
    data[idx_tendon_gast, :] = tendon_length_gast

    print(f"✓ Tendon lengths calculées et ajoutées")
    print(f"  Tibialis   : min={np.nanmin(tendon_length_ta):.4f}, max={np.nanmax(tendon_length_ta):.4f}")
    print(f"  Soleus     : min={np.nanmin(tendon_length_sol):.4f}, max={np.nanmax(tendon_length_sol):.4f}")
    print(f"  Gastrocnemius : min={np.nanmin(tendon_length_gast):.4f}, max={np.nanmax(tendon_length_gast):.4f}")

    return data

def get_initial_guess(muscle_tendon_parameters_num, data):
    """
    Génère initial_guess, lower_band et upper_band à partir des données de test
    et des paramètres musculo-tendineux.

    Args:
        muscle_tendon_parameters_num (np.ndarray): Shape (12,) contenant :
            [lom_ta, lom_sol, lom_gast,
             phi0_ta, phi0_sol, phi0_gast,
             Fom_ta, Fom_sol, Fom_gast,
             lst_ta, lst_sol, lst_gast]

        data (np.ndarray): Shape (15, ntrials) contenant les données de test

    Returns:
        tuple: (initial_guess, upper_band, lower_band)
            - initial_guess (np.ndarray): Shape (12,) - valeurs initiales
            - upper_band (np.ndarray): Shape (12,) - bornes supérieures
            - lower_band (np.ndarray): Shape (12,) - bornes inférieures
    """

    # Indices des variables dans data
    # [ankle_torque, q_knee, q_ankle,
    #  a_tibialis, a_soleus, a_gastrocnemius,
    #  fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
    #  pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius,
    #  tendon_length_tibialis, tendon_length_soleus, tendon_length_gastrocnemius]
    range_band = .8 # % de variation

    idx_fiber_length_ta = 6
    idx_fiber_length_sol = 7
    idx_fiber_length_gast = 8
    idx_pennation_ta = 9
    idx_pennation_sol = 10
    idx_pennation_gast = 11
    idx_tendon_ta = 12
    idx_tendon_sol = 13
    idx_tendon_gast = 14

    # Extraire les données pertinentes
    fiber_length_ta_test = data[idx_fiber_length_ta, :]
    fiber_length_sol_test = data[idx_fiber_length_sol, :]
    fiber_length_gast_test = data[idx_fiber_length_gast, :]
    pennation_ta_test = data[idx_pennation_ta, :]
    pennation_sol_test = data[idx_pennation_sol, :]
    pennation_gast_test = data[idx_pennation_gast, :]
    tendon_ta_test = data[idx_tendon_ta, :]
    tendon_sol_test = data[idx_tendon_sol, :]
    tendon_gast_test = data[idx_tendon_gast, :]

    # Paramètres musculo-tendineux (indices)
    idx_lom_ta = 0
    idx_lom_sol = 1
    idx_lom_gast = 2
    idx_phi0_ta = 3
    idx_phi0_sol = 4
    idx_phi0_gast = 5
    idx_Fom_ta = 6
    idx_Fom_sol = 7
    idx_Fom_gast = 8
    idx_lst_ta = 9
    idx_lst_sol = 10
    idx_lst_gast = 11

    # Initialiser les arrays
    initial_guess = np.zeros(12)
    lower_band = np.zeros(12)
    upper_band = np.zeros(12)

    # Pour chaque paramètre
    # lom_ta
    min_val = np.nanmin(fiber_length_ta_test)
    max_val = np.nanmax(fiber_length_ta_test)
    lower_band[idx_lom_ta] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_lom_ta] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_lom_ta] = (lower_band[idx_lom_ta] + upper_band[idx_lom_ta]) / 2

    # lom_sol
    min_val = np.nanmin(fiber_length_sol_test)
    max_val = np.nanmax(fiber_length_sol_test)
    lower_band[idx_lom_sol] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_lom_sol] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_lom_sol] = (lower_band[idx_lom_sol] + upper_band[idx_lom_sol]) / 2

    # lom_gast
    min_val = np.nanmin(fiber_length_gast_test)
    max_val = np.nanmax(fiber_length_gast_test)
    lower_band[idx_lom_gast] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_lom_gast] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_lom_gast] = (lower_band[idx_lom_gast] + upper_band[idx_lom_gast]) / 2

    # phi0_ta
    min_val = np.nanmin(pennation_ta_test)
    max_val = np.nanmax(pennation_ta_test)
    lower_band[idx_phi0_ta] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_phi0_ta] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_phi0_ta] = (lower_band[idx_phi0_ta] + upper_band[idx_phi0_ta]) / 2

    # phi0_sol
    min_val = np.nanmin(pennation_sol_test)
    max_val = np.nanmax(pennation_sol_test)
    lower_band[idx_phi0_sol] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_phi0_sol] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_phi0_sol] = (lower_band[idx_phi0_sol] + upper_band[idx_phi0_sol]) / 2

    # phi0_gast
    min_val = np.nanmin(pennation_gast_test)
    max_val = np.nanmax(pennation_gast_test)
    lower_band[idx_phi0_gast] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_phi0_gast] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_phi0_gast] = (lower_band[idx_phi0_gast] + upper_band[idx_phi0_gast]) / 2

    # Fom : conserver directement les valeurs de muscle_tendon_parameters_num
    initial_guess[idx_Fom_ta] = muscle_tendon_parameters_num[idx_Fom_ta]
    lower_band[idx_Fom_ta] = muscle_tendon_parameters_num[idx_Fom_ta] * (1-range_band)
    upper_band[idx_Fom_ta] = muscle_tendon_parameters_num[idx_Fom_ta] * (1+range_band)

    initial_guess[idx_Fom_sol] = muscle_tendon_parameters_num[idx_Fom_sol]
    lower_band[idx_Fom_sol] = muscle_tendon_parameters_num[idx_Fom_sol] * (1-range_band)
    upper_band[idx_Fom_sol] = muscle_tendon_parameters_num[idx_Fom_sol] * (1+range_band)

    initial_guess[idx_Fom_gast] = muscle_tendon_parameters_num[idx_Fom_gast]
    lower_band[idx_Fom_gast] = muscle_tendon_parameters_num[idx_Fom_gast] * (1-range_band)
    upper_band[idx_Fom_gast] = muscle_tendon_parameters_num[idx_Fom_gast] * (1+range_band)

    # lst : à partir des données de test
    min_val = np.nanmin(tendon_ta_test)
    max_val = np.nanmax(tendon_ta_test)
    lower_band[idx_lst_ta] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_lst_ta] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_lst_ta] = (lower_band[idx_lst_ta] + upper_band[idx_lst_ta]) / 2

    min_val = np.nanmin(tendon_sol_test)
    max_val = np.nanmax(tendon_sol_test)
    lower_band[idx_lst_sol] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_lst_sol] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_lst_sol] = (lower_band[idx_lst_sol] + upper_band[idx_lst_sol]) / 2

    min_val = np.nanmin(tendon_gast_test)
    max_val = np.nanmax(tendon_gast_test)
    lower_band[idx_lst_gast] = min_val - 0.1 * np.abs(min_val)
    upper_band[idx_lst_gast] = max_val + 0.1 * np.abs(max_val)
    initial_guess[idx_lst_gast] = (lower_band[idx_lst_gast] + upper_band[idx_lst_gast]) / 2

    print(f"✓ Initial guess generé")
    print(f"  Shape : {initial_guess.shape}")

    return initial_guess, upper_band, lower_band

def plot_data(data, muscle_names=None, trial_indices=None, save_path=None):
    """
    Visualise les données expérimentales utilisées pour l'identification MTU.

    Args:
        data (np.ndarray): Shape (15, n_trials).
            Lignes :
              [0]    : torque mesuré (N.m)
              [1:3]  : angles articulaires q (rad)
              [3:6]  : activations musculaires normalisées [0, 1]
              [6:9]  : longueurs de fibres mesurées (m)
              [9:12] : angles de pennation mesurés (rad)
              [12:15]: longueurs de tendon mesurées (m)
        muscle_names (list[str], optional): noms des 3 muscles. Défaut : ['M1','M2','M3'].
        trial_indices (array-like, optional): indices x à afficher. Défaut : range(n_trials).
        save_path (str, optional): si fourni, sauvegarde la figure.

    Returns:
        fig, axes : objets matplotlib.
    """
    n_trials = data.shape[1]
    if muscle_names is None:
        muscle_names = ['M1', 'M2', 'M3']
    if trial_indices is None:
        trial_indices = np.arange(n_trials)

    # Extraction (data en shape (15, n_trials))
    torque   = data[0, :]
    q        = data[1:3, :]    # (2, n_trials)
    a        = data[3:6, :]    # (3, n_trials)
    fiber_l  = data[6:9, :]    # (3, n_trials)
    penn     = data[9:12, :]   # (3, n_trials)
    tendon_l = data[12:15, :]  # (3, n_trials)

    # Couleurs cohérentes par muscle
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, 3))

    fig, axes = plt.subplots(3, 2, figsize=(13, 10), constrained_layout=True)
    fig.suptitle(f'Données expérimentales — {n_trials} trials', fontsize=13, fontweight='bold')

    # --- Torque ---
    ax = axes[0, 0]
    ax.plot(trial_indices, torque, 'o-', color='black', markersize=4, linewidth=1)
    ax.set_ylabel('Torque (N.m)')
    ax.set_title('Couple articulaire mesuré')
    ax.grid(alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5)

    # --- Angles articulaires q ---
    ax = axes[0, 1]
    ax.plot(trial_indices, np.rad2deg(q[0, :]), 'o-', label='q1', markersize=4, linewidth=1)
    ax.plot(trial_indices, np.rad2deg(q[1, :]), 's-', label='q2', markersize=4, linewidth=1)
    ax.set_ylabel('Angle (°)')
    ax.set_title('Angles articulaires')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)

    # --- Activations ---
    ax = axes[1, 0]
    for m in range(3):
        ax.plot(trial_indices, a[m, :], 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Activation [0–1]')
    ax.set_title('Activations musculaires')
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)

    # --- Longueurs de fibres ---
    ax = axes[1, 1]
    for m in range(3):
        ax.plot(trial_indices, fiber_l[m, :] * 1000, 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Longueur fibre (mm)')
    ax.set_title('Longueurs de fibres mesurées')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)
    if (fiber_l < 0).any():
        ax.axhline(0, color='red', linewidth=0.8, linestyle='--', alpha=0.6)

    # --- Angles de pennation ---
    ax = axes[2, 0]
    for m in range(3):
        ax.plot(trial_indices, np.rad2deg(penn[m, :]), 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Pennation (°)')
    ax.set_xlabel('Trial')
    ax.set_title('Angles de pennation mesurés')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5)

    # --- Longueurs de tendon ---
    ax = axes[2, 1]
    for m in range(3):
        ax.plot(trial_indices, tendon_l[m, :] * 1000, 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Longueur tendon (mm)')
    ax.set_xlabel('Trial')
    ax.set_title('Longueurs de tendon mesurées')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)
    if (tendon_l < 0).any():
        ax.axhline(0, color='red', linewidth=0.8, linestyle='--', alpha=0.6)

    # --- Diagnostic console : valeurs anormales ---
    print("\n=== Diagnostic data ===")
    print(f"n_trials : {n_trials}")
    print(f"Torque         : [{torque.min():.2f}, {torque.max():.2f}] N.m")
    print(f"Fiber length   : [{fiber_l.min()*1000:.2f}, {fiber_l.max()*1000:.2f}] mm")
    print(f"Pennation      : [{np.rad2deg(penn.min()):.2f}, {np.rad2deg(penn.max()):.2f}] °")
    print(f"Tendon length  : [{tendon_l.min()*1000:.2f}, {tendon_l.max()*1000:.2f}] mm")

    neg_fl = np.argwhere(fiber_l < 0)   # (muscle, trial)
    neg_tl = np.argwhere(tendon_l < 0)
    if len(neg_fl):
        print(f"⚠ Fiber length négatif en (muscle, trial) : {neg_fl.tolist()}")
    if len(neg_tl):
        print(f"⚠ Tendon length négatif en (muscle, trial) : {neg_tl.tolist()}")
    if not len(neg_fl) and not len(neg_tl):
        print("✓ Toutes les longueurs sont positives")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure sauvegardée : {save_path}")

    plt.show()
    return fig, axes

def save_muscle_tendon_parameters(muscle_tendon_parameters,
                                  osim_folder,
                                  filename='muscle_tendon_parameters_saved',
                                  muscle_names=('tibialis', 'soleus', 'gastrocnemius'),
                                  save_csv=True,
                                  save_npy=True,
                                  verbose=True):
    """
    Sauvegarde les paramètres muscle-tendon optimisés dans `osim_folder`.

    L'ordre attendu pour `muscle_tendon_parameters` (taille 12) est :
        [ℓom_TA, ℓom_SOL, ℓom_GAST,
         φo_TA,  φo_SOL,  φo_GAST,
         Fom_TA, Fom_SOL, Fom_GAST,
         ℓst_TA, ℓst_SOL, ℓst_GAST]

    Args:
        muscle_tendon_parameters (array-like, 12): paramètres optimisés.
        osim_folder (str): dossier de destination. Créé si inexistant.
        filename (str): nom de base (sans extension).
        muscle_names (tuple[str]): noms des 3 muscles dans l'ordre [TA, SOL, GAST].
        save_csv, save_npy (bool): formats à écrire.
        verbose (bool): logs.

    Returns:
        dict: chemins des fichiers écrits {'csv': ..., 'npy': ...}.
    """

    # --- Validation ---
    params = np.asarray(muscle_tendon_parameters, dtype=float).flatten()
    assert params.size == 12, f"Attendu 12 paramètres, reçu {params.size}"
    assert len(muscle_names) == 3, "3 noms de muscles attendus"

    # --- Structuration par muscle ---
    param_labels = ['optimal_fiber_length_m',
                    'pennation_angle_at_optimal_rad',
                    'max_isometric_force_N',
                    'tendon_slack_length_m']

    # Reshape : (4 paramètres, 3 muscles)
    params_matrix = params.reshape(4, 3)

    df = pd.DataFrame(params_matrix,
                      index=param_labels,
                      columns=list(muscle_names))

    # --- Préparation sortie ---
    os.makedirs(osim_folder, exist_ok=True)
    base = os.path.splitext(filename)[0]
    paths = {}

    # --- CSV : lisible pour inspection / annexes ---
    if save_csv:
        path_csv = os.path.join(osim_folder, base + '.csv')
        df_export = df.copy()
        # Ligne bonus en degrés pour lecture humaine
        df_export.loc['pennation_angle_at_optimal_deg'] = np.rad2deg(df.loc['pennation_angle_at_optimal_rad'])
        df_export.to_csv(path_csv)
        paths['csv'] = path_csv
        if verbose:
            print(f"  → CSV : {path_csv}")

    # --- NPY : recharge rapide dans le pipeline ---
    if save_npy:
        path_npy = os.path.join(osim_folder, base + '.npy')
        np.save(path_npy, params)
        paths['npy'] = path_npy
        if verbose:
            print(f"  → NPY : {path_npy}  (vecteur plat taille 12, unités SI)")

    # --- Console summary ---
    if verbose:
        print(f"\n  Timestamp : {datetime.now().isoformat(timespec='seconds')}")
        print('  Muscle-tendon parameters saved:')
        print(df.round(5).to_string())
        print()

    return paths
"""
Variante : paramètres MTU SCALÉS (ratios O(1)) + coût PONDÉRÉ par 1/sigma**2.
Fonction autonome — comportement identique à l'original.
"""

import numpy as np
from casadi import SX, vertcat, sum1, nlpsol
epsilon = 10-16


def optimization_nlp_scaled_weighted(data, initial_guess, lower_band, upper_band,
                                     skeleton_num, muscle_tendon_parameters_num,
                                     unknown_parameters, casadi_function, param_index):
    """
    Identify the muscle-tendon parameters declared 'sym' in param_config, via an NLP.

    Paramètres optimisés : SCALÉS (on optimise up = valeur_physique / initial_guess).
    Fonction de coût      : PONDÉRÉE par 1/sigma**2.

    Args:
        data (np.ndarray): Shape (15, n_trials).
            Rows:
              [0]     : measured joint torque (N.m)
              [1:3]   : joint angles q (rad)
              [3:6]   : muscle activations [0, 1]
              [6:9]   : measured fiber lengths (m)
              [9:12]  : measured pennation angles (rad)
              [12:15] : measured tendon lengths (m)
        initial_guess, lower_band, upper_band : optimized MTU parameters (size n_opt).
        skeleton_num : skeleton geometry (.osim).
        muscle_tendon_parameters_num : reference values (size n_opt).
        unknown_parameters : SX vector of the n_opt unknowns.
        casadi_function : dict of CasADi functions.
        param_index (dict): {param_name: {muscle: idx}} pour le reporting.

    Returns:
        param_opt (np.ndarray): Shape (n_opt,) — estimated parameters.
    """

    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data must be (15, n_trials), got {data.shape}"

    initial_guess = np.asarray(initial_guess, dtype=float)
    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)

    n_opt = unknown_parameters.shape[0]
    assert initial_guess.shape[0] == n_opt, \
        f"initial_guess ({initial_guess.shape[0]}) != unknown_parameters ({n_opt})"
    assert lower_band.shape[0] == n_opt and upper_band.shape[0] == n_opt, \
        "lower_band/upper_band incohérents avec unknown_parameters"
    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band at: {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ Scaling des paramètres MTU ============ #
    scale = initial_guess.copy()
    assert np.all(scale > 0), \
        "initial_guess doit être strictement positif pour servir d'échelle"

    up = unknown_parameters            # variable IPOPT (sans dimension, ~ O(1))
    up_phys = up * scale               # valeurs physiques

    # ============ NLP set up ============ #
    w, w0, lbw, ubw = [], [], [], []
    g, lbg, ubg = [], [], []
    j = SX(0)

    e_torque, e_fiber, e_pennation = [], [], []

    # Unknown MTU parameters (n_opt values), EN UNITÉS SCALÉES
    w += [up]
    w0 += [1.0] * n_opt                       # initial_guess / scale = 1
    lbw += list(lower_band / scale)
    ubw += list(upper_band / scale)

    # ============ Poids du coût = 1 / sigma**2 ============ #
    sigma_torque = 0.02                                 # N.m
    sigma_fiber = 0.005                                # ~5 mm en m
    sigma_penn = np.deg2rad(5)                         # ~5° en rad

    for trial in range(n_trials):
        data_trial = data[:, trial]

        measured_torque = data_trial[0]                 # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]                        # [0, 1]
        measured_fiber_length = data_trial[6:9]          # m
        measured_pennation_angle = data_trial[9:12]      # rad
        measured_tendon_length = data_trial[12:15]       # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        tendon_length_k = SX.sym(f"Tendon_Length_{trial}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial}", n_muscle)

        fl_meas = np.abs(measured_fiber_length)
        tl_meas = np.abs(measured_tendon_length)
        pa_meas = np.abs(measured_pennation_angle)

        lb_fl = fl_meas - fl_meas * 0.1
        ub_fl = fl_meas + fl_meas * 0.1

        lb_pa = pa_meas - np.deg2rad(2)
        ub_pa = pa_meas + np.deg2rad(2)
        lb_pa = np.clip(lb_pa, 0 + epsilon, np.pi - epsilon)
        ub_pa = np.clip(ub_pa, 0 + epsilon, np.pi - epsilon)

        lb_tl = tl_meas - tl_meas * 0.005
        ub_tl = tl_meas + tl_meas * 0.005

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        w0_k = np.concatenate([fl_meas, pa_meas, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        assert np.all(lb_block <= ub_block), \
            f"Inverted bounds at trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 out of bounds at trial {trial}: " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- Equilibrium constraints (9 per trial) --- #
        k = vertcat(a_trial, mtu_length, up_phys)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Torque simulation --- #
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, up_phys)

        # --- Residuals --- #
        e_torque_trials = measured_torque - torque_simulated
        e_fiber_trials = measured_fiber_length - fiber_length_k
        e_pennation_trials = measured_pennation_angle - pennation_angle_k

        j += (e_torque_trials ** 2) / sigma_torque ** 2
        j += sum1((e_fiber_trials ** 2) / sigma_fiber ** 2)
        j += sum1((e_pennation_trials ** 2) / sigma_penn ** 2)

        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    # ============ Assembly ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"[SCALED | WEIGHTED] n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    if np.any(np.isnan(w0)):
        print('NaNs in w0 at:', np.where(np.isnan(w0)))
        return None
    if np.any(np.isinf(w0)):
        print('Infs in w0 at:', np.where(np.isinf(w0)))
        return None
    bad = np.where(lbw > ubw)[0]
    if len(bad):
        print(f'lbw > ubw at indices: {bad[:20]}...')
        return None

    print('w0 is valid')

    # ============ NLP solver ============ #
    opts_ipopt = {
        "ipopt.max_iter": 2500,
        "ipopt.tol": 1e-4,
        "ipopt.print_info_string": "yes",
        "ipopt.linear_solver": "mumps",
    }
    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp, opts_ipopt)

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    # ============ Extraction ============ #
    w_opt = sol['x'].full().flatten()
    cost = sol['f'].full().item()

    param_opt = w_opt[:n_opt] * scale          # dé-scaling

    err_param = np.abs(muscle_tendon_parameters_num - param_opt)

    # ============ Reporting ============ #
    _report_nlp_results(param_index, muscles=['ta', 'sol', 'gast'],
                        param_opt=param_opt,
                        ref=muscle_tendon_parameters_num,
                        err=err_param, cost=cost, n_trials=n_trials)

    return param_opt



"""
Variante : paramètres MTU SCALÉS (ratios O(1)) + coût NON NORMALISÉ
(résidus bruts au carré, aucune pondération).
Fonction autonome.
"""



def optimization_nlp_scaled_plain(data, initial_guess, lower_band, upper_band,
                                  skeleton_num, muscle_tendon_parameters_num,
                                  unknown_parameters, casadi_function, param_index):
    """
    Identify the muscle-tendon parameters declared 'sym' in param_config, via an NLP.

    Paramètres optimisés : SCALÉS (on optimise up = valeur_physique / initial_guess).
    Fonction de coût      : NON NORMALISÉE (somme des résidus au carré, poids = 1).

    Args / Returns : voir docstring de la version d'origine.
    """

    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data must be (15, n_trials), got {data.shape}"

    initial_guess = np.asarray(initial_guess, dtype=float)
    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)

    n_opt = unknown_parameters.shape[0]
    assert initial_guess.shape[0] == n_opt, \
        f"initial_guess ({initial_guess.shape[0]}) != unknown_parameters ({n_opt})"
    assert lower_band.shape[0] == n_opt and upper_band.shape[0] == n_opt, \
        "lower_band/upper_band incohérents avec unknown_parameters"
    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band at: {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ Scaling des paramètres MTU ============ #
    scale = initial_guess.copy()
    assert np.all(scale > 0), \
        "initial_guess doit être strictement positif pour servir d'échelle"

    up = unknown_parameters            # variable IPOPT (sans dimension, ~ O(1))
    up_phys = up * scale               # valeurs physiques

    # ============ NLP set up ============ #
    w, w0, lbw, ubw = [], [], [], []
    g, lbg, ubg = [], [], []
    j = SX(0)

    e_torque, e_fiber, e_pennation = [], [], []

    # Unknown MTU parameters (n_opt values), EN UNITÉS SCALÉES
    w += [up]
    w0 += [1.0] * n_opt                       # initial_guess / scale = 1
    lbw += list(lower_band / scale)
    ubw += list(upper_band / scale)

    for trial in range(n_trials):
        data_trial = data[:, trial]

        measured_torque = data_trial[0]                 # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]                        # [0, 1]
        measured_fiber_length = data_trial[6:9]          # m
        measured_pennation_angle = data_trial[9:12]      # rad
        measured_tendon_length = data_trial[12:15]       # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        tendon_length_k = SX.sym(f"Tendon_Length_{trial}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial}", n_muscle)

        fl_meas = np.abs(measured_fiber_length)
        tl_meas = np.abs(measured_tendon_length)
        pa_meas = np.abs(measured_pennation_angle)

        lb_fl = fl_meas - fl_meas * 0.1
        ub_fl = fl_meas + fl_meas * 0.1

        lb_pa = pa_meas - np.deg2rad(2)
        ub_pa = pa_meas + np.deg2rad(2)
        lb_pa = np.clip(lb_pa, 0 + epsilon, np.pi - epsilon)
        ub_pa = np.clip(ub_pa, 0 + epsilon, np.pi - epsilon)

        lb_tl = tl_meas - tl_meas * 0.005
        ub_tl = tl_meas + tl_meas * 0.005

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        w0_k = np.concatenate([fl_meas, pa_meas, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        assert np.all(lb_block <= ub_block), \
            f"Inverted bounds at trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 out of bounds at trial {trial}: " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- Equilibrium constraints (9 per trial) --- #
        k = vertcat(a_trial, mtu_length, up_phys)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Torque simulation --- #
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, up_phys)

        # --- Residuals --- #
        e_torque_trials = measured_torque - torque_simulated
        e_fiber_trials = measured_fiber_length - fiber_length_k
        e_pennation_trials = measured_pennation_angle - pennation_angle_k

        # ===== Coût NON NORMALISÉ : résidus bruts au carré ===== #
        j += (e_torque_trials ** 2)
        j += sum1(e_fiber_trials ** 2)
        j += sum1(e_pennation_trials ** 2)

        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    # ============ Assembly ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"[SCALED | PLAIN] n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    if np.any(np.isnan(w0)):
        print('NaNs in w0 at:', np.where(np.isnan(w0)))
        return None
    if np.any(np.isinf(w0)):
        print('Infs in w0 at:', np.where(np.isinf(w0)))
        return None
    bad = np.where(lbw > ubw)[0]
    if len(bad):
        print(f'lbw > ubw at indices: {bad[:20]}...')
        return None

    print('w0 is valid')

    # ============ NLP solver ============ #
    opts_ipopt = {
        "ipopt.max_iter": 2500,
        "ipopt.tol": 1e-4,
        "ipopt.print_info_string": "yes",
        "ipopt.linear_solver": "mumps",
    }
    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp, opts_ipopt)

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    # ============ Extraction ============ #
    w_opt = sol['x'].full().flatten()
    cost = sol['f'].full().item()

    param_opt = w_opt[:n_opt] * scale          # dé-scaling

    err_param = np.abs(muscle_tendon_parameters_num - param_opt)

    # ============ Reporting ============ #
    _report_nlp_results(param_index, muscles=['ta', 'sol', 'gast'],
                        param_opt=param_opt,
                        ref=muscle_tendon_parameters_num,
                        err=err_param, cost=cost, n_trials=n_trials)

    return param_opt



"""
Variante : paramètres MTU BRUTS (valeurs physiques optimisées directement)
+ coût PONDÉRÉ par 1/sigma**2.
Fonction autonome.
"""




def optimization_nlp_raw_weighted(data, initial_guess, lower_band, upper_band,
                                  skeleton_num, muscle_tendon_parameters_num,
                                  unknown_parameters, casadi_function, param_index):
    """
    Identify the muscle-tendon parameters declared 'sym' in param_config, via an NLP.

    Paramètres optimisés : BRUTS (la variable IPOPT EST la valeur physique).
    Fonction de coût      : PONDÉRÉE par 1/sigma**2.

    Args / Returns : voir docstring de la version d'origine.
    """

    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data must be (15, n_trials), got {data.shape}"

    initial_guess = np.asarray(initial_guess, dtype=float)
    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)

    n_opt = unknown_parameters.shape[0]
    assert initial_guess.shape[0] == n_opt, \
        f"initial_guess ({initial_guess.shape[0]}) != unknown_parameters ({n_opt})"
    assert lower_band.shape[0] == n_opt and upper_band.shape[0] == n_opt, \
        "lower_band/upper_band incohérents avec unknown_parameters"
    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band at: {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ Paramètres MTU BRUTS (pas de scaling) ============ #
    up = unknown_parameters            # variable IPOPT = valeur physique
    up_phys = up                       # identité

    # ============ NLP set up ============ #
    w, w0, lbw, ubw = [], [], [], []
    g, lbg, ubg = [], [], []
    j = SX(0)

    e_torque, e_fiber, e_pennation = [], [], []

    # Unknown MTU parameters (n_opt values), EN UNITÉS PHYSIQUES
    w += [up]
    w0 += list(initial_guess)
    lbw += list(lower_band)
    ubw += list(upper_band)

    # ============ Poids du coût = 1 / sigma**2 ============ #
    sigma_torque = 0.02                                 # N.m
    sigma_fiber = 0.005                                # ~5 mm en m
    sigma_penn = np.deg2rad(5)                         # ~5° en rad

    for trial in range(n_trials):
        data_trial = data[:, trial]

        measured_torque = data_trial[0]                 # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]                        # [0, 1]
        measured_fiber_length = data_trial[6:9]          # m
        measured_pennation_angle = data_trial[9:12]      # rad
        measured_tendon_length = data_trial[12:15]       # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        tendon_length_k = SX.sym(f"Tendon_Length_{trial}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial}", n_muscle)

        fl_meas = np.abs(measured_fiber_length)
        tl_meas = np.abs(measured_tendon_length)
        pa_meas = np.abs(measured_pennation_angle)

        lb_fl = fl_meas - fl_meas * 0.1
        ub_fl = fl_meas + fl_meas * 0.1

        lb_pa = pa_meas - np.deg2rad(2)
        ub_pa = pa_meas + np.deg2rad(2)
        lb_pa = np.clip(lb_pa, 0 + epsilon, np.pi - epsilon)
        ub_pa = np.clip(ub_pa, 0 + epsilon, np.pi - epsilon)

        lb_tl = tl_meas - tl_meas * 0.005
        ub_tl = tl_meas + tl_meas * 0.005

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        w0_k = np.concatenate([fl_meas, pa_meas, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        assert np.all(lb_block <= ub_block), \
            f"Inverted bounds at trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 out of bounds at trial {trial}: " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- Equilibrium constraints (9 per trial) --- #
        k = vertcat(a_trial, mtu_length, up_phys)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Torque simulation --- #
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, up_phys)

        # --- Residuals --- #
        e_torque_trials = measured_torque - torque_simulated
        e_fiber_trials = measured_fiber_length - fiber_length_k
        e_pennation_trials = measured_pennation_angle - pennation_angle_k

        j += (e_torque_trials ** 2) / sigma_torque ** 2
        j += sum1((e_fiber_trials ** 2) / sigma_fiber ** 2)
        j += sum1((e_pennation_trials ** 2) / sigma_penn ** 2)

        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    # ============ Assembly ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"[RAW | WEIGHTED] n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    if np.any(np.isnan(w0)):
        print('NaNs in w0 at:', np.where(np.isnan(w0)))
        return None
    if np.any(np.isinf(w0)):
        print('Infs in w0 at:', np.where(np.isinf(w0)))
        return None
    bad = np.where(lbw > ubw)[0]
    if len(bad):
        print(f'lbw > ubw at indices: {bad[:20]}...')
        return None

    print('w0 is valid')

    # ============ NLP solver ============ #
    opts_ipopt = {
        "ipopt.max_iter": 2500,
        "ipopt.tol": 1e-4,
        "ipopt.print_info_string": "yes",
        "ipopt.linear_solver": "mumps",
    }
    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp, opts_ipopt)

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    # ============ Extraction ============ #
    w_opt = sol['x'].full().flatten()
    cost = sol['f'].full().item()

    param_opt = w_opt[:n_opt]                  # pas de dé-scaling (params bruts)

    err_param = np.abs(muscle_tendon_parameters_num - param_opt)

    # ============ Reporting ============ #
    _report_nlp_results(param_index, muscles=['ta', 'sol', 'gast'],
                        param_opt=param_opt,
                        ref=muscle_tendon_parameters_num,
                        err=err_param, cost=cost, n_trials=n_trials)

    return param_opt



"""
Variante : paramètres MTU BRUTS (valeurs physiques optimisées directement)
+ coût NON NORMALISÉ (résidus bruts au carré, aucune pondération).
Fonction autonome — version la plus 'nue' (ni scaling, ni pondération).
"""


def optimization_nlp_raw_plain(data, initial_guess, lower_band, upper_band,
                               skeleton_num, muscle_tendon_parameters_num,
                               unknown_parameters, casadi_function, param_index):
    """
    Identify the muscle-tendon parameters declared 'sym' in param_config, via an NLP.

    Paramètres optimisés : BRUTS (la variable IPOPT EST la valeur physique).
    Fonction de coût      : NON NORMALISÉE (somme des résidus au carré, poids = 1).

    Args / Returns : voir docstring de la version d'origine.
    """

    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data must be (15, n_trials), got {data.shape}"

    initial_guess = np.asarray(initial_guess, dtype=float)
    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)

    n_opt = unknown_parameters.shape[0]
    assert initial_guess.shape[0] == n_opt, \
        f"initial_guess ({initial_guess.shape[0]}) != unknown_parameters ({n_opt})"
    assert lower_band.shape[0] == n_opt and upper_band.shape[0] == n_opt, \
        "lower_band/upper_band incohérents avec unknown_parameters"
    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band at: {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ Paramètres MTU BRUTS (pas de scaling) ============ #
    up = unknown_parameters            # variable IPOPT = valeur physique
    up_phys = up                       # identité

    # ============ NLP set up ============ #
    w, w0, lbw, ubw = [], [], [], []
    g, lbg, ubg = [], [], []
    j = SX(0)

    e_torque, e_fiber, e_pennation = [], [], []

    # Unknown MTU parameters (n_opt values), EN UNITÉS PHYSIQUES
    w += [up]
    w0 += list(initial_guess)
    lbw += list(lower_band)
    ubw += list(upper_band)

    for trial in range(n_trials):
        data_trial = data[:, trial]

        measured_torque = data_trial[0]                 # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]                        # [0, 1]
        measured_fiber_length = data_trial[6:9]          # m
        measured_pennation_angle = data_trial[9:12]      # rad
        measured_tendon_length = data_trial[12:15]       # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        tendon_length_k = SX.sym(f"Tendon_Length_{trial}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial}", n_muscle)

        fl_meas = np.abs(measured_fiber_length)
        tl_meas = np.abs(measured_tendon_length)
        pa_meas = np.abs(measured_pennation_angle)

        lb_fl = fl_meas - fl_meas * 0.1
        ub_fl = fl_meas + fl_meas * 0.1

        lb_pa = pa_meas - np.deg2rad(2)
        ub_pa = pa_meas + np.deg2rad(2)
        lb_pa = np.clip(lb_pa, 0 + epsilon, np.pi - epsilon)
        ub_pa = np.clip(ub_pa, 0 + epsilon, np.pi - epsilon)

        lb_tl = tl_meas - tl_meas * 0.005
        ub_tl = tl_meas + tl_meas * 0.005

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        w0_k = np.concatenate([fl_meas, pa_meas, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        assert np.all(lb_block <= ub_block), \
            f"Inverted bounds at trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 out of bounds at trial {trial}: " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- Equilibrium constraints (9 per trial) --- #
        k = vertcat(a_trial, mtu_length, up_phys)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Torque simulation --- #
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, up_phys)

        # --- Residuals --- #
        e_torque_trials = measured_torque - torque_simulated
        e_fiber_trials = measured_fiber_length - fiber_length_k
        e_pennation_trials = measured_pennation_angle - pennation_angle_k

        # ===== Coût NON NORMALISÉ : résidus bruts au carré ===== #
        j += (e_torque_trials ** 2)
        j += sum1(e_fiber_trials ** 2)
        j += sum1(e_pennation_trials ** 2)

        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    # ============ Assembly ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"[RAW | PLAIN] n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    if np.any(np.isnan(w0)):
        print('NaNs in w0 at:', np.where(np.isnan(w0)))
        return None
    if np.any(np.isinf(w0)):
        print('Infs in w0 at:', np.where(np.isinf(w0)))
        return None
    bad = np.where(lbw > ubw)[0]
    if len(bad):
        print(f'lbw > ubw at indices: {bad[:20]}...')
        return None

    print('w0 is valid')

    # ============ NLP solver ============ #
    opts_ipopt = {
        "ipopt.max_iter": 2500,
        "ipopt.tol": 1e-4,
        "ipopt.print_info_string": "yes",
        "ipopt.linear_solver": "mumps",
    }
    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp, opts_ipopt)

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    # ============ Extraction ============ #
    w_opt = sol['x'].full().flatten()
    cost = sol['f'].full().item()

    param_opt = w_opt[:n_opt]                  # pas de dé-scaling (params bruts)

    err_param = np.abs(muscle_tendon_parameters_num - param_opt)

    # ============ Reporting ============ #
    _report_nlp_results(param_index, muscles=['ta', 'sol', 'gast'],
                        param_opt=param_opt,
                        ref=muscle_tendon_parameters_num,
                        err=err_param, cost=cost, n_trials=n_trials)

    return param_opt


















def _report_nlp_results(param_index, muscles, param_opt, ref, err, cost,
                        n_trials):
    """Affichage des résultats, ordre et contenu dérivés de param_index."""
    # Unités lisibles selon le paramètre (affichage uniquement)
    UNIT = {
        'l0m':  ('m',   1.0),
        'lst':  ('m',   1.0),
        'phi0': ('deg', None),   # converti via rad2deg
        'f0m':  ('N',   1.0),
        'km':   ('-',   1.0),
        'kt':   ('-',   1.0),
    }

    def _fmt(p, idx):
        val_opt = param_opt[idx]
        val_ref = ref[idx]
        val_err = err[idx]
        unit, _ = UNIT.get(p, ('', 1.0))
        if p == 'phi0':
            return (f"{np.rad2deg(val_ref):8.2f} | {np.rad2deg(val_opt):8.2f} | "
                    f"{np.rad2deg(val_err):7.2f}  [deg]")
        return f"{val_ref:8.4f} | {val_opt:8.4f} | {val_err:7.4f}  [{unit}]"

    print("\n" + "=" * 70)
    print(f"Résultats NLP — n_trials = {n_trials} | coût = {cost:.4e}")
    print("=" * 70)
    print(f"{'Param':<10} {'Ref':>8} | {'Estimé':>8} | {'|Diff|':>7}")
    print("-" * 70)
    for p in param_index:                      # ordre = ordre d'optimisation
        for m in muscles:
            idx = param_index[p][m]
            print(f"{p+'_'+m:<10} {_fmt(p, idx)}")
    print("=" * 70)

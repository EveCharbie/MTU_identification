import os
import numpy as np
import pandas as pd
from casadi import SX, DM, vertcat, horzcat, Function, sqrt, exp, if_else, sum1, jacobian, rootfinder, nlpsol, cos, sin
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib

matplotlib.use('TkAgg')  # or 'Qt5Agg', depending on your environment


def get_skeleton():
    """
                              # ========= 1. musculoskeletal Geometry  ========= #
    # 1.1 Topology
    B : Bones
    J : Joins
    df : degree of freedom

    Total df = 6 df
    --> J : Hip (4 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_hip ; X = x ; Y = y; Z = z )
        --> B : Thigh
            --> J : knee (1 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_knee ; X = 0 ; Y = 0 ; Z = 0 )
                --> B : Leg
                    --> J : Ankle (1 dF: Rotx= 0 ; Roty = 0 ; Rotz = theta_ankle ; X = 0 ; Y = 0 ; Z = 0 )
                        --> B : Talus
                            --> J : Subtalar (0 dF: Rotx= 0 ; Roty = 0 ; Rotz = 0 ; X = 0 ; Y = 0 ; Z = 0 )
                                --> B : Calc


    n_muscles = 3 ;
    NameMuscles = {'TibialisAnterior','Soleus','Gastrocnemius'} ;
    nBones = 4 ; NameBones = {'Thigh','Leg','Talus','Calcaneus'} ;
    nJoints = 4 ; NameJoints = {'Hip','Knee','Ankle','Subtalar'} ;
    """

    # ========= useful functions  ========= #
    def rototranslation_rz(translation, theta):
        # Build 4x4 homogeneous transformation matrix for rotation about Z and translation
        from casadi import cos, sin
        tx, ty, tz = translation[0], translation[1], translation[2]
        trans_matrix = SX.zeros(4, 4)
        trans_matrix[0, 0] = cos(theta);
        trans_matrix[0, 1] = -sin(theta)
        trans_matrix[1, 0] = sin(theta);
        trans_matrix[1, 1] = cos(theta)
        trans_matrix[0, 3] = tx
        trans_matrix[1, 3] = ty
        trans_matrix[2, 3] = tz
        trans_matrix[2, 2] = 1
        trans_matrix[3, 3] = 1
        return trans_matrix

    def to_global(trans_matrix, p_local):
        # Apply rotation and translation from 4x4 matrix
        return (trans_matrix @ vertcat(p_local, SX(1)))[:3]

    # Constants
    n_muscles = 3

    # ── generalised coordinates  q ∈ ℝ⁶ ───────────────────────────────────────
    x, y, z = SX.sym('x'), SX.sym('y'), SX.sym('z')
    theta_hip = SX.sym('theta_hip')
    theta_knee = SX.sym('theta_knee')
    theta_ankle = SX.sym('theta_ankle')
    q = vertcat(x, y, z, theta_hip, theta_knee, theta_ankle)

    # ── segment geometry — distal joint position in parent frame ──────────────
    thigh_geom = SX.sym('thigh_geom', 3)  # knee position in thigh frame
    leg_geom = SX.sym('leg_geom', 3)  # ankle position in leg frame
    talus_geom = SX.sym('talus_geom', 3)  # calcaneus position in talus frame
    calc_geom = SX.sym('calc_geom', 3)  # toe position in calcaneus frame
    segment_geometry = vertcat(thigh_geom, leg_geom, talus_geom, calc_geom)

    # ── muscle attachment points (local frames) ───────────────────────────────
    # tibialis anterior : leg → calcaneus, with one via point on the leg
    origin_tibialis = SX.sym('origin_tibialis', 3)
    via_tibialis = SX.sym('via_tibialis', 3)
    insertion_tibialis = SX.sym('insertion_tibialis', 3)

    # soleus : leg → calcaneus
    origin_soleus = SX.sym('origin_soleus', 3)
    insertion_soleus = SX.sym('insertion_soleus', 3)

    # gastrocnemius : thigh → calcaneus (bi-articular)
    origin_gastrocnemius = SX.sym('origin_gastrocnemius', 3)
    insertion_gastrocnemius = SX.sym('insertion_gastrocnemius', 3)

    muscle_geom = vertcat(
        origin_tibialis, origin_soleus, origin_gastrocnemius,
        insertion_tibialis, insertion_soleus, insertion_gastrocnemius,
        via_tibialis
    )

    musculoskeletal_params = vertcat(segment_geometry, muscle_geom)

    # ── forward kinematics : segment frames in global ─────────────────────────
    r_0_thigh = rototranslation_rz(vertcat(x, y, z), theta_hip)
    r_thigh_leg = rototranslation_rz(thigh_geom, -theta_knee)
    r_leg_talus = rototranslation_rz(leg_geom, theta_ankle)
    r_talus_calc = rototranslation_rz(talus_geom, 0)

    # ── express in R0 ─────────────────────────────────────────────────────────
    r_0_leg = r_0_thigh @ r_thigh_leg
    r_0_talus = r_0_leg @ r_leg_talus
    r_0_calc = r_0_talus @ r_talus_calc

    # ── muscle attachments in global frame ────────────────────────────────────
    origin_tibialis_global = to_global(r_0_leg, origin_tibialis)
    via_tibialis_global = to_global(r_0_leg, via_tibialis)
    insertion_tibialis_global = to_global(r_0_calc, insertion_tibialis)

    origin_soleus_global = to_global(r_0_leg, origin_soleus)
    insertion_soleus_global = to_global(r_0_calc, insertion_soleus)

    origin_gastrocnemius_global = to_global(r_0_thigh, origin_gastrocnemius)
    insertion_gastrocnemius_global = to_global(r_0_calc, insertion_gastrocnemius)

    # ── joint and markers ────────────────────────────────────────────────────
    origin_global = horzcat(origin_tibialis_global, origin_soleus_global, origin_gastrocnemius_global)
    insertion_global = horzcat(insertion_tibialis_global, insertion_soleus_global, insertion_gastrocnemius_global)
    via_global = horzcat(via_tibialis_global)
    hjc = r_0_thigh[0:3, 3]
    kjc = r_0_leg[0:3, 3]
    ajc = r_0_talus[0:3, 3]
    tjc = to_global(r_0_calc, calc_geom)
    calc = r_0_calc[0:3, 3]
    markers = horzcat(hjc, kjc, ajc, tjc, calc)

    musculoskeletal_states = vertcat(q, musculoskeletal_params)

    # ── MTU lengths and moment arms ───────────────────────────────────────────
    umt_length = vertcat(
        sqrt(sum1((insertion_tibialis_global - via_tibialis_global) ** 2)) +
        sqrt(sum1((via_tibialis_global - origin_tibialis_global) ** 2)),  # tibialis
        sqrt(sum1((insertion_soleus_global - origin_soleus_global) ** 2)),  # soleus
        sqrt(sum1((insertion_gastrocnemius_global - origin_gastrocnemius_global) ** 2))  # gastroc
    )

    moment_arm = jacobian(umt_length, q)

    # ── CasADi exported functions ─────────────────────────────────────────────
    forward_kinematics = Function("forward_kinematics",
                                  [musculoskeletal_states],
                                  [origin_global, insertion_global, via_global, markers],
                                  ["musculoskeletal_states"],
                                  ["origin_global", "insertion_global", "via_global", "markers"])

    get_mtu_length = Function("get_mtu_length",
                              [musculoskeletal_states],
                              [umt_length],
                              ["musculoskeletal_states"],
                              ["umt_length (tibialis, soleus, gastrocnemius)"])

    get_moment_arm = Function("get_moment_arm",
                              [musculoskeletal_states],
                              [moment_arm],
                              ["musculoskeletal_states"],
                              ["moment Arm (tibialis, soleus, gastrocnemius)"])

    return (q, moment_arm, musculoskeletal_params, musculoskeletal_states,
            forward_kinematics, get_mtu_length, get_moment_arm)

def get_fiber_active_force_length(a, normalized_fiber_length, maximal_isometric_force):
    # === Active Force-Length (S2) ===
    # First Gaussian coefficients
    b11, b21, b31, b41 = 0.814483478343008, 1.055033428970575, 0.162384573599574, 0.063303448465465
    # Second Gaussian coefficients
    b12, b22, b32, b42 = 0.433004984392647, 0.716775413397760, -0.029947116970696, 0.200356847296188
    # Third Gaussian coefficients
    b13, b23, b33, b43 = 0.100, 1.000, 0.5 * np.sqrt(0.5), 0.000

    # Assume these are defined: normalizedFiberLength, a, maximal_isometric_force
    # normalizedFiberLength should be a float or a NumPy array

    # Gaussian 1
    num1 = normalized_fiber_length - b21
    den1 = b31 + b41 * normalized_fiber_length
    fm_tilde1 = b11 * np.exp(-0.5 * (num1 ** 2) / (den1 ** 2))

    # Gaussian 2
    num2 = normalized_fiber_length - b22
    den2 = b32 + b42 * normalized_fiber_length
    fm_tilde2 = b12 * np.exp(-0.5 * (num2 ** 2) / (den2 ** 2))

    # Gaussian 3
    num3 = normalized_fiber_length - b23
    den3 = b33 + b43 * normalized_fiber_length
    fm_tilde3 = b13 * np.exp(-0.5 * (num3 ** 2) / (den3 ** 2))

    # Total normalized active force-length
    normalized_fiber_active_force_length = fm_tilde1 + fm_tilde2 + fm_tilde3

    # Non-normalized active force
    fiber_active_force_length = a * normalized_fiber_active_force_length * maximal_isometric_force
    return fiber_active_force_length

def get_fiber_passive_force_length(normalized_fiber_length, k_fiber, maximal_isometric_force):
    # === Passive Force-Length (S3) ===
    e0 = 0.6

    normalized_fiber_passive_force_part_1 = 0
    normalized_fiber_passive_force_part_2 = (exp(((k_fiber * (
            normalized_fiber_length - 1)) / e0)) - 1) / (exp(k_fiber) - 1)

    normalized_fiber_passive_force = if_else(normalized_fiber_length < 1,
                                             normalized_fiber_passive_force_part_1,
                                             normalized_fiber_passive_force_part_2)  # if normalized length under 0 the force = 0 %Normalized equation

    fiber_passive_force = normalized_fiber_passive_force * maximal_isometric_force  # Non - normalized equation
    return fiber_passive_force

def get_tendon_force_length(normalized_tendon_length, k_tendon, maximal_isometric_force):
    # Tendon force-length (S1)
    c1 = 0.200;
    c2 = 0.995;
    c3 = 0.250  # tendon parameters

    normalized_tendon_force_part_1 = 0
    normalized_tendon_force_part_2 = c1 * exp(
        k_tendon * (normalized_tendon_length - c2)) - c3  # Normalized equation

    normalized_tendon_force = if_else(normalized_tendon_length < 1,
                                      normalized_tendon_force_part_1,
                                      normalized_tendon_force_part_2)  # if normalized length under 0 the force = 0   # Normalized equation

    tendon_force = normalized_tendon_force * maximal_isometric_force  # Non-normalized equation
    return tendon_force

def get_muscle_total_force(fiber_active_force_length, normalized_fiber_force_velocity, fiber_passive_force):
    return fiber_active_force_length * normalized_fiber_force_velocity + fiber_passive_force

def get_muscle_dynamic(q, moment_arm, musculoskeletal, n_muscles,
                       param_config=None, fixed_values=None):
    """
    Build CasADi functions for muscle-tendon dynamics (De Groote 2016 model),
    including equilibrium rootfinders (per-muscle and all-muscles).

    Muscle-tendon parameters (ℓom, φo, Fom, ℓst, kt, km) can be declared either
    as symbolic variables (for optimization / identification) or hardcoded as
    numerical constants, independently for each parameter.

    Parameters
    ----------
    q, musculoskeletal : CasADi SX
    n_muscles : int
    param_config : dict {param_name: 'sym' | 'fixed'}, optional
        Default : all 'sym'.
    fixed_values : dict {param_name: scalar or array-like (n_muscles,)}, optional
        Required for every parameter marked 'fixed'.

    Returns
    -------
    dict with CasADi Functions and rootfinders.

    References
    ----------
    De Groote F, Kinney AL, Rao A V, Fregly BJ. Ann Biomed Eng. 2016;
    44: 2922–2936. doi:10.1007/s10439-016-1591-9
    """

    # ========= 1. Résolution de la configuration des paramètres ========= #
    MTU_PARAM_SPEC = [
        ('l0m', 'ℓom'),
        ('phi0', 'φo'),
        ('f0m', 'Fom'),
        ('km', 'km'),
        ('lst', 'ℓst'),
        ('kt', 'kt'),
    ]
    valid_keys = [k for k, _ in MTU_PARAM_SPEC]

    if param_config is None:
        param_config = {k: 'sym' for k in valid_keys}
    fixed_values = fixed_values or {}

    for k in valid_keys:
        if k not in param_config:
            raise ValueError(
                f"param_config missing entry for '{k}'. Required keys: {valid_keys}"
            )
        if param_config[k] not in ('sym', 'fixed'):
            raise ValueError(
                f"param_config['{k}'] must be 'sym' or 'fixed', got {param_config[k]!r}"
            )

    for k, status in param_config.items():
        if status == 'fixed' and k not in fixed_values:
            raise ValueError(
                f"Parameter '{k}' is marked 'fixed' but missing from fixed_values."
            )

    # ========= 2. Construction des variables MTU ========= #
    def _build_mtu_param(key, display_name):
        if param_config[key] == 'sym':
            return SX.sym(display_name, n_muscles)
        val = np.asarray(fixed_values[key], dtype=float).reshape(-1)
        if val.size == 1:
            val = np.full(n_muscles, val.item())
        if val.size != n_muscles:
            raise ValueError(
                f"fixed_values['{key}'] has size {val.size}, expected 1 or {n_muscles}."
            )
        return SX(val)

    mtu_vars = {key: _build_mtu_param(key, name) for key, name in MTU_PARAM_SPEC}

    optimal_fiber_length = mtu_vars['l0m']
    phi0 = mtu_vars['phi0']
    maximal_isometric_force = mtu_vars['f0m']
    tendon_slack_length = mtu_vars['lst']
    k_tendon = mtu_vars['kt']
    k_fiber = mtu_vars['km']

    # Vecteur des paramètres SYMBOLIQUES (pour optim / rootfinder)
    sym_keys = [k for k, _ in MTU_PARAM_SPEC if param_config[k] == 'sym']
    sym_params = [mtu_vars[k] for k in sym_keys]
    if sym_params:
        muscle_tendon_parameters = vertcat(*sym_params)
    else:
        muscle_tendon_parameters = SX.sym('mtu_params_empty', 0)

    # Version "single muscle" : premier élément de chaque paramètre symbolique
    sym_params_single = [mtu_vars[k][0] for k in sym_keys]
    if sym_params_single:
        muscle_tendon_parameters_single = vertcat(*sym_params_single)
    else:
        muscle_tendon_parameters_single = SX.sym('mtu_params_single_empty', 0)

    # ========= 3. États muscle-tendon ========= #
    a = SX.sym('muscle_activation', n_muscles)
    fiber_length = SX.sym('fiber_length', n_muscles)
    tendon_length = SX.sym('tendon_length', n_muscles)
    pennation_angle = SX.sym('pennation_angle', n_muscles)

    rooted_variables = vertcat(fiber_length, pennation_angle, tendon_length)

    # ========= 4. Équations MTU (De Groote 2016) ========= #
    normalized_tendon_length = tendon_length / tendon_slack_length
    normalized_fiber_length = fiber_length / optimal_fiber_length

    tendon_force = get_tendon_force_length(
        normalized_tendon_length, k_tendon, maximal_isometric_force
    )
    fiber_active_force_length = get_fiber_active_force_length(
        a, normalized_fiber_length, maximal_isometric_force
    )
    fiber_passive_force = get_fiber_passive_force_length(
        normalized_fiber_length, k_fiber, maximal_isometric_force
    )
    normalized_fiber_force_velocity = 1.0  # isométrique
    fiber_force = get_muscle_total_force(
        fiber_active_force_length,
        normalized_fiber_force_velocity,
        fiber_passive_force,
    )

    # ========= 5. Wrap into CasADi functions ========= #
    neuromusculoskeletal_state = vertcat(a, q, musculoskeletal)
    all_states = vertcat(neuromusculoskeletal_state, rooted_variables)

    # ========= 6. Rootfinder : équilibre muscle-tendon ========= #
    # Input géométrique : longueur MTU (ℓmt) — donnée par la cinématique
    l_mtu = SX.sym('umt_length', n_muscles)

    # --- Contraintes d'équilibre (De Groote 2016, eqs S5-S7) --- #
    # g5 : ℓmt = cos(φ) · ℓm + ℓt          (longueur MTU)
    # g6 : ℓom · sin(φo) = ℓm · sin(φ)     (conservation du volume musculaire)
    # g7 : Ff · cos(φ) = Ft                (équilibre des forces axiales)
    g5 = l_mtu - (cos(pennation_angle) * fiber_length + tendon_length)
    g6 = (optimal_fiber_length * sin(phi0)) - (fiber_length * sin(pennation_angle))
    g7 = fiber_force * cos(pennation_angle) - tendon_force

    # --- Problème single-muscle : 3 inconnues, 3 équations --- #
    unknown_single = vertcat(fiber_length[0], pennation_angle[0], tendon_length[0])
    known_single = vertcat(a[0], l_mtu[0], muscle_tendon_parameters_single)

    equilibrium_error_single_muscle = Function(
        'equilibrium_error_single_muscle',
        [unknown_single, known_single],
        [vertcat(g5[0], g6[0], g7[0])],
        ['x', 'p'],
        ['residuals'],
    )

    opts_newton_single = {
        "abstol": 1e-8,
        "max_iter": 1000,
        "error_on_fail": False,
        "print_iteration": False,
    }
    equilibrate_muscle_tendon_single_muscle = rootfinder(
        'equilibrate_muscle_tendon_single_muscle',
        'newton',
        equilibrium_error_single_muscle,
        opts_newton_single,
    )

    # --- Problème all-muscles : 3·n_muscles inconnues, 3·n_muscles équations --- #
    unknown_all = vertcat(fiber_length, pennation_angle, tendon_length)
    known_all = vertcat(a, l_mtu, muscle_tendon_parameters)

    equilibrium_error_all_muscle = Function(
        'equilibrium_error_all_muscle',
        [unknown_all, known_all],
        [vertcat(g5, g6, g7)],
        ['x', 'p'],
        ['residuals'],
    )

    opts_newton_all = {
        "abstol": 1e-8,
        "max_iter": 1000,
        "error_on_fail": False,
        "print_iteration": False,
    }
    equilibrate_muscle_tendon_all = rootfinder(
        'equilibrate_muscle_tendon_all',
        'newton',
        equilibrium_error_all_muscle,
        opts_newton_all,
    )

    # ========= 4. Computing Joint Moments and Angles    ========= #
    joint_torque = moment_arm * tendon_force
    joint_torque = sum1(joint_torque[:, -1:])

    get_joint_moment = Function(
        'get_joint_moment',
        [all_states, muscle_tendon_parameters],
        [joint_torque],
        ['all_states', 'muscle_tendon_parameters'],
        ['joint_torque']
    )

    get_tendon_force_from_tendon_length = Function(
        'get_tendon_force',
        [tendon_length, muscle_tendon_parameters],
        [tendon_force],
        ['tendon_length', 'muscle_tendon_parameters'],
        ['tendon_force']
    )

    get_fiber_force_from_fiber_length = Function(
        'estimateFiberForce',
        [a, fiber_length, muscle_tendon_parameters],
        [fiber_force],
        ['muscle_activation', 'fiber_length', 'muscle_tendon_parameters'],
        ['fiber_force']
    )

    def _wrap(name, output, output_name):
        return Function(
            name,
            [all_states, muscle_tendon_parameters],
            [output],
            ['all_states', 'muscle_tendon_parameters'],
            [output_name],
        )

    get_tendon_force_fn = _wrap('get_tendon_force', tendon_force, 'tendon_force')
    get_fiber_force_fn = _wrap('get_fiber_force', fiber_force, 'fiber_force')
    get_fiber_active_force_fn = _wrap('get_fiber_active_force', fiber_active_force_length, 'fiber_active_force')
    get_fiber_passive_force_fn = _wrap('get_fiber_passive_force', fiber_passive_force, 'fiber_passive_force')
    get_joint_moment_fn = _wrap('get_joint_moment', joint_torque, 'joint_torque')

    # ========= 7. Assemblage du retour ========= #
    return {
        # Fonctions de force
        'get_tendon_force': get_tendon_force_fn,
        'get_fiber_force': get_fiber_force_fn,
        'get_fiber_active_force': get_fiber_active_force_fn,
        'get_fiber_passive_force': get_fiber_passive_force_fn,
        'get_joint_moment': get_joint_moment_fn,
        # Résidus d'équilibre (utiles pour diagnostic)
        'equilibrium_error_single_muscle': equilibrium_error_single_muscle,
        'equilibrium_error_all_muscle': equilibrium_error_all_muscle,
        # Rootfinders
        'equilibrate_muscle_tendon_single_muscle': equilibrate_muscle_tendon_single_muscle,
        'equilibrate_muscle_tendon_all': equilibrate_muscle_tendon_all,
        # Métadonnées
        'muscle_tendon_parameters': muscle_tendon_parameters,
        'muscle_tendon_parameters_single': muscle_tendon_parameters_single,
        'sym_param_order': sym_keys,
        'param_config': dict(param_config),
    }

def get_model_equation(param_config=None, fixed_values=None, n_muscles=3):
    """
    Build the full musculoskeletal model by merging skeleton kinematics
    and muscle-tendon dynamics into a unified set of CasADi functions.

    Parameters
    ----------
    param_config : dict, optional
        Per-parameter status ('sym' or 'fixed') for the 6 MTU parameters
        (l0m, phi0, f0m, lst, kt, km). Default : geometric params symbolic
        (l0m, phi0, f0m, lst), shape params fixed (kt, km).
    fixed_values : dict, optional
        Numerical values for parameters marked 'fixed'. Scalars are broadcast
        to all n_muscles. Default : kt=3.0, km=4.0.
    n_muscles : int, optional
        Number of muscles in the model. Default: 3 (tib_ant, soleus, gast).

    Returns
    -------
    casadi_functions : dict[str, casadi.Function]
        Unified collection of all CasADi functions describing the model :
        kinematics + muscle-tendon dynamics + equilibrium rootfinders.
    sym_param_order : list of str
        Ordered list of MTU parameter names that are symbolic, indicating
        the layout of the 'muscle_tendon_parameters' input vector.
    definition : dict
        Documentation of the symbolic variables used in the model.
    """

    # ========= 1. Default settings ========= #
    if param_config is None:
        param_config = {
            'l0m': 'sym',
            'phi0': 'sym',
            'f0m': 'sym',
            'km': 'fixed',
            'lst': 'sym',
            'kt': 'fixed',
        }
    if fixed_values is None:
        fixed_values = {
            'kt': 35,
            'km': 4.0,
        }

    # ========= 2. Skeleton (kinematics) ========= #
    (q, moment_arm, musculoskeletal_params, musculoskeletal_states,
     forward_kinematics, get_mtu_length, get_moment_arm) = get_skeleton()

    # ========= 3. Muscle-tendon dynamics ========= #
    funcs = get_muscle_dynamic(
        q, moment_arm, musculoskeletal_params,
        n_muscles=n_muscles,
        param_config=param_config,
        fixed_values=fixed_values,
    )

    # ========= 4. Merge all CasADi functions ========= #
    casadi_functions = {
        name: obj for name, obj in funcs.items()
        if isinstance(obj, Function)
    }

    extras = {
        'forward_kinematics': forward_kinematics,
        'get_mtu_length': get_mtu_length,
        'get_moment_arm': get_moment_arm,
        'mtu_parameters_sym': funcs['muscle_tendon_parameters'],
        'mtu_parameters_name': funcs['sym_param_order'],
    }

    name_clash = set(casadi_functions) & set(extras)
    if name_clash:
        raise ValueError(f"Key collision in casadi_functions: {name_clash}")
    casadi_functions.update(extras)

    # ========= 5. Documentation des variables symboliques ========= #
    definition = {
        'a': "Neuromuscular activation (between 0 and 1)",
        'q': "Spatial skeleton configuration "
             "(x, y, z, theta_hip, theta_knee, theta_ankle)",
        'musculoskeletal': "Musculoskeletal parameters "
                           "(segment_geometry, muscle_insertion, "
                           "Local_ViaPoint_tibialis_anterior)",
        'neuromusculoskeletal_state': "Concatenation (a, q, musculoskeletal)",
        'muscle_tendon_parameters': f"Symbolic MTU parameters in order: "
                                    f"{funcs['sym_param_order']}",
        'rooted_variables': "Hidden states resolved by rootfinder "
                            "(fiber_length, pennation_angle, tendon_length)",
        'l_mtu': "Total muscle-tendon unit length (input to equilibrium)",
        'p': "Full input vector for force functions: "
             "(neuromusculoskeletal_state, muscle_tendon_parameters)",
    }

    return casadi_functions, funcs['muscle_tendon_parameters'], definition

def test_model(skeleton_num, muscle_tendon_parameters_num, casadi_function, a_num, q_num):
    # === plot musculosquelton system ===
    origin_global_num, insertion_global_num, via_global, markers_num = casadi_function['forward_kinematics'](
        np.concatenate((q_num, skeleton_num)).reshape(-1, 1))

    # === get muscle length ===
    mtu_length = casadi_function['get_mtu_length'](np.concatenate((q_num, skeleton_num)).reshape(-1, 1))

    tibialis_length = float(mtu_length[0])
    soleus_length = float(mtu_length[1])
    gastrocnemius_length = float(mtu_length[2])

    print('\n', '================== mtu length length ==================')
    print('tibialis: ', tibialis_length, 'm')
    print('soleus: ', soleus_length, 'm')
    print('gastrocnemius: ', gastrocnemius_length, 'm')

    musculoskeletal_states_num = np.concatenate((q_num, skeleton_num))
    neuromusculoskeletal_state_num = np.concatenate((a_num, musculoskeletal_states_num))

    # === get muscle moment arm ===
    mtu_moment_arm = casadi_function['get_moment_arm'](musculoskeletal_states_num)

    tibialis_moment_arm = float(mtu_moment_arm[0, -1])
    soleus_moment_arm = float(mtu_moment_arm[1, -1])
    gastrocnemius_moment_arm = float(mtu_moment_arm[2, -1])

    print('\n', '================== mtu moment arm ==================')
    print('tibialis: ', tibialis_moment_arm, 'm-1')
    print('soleus: ', soleus_moment_arm, 'm-1')
    print('gastrocnemius: ', gastrocnemius_moment_arm, 'm-1')

    # === muscle equilibrium ===

    """
    muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num)

    n_muscle = 3
    n_param = len(muscle_tendon_parameters_num)/n_muscle
    muscle_tendon_parameters_num_tibialis = muscle_tendon_parameters_num[[0, 3, 6, 9]]
    x_opt_tibialis,state_tibialis = root_muscle_dynamics(a_num[0], float(mtu_length[0]),muscle_tendon_parameters_num[[0, 3, 6, 9]], 'tibialis', casadi_function)
    x_opt_soleus,state_soleus = root_muscle_dynamics(a_num[1], float(mtu_length[1]), muscle_tendon_parameters_num[[1, 4, 7, 10]], 'soleus', casadi_function)
    x_opt_gastrocnemius,state_gastrocnemius = root_muscle_dynamics(a_num[2], float(mtu_length[2]),muscle_tendon_parameters_num[[3, 5, 8, 11]], 'gastrocnemius',casadi_function)

    # rooted_variables = [fiber length, pennation angle and tendon length]
    rooted_variables = np.array([x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0],
                                 x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0],
                                 x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()
    """

    results = solve_all_muscles(a_num, mtu_length, muscle_tendon_parameters_num, casadi_function)

    x_opt_tibialis = results['tibialis']['x_opt']
    x_opt_soleus = results['soleus']['x_opt']
    x_opt_gastrocnemius = results['gastrocnemius']['x_opt']
    # rooted_variables = [fiber length, pennation angle and tendon length]

    rooted_variables = np.array([x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0],
                                 x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0],
                                 x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()

    all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

    tendon_force = casadi_function['get_tendon_force'](all_state, muscle_tendon_parameters_num)
    fiber_force = casadi_function['get_fiber_force'](all_state, muscle_tendon_parameters_num)
    fiber_passive_force = casadi_function['get_fiber_passive_force'](all_state, muscle_tendon_parameters_num)
    fiber_active_force = casadi_function['get_fiber_active_force'](all_state, muscle_tendon_parameters_num)

    ankle_torque = casadi_function['get_joint_moment'](all_state, muscle_tendon_parameters_num)

    print('\n================== simulatied ==================')
    print('\n   mtu architecture:')
    print('     - fiber length:')
    print(f'tibialis: {float(x_opt_tibialis[0, 0]):.4f} m')
    print(f'soleus: {float(x_opt_soleus[0, 0]):.4f} m')
    print(f'gastrocnemius: {float(x_opt_gastrocnemius[0, 0]):.4f} m')
    print('     - pennation angle:')
    print(f'tibialis: {float(np.rad2deg(float(np.array(x_opt_tibialis[1, 0]).flatten()[0]))):.4f} deg')
    print(f'soleus: {float(np.rad2deg(float(np.array(x_opt_soleus[1, 0]).flatten()[0]))):.4f} deg')
    print(f'gastrocnemius: {float(np.rad2deg(float(np.array(x_opt_gastrocnemius[1, 0]).flatten()[0]))):.4f} deg')
    print('     - tendon length:')
    print(f'tibialis: {float(x_opt_tibialis[2, 0]):.4f} m')
    print(f'soleus: {float(x_opt_soleus[2, 0]):.4f} m')
    print(f'gastrocnemius: {float(x_opt_gastrocnemius[2, 0]):.4f} m')

    print('\n   tendon force:')
    print(f'tibialis: {float(tendon_force[0]):.4f} N')
    print(f'soleus: {float(tendon_force[1]):.4f} N')
    print(f'gastrocnemius: {float(tendon_force[2]):.4f} N')

    print('\n   fiber force:')
    print('     - fiber force:')
    print(f'tibialis: {float(fiber_force[0]):.4f} N')
    print(f'soleus: {float(fiber_force[1]):.4f} N')
    print(f'gastrocnemius: {float(fiber_force[2]):.4f} N')

    print('     - fiber passive force:')
    print(f'tibialis: {float(fiber_passive_force[0]):.4f} N')
    print(f'soleus: {float(fiber_passive_force[1]):.4f} N')
    print(f'gastrocnemius: {float(fiber_passive_force[2]):.4f} N')

    print('     - fiber active force:')
    print(f'tibialis: {float(fiber_active_force[0]):.4f} N')
    print(f'soleus: {float(fiber_active_force[1]):.4f} N')
    print(f'gastrocnemius: {float(fiber_active_force[2]):.4f} N')

    print('\n   joint torque:')
    print(f'ankle moment: {float(ankle_torque[0]):.4f} N.m')

def solve_all_muscles(a_num, mtu_length, mtu_params, casadi_function):
    """
    Solve the muscle-tendon equilibrium for every muscle.

    mtu_params is assumed to be in PARAM-MAJOR layout :
        [p1_m1, p1_m2, ..., p2_m1, p2_m2, ...]
    For muscle index m, parameters are at indices m, m+n_muscles, m+2*n_muscles...

    Returns
    -------
    dict {muscle_name: {'x_opt': ..., 'state': ..., 'params': ...}}
    """
    muscle_names = 'tibialis', 'soleus', 'gastrocnemius'
    mtu_params = np.asarray(mtu_params)
    n_muscles = len(muscle_names)

    results = {}
    for m_idx, m_name in enumerate(muscle_names):
        params_m = mtu_params[m_idx::n_muscles]
        x_opt, state = root_muscle_dynamics(
            a_num[m_idx],
            float(mtu_length[m_idx]),
            params_m,
            m_name,
            casadi_function,
        )
        results[m_name] = {'x_opt': x_opt, 'state': state, 'params': params_m}
    return results

def root_muscle_dynamics(a, lmtu, parameters, muscle_name, casadi_function):
    """
    Solves for the muscle-tendon equilibrium state using root-finding with CasADi.

    Parameters:
    - a (float): Muscle activation level.
    - lmtu (float): Muscle-tendon unit length.
    - parameters (np.ndarray): Array of muscle-specific parameters required by the model.
    - muscle_name (str): Name of the muscle (used for print/logging purposes).
    - casadi_function (dict): Dictionary containing CasADi functions:
        - 'equilibrate_muscle_tendon_single_muscle': Root-finding function to solve muscle-tendon equilibrium.
        - 'equilibrium_error_single_muscle': Function to evaluate equilibrium error (residuals).

    Returns:
    - x_opt (np.ndarray): Optimized muscle-tendon state vector satisfying equilibrium conditions.
    - equilibrium_status (str): Status indicating whether equilibrium was successfully achieved ('success' or 'fail').

    Description:
    This function attempts to solve the equilibrium condition of a single muscle-tendon unit
    using a root-finding algorithm. It allows up to `max_attempts` (default: 100) to find a
    solution where all residuals fall below a defined tolerance (`lim_residuals`). If successful,
    it returns the solution and the status. The process and results are printed for each attempt.
    """
    lim_residuals = 1e-8
    equilibrium_status = 'fail'
    max_attempts = 100
    attempt = 0

    parameters_name = casadi_function['mtu_parameters_name']
    required = 'l0m', 'phi0', 'f0m', 'lst' #  extract (ℓom, φo, Fom, ℓst)
    params = get_named_params(parameters, parameters_name, required)

    while equilibrium_status == 'fail' and attempt < max_attempts:
        attempt += 1
        x_start = x_start_equi(a, lmtu, params)
        print(x_start)
        try:
            x_opt = casadi_function['equilibrate_muscle_tendon_single_muscle'](
                x_start,
                np.array([a, lmtu] + parameters.tolist())
            )
            if x_opt[1] < 0 or x_opt[1] > math.pi / 2:
                x_opt[1] = DM(float(x_opt[1]) % math.pi / 2)
                x_opt[2] = DM(abs(float(x_opt[2])))

            residuals_print = casadi_function['equilibrium_error_single_muscle'](
                x_opt,
                np.array([a, lmtu] + parameters.tolist())
            )
            residuals = np.array(residuals_print).flatten()
        except Exception as e:
            print("Something went wrong:", e)

        if np.any(residuals > lim_residuals):
            print(f"[Attempt {attempt}] Root-finding failed: residual too high ({residuals})")
        else:
            print(f"[Attempt {attempt}] Success: residual within tolerance")
            equilibrium_status = 'success'

    print('\n', '=============== dynamics: ', muscle_name, ' ===============')
    print('residuals: length equilibrium, architecture equilibrium and force equilibrium')
    print('residuals ', residuals_print)
    print('x_opt: fiber length, pennation angle and tendon lenght')
    print('x_opt ', x_opt)
    print('equilibrium status: ', equilibrium_status)
    print('=============== ========================= ===============', '\n', '\n', '\n', )
    return x_opt, equilibrium_status


def x_start_equi(a, lmtu, parameters, rng=None):
    """ optimal start of the fonction equilibrium: it permit to start in respect for force equilibrium and geometry
   equilibrium.

   equilibrium function:
    g5 = LUMT - (np.cos(pennationAngle) * fiberLength + tendonLength)
    g6 = (optimal_fiber_length * np.sin(phi0)) - (fiberLength * np.sin(pennationAngle))
    g7 = FM * np.cos(pennationAngle) - FT

    Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)

           Returns:
        X_start (np.ndarray): Shape (3,) fiber length, pennation angle and tendon length,

   """
    rng = rng if rng is not None else np.random

    l0m = parameters[0]; phi0 = parameters[1]; lst = parameters[3]

    tendon_length = lst * (1.0 + a * 0.05)
    pennation_angle = phi0 * (1.0 + rng.uniform(-0.2, 0.2))
    fiber_length = (lmtu - tendon_length) / np.cos(pennation_angle)

    x_start = np.array([fiber_length, pennation_angle, tendon_length])

    if np.any(x_start < 0.0) or np.any(np.isinf(x_start)):
        x_start = np.array([
            l0m * rng.uniform(0.8, 1.2),
            phi0,
            lst * rng.uniform(0.8, 1.02),
        ])

    return np.round(x_start, 4)



def get_named_params(parameters, param_labels, required):
    """
    Extract specific parameters by name from a labelled vector.

    Parameters
    ----------
    parameters : array-like
    param_labels : list of str
    required : list of str
        Names of the parameters to extract, in desired order.

    Returns
    -------
    tuple of float
        Values in the same order as `required`.

    Raises
    ------
    KeyError if any required label is missing from param_labels.
    """
    params = dict(zip(param_labels, parameters))
    missing = [k for k in required if k not in params]
    if missing:
        raise KeyError(f"Missing parameters: {missing}. "
                       f"Available: {list(params.keys())}")
    return tuple(params[k] for k in required)

def plot_force_length(a, fiber_length, param, casadi_function):
    # parameter
    l0m = param[0]
    phi0 = param[1]
    f0m = param[2]
    lst = param[3]

    passive_force_ = float(casadi_function['representationMusclePassiveForce'](fiber_length, l0m, f0m))
    active_force_ = float(casadi_function['representationMuscleActiveForceLength'](a, fiber_length, l0m, f0m))
    total_force_ = passive_force_ + active_force_

    fiber_length_range = np.linspace(0.5 * l0m, 1.5 * l0m, 50)
    a_range = np.linspace(0.0, 1.0, 50)

    A, L = np.meshgrid(a_range, fiber_length_range)

    # Evaluate over the grid
    passive_force = np.zeros_like(L)
    active_force = np.zeros_like(L)
    total_force = np.zeros_like(L)

    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            l = L[i, j]
            act = A[i, j]
            try:
                passive_force[i, j] = float(casadi_function['representationMusclePassiveForce'](l, l0m, f0m))
                active_force[i, j] = float(casadi_function['representationMuscleActiveForceLength'](act, l, l0m, f0m))
                total_force[i, j] = passive_force[i, j] + active_force[i, j]
            except Exception as e:
                print(f"Error at i={i}, j={j}, l={l}, a={a}: {e}")
                passive_force[i, j] = np.nan
                active_force[i, j] = np.nan
                total_force[i, j] = np.nan

    # Plotting all three force surfaces
    fig = plt.figure(figsize=(18, 5))

    # Passive Force
    ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    ax1.plot_surface(A, L, passive_force, cmap='plasma')
    ax1.plot(a, fiber_length, passive_force_, color='k', marker='*', markersize=5, markeredgecolor='k',
             markeredgewidth=20)
    ax1.set_title('Passive Muscle Force')
    ax1.set_xlabel('Activation (a)')
    ax1.set_ylabel('Fiber Length (m)')
    ax1.set_zlabel('Passive Force (N)')
    ax1.view_init(elev=10, azim=200)

    # Active Force
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax2.plot_surface(A, L, active_force, cmap='plasma')
    ax2.plot(a, fiber_length, active_force_, color='k', marker='*', markersize=5, markeredgecolor='k',
             markeredgewidth=20)
    ax2.set_title('Active Muscle Force')
    ax2.set_xlabel('Activation (a)')
    ax2.set_ylabel('Fiber Length (m)')
    ax2.set_zlabel('Active Force (N)')
    ax2.view_init(elev=10, azim=200)

    # Total Force
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    ax3.plot_surface(A, L, total_force, cmap='plasma')
    ax3.plot(a, fiber_length, total_force_, color='k', marker='*', markersize=5, markeredgecolor='k',
             markeredgewidth=20)
    ax3.set_title('Total Muscle Force')
    ax3.set_xlabel('Activation (a)')
    ax3.set_ylabel('Fiber Length (m)')
    ax3.set_zlabel('Total Force (N)')
    ax3.view_init(elev=10, azim=200)

    plt.tight_layout()
    plt.show()


def plotmodel(ax, muscle_origins, muscle_insertions, joint_centers, via_point):
    ax.clear()
    ax.set_title("Model", fontsize=14)

    origins = np.array(muscle_origins)
    insertions = np.array(muscle_insertions)
    markers = np.array(joint_centers)
    via_point = np.array(via_point)

    ax.scatter(origins[0], origins[1], origins[2], c='r', label='muscle origin', marker='o')
    ax.scatter(insertions[0], insertions[1], insertions[2], c='b', label='muscle insertion', marker='o')
    ax.scatter(markers[0], markers[1], markers[2], c='k', label='joint center', marker='o')
    ax.scatter(via_point[0], via_point[1], via_point[2], c='y', label='tibialis via point', marker='o')

    for i in range(markers.shape[1] - 1):
        ax.plot([markers[0, i], markers[0, i + 1]],
                [markers[1, i], markers[1, i + 1]],
                [markers[2, i], markers[2, i + 1]], color='black')

    for i in range(1, origins.shape[1]):
        ax.plot([origins[0, i], insertions[0, i]],
                [origins[1, i], insertions[1, i]],
                [origins[2, i], insertions[2, i]], color='red')

    ax.plot([origins[0, 0], via_point[0, 0]],
            [origins[1, 0], via_point[1, 0]],
            [origins[2, 0], via_point[2, 0]], color='red')

    ax.plot([via_point[0, 0], insertions[0, 0]],
            [via_point[1, 0], insertions[1, 0]],
            [via_point[2, 0], insertions[2, 0]], color='red')

    ax.legend()
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')


def interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function):
    q_init = np.zeros(6)
    muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Initial plot
    def update_plot(val=None):
        a = np.array([muscle_sliders[i].val for i in range(3)])
        q = np.array([slider_q[i].val for i in range(6)])
        q[3:] = np.deg2rad(q[3:])  # Convert q4, q5, q6 from degrees to radians

        musculoskeletal_states_num = np.concatenate((q, skeleton_num))
        neuromusculoskeletal_state_num = np.concatenate([a, musculoskeletal_states_num])

        origins, insertions, via_point, markers = casadi_function['forward_kinematics'](musculoskeletal_states_num)
        plotmodel(ax, origins, insertions, markers, via_point)

        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_num)
        mtu_moment_arm = casadi_function['get_moment_arm'](musculoskeletal_states_num)

        print('\n================== MTU Length ==================')
        print(f'tibialis: {float(mtu_length[0]):.4f} m')
        print(f'soleus: {float(mtu_length[1]):.4f} m')
        print(f'gastrocnemius: {float(mtu_length[2]):.4f} m')

        print('\n================== Moment Arms ==================')
        print(f'tibialis: {float(mtu_moment_arm[0, -1]):.4f} m^-1')
        print(f'soleus: {float(mtu_moment_arm[1, -1]):.4f} m^-1')
        print(f'gastrocnemius: {float(mtu_moment_arm[2, -1]):.4f} m^-1')

        # === muscle equilibrium ===
        results = solve_all_muscles(a, mtu_length, muscle_tendon_parameters_num, casadi_function)

        x_opt_tibialis = results['tibialis']['x_opt']
        x_opt_soleus = results['soleus']['x_opt']
        x_opt_gastrocnemius = results['gastrocnemius']['x_opt']
        # rooted_variables = [fiber length, pennation angle and tendon length]

        rooted_variables = np.array([x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0],
                                     x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0],
                                     x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()

        all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

        tendon_force = casadi_function['get_tendon_force'](all_state, muscle_tendon_parameters_num)
        fiber_force = casadi_function['get_fiber_force'](all_state, muscle_tendon_parameters_num)
        fiber_passive_force = casadi_function['get_fiber_passive_force'](all_state, muscle_tendon_parameters_num)
        fiber_active_force = casadi_function['get_fiber_active_force'](all_state, muscle_tendon_parameters_num)

        ankle_torque = casadi_function['get_joint_moment'](all_state, muscle_tendon_parameters_num)

        print('\n================== simulated ==================')
        print('\n   mtu architecture:')
        print('     - fiber length:')
        print(f'tibialis: {float(x_opt_tibialis[0, 0]):.4f} m')
        print(f'soleus: {float(x_opt_soleus[0, 0]):.4f} m')
        print(f'gastrocnemius: {float(x_opt_gastrocnemius[0, 0]):.4f} m')
        print('     - pennation angle:')
        print(f'tibialis: {float(np.rad2deg(float(np.array(x_opt_tibialis[1, 0]).flatten()[0]))):.4f} deg')
        print(f'soleus: {float(np.rad2deg(float(np.array(x_opt_soleus[1, 0]).flatten()[0]))):.4f} deg')
        print(f'gastrocnemius: {float(np.rad2deg(float(np.array(x_opt_gastrocnemius[1, 0]).flatten()[0]))):.4f} deg')
        print('     - tendon length:')
        print(f'tibialis: {float(x_opt_tibialis[2, 0]):.4f} m')
        print(f'soleus: {float(x_opt_soleus[2, 0]):.4f} m')
        print(f'gastrocnemius: {float(x_opt_gastrocnemius[2, 0]):.4f} m')

        print('\n   tendon force:')
        print(f'tibialis: {float(tendon_force[0]):.4f} N')
        print(f'soleus: {float(tendon_force[1]):.4f} N')
        print(f'gastrocnemius: {float(tendon_force[2]):.4f} N')

        print('\n   fiber force:')
        print('     - fiber force:')
        print(f'tibialis: {float(fiber_force[0]):.4f} N')
        print(f'soleus: {float(fiber_force[1]):.4f} N')
        print(f'gastrocnemius: {float(fiber_force[2]):.4f} N')

        print('     - fiber passive force:')
        print(f'tibialis: {float(fiber_passive_force[0]):.4f} N')
        print(f'soleus: {float(fiber_passive_force[1]):.4f} N')
        print(f'gastrocnemius: {float(fiber_passive_force[2]):.4f} N')

        print('     - fiber active force:')
        print(f'tibialis: {float(fiber_active_force[0]):.4f} N')
        print(f'soleus: {float(fiber_active_force[1]):.4f} N')
        print(f'gastrocnemius: {float(fiber_active_force[2]):.4f} N')

        print('\n   joint torque:')
        print(f'ankle moment: {float(ankle_torque[0]):.4f} N.m')

        fig.canvas.draw_idle()

    # Add sliders
    slider_limits = [
        (-1, 1),  # q1
        (-1, 1),  # q2
        (-1, 1),  # q3
        (-90, 90),  # q4
        (0, 90),  # q5
        (-40, 40)  # q6
    ]
    axcolor = 'lightgoldenrodyellow'
    slider_q = []
    ax_slider = []

    # skeleton configuration
    for i, (min_val, max_val) in enumerate(slider_limits):
        ax_slider = plt.axes([0.15, 0.02 + i * 0.035, 0.65, 0.02], facecolor=axcolor)
        slider = Slider(ax_slider, f'q{i + 1}', min_val, max_val, valinit=0.0)
        slider.on_changed(update_plot)
        slider_q.append(slider)

    # Muscle activity sliders: tibialis, soleus, gastrocnemius
    muscle_names = ['tibialis', 'soleus', 'gastroc']
    muscle_sliders = []
    slider_width = 0.12
    slider_height = 0.02
    top = 0.95

    for i in range(3):
        ax_muscle = plt.axes(
            [0.83, top - i * (slider_height + 0.03), slider_width, slider_height],
            facecolor='mistyrose'
        )
        muscle_slider = Slider(ax_muscle, muscle_names[i], 0.0, 1.0, valinit=0.0, color='red')
        muscle_slider.on_changed(update_plot)
        muscle_sliders.append(muscle_slider)

    update_plot()  # initial plot
    ax.view_init(elev=90, azim=-90)
    plt.show()

def hypothetical_data_generator(skeleton_num, muscle_tendon_parameters_num, casadi_function):
    """ generate hypothetical data to make the NLP
   .
    rooted = vertcat(tendon_force,muscle_force,tendon_length, fiber_length, pennation_angle)


           Returns:
        hypothetical_data (np.ndarray): Shape (15,ntrials) including
        [ankle_torque, q_knee, q_ankle
        a_tibialis, a_soleus, a_gastrocnemius,
        fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
        pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius,
        tendon_length_tibialis,tendon_length_soleus,tendon_length_gastrocnemius]
   """

    header = [
        'ankle_torque',
        'q_knee', 'q_ankle',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius',
        'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
        'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
        'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
    ]

    n_trials_fail = 0  # compteur of trials non - optimized(p > 10e-5)
    n_trials_succeeds = 0  # compteur of trials optimized(p > 10e-5)

    # ========= 1 Muscle-Tendon Architecture Equations   ========= #
    # ========= 1.1 Musculo skeletal configuration during trial(input)   ========= #
    qknee = np.linspace(0, 80, 2)  # example knee angles
    qankle = np.linspace(-25, 20, 5)  # example ankle angles

    qknee = np.deg2rad(qknee)
    qankle = np.deg2rad(qankle)

    # ========= 1.2 Neuronal activation(input)   ========= #
    # Muscle activation [Tibialis Anterior, Soleus, Gastrocnemius]
    a_num = np.array([
        [0, 0, 0],
        [0, 0.1, 0.1],
        [0, 0.2, 0.2],
        [0, 0.3, 0.3],
        [0, 0.4, 0.4],
        [0, 0.5, 0.5],
        [0, 0.6, 0.6],
        [0.1, 0, 0],
        [0.2, 0, 0],
        [0.3, 0, 0],
        [0.4, 0, 0],
        [0.5, 0, 0],
        [0.6, 0, 0]
    ])


    # ========= 2.1 data generator   ========= #
    n_trials = len(a_num) * len(qknee) * len(qankle)
    hypothetical_data = np.zeros((n_trials, 15))

    trials = 0

    for i in range(len(a_num)):  # for each muscle activation
        for ii in range(len(qknee)):  # for each knee angle
            for iii in range(len(qankle)):  # for each ankle angle
                # 2.1 Progression
                trials += 1
                percent = round(trials / n_trials, 2)
                print(f"Progressing: {percent * 100:.0f} %")

                # 2.2 Neuromusculoskeletal configuration
                q_num = [0, 0, 0, 0, qknee[ii], qankle[iii]]
                musculoskeletal_states_num = q_num + list(skeleton_num)
                neuromusculoskeletal_state_num = np.concatenate([a_num[i], musculoskeletal_states_num])

                # 2.3 UMT length
                mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_num)


                results = solve_all_muscles(a_num[i], mtu_length, muscle_tendon_parameters_num, casadi_function)

                x_opt_tibialis = results['tibialis']['x_opt']
                x_opt_soleus = results['soleus']['x_opt']
                x_opt_gastrocnemius = results['gastrocnemius']['x_opt']
                # rooted_variables = [fiber length, pennation angle and tendon length]

                state_tibialis = results['tibialis']['state']
                state_soleus = results['soleus']['state']
                state_gastrocnemius = results['gastrocnemius']['state']



                if state_tibialis == 'fail' or state_soleus == 'fail' or state_tibialis == 'fail' or state_gastrocnemius == 'fail':
                    n_trials_fail += 1
                else:
                    n_trials_succeeds += 1

                    rooted_fiber_length = np.array(
                        [x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0]]).flatten()
                    rooted_pennation_angle = np.array(
                        [x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0]]).flatten()
                    rooted_tendon_length = np.array(
                        [x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()


                    rooted_variables = np.array([x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0],
                                                 x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0],
                                                 x_opt_tibialis[2, 0], x_opt_soleus[2, 0],
                                                 x_opt_gastrocnemius[2, 0]]).flatten()

                    all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

                    ankle_torque = casadi_function['get_joint_moment'](all_state, muscle_tendon_parameters_num)

                    # 2.7 extract variable
                    # rooted = vertcat(tendon length,fiber lenght,,pennationAngle)
                    tibialis_fiber_length = float(x_opt_tibialis[0])
                    tibialis_pennation_angle = float(x_opt_tibialis[1])
                    tibialis_tendon_length = float(x_opt_tibialis[2])

                    soleus_fiber_length = float(x_opt_soleus[0])
                    soleus_pennation_angle = float(x_opt_soleus[1])
                    soleus_tendon_length = float(x_opt_soleus[2])

                    gastrocnemius_fiber_length = float(x_opt_gastrocnemius[0])
                    gastrocnemius_pennation_angle = float(x_opt_gastrocnemius[1])
                    gastrocnemius_tendon_length = float(x_opt_gastrocnemius[2])

                    hypothetical_data[n_trials_succeeds - 1] = [
                        float(ankle_torque),
                        qknee[ii], qankle[iii],
                        float(a_num[i, 0]), float(a_num[i, 1]), float(a_num[i, 2]),
                        tibialis_fiber_length, soleus_fiber_length, gastrocnemius_fiber_length,
                        tibialis_pennation_angle, soleus_pennation_angle, gastrocnemius_pennation_angle,
                        tibialis_tendon_length, soleus_tendon_length, gastrocnemius_tendon_length
                    ]
    hypothetical_data = hypothetical_data[0:n_trials_succeeds]

    hypothetical_data = hypothetical_data.T
    print('\n================== hypotetical data generator state ==================')
    print(f"fail trials: {(n_trials_fail / n_trials) * 100:.2f} %")

    return header, hypothetical_data

def nlp_identification(skeleton_num, muscle_tendon_parameters_num, unknown_parameters, casadi_function, data, opts,
                       initial_guess):
    """ root muscle tendon parameter (ℓom, φo, Fom, ℓst)
   .
    skeleton_num: .osim scale squeleton
    muscle_tendon_parameters_num: 
    unknown_parameters
    casadi_function
    data
    opts
    initial_guess

           Returns:
               xopt (np.ndarray): Shape (12,ntrials) (ℓom, φo, Fom, ℓst)

   """

    xopt = []
    n_muscle = 3  # Number of muscle in our model[TibialisAnterior, Soleus, Gastrocnemius]
    n_trials = 20  # Number of trials selected for the estimation of muscle tendon parameters
    muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num)
    initial_guess = np.array(initial_guess)

    # ============ 2. Selection of trials  ============ #
    # ============ 2.1 Interst  ============ #
    if opts == 'random':  # (Random seletion)
        selection = np.random.randint(1, data.shape[0], n_trials)
        selection.sort()
        data = data[selection]

    elif opts == 'chosen':
        n_trials = data.shape[0]

    # ============ 3. NLP configuration  ============ #
    # ============ 3.1 NLP set Up ============ #
    # Initialize variables
    w = []  # decision variables
    w0 = []  # initial guess
    lbw = []  # lower bounds
    ubw = []  # upper bounds
    j = 0  # cost function
    g = []  # constraints
    lbg = []  # lower bounds for constraints
    ubg = []  # upper bounds for constraints

    # Unknown parameter vector concatenation
    up = unknown_parameters

    e_torque = []
    e_fiber = []
    e_pennation = []

    # Create muscle parameters variables
    w += [up]
    w0 += list(initial_guess)  # transpose to match shapes
    lbw += list(initial_guess * 0.05)
    ubw += list(initial_guess * 3)

    # Weightings in cost function
    w_torque = 1  # Nm
    w_length = 0.005  # Nm/mm
    w_angle = (3 / 180) * np.pi  # radians

    # NLP formulation
    for trial in range(n_trials):
        # Extract measured data
        data_trials = data[trial, :]
        a_trial = data_trials[3:6]
        q_trial = [0, 0, 0, 0] + list(np.deg2rad(data_trials[1:3]))

        mesured_torque = data_trials[0]
        mesured_fiber_length = data_trials[6:9]
        mesured_pennation_angle = np.deg2rad(data_trials[9:12])  # mesured in deg but in rad in NLP
        mesured_tendon_length = data_trials[12:15]

        estimated_tendon_force = np.array(
            casadi_function['get_tendon_force_from_tendon_length'](mesured_tendon_length, initial_guess)).flatten()
        estimated_fiber_force = np.array(
            casadi_function['get_fiber_force_from_fiber_length'](a_trial, mesured_fiber_length,
                                                                 initial_guess)).flatten()

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate([a_trial, musculoskeletal_states_trial])

        # Compute UMT length
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        """
        ################################################################################################################
        #                                   equilibrium with 5 constraints function as
        #     g3 = tendon_force_SX - tendon_force
        #     g4 = fiber_force_SX - fiber_force
        #     g5 = l_mtu - (np.cos(pennation_angle) * fiber_length + tendon_length)
        #     g6 = (optimal_fiber_length * np.sin(phi0)) - (fiber_length * np.sin(pennation_angle))
        #     g7 = fiber_force_SX * np.cos(pennation_angle) - tendon_force_SX
        # 15 constraints at total
        ################################################################################################################

        # Define decision variables for the trial: rooted_variables = (tendon force, muscle force, tendon length, fiber length, pennation angle)
        tendon_length_k = SX.sym(f"Tendon_Length_{trial + 1}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial + 1}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial + 1}", n_muscle)
        tendon_force_k = SX.sym(f"Tendon_Force_{trial + 1}", n_muscle)
        fiber_force_k = SX.sym(f"Fiber_Force_{trial + 1}", n_muscle)

        w0_k = np.concatenate([estimated_tendon_force,estimated_fiber_force,mesured_fiber_length,mesured_pennation_angle,mesured_tendon_length]) # zero-based indexing
        w_k = vertcat( tendon_force_k,fiber_force_k,fiber_length_k, pennation_angle_k,tendon_length_k)

        # Append to global variables
        w += [w_k]
        w0 += list(w0_k)
        lbw += list(w0_k * 0.1)
        ubw += list(w0_k * 3)

        # Form the input vector K
        k = vertcat(a_trial, mtu_length, up)

        # Muscle-tendon equilibrium constraints
        constraints = casadi_function['equilibriumError'](w_k, k)
        g += [constraints]
        lbg += [0] * 15
        ubg += [0] * 15

        all_states_nlp = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)

        # Simulate torque and compute errors
        torque_simulated = casadi_function['get_joint_momentNLP'](all_states_nlp, unknown_parameters)
        ################################################################################################################
        """

        ################################################################################################################
        #                                   equilibrium with 3 constraints function as
        #     g5 = l_mtu - (np.cos(pennation_angle) * fiber_length + tendon_length)
        #     g6 = (optimal_fiber_length * np.sin(phi0)) - (fiber_length * np.sin(pennation_angle))
        #     g7 = fiber_force * np.cos(pennation_angle) - tendon_force
        # 9 constraints at total
        ################################################################################################################
        # Define decision variables for the trial: rooted_variables = (tendon force, muscle force, tendon length, fiber length, pennation angle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial + 1}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial + 1}", n_muscle)
        tendon_length_k = SX.sym(f"Tendon_Length_{trial + 1}", n_muscle)

        w0_k = np.concatenate(
            [mesured_fiber_length, mesured_pennation_angle, mesured_tendon_length])  # zero-based indexing
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        # Append to global variables
        w += [w_k]
        w0 += list(w0_k)
        lbw += list(w0_k * 0.8)
        ubw += list(w0_k * 1.2)

        # Form the input vector K
        k = vertcat(a_trial, mtu_length, up)

        # Muscle-tendon equilibrium constraints
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)  # changer le nom de la fonction

        # Simulate torque and compute errors
        torque_simulated = casadi_function['get_joint_moment'](all_states, unknown_parameters)
        ################################################################################################################

        e_torque_trials = mesured_torque - torque_simulated
        e_fiber_trials = mesured_fiber_length - fiber_length_k
        e_pennation_trials = mesured_pennation_angle - pennation_angle_k

        # Update cost
        # j_trials = w_torque * sqrt(e_torque_trials**2) + sum1(w_length * sqrt(e_fiber_trials**2)) + sum1(w_angle * sqrt(e_pennation_trials**2))
        # j += [j_trials]
        j += w_torque * e_torque_trials ** 2
        j += sum1(w_length * e_fiber_trials ** 2)
        j += sum1(w_angle * e_pennation_trials ** 2)

        # Store errors
        e_torque.append(e_torque_trials)  # err between measured and estimated torque (cost function)
        e_fiber.append(e_fiber_trials)  # err between measured and estimated fiber length (cost function)
        e_pennation.append(e_pennation_trials)  # err between measured and estimated pennation angle (cost function)

    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    J = jacobian(g, w)
    sparcity_dense = np.array(J.sparsity())
    zero_index = np.where(sparcity_dense == 0)

    J2 = jacobian(j, w)
    sparcity_dense2 = np.array(J2.sparsity())
    zero_index2 = np.where(sparcity_dense2 == 0)

    J2.is_zero()

    if np.any(np.isnan(w0)):
        print('NaNs found in w0 at indices:', np.where(np.isnan(w0)))
    if np.any(np.isinf(w0)):
        print('Infs found in w0 at indices:', np.where(np.isinf(w0)))
    if not (np.any(np.isnan(w0)) or np.any(np.isinf(w0))):
        print('w0 is valid')
        # ============ NLP solver ============ #
        # "x" opt parameters, 'f' function to minimized, 'g' contraint function
        nlp = {'x': w,
               'f': j,
               'g': g
               }

        solver = nlpsol('solver', 'ipopt', nlp)

        print(solver)

        # Solve the NLP
        sol = solver(
            x0=w0,
            lbx=lbw,
            ubx=ubw,
            lbg=lbg,
            ubg=ubg
        )
        # Extract solution
        w_opt = sol['x'].full().flatten()
        cost = sol['f'].full().item()

        param_opt = w_opt[:12]
        # temp_Function = Function('temp_Function', [w],[g])
        # temp_Function(w_opt)

        err_param = abs(muscle_tendon_parameters_num - param_opt)

        print("\n", "\n", "Error between Reel and Estimated Muscle-Tendon Parameters :")
        print(f"Number of trials : {n_trials}")
        print(f"err ℓom : {err_param[0:3]}")
        print(f"err φo : {err_param[3:6]}")
        print(f"err Fom : {err_param[6:9]}")
        print(f"err ℓst  : {err_param[9:12]}")

        print("Estimated Muscle-Tendon Parameters :")
        print(f"Cost : {cost}")
        print(f"ℓom : {param_opt[0:3]}")
        print(f"φo : {param_opt[3:6]}")
        print(f"Fom : {param_opt[6:9]}")
        print(f"ℓst  : {param_opt[9:12]}")

        xopt = param_opt[0:12]

    return xopt


def simulation_(skeleton_num, muscle_tendon_parameters_num, casadi_function, time_num, q5_num, q6_num, a_tibialis_num,
                a_soleus_num, a_gastrocnemius_num):
    """ simulation_
   .
    rooted = vertcat(tendon_force,muscle_force,tendon_length, fiber_length, pennation_angle)


           Returns:
               simd_data = (np.ndarray): Shape (16,n_trials) including
               [torque,
               fiber_force_tibialis, fiber_force_soleus, fiber_force_gastrocnemius,
               fiber_active_force_tibialis, fiber_active_force_soleus, fiber_force_active_gastrocnemius,
               fiber_passive_force_tibialis, fiber_passive_force_soleus, fiber_passive_force_gastrocnemius,
               fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
               pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius]
   """

    # ========= 0 set up & header & output  ========= #
    n_trials = len(q6_num)

    header_neuro_musculo_state = [
        'frame', 'time',
        'q5', 'q6',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius'
    ]
    simulated_data_neuro_musculo_state = np.zeros((n_trials, len(header_neuro_musculo_state)))

    header_mtu_architecture = [
        'frame', 'time',
        'tibialis_mtu_length', 'soleus_mtu_length', 'gastrocnemius_mtu_length',
        'tibialis_fiber_length', 'soleus_fiber_length', 'gastrocnemius_fiber_length',
        'tibialis_pennation_angle', 'soleus_pennation_angle', 'gastrocnemius_pennation_angle',
        'tibialis_tendon_length', 'soleus_tendon_length', 'gastrocnemius_tendon_length'
    ]
    simulated_data_mtu_architecture = np.zeros((n_trials, len(header_mtu_architecture)))

    header_mtu_forces = ['frame', 'time',
                         'torque_simulated',
                         'tendon_force_tibialis', 'tendon_force_soleus', 'tendon_force_gastrocnemius',
                         'fiber_force_tibialis', 'fiber_force_soleus', 'fiber_force_gastrocnemius',
                         'fiber_active_force_tibialis', 'fiber_active_force_soleus', 'fiber_active_force_gastrocnemius',
                         'fiber_passive_force_tibialis', 'fiber_passive_force_soleus',
                         'fiber_passive_force_gastrocnemius']
    simulated_data_mtu_forces = np.zeros((n_trials, len(header_mtu_forces)))

    n_trials_fail = 0  # compteur of trials non - optimized(p > 10e-5)
    n_trials_succeds = 0  # compteur of trials optimized(p > 10e-5)

    # ========= 1 Muscle-Tendon Architecture Equations   ========= #
    # ========= 1.1 Musculo skeletical configuration during trial(input)   ========= #

    q5_num = np.deg2rad(q5_num)
    q6_num = np.deg2rad(q6_num)

    # ========= 1.2 Neuronal activation(input)   ========= #
    # Muscle activation [a_tibialis_num, a_soleus_num, a_gastrocnemius_num]

    # ========= 1.3 Muscle Tendon Parameters(ℓom, φo, Fom, ℓst)(input)   ========= #
    muscle_tendon_parameters_ta = np.array(muscle_tendon_parameters_num)[[0, 3, 6, 9]]
    muscle_tendon_parameters_sol = np.array(muscle_tendon_parameters_num)[[1, 4, 7, 10]]
    muscle_tendon_parameters_gast = np.array(muscle_tendon_parameters_num)[[2, 5, 8, 11]]

    # ========= 2.1 data generator   ========= #

    trials = 0

    for i in range(n_trials):  # for each frame
        # 2.1 Progression
        trials += 1
        percent = round(trials / n_trials, 2)
        print(f"Progressing: {percent * 100:.0f} %")

        # 2.2 Neuromusculoskeletal configuration
        a_num = [a_tibialis_num[i], a_soleus_num[i], a_gastrocnemius_num[i]]
        q_num = [0, 0, 0, 0, q5_num[i], q6_num[i]]

        musculoskeletal_states_num = q_num + list(skeleton_num)
        neuromusculoskeletal_state_num = np.concatenate([a_num, musculoskeletal_states_num])

        # 2.3 UMT length
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_num)
        tibialis_mtu_length = float(mtu_length[0])
        soleus_mtu_length = float(mtu_length[1])
        gastrocnemius_mtu_length = float(mtu_length[2])

        # 2.4 muscle activities
        a_tibialis = a_num[0]
        a_soleus = a_num[1]
        a_gastrocnemius = a_num[2]

        # 2.5 Solver
        x_opt_tibialis, state_tibialis = root_muscle_dynamics(a_tibialis,
                                                              tibialis_mtu_length,
                                                              muscle_tendon_parameters_ta,
                                                              'tibialis',
                                                              casadi_function)

        x_opt_soleus, state_soleus = root_muscle_dynamics(a_soleus,
                                                          soleus_mtu_length,
                                                          muscle_tendon_parameters_sol,
                                                          'soleus',
                                                          casadi_function)

        x_opt_gastrocnemius, state_gastrocnemius = root_muscle_dynamics(a_gastrocnemius,
                                                                        gastrocnemius_mtu_length,
                                                                        muscle_tendon_parameters_gast,
                                                                        'gastrocnemius',
                                                                        casadi_function)

        if state_tibialis == 'fail' or state_soleus == 'fail' or state_tibialis == 'fail' or state_gastrocnemius == 'fail':
            n_trials_fail += 1
        else:
            n_trials_succeds += 1

            rooted_fiber_length = np.array(
                [x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0]]).flatten()
            rooted_pennation_angle = np.array(
                [x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0]]).flatten()
            rooted_tendon_length = np.array(
                [x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()

            rooted_variables = np.concatenate([rooted_fiber_length, rooted_pennation_angle, rooted_tendon_length])

            all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

            # 2.7 extract variable
            # about architecture
            # rooted = vertcat(tendon length,fiber lenght,,pennationAngle)
            tibialis_fiber_length = float(x_opt_tibialis[0])
            tibialis_pennation_angle = x_opt_tibialis[1]
            tibialis_pennation_angle = float(np.rad2deg(tibialis_pennation_angle))
            tibialis_tendon_length = float(x_opt_tibialis[2])

            soleus_fiber_length = float(x_opt_soleus[0])
            soleus_pennation_angle = x_opt_soleus[1]
            soleus_pennation_angle = float(np.rad2deg(soleus_pennation_angle))
            soleus_tendon_length = float(x_opt_soleus[2])

            gastrocnemius_fiber_length = float(x_opt_gastrocnemius[0])
            gastrocnemius_pennation_angle = x_opt_gastrocnemius[1]
            gastrocnemius_pennation_angle = float(np.rad2deg(gastrocnemius_pennation_angle))
            gastrocnemius_tendon_length = float(x_opt_gastrocnemius[2])

            # about forces
            ankle_torque = casadi_function['get_joint_moment'](all_state, muscle_tendon_parameters_num)
            torque_simulated = float(ankle_torque)

            tendon_force = casadi_function['get_tendon_force'](all_state, muscle_tendon_parameters_num)
            tendon_force_tibialis = float(tendon_force[0])
            tendon_force_soleus = float(tendon_force[1])
            tendon_force_gastrocnemius = float(tendon_force[2])

            fiber_force = casadi_function['get_fiber_force'](all_state, muscle_tendon_parameters_num)
            fiber_force_tibialis = float(fiber_force[0])
            fiber_force_soleus = float(fiber_force[1])
            fiber_force_gastrocnemius = float(fiber_force[2])

            fiber_active_force = casadi_function['get_fiber_active_force'](all_state, muscle_tendon_parameters_num)
            fiber_active_force_tibialis = float(fiber_active_force[0])
            fiber_active_force_soleus = float(fiber_active_force[1])
            fiber_active_force_gastrocnemius = float(fiber_active_force[2])

            fiber_passive_force = casadi_function['get_fiber_passive_force'](all_state, muscle_tendon_parameters_num)
            fiber_passive_force_tibialis = float(fiber_passive_force[0])
            fiber_passive_force_soleus = float(fiber_passive_force[1])
            fiber_passive_force_gastrocnemius = float(fiber_passive_force[2])

            # 2.output
            # header_neuro_musculo_state = [
            #     'frame', 'time',
            #     'q5', 'q6',
            #     'a_tibialis', 'a_soleus', 'a_gastrocnemius'
            # ]

            simulated_data_neuro_musculo_state[i, :] = [
                trials, time_num[i],
                np.rad2deg(q_num[4]), np.rad2deg(q_num[5]),
                a_num[0] * 100, a_num[1] * 100, a_num[2] * 100
            ]
            # header_mtu_architecture = [
            #     'frame', 'time',
            #     'tibialis_mtu_length', 'soleus_mtu_length', 'gastrocnemius_mtu_length',
            #     'tibialis_fiber_length', 'soleus_fiber_length', 'gastrocnemius_fiber_length',
            #     'tibialis_pennation_angle', 'soleus_pennation_angle','gastrocnemius_pennation_angle',
            #     'tibialis_tendon_length', 'soleus_tendon_length', 'gastrocnemius_tendon_length'
            # ]
            simulated_data_mtu_architecture[i, :] = [
                trials, time_num[i],
                tibialis_mtu_length, soleus_mtu_length, gastrocnemius_mtu_length,
                tibialis_fiber_length, soleus_fiber_length, gastrocnemius_fiber_length,
                tibialis_pennation_angle, soleus_pennation_angle, gastrocnemius_pennation_angle,
                tibialis_tendon_length, soleus_tendon_length, gastrocnemius_tendon_length
            ]

            # header_mtu_forces = [
            #     'frame', 'time',
            #     'torque_simulated',
            #     'tendon_force_tibialis', 'tendon_force_soleus', 'tendon_force_gastrocnemius',
            #     'fiber_force_tibialis', 'fiber_force_soleus', 'fiber_force_gastrocnemius',
            #     'fiber_active_force_tibialis', 'fiber_active_force_soleus', 'fiber_active_force_gastrocnemius',
            #     'fiber_passive_force_tibialis', 'fiber_passive_force_soleus', 'fiber_passive_force_gastrocnemius'
            # ]
            simulated_data_mtu_forces[i, :] = [
                trials, time_num[i],
                torque_simulated,
                tendon_force_tibialis, tendon_force_soleus, tendon_force_gastrocnemius,
                fiber_force_tibialis, fiber_force_soleus, fiber_force_gastrocnemius,
                fiber_active_force_tibialis, fiber_active_force_soleus, fiber_active_force_gastrocnemius,
                fiber_passive_force_tibialis, fiber_passive_force_soleus, fiber_passive_force_gastrocnemius
            ]

    print('\n================== simulation exit ==================')
    print(f"fail trials: {(n_trials_fail / n_trials) * 100:.2f} %")

    return header_neuro_musculo_state, simulated_data_neuro_musculo_state, header_mtu_architecture, simulated_data_mtu_architecture, header_mtu_forces, simulated_data_mtu_forces

def add_noise(data, sigma_torque=0.1, sigma_length=0.002, sigma_angle_deg=3.0,
              add_internal_noise=False, rng=None):
    """
    Add Gaussian sensor noise to hypothetical_data IN PLACE.

    Noise is applied to measurable quantities only (torque, joint angles).
    Internal Hill-model variables (fiber length, pennation, tendon length,
    activations) are left untouched by default — set add_internal_noise=True
    to also perturb them, e.g. to simulate ultrasound measurement noise on
    fiber length and pennation.

    Parameters
    ----------
    data : np.ndarray, shape (15, n_trials)
        Hypothetical data matrix. See module docstring for row layout.
        Modified in place.
    sigma_torque : float, default 0.1
        Std of torque noise [N.m].
    sigma_length : float, default 0.002
        Std of length noise [m] (= 0.2 cm).
    sigma_angle_deg : float, default 3.0
        Std of angle noise [degrees]. Internally converted to radians.
    add_internal_noise : bool, default False
        If True, also adds noise to fiber length, pennation angle, and
        tendon length (typical for ultrasound-based estimates).
    rng : np.random.Generator, optional
        For reproducibility. If None, uses np.random.default_rng().

    Returns
    -------
    np.ndarray
        The same data array (modified in place), returned for chaining.
    """
    # Indices des lignes dans la matrice (15, n_trials)
    ROW_TORQUE = 0  # 1 ligne   : couple [N.m]
    ROW_Q = slice(1, 3)  # 2 lignes  : angles articulaires [rad]
    ROW_ACTIVATION = slice(3, 6)  # 3 lignes  : activations [0, 1]
    ROW_FIBER_LENGTH = slice(6, 9)  # 3 lignes  : longueurs fibres [m]
    ROW_PENNATION = slice(9, 12)  # 3 lignes : angles pennation [rad]
    ROW_TENDON_LENGTH = slice(12, 15)  # 3 lignes : longueurs tendon [m]


    rng = rng if rng is not None else np.random.default_rng()
    sigma_angle_rad = np.deg2rad(sigma_angle_deg)
    n_trials = data.shape[1]

    # === Bruit sur les mesures (capteur) === #
    # Couple : 1 ligne
    data[ROW_TORQUE, :] += rng.normal(0.0, sigma_torque, size=n_trials)

    # Angles articulaires : 2 lignes
    data[ROW_Q, :] += rng.normal(0.0, sigma_angle_rad, size=(2, n_trials))

    # === Bruit sur les estimations échographiques (optionnel) === #
    if add_internal_noise:
        # Longueurs de fibres : 3 lignes
        data[ROW_FIBER_LENGTH, :] += rng.normal(0.0, sigma_length, size=(3, n_trials))

        # Angles de pennation : 3 lignes
        data[ROW_PENNATION, :] += rng.normal(0.0, sigma_angle_rad, size=(3, n_trials))

        # Longueurs de tendon : 3 lignes
        data[ROW_TENDON_LENGTH, :] += rng.normal(0.0, sigma_length, size=(3, n_trials))

    return data

def optimization_nlp(data, initial_guess, lower_band, upper_band, skeleton_num,
                     muscle_tendon_parameters_num, unknown_parameters, casadi_function):
    """
    Identification des paramètres muscle-tendon (ℓom, φo, Fom, ℓst) via NLP.

    Args:
        data (np.ndarray): Shape (15, n_trials).
            Lignes :
              [0]     : torque mesuré (N.m)
              [1:3]   : angles articulaires q (rad)
              [3:6]   : activations musculaires [0, 1]
              [6:9]   : longueurs de fibres mesurées (m)
              [9:12]  : angles de pennation mesurés (rad)
              [12:15] : longueurs de tendon mesurées (m)
        initial_guess, lower_band, upper_band : paramètres muscle-tendon (taille 12)
        skeleton_num : géométrie squelette (.osim)
        muscle_tendon_parameters_num : valeurs de référence pour comparaison
        unknown_parameters : vecteur SX des 12 inconnues
        casadi_function : dict de fonctions CasADi

    Returns:
        xopt (np.ndarray): Shape (12,) — (ℓom, φo, Fom, ℓst) estimés.
    """

    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data doit être (15, n_trials), reçu {data.shape}"

    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)
    initial_guess = np.asarray(initial_guess, dtype=float)
    lower_band = np.asarray(initial_guess, dtype=float)
    upper_band = np.asarray(initial_guess, dtype=float)

    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band à : {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ NLP set up ============ #
    w, w0, lbw, ubw = [], [], [], []
    j = 0
    g, lbg, ubg = [], [], []

    up = unknown_parameters

    e_torque, e_fiber, e_pennation = [], [], []

    # Paramètres muscle-tendon inconnus (12)
    w += [up]
    w0 += list(initial_guess)
    lbw += list(lower_band)
    ubw += list(upper_band)

    # Poids dans la fonction coût
    w_torque = 1  # N.m
    w_length = 0.005  # mm
    w_angle = (3 / 180) * np.pi  # rad

    for trial in range(n_trials):
        # --- Extraction (colonne = trial) ---
        data_trial = data[:, trial]

        mesured_torque = data_trial[0]  # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]  # [0, 1]
        mesured_fiber_length = data_trial[6:9]  # m
        mesured_pennation_angle = data_trial[9:12]  # rad
        mesured_tendon_length = data_trial[12:15]  # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        # Longueur MTU
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        # --- Variables de décision du trial ---
        tendon_length_k = SX.sym(f"Tendon_Length_{trial + 1}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial + 1}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial + 1}", n_muscle)

        # --- Bornes physiologiques (robustes aux signes) ---
        fl_meas = np.abs(mesured_fiber_length)
        tl_meas = np.abs(mesured_tendon_length)

        # fiber length : positive, autour de la mesure
        lb_fl = fl_meas * 0.3
        ub_fl = fl_meas * 2.0

        # pennation : bornes physiologiques, signe libre
        lb_pa = np.full(n_muscle, -np.pi / 2 + 1e-3)
        ub_pa = np.full(n_muscle, np.pi / 2 - 1e-3)

        # tendon length : positive, range resserré
        lb_tl = tl_meas * 0.05
        ub_tl = tl_meas * 1.05

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        # Initial guess avec longueurs rendues positives
        w0_k = np.concatenate([fl_meas, mesured_pennation_angle, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        # Sécurité : bornes valides et w0 dans les bornes
        assert np.all(lb_block <= ub_block), \
            f"Bornes inversées au trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 hors bornes au trial {trial} : " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- equilibrium constraints (9 par trial) ---
        k = vertcat(a_trial, mtu_length, up)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Simulation du torque ---
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, unknown_parameters)

        # --- errors ---
        e_torque_trials = mesured_torque - torque_simulated
        e_fiber_trials = mesured_fiber_length - fiber_length_k
        e_pennation_trials = mesured_pennation_angle - pennation_angle_k

        j += w_torque * e_torque_trials ** 2
        j += sum1(w_length * e_fiber_trials ** 2)
        j += sum1(w_angle * e_pennation_trials ** 2)

        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    # ============ Assemblage ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    # Vérifications globales
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
    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp)
    print(solver)

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    # ============ Extraction ============ #
    w_opt = sol['x'].full().flatten()
    cost = sol['f'].full().item()

    param_opt = w_opt[:12]
    err_param = np.abs(muscle_tendon_parameters_num - param_opt)

    print("\n\ndifference between input and Estimated Muscle-Tendon Parameters :")
    print(f"Number of trials : {n_trials}")
    print(f"diff ℓom : {err_param[0:3]}")
    print(f"diff φo  : {err_param[3:6]}")
    print(f"diff Fom : {err_param[6:9]}")
    print(f"diff ℓst : {err_param[9:12]}")

    print("Estimated Muscle-Tendon Parameters :")
    print(f"Cost : {cost}")
    print(f"ℓom : {param_opt[0:3]}")
    print(f"φo  : {param_opt[3:6]}")
    print(f"Fom : {param_opt[6:9]}")
    print(f"ℓst : {param_opt[9:12]}")

    return param_opt


def generate_estimated_data(data, skeleton_num, muscle_tendon_parameters_num,
                            casadi_function,
                            output_dir='data_estime', filename='data_estime.csv',
                            save_npy=True, save_csv=True, verbose=True):
    """
    Génère un dataset "estimé" par le modèle direct à partir d'un dataset d'entrée.

    Pour chaque trial de `data`, on réutilise les angles articulaires et les
    activations mesurées, puis on résout l'équilibre muscle-tendon pour
    recalculer :
        - les longueurs de fibres
        - les angles de pennation
        - les longueurs de tendon
        - le torque articulaire
    avec les paramètres muscle-tendon fournis.

    Utile pour :
        - valider la calibration (mesure vs modèle sur mêmes q, a)
        - générer un jeu de référence pour évaluer l'optim

    Args:
        data (np.ndarray): Shape (15, n_trials). Seules les lignes [1:3] (q)
            et [3:6] (a) sont utilisées ; les autres sont recalculées.
            Unités : q en rad, a ∈ [0, 1].
        skeleton_num : géométrie squelette (.osim scalé).
        muscle_tendon_parameters_num (array-like, 12): (ℓom, φo, Fom, ℓst) × 3 muscles.
            Ordre : [TA, SOL, GAST] pour chaque paramètre.
        casadi_function (dict): fonctions CasADi.
        output_dir (str): dossier de sortie. Créé si inexistant.
        filename (str): nom de base du fichier.
        save_npy, save_csv (bool): formats de sauvegarde.
        verbose (bool): logs de progression.

    Returns:
        header (list[str]): noms des lignes.
        data_est (np.ndarray): Shape (15, n_success_trials), mêmes unités que `data`
            (angles en rad, longueurs en m, torque en N.m).
        path_csv (str | None): chemin du CSV écrit.
    """

    header = [
        'ankle_torque', 'q_knee', 'q_ankle',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius',
        'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
        'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
        'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
    ]

    # --- Validation input ---
    assert data.shape[0] == 15, f"data doit être (15, n_trials), reçu {data.shape}"
    n_trials = data.shape[1]

    # Paramètres muscle par muscle
    mtp = np.array(muscle_tendon_parameters_num)
    mtp_ta = mtp[[0, 3, 6, 9]]
    mtp_sol = mtp[[1, 4, 7, 10]]
    mtp_gast = mtp[[2, 5, 8, 11]]

    data_est = np.zeros((15, n_trials))
    n_fail, n_ok = 0, 0

    if verbose:
        print(f'\n================== Generating estimated data ==================')
        print(f"  n_trials input : {n_trials}")

    for trial in range(n_trials):
        if verbose and (trial + 1) % max(1, n_trials // 20) == 0:
            print(f"  progress : {100 * (trial + 1) / n_trials:.0f} %")

        col = data[:, trial]

        # --- Récupération des inputs (q, a) ---
        q_trial = col[1:3]  # [q_knee, q_ankle] en rad
        a_trial = col[3:6]  # [a_TA, a_SOL, a_GAST]
        a_ta, a_sol, a_gast = a_trial

        # État musculosquelettique complet
        q_num = [0, 0, 0, 0] + list(q_trial)
        musc_states = q_num + list(skeleton_num)
        neuromusc_state = np.concatenate([a_trial, musc_states])

        # Longueur MTU
        mtu_length = casadi_function['get_mtu_length'](musc_states)
        l_ta = float(mtu_length[0])
        l_sol = float(mtu_length[1])
        l_gast = float(mtu_length[2])

        # --- Équilibre muscle-tendon par muscle ---
        x_ta, s_ta = root_muscle_dynamics(a_ta, l_ta, mtp_ta, 'tibialis', casadi_function)
        x_sol, s_sol = root_muscle_dynamics(a_sol, l_sol, mtp_sol, 'soleus', casadi_function)
        x_gast, s_gast = root_muscle_dynamics(a_gast, l_gast, mtp_gast, 'gastrocnemius', casadi_function)

        if 'fail' in (s_ta, s_sol, s_gast):
            n_fail += 1
            continue

        # Variables rootées
        fl = np.array([x_ta[0, 0], x_sol[0, 0], x_gast[0, 0]]).flatten()
        penn = np.array([x_ta[1, 0], x_sol[1, 0], x_gast[1, 0]]).flatten()
        tl = np.array([x_ta[2, 0], x_sol[2, 0], x_gast[2, 0]]).flatten()

        rooted = np.concatenate([fl, penn, tl])
        all_state = np.concatenate([neuromusc_state, rooted])

        # Torque estimé
        ankle_torque = float(casadi_function['get_joint_moment'](all_state, mtp))

        # --- Remplissage (mêmes unités que `data`) ---
        data_est[:, n_ok] = [
            ankle_torque,
            q_trial[0], q_trial[1],  # rad
            a_ta, a_sol, a_gast,
            fl[0], fl[1], fl[2],  # m
            penn[0], penn[1], penn[2],  # rad
            tl[0], tl[1], tl[2],  # m
        ]
        n_ok += 1

    data_est = data_est[:, :n_ok]

    if verbose:
        print(f"  → {n_ok} / {n_trials} trials OK ({100 * n_fail / n_trials:.1f} % fail)")

    # --- Sauvegarde ---
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.splitext(filename)[0]
    path_csv = None

    if save_csv:
        # Pour le CSV : transposer en (n_trials, 15) + angles en degrés pour lisibilité
        path_csv = os.path.join(output_dir, base + '.csv')
        df = pd.DataFrame(data_est.T, columns=header)
        df['q_knee'] = np.rad2deg(df['q_knee'])
        df['q_ankle'] = np.rad2deg(df['q_ankle'])
        for m in ['tibialis', 'soleus', 'gastrocnemius']:
            df[f'pennation_angle_{m}'] = np.rad2deg(df[f'pennation_angle_{m}'])
        df.to_csv(path_csv, index=False)
        if verbose:
            print(f"  → CSV (deg, n_trials × 15) : {path_csv}")

    if save_npy:
        # NPY : shape (15, n_trials) en unités SI (rad, m, N.m) — directement utilisable
        path_npy = os.path.join(output_dir, base + '.npy')
        np.save(path_npy, data_est)
        if verbose:
            print(f"  → NPY (rad, 15 × n_trials) : {path_npy}")

    return header, data_est, path_csv

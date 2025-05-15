import os
import xml.etree.ElementTree as ET
import numpy as np
from casadi import SX, DM, vertcat, horzcat, Function, sqrt, exp, if_else, sum1, jacobian, rootfinder
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', depending on your environment

def model_osim2mat(oism_path, oism_file):
    """
    Extract muscle-tendon and joint parameters from an OpenSim .osim/.oism file.
    this function is for extract from opensim scale subjects -->
    - known_parameters_num :
       - segment_geometry_num
       - muscle_origin_num
       - muscle_insersion_num
    output type (30,1) : known_parameters_num = [segment_geometry_num; muscle_origin_num; muscle_insersion_num];
    output type (33,1) : known_parameters_num = [segment_geometry_num; muscle_origin_num; muscle_insersion_num, muscle_viaPoint];

    - muscle_tendon_parameters_num :
       - l0m_num
       - phi0_num
       - f0m_num
       - lst_num
    output type (12,1) : muscle_tendon_parameters_num =  [l0m_num , phi0_num , f0m_num , lst_num];

    Returns:
        known_parameters_num (np.ndarray): Shape (33,) including segment geometry and muscle attachment points.
        muscle_tendon_parameters_num (np.ndarray): Shape (12,) including l0m, phi0, f0m, lst values.
    """

    # Target bones and muscles
    interest_muscle_tendon_parameters = ['l0m','phi0','f0m','lst']
    interst_body = ['tibia_r','talus_r','calcn_r','toes_r']
    interest_joints = ['knee_r', 'ankle_r', 'subtalar_r', 'mtp_r']
    interest_muscles = ['tib_ant_r', 'soleus_r', 'gast_r']

    try:
        # Load XML file
        full_path = os.path.join(oism_path, oism_file)
        tree = ET.parse(full_path)
        root = tree.getroot()
        model = root.find('Model')

        # get in the model
        if model is not None:
            # === segment_geometry_num ===
            bodies = {}  # get information for all bodies (segment_geometry_num)
            body_set = model.find('BodySet')
            if body_set is not None:
                body_object =  body_set.find('objects')
                if body_object is not None:
                    for body in body_object:
                        body_name = body.attrib.get('name')
                        joint = body.find('Joint')
                        if joint is not None:
                            customJoint = joint.find('CustomJoint')
                            if customJoint is not None:
                                loc_elem = customJoint.find('location_in_parent')
                                if loc_elem is not None:
                                    x, y, z = map(float, loc_elem.text.strip().split())
                                    bodies[body_name] = {'X': x, 'Y': y, 'Z': z}

            # === muscle_origin_num & muscle_insersion_num & muscle_tendon_parameters_num ===
            muscle_pass_point = {}  # get information for all muscles
            muscle_tendon_parameters = {} # get information muscle_tendon_parameters_num

            muscle_set = model.find('ForceSet')
            if muscle_set is not None:
                muscle_object = muscle_set.find('objects')
                for muscle in muscle_object:
                    # === geometry === #
                    muscle_name = muscle.attrib.get('name')
                    muscle_pass_point[muscle_name] = {}
                    muscle_geometry = muscle.find('GeometryPath')
                    if muscle_geometry is not None:
                        muscle_path_point_set = muscle_geometry.find('PathPointSet')
                        if muscle_path_point_set is not None:
                            muscle_path_point_objects = muscle_path_point_set.find('objects')
                            for pass_point in muscle_path_point_objects:
                                pass_point_name = pass_point.attrib.get('name')
                                loc_elem = pass_point.find('location')
                                if loc_elem is not None:
                                    x, y, z = map(float, loc_elem.text.strip().split())
                                    muscle_pass_point[muscle_name][pass_point_name] = {'X': x, 'Y': y, 'Z': z}



                    # === Parameters === #
                    def get_float_value(tag):
                        elem = muscle.find(tag)
                        return float(elem.text) if elem is not None and elem.text else None

                    muscle_tendon_parameters[muscle_name] = {
                        'l0m': get_float_value('optimal_fiber_length'),
                        'phi0': get_float_value('pennation_angle_at_optimal'),
                        'f0m': get_float_value('max_isometric_force'),
                        'lst': get_float_value('tendon_slack_length'),
                    }
                    # output muscle_tendon_parameters_num[muscle_name] = {'l0m_num': optimal_fiber_length,
                    # 'phi0_num': pennation_angle_at_optimal,
                    # 'f0m_num': max_isometric_force,
                    # 'lst_num': tendon_slack_length}
                    # [l0m_num , phi0_num , f0m_num , lst_num];

        # === Combine gastrocnemius ===
            # === Parameters === #
        gast = 'gast_r'
        muscle_tendon_parameters[gast] = {'l0m': np.mean([
            muscle_tendon_parameters['lat_gas_r']['l0m'],
                muscle_tendon_parameters['med_gas_r']['l0m']
            ]),
            'phi0': np.mean([
                muscle_tendon_parameters['lat_gas_r']['phi0'],
                muscle_tendon_parameters['med_gas_r']['phi0']
            ]),
            'f0m': np.sum([
                muscle_tendon_parameters['lat_gas_r']['f0m'],
                muscle_tendon_parameters['med_gas_r']['f0m']
            ]),
            'lst': np.mean([
                muscle_tendon_parameters['lat_gas_r']['lst'],
                muscle_tendon_parameters['med_gas_r']['lst']
            ])
        }

            # === Geometry === #
        muscle_pass_point[gast] = {}
        muscle_pass_point[gast]['gast-P1']= {}
        muscle_pass_point[gast]['gast-P1'] = {'X': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P1']['X'],
            muscle_pass_point['med_gas_r']['med_gas_r-P1']['X']
            ]),
            'Y': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P1']['Y'],
            muscle_pass_point['med_gas_r']['med_gas_r-P1']['Y']
            ]),
            'Z': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P1']['Z'],
            muscle_pass_point['med_gas_r']['med_gas_r-P1']['Z']
            ])
            }

        muscle_pass_point[gast]['gast-P2']= {}
        muscle_pass_point[gast]['gast-P2'] = {'X': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P2']['X'],
            muscle_pass_point['med_gas_r']['med_gas_r-P2']['X']
            ]),
            'Y': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P2']['Y'],
            muscle_pass_point['med_gas_r']['med_gas_r-P2']['Y']
            ]),
            'Z': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P2']['Z'],
            muscle_pass_point['med_gas_r']['med_gas_r-P2']['Z']
            ])
            }

        muscle_pass_point[gast]['gast-P3']= {}
        muscle_pass_point[gast]['gast-P3'] = {'X': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P3']['X'],
            muscle_pass_point['med_gas_r']['med_gas_r-P3']['X']
            ]),
            'Y': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P3']['Y'],
            muscle_pass_point['med_gas_r']['med_gas_r-P3']['Y']
            ]),
            'Z': np.mean([
            muscle_pass_point['lat_gas_r']['lat_gas_r-P3']['Z'],
            muscle_pass_point['med_gas_r']['med_gas_r-P3']['Z']
            ])
            }

        # === output mtu parameters === (12 values)#
        muscle_tendon_parameters_num = []
        for param in interest_muscle_tendon_parameters:
            for muscle in interest_muscles:
                value = round(muscle_tendon_parameters[muscle][param], 3)
                muscle_tendon_parameters_num.append(float(value))

        # === Segment geometry (12 values) === #
        segment_geometry_num = []
        for body_name in interst_body:
            coord = bodies[body_name]
            segment_geometry_num.extend([coord['X'], coord['Y'], 0.0])

        # === Muscle origins (3 muscles × 3 values) ===
        muscle_origin = []
        for m in ['tib_ant_r', 'soleus_r', gast]:
            first_point_name = next(iter(muscle_pass_point[m]))
            first_point_coords = muscle_pass_point[m][first_point_name]
            muscle_origin.extend([first_point_coords['X'], first_point_coords['Y'], 0.0])

        # === Muscle insertions (3 muscles × 3 values) ===
        muscle_insertion = []
        for m in ['tib_ant_r', 'soleus_r', gast]:
            last_point_name = next(reversed(muscle_pass_point[m]))
            last_point_coords = muscle_pass_point[m][last_point_name]
            muscle_insertion.extend([last_point_coords['X'], last_point_coords['Y'], 0.0])

        # === Muscle via point (only tib_ant_r) ===
        muscle_via = [muscle_pass_point['tib_ant_r']['tib_ant_r-P2']['X'],
                      muscle_pass_point['tib_ant_r']['tib_ant_r-P2']['Y'],
                      0.0]

        known_parameters_num = np.array(segment_geometry_num + muscle_origin + muscle_insertion + muscle_via)

        return known_parameters_num, muscle_tendon_parameters_num

    except Exception as e:
        print(f"Error extracting parameters: {e}")
        return np.zeros(33), np.zeros(12)

def DeGrooteFunction():
    """
        Casadi funtion for neuromusculo-model
    # ========= variables  ========= #
    "a" is neuromuscular activation (between 0 and 1) [3]
    "q" is the spatial skeleton configuration (x, y, z, theta_hip,theta_knee, theta_ankle) [6]
    "knownparameter" is subject segment length and muscle insertion position. (segment_geometry, muscle_insertion, Local_ViaPoint_tibialis_anterior)[33]
    "musculoskeletal_states" = [q, known_parameter]
    "neuromusculoskeletal_state" = [a, q, known_parameter]
    "muscleTendonParameters" : musculoskeletal parameters that are generally assumed not to change: maximal isometric muscle force (Fom), optimal fiber length (ℓom), tendon slack length (ℓst), and pennation angle at optimal fiber length (φo).
    "rootedvariables" = (tendonLengthening, fiberLength)
    "all_state" = (neuromusculoskeletal_state, rootedvariables)

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


    nMuscles = 3 ;
    NameMuscles = {'TibialisAnterior','Soleus','Gastrocnemius'} ;
    nBones = 4 ; NameBones = {'Thigh','Leg','Talus','Calcaneus'} ;
    nJoints = 4 ; NameJoints = {'Hip','Knee','Ankle','Subtalar'} ;


        Returns:
            known_parameters_num (np.ndarray): Shape (33,) including segment geometry and muscle attachment points.
            muscle_tendon_parameters_num (np.ndarray): Shape (12,) including l0m, phi0, f0m, lst values.
        """
    # ========= usful functions  ========= #
    def Rototranslation_Rz(x, y, z, theta):
        # Build 4x4 homogeneous transformation matrix for rotation about Z and translation
        from casadi import cos, sin, DM
        R = SX.zeros(4, 4)
        R[0, 0] = cos(theta)
        R[0, 1] = -sin(theta)
        R[1, 0] = sin(theta)
        R[1, 1] = cos(theta)
        R[0, 3] = x
        R[1, 3] = y
        R[2, 3] = z
        R[2, 2] = 1
        R[3, 3] = 1
        return R

    def Rototranslate(R, p_local):
        # Apply rotation and translation from 4x4 matrix
        p_local_h = vertcat(p_local, SX(1))
        p_global_h = R @ p_local_h
        return p_global_h[:3]

    # Constants
    nMuscles = 3

    # Define states
    x = SX.sym('x')
    y = SX.sym('y')
    z = SX.sym('z')
    theta_hip = SX.sym('theta_hip')
    theta_knee = SX.sym('theta_knee')
    theta_ankle = SX.sym('theta_ankle')
    q = vertcat(x, y, z, theta_hip, theta_knee, theta_ankle)
    # Segment geometry parameters
    thigh = SX.sym('thigh', 3)
    leg = SX.sym('leg', 3)
    talus = SX.sym('talus', 3)
    clac = SX.sym('foot', 3)
    segment_geometry = vertcat(thigh, leg, talus, clac)

    # Muscle insertions and origins
    Local_Insertion_ta = SX.sym('Local_Insertion_tibialis_anterior', 3)
    Local_Origin_ta = SX.sym('Local_Origin_tibialis_anterior', 3)
    Local_ViaPoint_ta = SX.sym('Local_ViaPoint_tibialis_anterior', 3)

    Local_Insertion_sol = SX.sym('Local_Insertion_soleus', 3)
    Local_Origin_sol = SX.sym('Local_Origin_soleus', 3)
    Local_Insertion_gas = SX.sym('Local_Insertion_gastrocnemius', 3)
    Local_Origin_gas = SX.sym('Local_Origin_gastrocnemius', 3)

    muscle_insertion = vertcat(Local_Origin_ta, Local_Origin_sol, Local_Origin_gas,
                               Local_Insertion_ta, Local_Insertion_sol, Local_Insertion_gas)

    known_parameters = vertcat(segment_geometry, muscle_insertion, Local_ViaPoint_ta)

    # Rototranslation
    R_0_thigh = Rototranslation_Rz(x, y, z, 0)
    R_thigh_leg = Rototranslation_Rz(thigh[0], thigh[1], thigh[2], -theta_knee)
    R_leg_talus = Rototranslation_Rz(leg[0], leg[1], leg[2], -theta_ankle)
    R_talus_calc = Rototranslation_Rz(talus[0], talus[1], talus[2], 0)

    # express in thigh
    R_thigh_talus = R_thigh_leg @ R_leg_talus
    R_thigh_calc = R_thigh_talus @ R_talus_calc


    # Origins and insertions
    Origin_ta = Rototranslate(R_thigh_leg, Local_Origin_ta)
    Insertion_ta = Rototranslate(R_thigh_calc, Local_Insertion_ta)
    ViaPoint_ta = Rototranslate(R_thigh_leg, Local_ViaPoint_ta)

    Origin_sol = Rototranslate(R_thigh_leg, Local_Origin_sol)
    Insertion_sol = Rototranslate(R_thigh_calc, Local_Insertion_sol)

    Origin_gas = Local_Origin_gas
    Insertion_gas = Rototranslate(R_thigh_calc, Local_Insertion_gas)

    # Forward kinematics
    Origin = horzcat(Origin_ta, Origin_sol, Origin_gas)
    Insertion = horzcat(Insertion_ta, Insertion_sol, Insertion_gas)
    ViaPoint = horzcat(ViaPoint_ta)
    HJC = [0,0,0]
    KJC = R_thigh_leg[0:3, 3]
    AJC = R_thigh_talus[0:3, 3]
    TJC = Rototranslate(R_thigh_calc, clac)
    CALC = R_thigh_calc[0:3, 3]
    Markers = horzcat(HJC, KJC, AJC, TJC, CALC)

    musculoskeletal_states = vertcat(q, known_parameters)

    # ============ forward kinematics function ============ #
    forwardKinematics = Function("forwardKinematics",
                                [musculoskeletal_states],
                                [Origin, Insertion, ViaPoint, Markers],
                                ["musculoskeletal_states"],
                                ["Origin", "Insertion", "ViaPoint", "Markers"])

    # ============ umt length: function ============ #
    umtLength = vertcat(
        sqrt(sum1((Insertion_ta - ViaPoint_ta) ** 2)) +
        sqrt(sum1((ViaPoint_ta - Origin_ta) ** 2)),  # tibialis
        sqrt(sum1((Insertion_sol - Origin_sol) ** 2)),  # soleus
        sqrt(sum1((Insertion_gas - Origin_gas) ** 2))  # gastroc
    )


    getMTULength = Function("getMTULength",
                                [musculoskeletal_states],
                                [umtLength],
                                ["musculoskeletal_states"],
                                ["umtLength (tibialis, soleus, gastrocnemius)"])

    # ============ moment arm: function  ============ #
    momentArm = jacobian(umtLength, q)

    getMomentArm = Function("getMomentArm",
                            [musculoskeletal_states],
                            [momentArm],
                            ["musculoskeletal_states"],
                            ["moment Arm (tibialis, soleus, gastrocnemius)"])

                        # ========= 2. Dynamics equations  ========= #
    # ========= ref  ========= #
    # Groote F De, Kinney AL, Rao A V, Fregly BJ. Evaluation of Direct Collocation
    # Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem.
    # Ann Biomed Eng. 2016; 44: 2922–2936. https://doi.org/10.1007/s10439-016-1591-9
    # PMID: 27001399
    nMuscles = 3

    # ========= 2.1 Muscle Tendon Parameters(ℓom, φo, Fom, ℓst)  ========= #
    optimalFiberLength = SX.sym('Optimal_fiber_length', nMuscles)
    phi0 = SX.sym('Pennation_angle_at_muscle_optimal_fiber_length', nMuscles)
    maximalIsometricForce = SX.sym('Maximal_isometric_muscle_force', nMuscles)
    tendonSlackLength = SX.sym('Tendon_slack_length', nMuscles)

    muscleTendonParameters = vertcat(optimalFiberLength, phi0, maximalIsometricForce, tendonSlackLength)

    # ========= Muscle tendon states input of f(a,l,l.)
    fiberLength = SX.sym('Fiber_length', nMuscles)
    tendonLengthening = SX.sym('Tendon_Lengthening', nMuscles)
    a = SX.sym('Muscle_Activation', nMuscles)

    rootedvariables = vertcat(fiberLength, tendonLengthening)

    # ========= 2.2 Activation Dynamics (Not use in our modelisation)   ========= #
    # muscle activation is described by two nonlinear, first ODE
    # parameter value
    # taua = 0.015
    # taud = 0.060
    # b = 0.1
    # taua = 0.015; %SX.sym('Activation_time_constant',nMuscles);
    # taud = 0.060; %SX.sym('Desactivation_time_constant',nMuscles);
    # b = 0.1; %SX.sym('Transition_smoothness', nMuscles);
    #
    # fa = 0.5 .* tanh(b .* (e-a)) ;
    # da_dt = ((1 ./ taua .* (0.5 + 1.5 .* a)) .* (fa + 0.5) + ...
    #     ((0.5 + 1.5 .* a) ./ (taud)) .* (-fa + 0.5)) .* (e - a);
    # ode = struct('x',a, 'u', e, 'ode', da_dt);
    # tgrid = linspace(0, 1, 100);
    # activation_dynamics = integrator('activation_dynamics', 'cvodes', ode, tgrid(1), tgrid(1:end));
    # # example:
    # next_a = activation_dynamics('x0',[0;0;0],'u',[0.5; 0.7; 1]);

    # ========= 2.3 Muscle-Tendon Architecture Equations   ========= #
    tendonLength = tendonSlackLength + tendonLengthening #umtLength - muscleLength;

    normalizedTendonLength = tendonLength/ tendonSlackLength
    normalizedFiberLength = fiberLength/ optimalFiberLength

    # ========= 2.4  Muscle-Tendon Forces Equations   ========= #
        # Tendon force-length (S1)
    kT = 35; c1 = 0.200; c2 = 0.995; c3 = 0.250 # tendon parameters

    normalizedTendonForcePart1 =  0
    normalizedTendonForcePart2 = c1 * exp(kT * (normalizedTendonLength - c2)) - c3 # Normalized equation

    normalizedTendonForce = if_else(normalizedTendonLength < 1,
        normalizedTendonForcePart1,
        normalizedTendonForcePart2) # if normalized length under 0 the force = 0   # Normalized equation

    tendonForce= normalizedTendonForce * maximalIsometricForce # Non-normalized equation

        # === Active Force-Length (S2) ===
    # First Gaussian coefficients
    b11, b21, b31, b41 = 0.814483478343008, 1.055033428970575, 0.162384573599574, 0.063303448465465
    # Second Gaussian coefficients
    b12, b22, b32, b42 = 0.433004984392647, 0.716775413397760, -0.029947116970696, 0.200356847296188
    # Third Gaussian coefficients
    b13, b23, b33, b43 = 0.100, 1.000, 0.5 * np.sqrt(0.5), 0.000

    # Assume these are defined: normalizedFiberLength, a, maximalIsometricForce
    # normalizedFiberLength should be a float or a NumPy array

    # Gaussian 1
    num1 = normalizedFiberLength - b21
    den1 = b31 + b41 * normalizedFiberLength
    FMtilde1 = b11 * np.exp(-0.5 * (num1 ** 2) / (den1 ** 2))

    # Gaussian 2
    num2 = normalizedFiberLength - b22
    den2 = b32 + b42 * normalizedFiberLength
    FMtilde2 = b12 * np.exp(-0.5 * (num2 ** 2) / (den2 ** 2))

    # Gaussian 3
    num3 = normalizedFiberLength - b23
    den3 = b33 + b43 * normalizedFiberLength
    FMtilde3 = b13 * np.exp(-0.5 * (num3 ** 2) / (den3 ** 2))

    # Total normalized active force-length
    normalizedMuscleActiveForceLength = FMtilde1 + FMtilde2 + FMtilde3

    # Non-normalized active force
    MuscleActiveForceLength = a * normalizedMuscleActiveForceLength * maximalIsometricForce

        # === Passive Force-Length (S3) ===
    kpe = 4.0
    e0 = 0.6

    normalizedMusclePassiveForce = if_else(
        normalizedFiberLength < 1,
        0.0,
        (np.exp((kpe * (normalizedFiberLength - 1)) / e0) - 1) / (np.exp(kpe) - 1)
    )
    normalizedMusclePassiveForcePart1 = 0
    normalizedMusclePassiveForcePart2 = (exp(((kpe * (normalizedFiberLength - 1)) / e0)) - 1) / (exp(kpe) - 1)
    normalizedMusclePassiveForce = if_else(normalizedFiberLength < 1,
    normalizedMusclePassiveForcePart1,
    normalizedMusclePassiveForcePart2) # if normalized length under 0 the force = 0 %Normalized equation
    musclePassiveForce = normalizedMusclePassiveForce * maximalIsometricForce #Non - normalized equation

    musclePassiveForce = normalizedMusclePassiveForce * maximalIsometricForce

        # === Force-Velocity (S4) ===
    normalizedMuscleForceVelocity = 1.0  # Assuming velocity = 0

        # === Total Force ===
    normalizedMuscleForce = (
            a * normalizedMuscleActiveForceLength * normalizedMuscleForceVelocity
            + normalizedMusclePassiveForce
    )
    muscleForce = normalizedMuscleForce * maximalIsometricForce

    # ========= 2.5  Casadi functions about model Muscle-Tendon Forces   ========= #
    neuromusculoskeletal_state = vertcat(a, q, known_parameters)
    all_states = vertcat(neuromusculoskeletal_state, rootedvariables)

    getTendonForce = Function(
        'getTendonForce',
        [all_states, muscleTendonParameters],
        [tendonForce],
        ['all_states', 'muscle_tendon_parameters'],
        ['tendon_force']
    )

    getMuscleForce = Function(
        'getMuscleForce',
        [all_states, muscleTendonParameters],
        [muscleForce],
        ['all_states', 'muscle_tendon_parameters'],
        ['muscle_force']
    )

    normalizeTendonForce = Function(
        'normalizeTendonForce',
        [all_states, muscleTendonParameters],
        [normalizedTendonForce],
        ['all_states', 'muscleTendonParameters'],
        ['NormalizedTendonForce']
    )

    getMusclePassiveForce = Function(
        'getMusclePassiveForce',
        [all_states, muscleTendonParameters],
        [musclePassiveForce],
        ['all_states', 'muscle_tendon_parameters'],
        ['MusclePassiveForce']
    )

    getMuscleActiveForce = Function(
        'getMuscleActiveForce',
        [all_states, muscleTendonParameters],
        [MuscleActiveForceLength],
        ['all_states', 'muscle_tendon_parameters'],
        ['MuscleActiveForce']
    )

    representationMusclePassiveForce = Function(
        'representationMusclePassiveForce',
        [fiberLength[0], optimalFiberLength[0], maximalIsometricForce[0]],
        [musclePassiveForce[0]],
        ['fiberLength', 'optimalFiberLength', 'maximalIsometricForce'],
        ['MusclePassiveForce']
    )

    representationMuscleActiveForceLength = Function(
        'representationMuscleActiveForceLength',
        [a[0], fiberLength[0], optimalFiberLength[0], maximalIsometricForce[0]],
        [MuscleActiveForceLength[0]],
        ['a', 'fiberLength', 'optimalFiberLength', 'maximalIsometricForce'],
        ['MuscleActiveForceLength']
    )

    representationTendonForce = Function(
        'representationTendonForce',
        [tendonSlackLength[0], tendonLengthening[0], maximalIsometricForce[0]],
        [tendonForce[0]],
        ['tendonSlackLength', 'tendonLengthening', 'maximalIsometricForce'],
        ['tendonForce']
    )

    normalizeTendonLength = Function(
        'normalizeTendonLength',
        [all_states, muscleTendonParameters],
        [normalizedTendonLength],
        ['all_states', 'muscle_tendon_parameters'],
        ['length_tendon']
    )

    normalizeFiberLength = Function(
        'normalizeFiberLength',
        [all_states, muscleTendonParameters],
        [normalizedFiberLength],
        ['all_states', 'muscle_tendon_parameters'],
        ['normalizeFiberLength']
    )

                        # ========= 3. equilibrium functions   ========= #
    # ========= 3.1 inputed and rooted variables of equilibrium functions   ========= #
        # === input ==
    LUMT = SX.sym('UMT_length', nMuscles)

        # === unknown ==
    FT = SX.sym('Tendon_force', nMuscles)
    FM = SX.sym('Muscle_force', nMuscles)
    pennationAngle = SX.sym('Pennation_angle', nMuscles)

    # ========= 3.2 constraint functions   ========= #
    g3 = FT - tendonForce
    g4 = FM - muscleForce
    g5 = LUMT - (np.cos(pennationAngle) * fiberLength + tendonLength)
    g6 = (optimalFiberLength * np.sin(phi0)) - (fiberLength * np.sin(pennationAngle))
    g7 = FM * np.cos(pennationAngle) - FT

    # ========= 3.2 Muscle - tendon equilibrium   ========= #
        #  === all muscle  ===
    unknown = vertcat(FT, FM, tendonLengthening, fiberLength, pennationAngle)
    known = vertcat(a, LUMT, muscleTendonParameters)

    # === Options dictionary ===
    scale_const = 1/DM([1000, 1000, 1000,
            1000, 1000, 1000,
            0.1, 0.1, 0.1,
            0.01, 0.01, 0.01,
            0.1, 0.1, 0.1])

    opts_kinsol = {
        "constraints": tuple(1 for _ in range(len(scale_const.elements()))),  # [1,1,1,1,1], #SX.ones(5, 1),
        "abstol": 1e-15,
        "u_scale": scale_const.elements(),
        "error_on_fail": False,
        "max_iter": 1000,
        "iterative_solver": 'bcgstab',
        "print_level": 3,
        "disable_internal_warnings": True
    }

    # Define the residual function
    equilibriumError = Function(
        'equilibriumError',
        [unknown, known],
        [vertcat(g3, g4, g5, g6, g7)],
        ['x', 'p'],
        ['residuals']
    )

    # Define the rootfinder
    equilibrateMuscleTendon = rootfinder(
        'equilibrateMuscleTendon',  # name
        'kinsol',  # solver type
        equilibriumError,  # function
        opts_kinsol  # options dictionary
    )
    # ========= 3.3.1 single muscle   ========= #
    unknown = []
    known = []
    unknown = vertcat(FT[0], FM[0], tendonLengthening[0], fiberLength[0], pennationAngle[0])
    known = vertcat(a[0], LUMT[0], optimalFiberLength[0], phi0[0], maximalIsometricForce[0], tendonSlackLength[0])

    equilibriumErrorSingleMuscle = Function(
        'equilibriumErrorSingleMuscle',
        [unknown, known],
        [vertcat(g3[0], g4[0], g5[0], g6[0], g7[0])],
        ['x', 'p'],
        ['residuals'],
    )

    # solver_alorithm = 'kinsol'
    solver_alorithm = 'newton'

    if solver_alorithm == 'kinsol':
        # kinsol
        scale_const = 1/DM([500, 500, 0.1, 0.001, 0.1])
        opts_kinsolSM = {
            "constraints" : tuple(1 for _ in range(len(scale_const.elements())))  , #[1,1,1,1,1], #SX.ones(5, 1),
            "abstol" : 1e-8,
            "u_scale" : scale_const.elements(),
            "fd_method" :'central',
            "error_on_fail" : False,
            "max_iter" : 1000,
            "iterative_solver" : 'bcgstab',
            "print_level" : 3,
            "disable_internal_warnings" : True
        }

        equilibrateMuscleTendonSingleMuscle = rootfinder(
              'equilibrateMuscleTendonSingleMuscle',
              'kinsol',
              equilibriumErrorSingleMuscle,
              opts_kinsolSM
        )

    elif solver_alorithm == 'newton':
        # newton
        opts_newtonSM = {
            "abstol" : 1e-8,
            "max_iter" : 1000,
            "print_iteration" : True,
        }

        equilibrateMuscleTendonSingleMuscle = rootfinder(
            'equilibrateMuscleTendonSingleMuscle',
            'newton',
            equilibriumErrorSingleMuscle,
            opts_newtonSM
        )



                                # ========= 4. Computing Joint Moments and Angles    ========= #
    jointMoment = momentArm * tendonForce
    jointMoment = sum1(jointMoment[:,-1:])

    getJointMoment = Function(
        'getJointMoment',
        [all_states, muscleTendonParameters],
        [jointMoment],
        ['all_states', 'muscle_tendon_parameters'],
        ['jointMoment']
    )

    jointMoment = momentArm * FT
    jointMoment = sum1(jointMoment[:,-1:])

    getJointMoment2 = Function(
        'getJointMoment2',
        [musculoskeletal_states, FT],
        [jointMoment],
        ['musculoskeletal_states', 'tendon_force'],
        ['jointMoment']
    )

                            # ========= 5. dictionnary of casadi function     ========= #

    casadi_function = {
        "forwardKinematics": forwardKinematics ,
        "getMTULength": getMTULength,
        "getMomentArm": getMomentArm,
        "getTendonForce": getTendonForce,
        "getMuscleForce": getMuscleForce,
        "getMusclePassiveForce": getMusclePassiveForce,
        "getMuscleActiveForce": getMuscleActiveForce,
        "normalizeTendonForce": normalizeTendonForce,
        "normalizeTendonLength": normalizeTendonLength,
        "normalizeFiberLength": normalizeFiberLength,
        "representationMusclePassiveForce": representationMusclePassiveForce,
        "representationMuscleActiveForceLength": representationMuscleActiveForceLength,
        "representationTendonForce": representationTendonForce,
        "equilibriumError": equilibriumError,
        "equilibrateMuscleTendon": equilibrateMuscleTendon,
        "equilibriumErrorSingleMuscle": equilibriumErrorSingleMuscle,
        "equilibrateMuscleTendonSingleMuscle": equilibrateMuscleTendonSingleMuscle,
        "getJointMoment": getJointMoment,
        "getJointMoment2": getJointMoment2
    }

    definition = ["a : Neuromuscular activation (between 0 and 1)",
                  "q : Spatial skeleton configuration (x, y, z, theta_hip, theta_knee, theta_ankle)",
                  "known_parameters : musculoskeletal parameters (segment_geometry, muscle_insertion, Local_ViaPoint_tibialis_anterior)",
                  "neuromusculoskeletal_state : neuromusculoskeletal parameters(a, q, known_parameters)",
                  "muscleTendonParameters : Muscle Tendon Parameters (ℓom, φo, Fom, ℓst) ",
                  "p : neuromusculoskeletal parameters + Muscle Tendon Parameters (neuromusculoskeletal_state, muscleTendonParameters) "]


    return casadi_function,muscleTendonParameters,definition

def test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num):
    # === plot musculosquelton system ===
    Origin_num, Insertion_num, ViaPoint, Markers_num = casadi_function['forwardKinematics'](np.concatenate((q_num, skeleton_num)).reshape(-1, 1))

    # === get muscle length ===
    mtu_length = casadi_function['getMTULength'](np.concatenate((q_num, skeleton_num)).reshape(-1, 1))

    tibialis_length = float(mtu_length[0])
    soleus_length = float(mtu_length[1])
    gastrocnemius_length = float(mtu_length[2])

    print('\n','================== mtu length length ==================')
    print('tibialis: ', tibialis_length, 'm')
    print('soleus: ', soleus_length, 'm')
    print('gastrocnemius: ', gastrocnemius_length, 'm')

    musculoskeletal_states_num = np.concatenate((q_num, skeleton_num)).reshape(-1, 1)


    # === get muscle moment arm ===
    mtu_moment_arm = casadi_function['getMomentArm'](musculoskeletal_states_num)

    tibialis_moment_arm = float(mtu_moment_arm[0,-1])
    soleus_moment_arm = float(mtu_moment_arm[1,-1])
    gastrocnemius_moment_arm = float(mtu_moment_arm[2,-1])

    print('\n','================== mtu moment arm ==================')
    print('tibialis: ', tibialis_moment_arm, 'm-1')
    print('soleus: ', soleus_moment_arm, 'm-1')
    print('gastrocnemius: ', gastrocnemius_moment_arm, 'm-1')

    # === muscle equilibrium ===
    muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num)
    x_opt_tibialis,state_tibialis = root_muscle_dynamics(a_num[0], float(mtu_length[0]),muscle_tendon_parameters_num[[0, 3, 6, 9]], 'tibialis', casadi_function)
    x_opt_soleus,state_soleus = root_muscle_dynamics(a_num[1], float(mtu_length[1]), muscle_tendon_parameters_num[[1, 4, 7, 10]], 'soleus', casadi_function)
    x_opt_gastrocnemius,state_gastrocnemius = root_muscle_dynamics(a_num[2], float(mtu_length[2]),muscle_tendon_parameters_num[[3, 5, 8, 11]], 'gastrocnemius',casadi_function)

    FT = np.array([float(x_opt_tibialis[0, 0]), float(x_opt_soleus[0, 0]), float(x_opt_gastrocnemius[0, 0])])

    ankle_moment = casadi_function['getJointMoment2'](musculoskeletal_states_num, FT)
    print('\n================== simulatied ==================')
    print(f'tibialis force: {float(FT[0]):.4f} N')
    print(f'soleus force: {float(FT[1]):.4f} N')
    print(f'gastrocnemius force: {float(FT[2]):.4f} N')
    print(f'ankle moment: {float(ankle_moment[0]):.4f} N.m')

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
        ax.plot([markers[0, i], markers[0, i+1]],
                [markers[1, i], markers[1, i+1]],
                [markers[2, i], markers[2, i+1]], color='black')

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
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

       # Initial plot
    def update_plot(val=None):
        q = np.array([slider_q[i].val for i in range(6)])
        q[3:] = np.deg2rad(q[3:])  # Convert q4, q5, q6 from degrees to radians
        musculoskeletal_states_num = np.concatenate((q, skeleton_num)).reshape(-1, 1)
        a = np.array([muscle_sliders[i].val for i in range(3)])

        origins, insertions, via_point, markers = casadi_function['forwardKinematics'](musculoskeletal_states_num)
        plotmodel(ax, origins, insertions, markers, via_point)

        mtu_length = casadi_function['getMTULength'](musculoskeletal_states_num)
        mtu_moment_arm = casadi_function['getMomentArm'](musculoskeletal_states_num)

        print('\n================== MTU Length ==================')
        print(f'tibialis: {float(mtu_length[0]):.4f} m')
        print(f'soleus: {float(mtu_length[1]):.4f} m')
        print(f'gastrocnemius: {float(mtu_length[2]):.4f} m')

        print('\n================== Moment Arms ==================')
        print(f'tibialis: {float(mtu_moment_arm[0,-1]):.4f} m^-1')
        print(f'soleus: {float(mtu_moment_arm[1,-1]):.4f} m^-1')
        print(f'gastrocnemius: {float(mtu_moment_arm[2,-1]):.4f} m^-1')

        x_opt_tibialis,state_tibialis = root_muscle_dynamics(a[0],float(mtu_length[0]),np.array(muscle_tendon_parameters_num)[[0, 3, 6, 9]],'tibialis',casadi_function)
        x_opt_soleus,state_soleus = root_muscle_dynamics(a[1],float(mtu_length[1]),np.array(muscle_tendon_parameters_num)[[1, 4, 7, 10]],'soleus',casadi_function)
        x_opt_gastrocnemius,state_gastrocnemius = root_muscle_dynamics(a[2],float(mtu_length[2]),np.array(muscle_tendon_parameters_num)[[3, 5, 8, 11]],'gastrocnemius',casadi_function)
        FT = np.array([float(x_opt_tibialis[0,0]),float(x_opt_soleus[0,0]),float(x_opt_gastrocnemius[0,0])])

        ankle_moment = casadi_function['getJointMoment2'](musculoskeletal_states_num,FT)
        print('\n================== simulatied ==================')
        print(f'tibialis force: {float(FT[0]):.4f} N')
        print(f'soleus force: {float(FT[1]):.4f} N')
        print(f'gastrocnemius force: {float(FT[2]):.4f} N')
        print(f'ankle moment: {float(ankle_moment[0]):.4f} N.m')

        fig.canvas.draw_idle()

   # Add sliders
    slider_limits = [
        (-1, 1),  # q1
        (-1, 1),  # q2
        (-1, 1),  # q3
        (-30, 30),  # q4
        (0, 90),  # q5
        (-30, 30)  # q6
    ]
    axcolor = 'lightgoldenrodyellow'
    slider_q = []
    slider_axes = []

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

def root_muscle_dynamics(a,lmtu,parameters,muscle_name,casadi_function):
    """
    Solves for the muscle-tendon equilibrium state using root-finding with CasADi.

    Parameters:
    - a (float): Muscle activation level.
    - lmtu (float): Muscle-tendon unit length.
    - parameters (np.ndarray): Array of muscle-specific parameters required by the model.
    - muscle_name (str): Name of the muscle (used for print/logging purposes).
    - casadi_function (dict): Dictionary containing CasADi functions:
        - 'equilibrateMuscleTendonSingleMuscle': Root-finding function to solve muscle-tendon equilibrium.
        - 'equilibriumErrorSingleMuscle': Function to evaluate equilibrium error (residuals).

    Returns:
    - x_opt (np.ndarray): Optimized muscle-tendon state vector satisfying equilibrium conditions.
    - equilibrium_status (str): Status indicating whether equilibrium was successfully achieved ('success' or 'fail').

    Description:
    This function attempts to solve the equilibrium condition of a single muscle-tendon unit
    using a root-finding algorithm. It allows up to `max_attempts` (default: 50) to find a
    solution where all residuals fall below a defined tolerance (`lim_residuals`). If successful,
    it returns the solution and the status. The process and results are printed for each attempt.
    """
    lim_residuals = 1e-10
    equilibrium_status = 'fail'
    max_attempts = 50
    attempt = 0

    while equilibrium_status == 'fail' and attempt < max_attempts:
        attempt += 1
        x_start = Xstart_equi(a, lmtu, parameters)

        x_opt = casadi_function['equilibrateMuscleTendonSingleMuscle'](
            x_start,
            np.array([a, lmtu] + parameters.tolist())
        )

        residuals_print = casadi_function['equilibriumErrorSingleMuscle'](
            x_opt,
            np.array([a, lmtu] + parameters.tolist())
        )
        residuals = np.array(residuals_print).flatten()

        if np.any(residuals > lim_residuals):
            print(f"[Attempt {attempt}] Root-finding failed: residual too high ({residuals})")
        else:
            print(f"[Attempt {attempt}] Success: residual within tolerance")
            equilibrium_status = 'success'

    print('\n','=============== dynamics: ', muscle_name,' ===============')
    print('residuals ', residuals_print)
    print('x_opt ', x_opt)
    print('equilibrium status: ', equilibrium_status)
    print('=============== ========================= ===============','\n','\n','\n',)
    return x_opt, equilibrium_status

def Xstart_equi(a, lmtu, parameters):
    """ optimal start of the fonction equilibrium: it permit to start in respect for force equilibrium and geometry
   equilibrium.

   equilibrium function:
    g3 = FT - tendonForce
    g4 = FM - muscleForce
    g5 = LUMT - (np.cos(pennationAngle) * fiberLength + tendonLength)
    g6 = (optimalFiberLength * np.sin(phi0)) - (fiberLength * np.sin(pennationAngle))
    g7 = FM * np.cos(pennationAngle) - FT

    Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)

           Returns:
        X_start (np.ndarray): Shape (5,) including tendon force, muscle force, tendon lengthening, fiber length and
        pennation angle
   """
    muscle_force = a * parameters[2]
    tendon_lengthening = a * 0.05
    tendon_length = parameters[3] + (parameters[3] * tendon_lengthening)
    pennation_angle = parameters[1] + parameters[1] * np.random.uniform(-0.2, 0.2)
    fiber_length = (lmtu - tendon_length) / np.cos(pennation_angle)
    tendon_force = a * parameters[2] * np.cos(pennation_angle)

    x_start = [tendon_force, muscle_force, tendon_lengthening, fiber_length, pennation_angle]
    x_start = np.round(x_start, 4)

    return x_start

def hypotetical_data_generator(skeleton_num, muscle_tendon_parameters_num, casadi_function):
    """ generate hypotetical data to make the NLP
   .
    rooted = vertcat(tendon_force,
    tendon_force,
    tendonLengthening,
    fiberLength,
    pennationAngle)


           Returns:
        hypotetical_data (np.ndarray): Shape (12,ntrials) including
        [ankle_torque, q_knee, q_ankle
        a_tibialis, a_soleus, a_gastrocnemius,
        fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
        pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius]
   """

    header = [
        'ankle_torque', 'q_knee', 'q_ankle',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius',
        'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
        'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius'
        ]

    ntrialsFail = 0 # compteur of trials non - optimized(p > 10e-5)
    ntrialsSucceds = 0 # compteur of trials optimized(p > 10e-5)

    # ========= 1 Muscle-Tendon Architecture Equations   ========= #
    # ========= 1.1 Musculo skeletical configuration during trial(input)   ========= #
    qknee = np.linspace(0, 90, 5)  # example knee angles
    qankle = np.linspace(-20, 30, 5)  # example ankle angles

    # ========= 1.2 Neuronal activation(input)   ========= #
    # Muscle activation [Tibialis Anterior, Soleus, Gastrocnemius]
    a_num = np.array([
        [0, 0, 0],
        [0, 0.2, 0.2],
        [0, 0.4, 0.4],
        [0, 0.6, 0.6],
        [0.2, 0, 0],
        [0.4, 0, 0],
        [0.6, 0, 0]
    ])

    # ========= 1.3 Muscle Tendon Parameters(ℓom, φo, Fom, ℓst)(input)   ========= #
    muscle_tendon_parameters_ta = np.array(muscle_tendon_parameters_num)[[0, 3, 6, 9]]
    muscle_tendon_parameters_sol = np.array(muscle_tendon_parameters_num)[[1, 4, 7, 10]]
    muscle_tendon_parameters_gast =np.array(muscle_tendon_parameters_num)[[2, 5, 8, 11]]

    # ========= 2.1 data generator   ========= #
    ntrials = len(a_num) * len(qknee) * len(qankle)
    hypotetical_data = np.zeros((ntrials, 12))

    trials = 0

    for i in range(len(a_num)):  # for each muscle activation
        for ii in range(len(qknee)):  # for each knee angle
            for iii in range(len(qankle)):  # for each ankle angle
                # 2.1 Progression
                trials += 1
                percent = round(trials / ntrials, 2)
                print(f"Progressing: {percent * 100:.0f} %")

                # 2.2 Neuromusculoskeletal configuration
                q_num = [0, 0, 0, 0, qknee[ii], qankle[iii]]
                musculoskeletal_states_num = q_num + list(skeleton_num)

                # 2.3 UMT length
                mtu_length = casadi_function['getMTULength'](musculoskeletal_states_num)
                mtu_length_tibialis = float(mtu_length[0])
                mtu_length_soleus = float(mtu_length[1])
                mtu_length_gastrocnemius = float(mtu_length[2])

                # 2.4 muscle activities
                a_tibialis = a_num[i,0]
                a_soleus = a_num[i,1]
                a_gastrocnemius = a_num[i,2]

                # 2.5 Solver
                x_opt_tibialis,state_tibialis = root_muscle_dynamics(a_tibialis,
                                                      mtu_length_tibialis,
                                                      muscle_tendon_parameters_ta,
                                                      'tibialis',
                                                      casadi_function)
                x_opt_soleus,state_soleus = root_muscle_dynamics(a_soleus,
                                                      mtu_length_soleus,
                                                      muscle_tendon_parameters_sol,
                                                      'soleus',
                                                      casadi_function)

                x_opt_gastrocnemius,state_gastrocnemius = root_muscle_dynamics(a_gastrocnemius,
                                                      mtu_length_gastrocnemius,
                                                      muscle_tendon_parameters_gast,
                                                      'gastrocnemius',
                                                      casadi_function)

                if state_tibialis == 'fail' or state_soleus == 'fail'or state_tibialis == 'fail' or state_gastrocnemius == 'fail':
                    ntrialsFail += 1
                else:
                    ntrialsSucceds += 1
                    # 2.6 ankle torque
                    FT = np.array(
                        [float(x_opt_tibialis[0, 0]),
                         float(x_opt_soleus[0, 0]),
                         float(x_opt_gastrocnemius[0, 0])])

                    ankle_torque = casadi_function['getJointMoment2'](musculoskeletal_states_num, FT)
                    print('\n================== simulatied ==================')
                    print(f'tibialis force: {float(FT[0]):.4f} N')
                    print(f'soleus force: {float(FT[1]):.4f} N')
                    print(f'gastrocnemius force: {float(FT[2]):.4f} N')
                    print(f'ankle moment: {float(ankle_torque[0]):.4f} N.m')

                    # 2.7 extract variable
                    # rooted = vertcat(tendon_force,tendon_force,tendonLengthening,fiberLength,pennationAngle)
                    tibialis_fiber_length = float(x_opt_tibialis[3])
                    tibialis_pennation_angle = x_opt_tibialis[4]
                    tibialis_pennation_angle = float(np.rad2deg(tibialis_pennation_angle))

                    soleus_fiber_length = float(x_opt_soleus[3])
                    soleus_pennation_angle = x_opt_soleus[4]
                    soleus_pennation_angle = float(np.rad2deg(soleus_pennation_angle))

                    gastrocnemius_fiber_length = float(x_opt_gastrocnemius[3])
                    gastrocnemius_pennation_angle = x_opt_gastrocnemius[4]
                    gastrocnemius_pennation_angle = float(np.rad2deg(gastrocnemius_pennation_angle))

                    hypotetical_data[ntrialsSucceds - 1] = [
                        float(ankle_torque),
                        qknee[ii], qankle[iii],
                        float(a_num[i, 0]), float(a_num[i, 1]), float(a_num[i, 2]),
                        tibialis_fiber_length, soleus_fiber_length, gastrocnemius_fiber_length,
                        tibialis_pennation_angle, soleus_pennation_angle, gastrocnemius_pennation_angle
                    ]

    print('\n================== hypotetical data generator state ==================')
    print(f"fail trials: {(ntrialsFail/ntrials)*100:.2f} %")

    return header, hypotetical_data
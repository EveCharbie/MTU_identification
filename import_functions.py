import os
import numpy as np
import xml.etree.ElementTree as ET
from scipy.optimize import brentq

# Coefficients par défaut de De Groote 2016 (eq. S1)
DEGROOTE_C1 = 0.200
DEGROOTE_C2 = 0.995
DEGROOTE_C3 = 0.200


def kT_degroote_to_kt_physical(kT, f0m, lst, eps0_t=0.04,
                               c1=DEGROOTE_C1, c2=DEGROOTE_C2):
    """
    Convert De Groote's normalized tendon shape factor kT (dimensionless)
    to a physical tendon stiffness k^T [N/m] at F0m.

    The physical stiffness is the slope of the dimensional force-length
    curve evaluated at fiber length = lst * (1 + eps0_t).

    Parameters
    ----------
    kT : float
        De Groote's tendon shape factor (~35 by default).
    f0m : float
        Maximum isometric force [N].
    lst : float
        Tendon slack length [m].
    eps0_t : float, optional
        Tendon strain at F0m. Default 0.04 (De Groote default).
    c1, c2 : float, optional
        Calibration coefficients of the S1 exponential.

    Returns
    -------
    float
        Physical tendon stiffness k^T [N/m].
    """
    normalized_slope = c1 * kT * np.exp(kT * (1 + eps0_t - c2))
    return (f0m / lst) * normalized_slope


def _kt_physical_to_kT_degroote(k_phys, f0m, lst, eps0_t,
                                c1=DEGROOTE_C1, c2=DEGROOTE_C2,
                                kT_bracket=(1.0, 200.0)):
    """Invert De Groote's S1 slope equation to find the shape factor kT."""
    target_normalized_slope = k_phys * lst / f0m

    def residual(kT):
        return c1 * kT * np.exp(kT * (1 + eps0_t - c2)) - target_normalized_slope

    return brentq(residual, *kT_bracket)


def _compute_kt_degroote(p):
    """Compute De Groote's tendon shape factor from scraped MTU parameters."""
    required = ('f0m', 'eps0_t', 'lst')
    if not all(p.get(k) not in (None, 0) for k in required):
        return None

    k_phys = p['f0m'] / (p['eps0_t'] * p['lst'])
    try:
        return _kt_physical_to_kT_degroote(k_phys, p['f0m'], p['lst'], p['eps0_t'])
    except ValueError:
        return None

MTU_PARAM_TAGS = {
    # Paramètres de Hill (géométrie MTU)
    'l0m':        'optimal_fiber_length',
    'phi0':       'pennation_angle_at_optimal',
    'f0m':        'max_isometric_force',
    'lst':        'tendon_slack_length',
    # Tendon
    'eps0_t':     'FmaxTendonStrain',
    # Muscle passif
    'eps0_m':     'FmaxMuscleStrain',
    # Forme des courbes F-L
    'kshape_active': 'KshapeActive',
    'km': 'KshapePassive',
    # Courbe F-v
    'af':         'Af',
    'flen':       'Flen',
    # Dynamique d'activation
    'tau_act':    'activation_time_constant',
    'tau_deact':  'deactivation_time_constant',
}

# Paramètres dérivés (calculés à partir des paramètres scrapés)
MTU_DERIVED_PARAMS = {
    'k_phys': (
        lambda p: p['f0m'] / (p['eps0_t'] * p['lst'])
                  if all(p.get(k) not in (None, 0) for k in ('f0m', 'eps0_t', 'lst'))
                  else None,
        ['f0m', 'eps0_t', 'lst'],
    ),
    'kt': (_compute_kt_degroote, ['f0m', 'eps0_t', 'lst']),
}

# Règles d'agrégation pour la fusion lat_gas_r + med_gas_r -> gast_r
MTU_AGGREGATION_RULES = {
    'f0m': 'sum',
    'kt': 'sum',
}

def get_model_osim_generic(oism_path, oism_file, mtu_params=None):
    """
    Extract muscle-tendon and joint parameters from an OpenSim .osim file
    (Thelen2003Muscle).

    Parameters
    ----------
    oism_path : str
        Directory containing the .osim file.
    oism_file : str
        Name of the .osim file.
    mtu_params : list of str, optional
        Ordered list of muscle-tendon parameters to extract.
        Valid keys : scrapés (MTU_PARAM_TAGS) + dérivés (MTU_DERIVED_PARAMS).
        If None, all parameters are extracted in declaration order.

    Returns
    -------
    musculoskeletal_num : np.ndarray, shape (33,)
        Segment geometry (12) + muscle origins (9) + insertions (9) + via point (3).
    muscle_tendon_parameters_num : np.ndarray
        Flattened param-major vector, shape (n_params * 3,).
    """

    all_valid = list(MTU_PARAM_TAGS.keys()) + list(MTU_DERIVED_PARAMS.keys())

    if mtu_params is None:
        mtu_params = all_valid

    unknown = [p for p in mtu_params if p not in all_valid]
    if unknown:
        raise ValueError(
            f"Unknown MTU parameter(s): {unknown}. Valid keys: {all_valid}"
        )

    # Dépendances des paramètres dérivés
    scraped_needed = set(p for p in mtu_params if p in MTU_PARAM_TAGS)
    for p in mtu_params:
        if p in MTU_DERIVED_PARAMS:
            _, deps = MTU_DERIVED_PARAMS[p]
            scraped_needed.update(deps)

    interst_body     = ['tibia_r', 'talus_r', 'calcn_r', 'toes_r']
    interest_muscles = ['tib_ant_r', 'soleus_r', 'gast_r']

    n_mtu_out = len(mtu_params) * len(interest_muscles)
    n_msk_out = 33

    try:
        full_path = os.path.join(oism_path, oism_file)
        tree = ET.parse(full_path)
        root = tree.getroot()
        model = root.find('Model')
        if model is None:
            raise ValueError("No <Model> element found in the .osim file.")

        # --- Segment geometry --- #
        bodies = {}
        body_set = model.find('BodySet')
        if body_set is not None:
            body_object = body_set.find('objects')
            if body_object is not None:
                for body in body_object:
                    body_name = body.attrib.get('name')
                    joint = body.find('Joint')
                    if joint is not None:
                        custom_joint = joint.find('CustomJoint')
                        if custom_joint is not None:
                            loc_elem = custom_joint.find('location_in_parent')
                            if loc_elem is not None:
                                x, y, z = map(float, loc_elem.text.strip().split())
                                bodies[body_name] = {'X': x, 'Y': y, 'Z': z}

        # --- Muscles : géométrie + paramètres scrapés --- #
        muscle_pass_point = {}
        muscle_tendon_parameters = {}

        muscle_set = model.find('ForceSet')
        if muscle_set is not None:
            muscle_object = muscle_set.find('objects')
            for muscle in muscle_object:
                muscle_name = muscle.attrib.get('name')
                muscle_pass_point[muscle_name] = {}

                # Path points
                muscle_geometry = muscle.find('GeometryPath')
                if muscle_geometry is not None:
                    path_point_set = muscle_geometry.find('PathPointSet')
                    if path_point_set is not None:
                        path_point_objects = path_point_set.find('objects')
                        for pass_point in path_point_objects:
                            pass_point_name = pass_point.attrib.get('name')
                            loc_elem = pass_point.find('location')
                            if loc_elem is not None:
                                x, y, z = map(float, loc_elem.text.strip().split())
                                muscle_pass_point[muscle_name][pass_point_name] = {
                                    'X': x, 'Y': y, 'Z': z
                                }

                # Paramètres MTU
                def get_float_value(tag, element=muscle):
                    elem = element.find(tag)
                    return float(elem.text) if elem is not None and elem.text else None

                muscle_tendon_parameters[muscle_name] = {
                    param: get_float_value(MTU_PARAM_TAGS[param])
                    for param in scraped_needed
                }

        # --- Fusion gastrocnémiens --- #
        gast = 'gast_r'

        def _aggregate(param, lat_val, med_val):
            if lat_val is None or med_val is None:
                return None
            rule = MTU_AGGREGATION_RULES.get(param, 'mean')
            if rule == 'sum':
                return float(np.sum([lat_val, med_val]))
            return float(np.mean([lat_val, med_val]))

        muscle_tendon_parameters[gast] = {
            param: _aggregate(
                param,
                muscle_tendon_parameters['lat_gas_r'][param],
                muscle_tendon_parameters['med_gas_r'][param],
            )
            for param in scraped_needed
        }

        muscle_pass_point[gast] = {}
        for i in range(1, 4):
            lat_pt = muscle_pass_point['lat_gas_r'][f'lat_gas_r-P{i}']
            med_pt = muscle_pass_point['med_gas_r'][f'med_gas_r-P{i}']
            muscle_pass_point[gast][f'gast-P{i}'] = {
                axis: float(np.mean([lat_pt[axis], med_pt[axis]]))
                for axis in ('X', 'Y', 'Z')
            }

        # --- Paramètres dérivés (calculés par muscle, après fusion gast) --- #
        for muscle in list(muscle_tendon_parameters.keys()):
            for dparam in mtu_params:
                if dparam in MTU_DERIVED_PARAMS:
                    func, _ = MTU_DERIVED_PARAMS[dparam]
                    muscle_tendon_parameters[muscle][dparam] = func(
                        muscle_tendon_parameters[muscle]
                    )

        muscle_tendon_parameters_num = []
        for param in mtu_params:
            for muscle in interest_muscles:
                value = muscle_tendon_parameters[muscle][param]
                muscle_tendon_parameters_num.append(round(float(value), 3))

        # --- Geometry outputs --- #
        segment_geometry_num = []
        for body_name in interst_body:
            coord = bodies[body_name]
            segment_geometry_num.extend([coord['X'], coord['Y'], 0.0])

        muscle_origin = []
        for m in interest_muscles:
            first_name = next(iter(muscle_pass_point[m]))
            pt = muscle_pass_point[m][first_name]
            muscle_origin.extend([pt['X'], pt['Y'], 0.0])

        muscle_insertion = []
        for m in interest_muscles:
            last_name = next(reversed(muscle_pass_point[m]))
            pt = muscle_pass_point[m][last_name]
            muscle_insertion.extend([pt['X'], pt['Y'], 0.0])

        via = muscle_pass_point['tib_ant_r']['tib_ant_r-P2']
        muscle_via = [via['X'], via['Y'], 0.0]

        musculoskeletal_num = np.array(
            segment_geometry_num + muscle_origin + muscle_insertion + muscle_via
        )
        muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num)

        return musculoskeletal_num, muscle_tendon_parameters_num

    except Exception as e:
        print(f"Error extracting parameters: {e}")
        return np.zeros(n_msk_out), np.zeros(n_mtu_out)

def get_model_osim_scaled(oism_path, oism_file, mtu_params=None):
    """
    Extract muscle-tendon and joint parameters from an OpenSim 4.x .osim file.

    Output vector is MUSCLE-MAJOR : for each muscle (in interest_muscles order),
    parameters are listed in the order given by mtu_params.

    Parameters
    ----------
    oism_path : str
    oism_file : str
    mtu_params : list of str, optional
        If None, all valid parameters (scrapés + dérivés) are returned.

    Returns
    -------
    musculoskeletal_num : np.ndarray, shape (33,)
        Segment geometry (12) + origins (9) + insertions (9) + via point (3).
    muscle_tendon_parameters_num : np.ndarray, shape (n_params * 3,)
    """

    all_valid = list(MTU_PARAM_TAGS.keys()) + list(MTU_DERIVED_PARAMS.keys())
    if mtu_params is None:
        mtu_params = all_valid

    unknown = [p for p in mtu_params if p not in all_valid]
    if unknown:
        raise ValueError(
            f"Unknown MTU parameter(s): {unknown}. Valid keys: {all_valid}"
        )

    scraped_needed = set(p for p in mtu_params if p in MTU_PARAM_TAGS)
    for p in mtu_params:
        if p in MTU_DERIVED_PARAMS:
            _, deps = MTU_DERIVED_PARAMS[p]
            scraped_needed.update(deps)

    interst_body     = ['tibia_r', 'talus_r', 'calcn_r', 'toes_r']
    interest_muscles = ['tib_ant_r', 'soleus_r', 'gast_r']

    n_mtu_out = len(mtu_params) * len(interest_muscles)
    n_msk_out = 33

    try:
        full_path = os.path.join(oism_path, oism_file)
        tree = ET.parse(full_path)
        root = tree.getroot()
        model = root.find('Model')
        if model is None:
            raise ValueError("No <Model> element found in the .osim file.")

        # --- Segment geometry : via JointSet (OpenSim 4.x) --- #
        # Pour chaque body d'intérêt, on récupère la position du joint qui le
        # connecte à son parent, exprimée dans le repère du parent.
        bodies = _build_joint_index_by_child_body(model)

        # --- Muscles : géométrie + paramètres scrapés --- #
        muscle_pass_point = {}
        muscle_tendon_parameters = {}

        muscle_set = model.find('ForceSet')
        if muscle_set is not None:
            for muscle in muscle_set.find('objects'):
                muscle_name = muscle.attrib.get('name')
                if muscle_name is None:
                    continue
                muscle_pass_point[muscle_name] = {}

                # Path points
                muscle_geometry = muscle.find('GeometryPath')
                if muscle_geometry is not None:
                    path_point_set = muscle_geometry.find('PathPointSet')
                    if path_point_set is not None:
                        for pass_point in path_point_set.find('objects'):
                            pp_name = pass_point.attrib.get('name')
                            if pp_name is None:
                                continue
                            loc = _extract_path_point_location(pass_point)
                            if loc is not None:
                                muscle_pass_point[muscle_name][pp_name] = loc

                # Paramètres MTU scrapés
                def get_float_value(tag, element=muscle):
                    elem = element.find(tag)
                    return float(elem.text) if elem is not None and elem.text else None

                muscle_tendon_parameters[muscle_name] = {
                    param: get_float_value(MTU_PARAM_TAGS[param])
                    for param in scraped_needed
                }

        # --- Fusion gastrocnémiens --- #
        gast = 'gast_r'
        if 'lat_gas_r' in muscle_tendon_parameters and 'med_gas_r' in muscle_tendon_parameters:

            def _aggregate(param, lat_val, med_val):
                if lat_val is None or med_val is None:
                    return None
                rule = MTU_AGGREGATION_RULES.get(param, 'mean')
                if rule == 'sum':
                    return float(np.sum([lat_val, med_val]))
                return float(np.mean([lat_val, med_val]))

            muscle_tendon_parameters[gast] = {
                param: _aggregate(
                    param,
                    muscle_tendon_parameters['lat_gas_r'][param],
                    muscle_tendon_parameters['med_gas_r'][param],
                )
                for param in scraped_needed
            }

            muscle_pass_point[gast] = {}
            for i in range(1, 4):
                lat_key = f'lat_gas_r-P{i}'
                med_key = f'med_gas_r-P{i}'
                if (lat_key in muscle_pass_point['lat_gas_r']
                        and med_key in muscle_pass_point['med_gas_r']):
                    lat_pt = muscle_pass_point['lat_gas_r'][lat_key]
                    med_pt = muscle_pass_point['med_gas_r'][med_key]
                    muscle_pass_point[gast][f'gast-P{i}'] = {
                        axis: float(np.mean([lat_pt[axis], med_pt[axis]]))
                        for axis in ('X', 'Y', 'Z')
                    }

        # --- Paramètres dérivés (calculés après fusion gast) --- #
        for muscle in list(muscle_tendon_parameters.keys()):
            for dparam in mtu_params:
                if dparam in MTU_DERIVED_PARAMS:
                    func, _ = MTU_DERIVED_PARAMS[dparam]
                    muscle_tendon_parameters[muscle][dparam] = func(
                        muscle_tendon_parameters[muscle]
                    )

        # --- Output MTU vector (MUSCLE-MAJOR) --- #
        muscle_tendon_parameters_num = []
        for param in mtu_params:
            for muscle in interest_muscles:
                value = muscle_tendon_parameters[muscle][param]
                muscle_tendon_parameters_num.append(round(float(value), 3))

        # --- Segment geometry (12) --- #
        segment_geometry_num = []
        for body_name in interst_body:
            coord = bodies[body_name]
            segment_geometry_num.extend([coord['X'], coord['Y'], 0.0])

        # --- Muscle origins (9) --- #
        muscle_origin = []
        for m in interest_muscles:
            first_name = next(iter(muscle_pass_point[m]))
            pt = muscle_pass_point[m][first_name]
            muscle_origin.extend([pt['X'], pt['Y'], 0.0])

        # --- Muscle insertions (9) --- #
        muscle_insertion = []
        for m in interest_muscles:
            last_name = next(reversed(muscle_pass_point[m]))
            pt = muscle_pass_point[m][last_name]
            muscle_insertion.extend([pt['X'], pt['Y'], 0.0])

        # --- Via point tib_ant_r (3) --- #
        via = muscle_pass_point['tib_ant_r']['tib_ant_r-P2']
        muscle_via = [via['X'], via['Y'], 0.0]

        musculoskeletal_num = np.array(
            segment_geometry_num + muscle_origin + muscle_insertion + muscle_via
        )
        muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num)

        return musculoskeletal_num, muscle_tendon_parameters_num

    except Exception as e:
        print(f"Error extracting parameters: {e}")
        return np.zeros(n_msk_out), np.zeros(n_mtu_out)

def _build_joint_index_by_child_body(model):
    """
    Build a mapping {child_body_name: translation_in_parent} from the JointSet.

    In OpenSim 4.x, each joint connects a parent_frame to a child_frame, both
    are PhysicalOffsetFrames. The joint "location in parent" (like the old
    OpenSim 3.x <location_in_parent>) corresponds to the <translation> of the
    PhysicalOffsetFrame referenced by socket_parent_frame.

    For each joint, we identify the child body (via the PhysicalOffsetFrame
    referenced by socket_child_frame, whose socket_parent gives us the body).

    Returns
    -------
    dict {child_body_name: {'X': float, 'Y': float, 'Z': float}}
    """
    index = {}
    joint_set = model.find('JointSet')
    if joint_set is None:
        return index

    for joint in joint_set.find('objects'):
        parent_frame_name = joint.find('socket_parent_frame')
        child_frame_name = joint.find('socket_child_frame')
        frames_elem = joint.find('frames')
        if parent_frame_name is None or child_frame_name is None or frames_elem is None:
            continue

        # Les noms de socket sont bruts (ex: "tibia_r_offset"), sans slash
        parent_fname = parent_frame_name.text.strip().split('/')[-1]
        child_fname = child_frame_name.text.strip().split('/')[-1]

        # Localiser les PhysicalOffsetFrame correspondants
        pof_by_name = {
            pof.attrib.get('name'): pof
            for pof in frames_elem.findall('PhysicalOffsetFrame')
        }
        parent_pof = pof_by_name.get(parent_fname)
        child_pof = pof_by_name.get(child_fname)
        if parent_pof is None or child_pof is None:
            continue

        # Body enfant (ex: '/bodyset/talus_r' -> 'talus_r')
        child_body_elem = child_pof.find('socket_parent')
        if child_body_elem is None or not child_body_elem.text:
            continue
        child_body = child_body_elem.text.strip().split('/')[-1]

        # Translation du parent offset frame = position du joint dans le parent body
        trans_elem = parent_pof.find('translation')
        if trans_elem is None or not trans_elem.text:
            continue
        x, y, z = map(float, trans_elem.text.strip().split())
        index[child_body] = {'X': x, 'Y': y, 'Z': z}

    return index

def _extract_path_point_location(pass_point):
    """
    Extract (X, Y, Z) location from a path point, handling OpenSim 4.x types.

    - PathPoint / ConditionalPathPoint : <location> tag (fixed position).
    - MovingPathPoint : position depends on joint coordinate via a SimmSpline.
      We evaluate the spline at coordinate = 0 (neutral pose).
    """
    tag = pass_point.tag

    if tag in ('PathPoint', 'ConditionalPathPoint'):
        loc_elem = pass_point.find('location')
        if loc_elem is not None and loc_elem.text:
            x, y, z = map(float, loc_elem.text.strip().split())
            return {'X': x, 'Y': y, 'Z': z}
        return None

    if tag == 'MovingPathPoint':
        coords = {}
        for axis, xml_tag in (('X', 'x_location'), ('Y', 'y_location'), ('Z', 'z_location')):
            func_elem = pass_point.find(xml_tag)
            coords[axis] = _evaluate_spline_at_zero(func_elem) if func_elem is not None else None
        if all(v is not None for v in coords.values()):
            return coords
        return None

    return None

def _evaluate_spline_at_zero(func_elem):
    """
    Evaluate a SimmSpline (optionally wrapped in MultiplierFunction) at x=0.
    Linear interpolation between the two knots surrounding 0.
    """
    scale_elem = func_elem.find('.//scale')
    scale = float(scale_elem.text) if scale_elem is not None and scale_elem.text else 1.0

    spline = func_elem.find('.//SimmSpline')
    if spline is None:
        return None

    x_elem = spline.find('x')
    y_elem = spline.find('y')
    if x_elem is None or y_elem is None:
        return None

    xs = np.array(list(map(float, x_elem.text.strip().split())))
    ys = np.array(list(map(float, y_elem.text.strip().split())))

    if len(xs) == 1:
        return float(ys[0]) * scale

    return float(np.interp(0.0, xs, ys)) * scale

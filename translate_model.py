from pathlib import Path
import numpy as np

from biobuddy import (
    BiomechanicalModelReal,
    MuscleStateType,
    MuscleType,
)


# Paths
participant_name = "TEJ_03"
current_path_file = Path(__file__).parent
biomod_filepath = f"{current_path_file}/num_data/{participant_name}/{participant_name}.bioMod"
osim_filepath = f"{current_path_file}/num_data/{participant_name}/{participant_name}.osim"
# geometry_path = f"{current_path_file}/../external/opensim-models/Geometry"
# geometry_cleaned_path = f"{current_path_file}/models/Geometry_cleaned"

# # Convert the vtp files
# mesh = MeshParser(geometry_folder=geometry_path)
# mesh.process_meshes(fail_on_error=False)
# mesh.write(geometry_cleaned_path, format=MeshFormat.VTP)

# --- Reading an .osim model and translating it to a .bioMod model --- #
# Read an .osim file
model = BiomechanicalModelReal().from_osim(
    filepath=osim_filepath,
    muscle_type=MuscleType.HILL_DE_GROOTE,
    muscle_state_type=MuscleStateType.DEGROOTE,
    mesh_dir="Geometry_cleaned",
)

# Fix the via points before translating to biomod as there are some conditional and moving via points
model.fix_via_points(q=np.zeros((model.nb_q,)))

# Remove everything unused in the model
segments_to_remove = [
    "pelvis_parent_offset",
    "pelvis_translation",
    "pelvis_rotation_transform",
    "pelvis_reset_axis",
    "pelvis_geom_2",
    "pelvis_geom_3",
    "femur_r_translation",
    "femur_r_rotation_transform",
    "femur_r_reset_axis",
    "tibia_r_translation",
    "talus_r_translation",
    "calcn_r_translation",
    "toes_r_parent_offset",
    "toes_r_translation",
    "toes_r_rotation_transform",
    "toes_r_reset_axis",
    "femur_l_parent_offset",
    "femur_l_translation",
    "femur_l_rotation_transform",
    "femur_l_reset_axis",
    "femur_l",
    "tibia_l_parent_offset",
    "tibia_l_translation",
    "tibia_l_rotation_transform",
    "tibia_l_reset_axis",
    "tibia_l_geom_2",
    "tibia_l",
    "talus_l_parent_offset",
    "talus_l_translation",
    "tibia_l_offset_ankle_angle_l",
    "talus_l_rotation_1",
    "talus_l_rotation_2",
    "talus_l_reset_axis",
    "talus_l",
    "calcn_l_parent_offset",
    "calcn_l_translation",
    "talus_l_offset_subtalar_angle_l",
    "calcn_l_rotation_1",
    "calcn_l_rotation_2",
    "calcn_l_reset_axis",
    "calcn_l",
    "toes_l_parent_offset",
    "toes_l_translation",
    "toes_l_rotation_transform",
    "toes_l_reset_axis",
    "toes_l",
    'torso_parent_offset',
    'torso_translation',
    'torso_rotation_transform',
    'torso_reset_axis',
    'torso_geom_2',
    'torso',
    'head_and_neck_parent_offset',
    'head_and_neck_translation',
    'head_and_neck_rotation_transform',
    'head_and_neck_reset_axis',
    'head_and_neck_geom_2',
    'head_and_neck',
    'humerus_r_parent_offset',
    'humerus_r_translation',
    'humerus_r_rotation_transform',
    'humerus_r_reset_axis',
    'humerus_r',
    'ulna_r_parent_offset',
    'ulna_r_translation',
    'humerus_r_offset_elbow_flex_r',
    'ulna_r_rotation_1',
    'ulna_r_rotation_2',
    'ulna_r_reset_axis',
    'ulna_r',
    'radius_r_parent_offset',
    'radius_r_translation',
    'ulna_r_offset_pro_sup_r',
    'radius_r_rotation_1',
    'radius_r_rotation_2',
    'radius_r_reset_axis',
    'radius_r',
    'lunate_r_parent_offset',
    'lunate_r_translation',
    'lunate_r_rotation_transform',
    'lunate_r_reset_axis',
    'lunate_r',
    'hand_r_parent_offset',
    'hand_r_translation',
    'hand_r_rotation_transform',
    'hand_r_reset_axis',
    'hand_r_geom_2',
    'hand_r_geom_3',
    'hand_r_geom_4',
    'hand_r_geom_5',
    'hand_r_geom_6',
    'hand_r_geom_7',
    'hand_r_geom_8',
    'hand_r_geom_9',
    'hand_r_geom_10',
    'hand_r_geom_11',
    'hand_r_geom_12',
    'hand_r_geom_13',
    'hand_r_geom_14',
    'hand_r',
    'fingers_r_parent_offset',
    'fingers_r_translation',
    'fingers_r_rotation_transform',
    'fingers_r_reset_axis',
    'fingers_r_geom_2',
    'fingers_r_geom_3',
    'fingers_r_geom_4',
    'fingers_r_geom_5',
    'fingers_r_geom_6',
    'fingers_r_geom_7',
    'fingers_r_geom_8',
    'fingers_r_geom_9',
    'fingers_r_geom_10',
    'fingers_r_geom_11',
    'fingers_r_geom_12',
    'fingers_r',
    'humerus_l_parent_offset',
    'humerus_l_translation',
    'humerus_l_rotation_transform',
    'humerus_l_reset_axis',
    'humerus_l',
    'ulna_l_parent_offset',
    'ulna_l_translation',
    'humerus_l_offset_elbow_flex_l',
    'ulna_l_rotation_1',
    'ulna_l_rotation_2',
    'ulna_l_reset_axis',
    'ulna_l',
    'radius_l_parent_offset',
    'radius_l_translation',
    'ulna_l_offset_pro_sup_l',
    'radius_l_rotation_1',
    'radius_l_rotation_2',
    'radius_l_reset_axis',
    'radius_l',
    'lunate_l_parent_offset',
    'lunate_l_translation',
    'lunate_l_rotation_transform',
    'lunate_l_reset_axis',
    'lunate_l',
    'hand_l_parent_offset',
    'hand_l_translation',
    'hand_l_rotation_transform',
    'hand_l_reset_axis',
    'hand_l_geom_2',
    'hand_l_geom_3',
    'hand_l_geom_4',
    'hand_l_geom_5',
    'hand_l_geom_6',
    'hand_l_geom_7',
    'hand_l_geom_8',
    'hand_l_geom_9',
    'hand_l_geom_10',
    'hand_l_geom_11',
    'hand_l_geom_12',
    'hand_l_geom_13',
    'hand_l_geom_14',
    'hand_l',
    'fingers_l_parent_offset',
    'fingers_l_translation',
    'fingers_l_rotation_transform',
    'fingers_l_reset_axis',
    'fingers_l_geom_2',
    'fingers_l_geom_3',
    'fingers_l_geom_4',
    'fingers_l_geom_5',
    'fingers_l_geom_6',
    'fingers_l_geom_7',
    'fingers_l_geom_8',
    'fingers_l_geom_9',
    'fingers_l_geom_10',
    'fingers_l_geom_11',
    'fingers_l_geom_12',
    'fingers_l',
]
for segment_name in segments_to_remove:
    model.remove_segment(segment_name)
model.update_segments()
print(model.segment_names)

model.remove_muscles_without_segment()

# Create a mean gastroc
model.muscle_groups["femur_r_to_calcn_r"].merge_muscles(["med_gas_r", "lat_gas_r"])

muscles_to_remove = [
    'glut_med1_r',
    'semiten_r',
    'bifemlh_r',
    'sar_r',
    'tfl_r',
    'tib_post_r',
    'per_long_r'
]
model.remove_muscles(muscles_to_remove)
print(model.muscle_names)

model.update_muscle_groups()
print(model.muscle_group_names)

# And convert it to a .bioMod file
model.to_biomod(biomod_filepath, with_mesh=False)

# Test that the model created is valid
try:
    import biorbd
except:
    raise ImportError("You must install biorbd to load the model with biorbd")
biorbd.Model(biomod_filepath)

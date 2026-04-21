import os
import numpy as np
import import_functions
import manipfun
import useful

#################################################
#     1. Organized path and files in a dictionary
#################################################
current_folder =  os.getcwd()
measured_data_folder = os.path.join(current_folder, "num_data")
subject_name = 'BOC_09'
condition_name = 'sPra' # 'sDyn' or 'sPra'
trial_name = subject_name + '_' + condition_name

folders = {
    "main":current_folder,
    "fun": os.path.join(current_folder, "pyFun"),
    "osim_model": os.path.join(current_folder, "osimModel"),
    "measured_data": os.path.join(measured_data_folder, subject_name)
}

#################################################
#     2. neuro-musculo-skeletal: generic model
#################################################
"""
        # 2.1 neuro-musculo-skeletal (ℓom, φo, Fom, ℓst)
# 2.1.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'lst']

# 2.1.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_generic(folders["osim_model"], "wholebody.osim",mtu_params = mtu_params)

# 2.1.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
casadi_function, unknown_parameters, definition = useful.get_model_equation()

# 2.1.4 Test the model
# skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 80 # phi angle knee
q6 = 30 # phi angle ankle

q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# muscle state (a)
a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast

# test
useful.test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)

# 2.1.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)
"""

"""
        # 2.2 neuro-musculo-skeletal (ℓom, φo, Fom, kpe, ℓst, kt)
# 2.2.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'km', 'lst','kt']

# 2.2.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_generic(folders["osim_model"], "wholebody.osim",mtu_params = mtu_params)

# 2.2.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; [ℓom, φo, Fom, kpe, ℓst, kt]).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'sym',
    'lst': 'sym',
    'kt': 'sym',
}
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config)

# 2.2.4 Test the model
# skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 80 # phi angle knee
q6 = 30 # phi angle ankle
q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# muscle state (a)
a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast

# test
useful.test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)

# 2.2.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)
"""
#################################################
#     3. neuro-musculo-skeletal: scaled model
#################################################
"""
        # 3.1 neuro-musculo-skeletal (ℓom, φo, Fom, ℓst)
# 3.1.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'lst']

# 3.1.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(folders["measured_data"], f"{subject_name}.osim",mtu_params = mtu_params)

# 3.1.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
casadi_function, unknown_parameters, definition = useful.get_model_equation()

# 3.1.4 Test the model
# skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 80 # phi angle knee
q6 = 30 # phi angle ankle

q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# muscle state (a)
a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast

# test
useful.test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)

# 3.1.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)

"""

"""
        # 3.2 neuro-musculo-skeletal (ℓom, φo, Fom, kpe, ℓst, kt)
# 3.2.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'km', 'lst','kt']

# 3.2.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(folders["measured_data"], f"{subject_name}.osim",mtu_params = mtu_params)

# 3.2.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; [ℓom, φo, Fom, kpe, ℓst, kt]).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'sym',
    'lst': 'sym',
    'kt': 'sym',
}
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config)

# 3.2.4 Test the model
# skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 80 # phi angle knee
q6 = 30 # phi angle ankle
q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# muscle state (a)
a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast

# test
useful.test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)

# 3.2.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)
"""

"""
useful.plot_force_length(0,1,[1,1,1,1],casadi_function)
"""

"""
#    4. hypothetical data generator
##########################################################################

header, hypothetical_data = useful.hypothetical_data_generator(skeleton_num, muscle_tendon_parameters_num, casadi_function)

# Save path
save_dir = os.path.join(folders['main'])
os.makedirs(save_dir, exist_ok=True)

# Convert to a DataFrame
df = pd.DataFrame(hypothetical_data, columns=header)

# Save to Excel
excel_path = os.path.join(save_dir, "hypothetical_data.xlsx")
df.to_excel(excel_path, index=False)

#    4. NLP  NonLinear Programming optimisation problem (ℓom, φo, Fom, ℓst)
##########################################################################
opts = 'chosen'
initial_guess = np.array(muscle_tendon_parameters_num) * np.random.uniform(0.91, 1.09)
useful.nlp_identification(skeleton_num,muscle_tendon_parameters_num,unknown_parameters,casadi_function,hypothetical_data,'chosen',initial_guess)
"""


#    5. train and test
##########################################################################
osim_path, train_path, test_path = manipfun.select_folder_and_get_files() #  get the interest folder

#  get the name of the osim files
osim_folder, osim_name = manipfun.split_path_name(osim_path)

# import the osim scaled as a python variable
mtu_params=['l0m', 'phi0', 'f0m', 'lst']

skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(folders["measured_data"], f"{subject_name}.osim",mtu_params = mtu_params)

casadi_function, unknown_parameters, definition = useful.get_model_equation()

# import training data set and import test data set
data_train = manipfun.import_data_from_excel(train_path)
data_test = manipfun.import_data_from_excel(test_path)

# add tendon length to data set and import test data set
data_train = manipfun.add_tendon_length_to_data(data_train,skeleton_num, casadi_function)
data_test = manipfun.add_tendon_length_to_data(data_test,skeleton_num, casadi_function)

#  data verification
manipfun.plot_data(data_train, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])

# set the initial guess according to the test data
initial_guess, upper_band, lower_band = manipfun.get_initial_guess(muscle_tendon_parameters_num, data_test)

# optimisation problem
muscle_tendon_parameters_opt = useful.optimization_nlp(data_train,initial_guess,lower_band,upper_band,skeleton_num,muscle_tendon_parameters_num,unknown_parameters,casadi_function)

# opt parameters save
manipfun.save_muscle_tendon_parameters(muscle_tendon_parameters_opt,
                                  osim_folder,
                                  filename='muscle_tendon_parameters_saved',
                                  muscle_names=('tibialis', 'soleus', 'gastrocnemius'),
                                  save_csv=True,
                                  save_npy=True,
                                  verbose=True)


# generate data with the optimized parameter
header, data_est_nopt, path_csv = useful.generate_estimated_data(data_test, skeleton_num, muscle_tendon_parameters_num,
                            casadi_function,
                            output_dir='data_generic', filename='data_estime_generic.csv',
                            save_npy=False, save_csv=False, verbose=True)

manipfun.plot_data(data_est_nopt, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])
manipfun.plot_data(data_test, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])

"""
header, data_est_nopt, path_csv = useful.generate_estimated_data(data_test, skeleton_num, muscle_tendon_parameters_opt,
                            casadi_function,
                            output_dir='data_optimized', filename='data_estime_optimized.csv',
                            save_npy=False, save_csv=False, verbose=True)

manipfun.plot_data(data_est_nopt, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])
"""
# manipfun.plot_data(data_test, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])




#header_neuro_musculo_state, simulated_data_neuro_musculo_state, header_mtu_architecture, simulated_data_mtu_architecture, header_mtu_forces, simulated_data_mtu_forces = useful.simulation_(skeleton_num, muscle_tendon_parameters_opt, casadi_function,time_num, q5_num, q6_num, a_tibialis_num, a_soleus_num, a_gastrocnemius_num)




#hypothetical_data
#initial_guess

#useful.nlp_identification(skeleton_num,muscle_tendon_parameters_num,unknown_parameters,casadi_function,hypothetical_data,'chosen',initial_guess)


"""
folders = {
    "main":current_folder,
    "fun": os.path.join(current_folder, "pyFun"),
    "osim_model": os.path.join(current_folder, "osimModel"),
    "measured_data": folder_path
}
"""



"""
#    6. Simulation based on (skeleton_num, muscle_tendon_parameters_num, casadi_function,time_num, q5_num, q6_num, aTibialis_num, aSoleus_num, aGastrocnemius_num)
##########################################################################

# find measured data
##########################
measured_data_folder_name = os.path.join(folders["measured_data"], trial_name + '.xlsx')
df = pd.read_excel(measured_data_folder_name) # import your file

# Assign columns to variables
time_series = df["time (s)"]      # pandas Series time
#q5_series = df["time (s)"]      # pandas Series knee
q6_series = df["ankle angle (deg)"]      # pandas Series ankle
a_tibialis_series = df["tibialis emg"]      # pandas Series ankle
a_soleus_series = df["soleus emg"]      # pandas Series ankle
a_gastrocnemius_series = df["gastrocnemius emg"]      # pandas Series ankle

# Example: Convert to Python lists
time_num = time_series.tolist()
q6_num = q6_series.tolist()
q6_num = [-x for x in q6_num] # change orientation frame
q5_num = [0] * len(q6_num)
a_tibialis_num = a_tibialis_series.tolist()
a_soleus_num = a_soleus_series.tolist()
a_gastrocnemius_num = a_gastrocnemius_series.tolist()


# change normalization of emg from % to ratio
# if any(x < 1 for x in a_tibialis_num):
    # divide every element by 100
a_tibialis_num = [x/100 for x in a_tibialis_num]

# if any(x < 1 for x in a_soleus_num):
    # divide every element by 100
a_soleus_num = [x/100 for x in a_soleus_num]

# if any(x < 1 for x in a_gastrocnemius_num):
    # divide every element by 100
a_gastrocnemius_num = [x / 100 for x in a_gastrocnemius_num]

# define activation limit

# simulation data
##########################
header_neuro_musculo_state, simulated_data_neuro_musculo_state, header_mtu_architecture, simulated_data_mtu_architecture, header_mtu_forces, simulated_data_mtu_forces = useful.simulation_(skeleton_num, muscle_tendon_parameters_num, casadi_function,time_num, q5_num, q6_num, a_tibialis_num, a_soleus_num, a_gastrocnemius_num)


# save simulated data
##########################
    # 1. Save neuro_musculo_state
    ##########################
# Save path
save_dir = os.path.join(folders["measured_data"])
os.makedirs(save_dir, exist_ok=True)

# Convert to a DataFrame
df = pd.DataFrame(simulated_data_neuro_musculo_state, columns=header_neuro_musculo_state)

# Save to Excel
excel_name = trial_name + '_neuro_musculo_state_simulated'
excel_path = os.path.join(save_dir, excel_name + '.xlsx')
df.to_excel(excel_path, index=False)

    # 2. Save mtu_architecture
    ##########################
# Convert to a DataFrame
df = pd.DataFrame(simulated_data_mtu_architecture, columns=header_mtu_architecture)

# Save to Excel
excel_name = trial_name + '_mtu_architecture'
excel_path = os.path.join(save_dir, excel_name + '.xlsx')
df.to_excel(excel_path, index=False)

    # 3. Save mtu_forces
    ##########################
# Convert to a DataFrame
df = pd.DataFrame(simulated_data_mtu_forces, columns=header_mtu_forces)

# Save to Excel
excel_name = trial_name + '_mtu_forces_simulated'
excel_path = os.path.join(save_dir, excel_name + '.xlsx')
df.to_excel(excel_path, index=False)
print('-------------------------------------------------------------------')
print('Process succeed')
print('simulated data saved in : ')
print(save_dir)

"""



import os
import numpy as np
import pandas as pd
from pyFun import useful

########################################################################################################################

#     1. Organized path and files in a dictionnary
##########################################################################
curent_folder =  os.getcwd()
folders = {
    "main":curent_folder,
    "fun": os.path.join(curent_folder, "pyFun"),
    "osim_model": os.path.join(curent_folder, "osimModel"),
}

#     2. Load Neuromusculoskeletal
##########################################################################
# Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)
known_parameters_num, muscle_tendon_parameters_num = useful.model_osim2mat(folders["osim_model"], "wholebody.osim")

# Import Muscle Contraction Dynamics (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
casadi_function,unknown_parameters,definition = useful.DeGrooteFunction()

#     2. test the neuromusculo model
##########################################################################
q1 = 0 # x coord pelvis [ to do ]
q2 = 0 # y coord pelvis [ to do ]
q3 = 0 # z coord pelvis [ to do ]
q4 = 0 # phi angle pelvis [ to do ]
q5 = 80# phi angle knee
q6 = 30 # phi angle ankle

# convert deg in rad
q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi


#concatenate
q_num = [q1,q2,q3,q4,q5,q6]


a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast

#useful.test_model(known_parameters_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)

#     2.bis test the neuromusculo model - interactive
##########################################################################
# useful.interactive_model(known_parameters_num, muscle_tendon_parameters_num, casadi_function)

#    3. hypotetical datavgenerator
##########################################################################

header, hypotetical_data = useful.hypotetical_data_generator(known_parameters_num, muscle_tendon_parameters_num, casadi_function)

# Save path
save_dir = os.path.join(folders['main'])
os.makedirs(save_dir, exist_ok=True)

# Convert to a DataFrame
df = pd.DataFrame(hypotetical_data, columns=header)

# Save to Excel
excel_path = os.path.join(save_dir, "hypothetical_data.xlsx")
# df.to_excel(excel_path, index=False)

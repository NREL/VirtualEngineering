import numpy as np

def run_model1(ve_params):
    A = ve_params.model1_in['freq_factor']     # frequancy factor [1/c]
    Ea = ve_params.model1_in['act_energy']    # activation energy [J/mol]
    T = ve_params.model1_in['temp']     # absolute temperature [K]
    
    R = 8.3145 # Universal gas constant [J/K*mol]

    # Reaction rate (Arrhenius equation)
    k = A*np.exp(-Ea/(R*T))

    return {"k": k}

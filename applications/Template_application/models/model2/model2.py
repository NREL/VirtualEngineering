import numpy as np

def run_model2(ve_params):
    k = ve_params.model1_out['k']     # concentration [mol]
    t = ve_params.model2_in['time']    # time [s]

    # concentration equation 
    C = np.exp(-k*t)

    return {"C": C}

t = 1


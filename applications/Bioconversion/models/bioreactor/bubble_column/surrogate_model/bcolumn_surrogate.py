import sys
import os
import numpy as np
from vebio.Utilities import dict_to_yaml, yaml_to_dict
from joblib import dump, load
import matplotlib as mpl
import warnings


def run_br_surrogate(ve_params, verbose=True):

    dt = 4

    our_base = 100.0 # units? mol/m^3/hr
    r = 0.75  # what is this parameter?

    rho_g = ve_params.eh_out['rho_g']
    rho_x = ve_params.eh_out['rho_x']
    rho_f = ve_params.eh_out['rho_f']
    # what is the source for this expression?
    our_max = our_base * np.exp(-rho_f/100.0) * (r + (1.0-r) * rho_g/(rho_g+rho_x))

    gas_velocity = ve_params.br_in['gas_velocity']
    column_height = ve_params.br_in['column_height']
    column_diameter = ve_params.br_in['column_diameter']
    bubble_diameter = ve_params.br_in['bubble_diameter']
    T = ve_params.br_in['t_final']

    if verbose:
        print('\nINPUTS')
        print(f'Gas velocity = {gas_velocity:.2f} m/s')
        print(f'Column height = {column_height:.2f} m')
        print(f'Column diameter = {column_diameter:.2f} m')
        print(f'Bubble diameter = {bubble_diameter:.4f} m')
        print(f'OUR_max = {our_max:.2f} mol/m^3/hr')
        print(f't_final = {T:.1f} s')

    lb, ub = np.array([0.01, 10., 1., 5., 0.003]), np.array([0.1, 50., 6., 100., 0.008])
    X = np.array([[gas_velocity, column_height, column_diameter, our_max, bubble_diameter]])
    X = 2.*(X - lb)/(ub - lb) - 1.

    W1 = np.load(os.path.join(os.path.dirname(__file__),'W1.npy'))

    idx = T/dt
    if (idx%1 == 0):
        idx, s = int(idx), -1
    else:
        idx, s = int(idx+1), idx%1

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        with open(os.path.join(os.path.dirname(__file__), 'gp_bub_col.pkl'), 'rb') as f_id:
            for i in range(idx+1):
                gp = load(f_id)
                if (s > 0) and (i == idx-1):
                    f_0 = gp.predict(X@W1, return_std=False)[0]

            if (s == -1):
                ff = np.power(10., gp.predict(X@W1, return_std=False)[0])
            else:
                f_1 = gp.predict(X@W1, return_std=False)[0]
                ff = np.power(10., f_0*(1-s) + f_1*s)

            ff = our_max * ff

    if verbose:
        print('\nFINAL OUTPUTS (at t = %.1f seconds)' % (T))
        print(F'OUR = {ff:.4f} mol/m^3/hr')

    # Save the outputs into a dictionary
    output_dict = {}
    output_dict['our'] = float(ff)
    
    return output_dict

if __name__ == '__main__':
    if len(sys.argv) > 1:
        # TODO: fix it
        # params_filename = sys.argv[1]
        # main(params_filename)
        pass
    else:
        raise Exception("VE parameters filename not provided. This model must be called with inputs specified via filename")
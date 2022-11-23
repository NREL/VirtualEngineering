import sys
import numpy as np
from vebio.Utilities import dict_to_yaml, yaml_to_dict
from joblib import dump, load
import matplotlib as mpl


def main(params_filename):

    ve_params = yaml_to_dict(params_filename)

    dt = 4

    our_base = 100.0 # units?
    r = 0.75  # what is this parameter?

    rho_g = ve_params['enzymatic_output']['rho_g']
    rho_x = ve_params['enzymatic_output']['rho_x']
    rho_f = ve_params['enzymatic_output']['rho_f']
    # what is the source for this expression?
    our_max = our_base * np.exp(-rho_f/100.0) * (r + (1.0-r) * rho_g/(rho_g+rho_x))

    gas_velocity = ve_params['bioreactor_input']['gas_velocity']
    column_height = ve_params['bioreactor_input']['column_height']
    column_diameter = ve_params['bioreactor_input']['column_diameter']
    bubble_diameter = ve_params['bioreactor_input']['bubble_diameter']

    T = ve_params['bioreactor_input']['t_final']

    # these outputs should show units 
    print('\nINPUTS')
    print(f'Gas velocity = {gas_velocity:.2f} m/s')
    print(f'Column height = {column_height:.2f} m')
    print(f'Column diameter = {column_diameter:.2f} m')
    print(f'Bubble diameter = {bubble_diameter:.4f} m')
    print(f'OUR_max = {our_max:.2f}' % (our_max))
    print(f't_final = {T:.1f} s')

    lb, ub = np.array([0.01, 10., 1., 5., 0.003]), np.array([0.1, 50., 6., 100., 0.008])
    X = np.array([[gas_velocity, column_height, column_diameter, our_max, bubble_diameter]])
    X = 2.*(X - lb)/(ub - lb) - 1.

    W1 = np.load('W1.npy')

    idx = T/dt
    if (idx%1 == 0):
        idx, s = int(idx), -1
    else:
        idx, s = int(idx+1), idx%1

    with open('gp_bub_col.pkl', 'rb') as f_id:
        for i in range(idx+1):
            gp = load(f_id)
            if (s > 0) and (i == idx-1):
                f_0 = gp.predict(X@W1, return_std=False)[0]

        if (s == -1):
            ff = np.power(10., gp.predict(X@W1, return_std=False)[0])
        else:
            f_1 = gp.predict(X@W1, return_std=False)[0]
            ff = np.power(10., f_0*(1-s) + f_1*s)
        
        ff = our_max * ff/(0.01 + ff)

    print('\nFINAL OUTPUTS (at t = %.1f seconds)' % (T))
    print('OUR = %.4f' % (ff))  # units???

    # Save the outputs into a dictionary
    output_dict = {}
    output_dict['our'] = float(ff)
    ve_params['bioreactor_output'] = output_dict
    dict_to_yaml(ve_params, params_filename)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        params_filename = sys.argv[1]
        print(params_filename)
        main(params_filename)
    else:
        raise Exception("VE parameters filename not provided. This model must be called with inputs specified via filename")
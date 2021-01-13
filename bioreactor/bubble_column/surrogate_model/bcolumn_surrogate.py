import sys
import numpy as np
from vebio.Utilities import dict_to_yaml, yaml_to_dict
from joblib import dump, load

if len(sys.argv) > 1:
    params_filename = sys.argv[1]
    ve_params = yaml_to_dict(params_filename)

else:
    ve_params = {}

font={'family':'Helvetica', 'size':'15'}
mpl.rc('font',**font)
mpl.rc('xtick',labelsize=14)
mpl.rc('ytick',labelsize=14)

our_base = 100.0
r = 0.75

gas_velocity = 0.08
column_height = 40.
column_diameter = 5.
bubble_diameter = 0.006

rho_g = ve_params['enzymatic_output']['rho_g']
rho_x = ve_params['enzymatic_output']['rho_x']
rho_f = ve_params['enzymatic_output']['rho_f']
our_max = our_base * np.exp(-rho_f/100.0) * (r + (1.0-r) * rho_g/(rho_g+rho_x));

T = ve_params['bioreactor_input']['t_final']

print('\nINPUTS')
print('Gas velocity = %.2f' % (gas_velocity))
print('Column height = %.2f' % (column_height))
print('Column diameter = %.2f' % (column_diameter))
print('Bubble diameter = %.4f' % (bubble_diameter))
print('OUR_max = %.2f' % (our_max))
print('t_final = %.1f' % (T))

lb, ub = np.array([0.01, 10., 1., 5., 0.003]), np.array([0.1, 50., 6., 100., 0.008])
X = np.array([[gas_velocity, column_height, column_diameter, our_max, bubble_diameter]])
X = 2.*(X - lb)/(ub - lb) - 1.

load_model = 'gp_model/gp_T{0:03d}.pkl'.format(int(T))
W1 = np.load('W1.npy')

gp = load(load_model)

ff = gp.predict(X@W1, return_std=False)
ff = np.power(10., ff[0, 0])

print('\nFINAL OUTPUTS (at t = %.1f seconds)' % (T))
print('OUR = %.4f' % (ff))

# Save the outputs into a dictionary
output_dict = {'bioreactor_output': {}}
output_dict['bioreactor_output']['our'] = float(ff)

dict_to_yaml([ve_params, output_dict], params_filename)
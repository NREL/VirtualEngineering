import sys
import numpy as np
from vebio.Utilities import dict_to_yaml, yaml_to_dict
from joblib import dump, load
import matplotlib as mpl

if len(sys.argv) > 1:
    params_filename = sys.argv[1]
    ve_params = yaml_to_dict(params_filename)

else:
    ve_params = {}

font={'family':'Helvetica', 'size':'15'}
mpl.rc('font',**font)
mpl.rc('xtick',labelsize=14)
mpl.rc('ytick',labelsize=14)

dt = 4

our_base = 100.0
r = 0.75

gas_velocity = 0.08
column_height = 40.
column_diameter = 5.
bubble_diameter = 0.006

fis0 = ve_params['enzymatic_input']['fis_0']
xG0 = ve_params['pretreatment_output']['X_G']
xX0 = ve_params['pretreatment_output']['X_X']
yF0 = 0.2 + 0.6*ve_params['pretreatment_output']['conv']
lmbdE = ve_params['enzymatic_input']['lambda_e']
dilution_factor = ve_params['enzymatic_input']['fis_0']/ve_params['pretreatment_output']['fis_0']
rhox0 = ve_params['pretreatment_output']['rho_x']*dilution_factor
omega = 3.14

rho_f = ve_params['pretreatment_output']['rho_f']

T = ve_params['enzymatic_input']['t_final']
'''
'insoluble_solids', 'solid_cellulose', 'solid_xylan', 'facile_glucan',
                    'enzyme_loading', 'xylose_conc', 'rot_vel'
print('\nINPUTS')
print('Gas velocity = %.2f' % (gas_velocity))
print('Column height = %.2f' % (column_height))
print('Column diameter = %.2f' % (column_diameter))
print('Bubble diameter = %.4f' % (bubble_diameter))
print('OUR_max = %.2f' % (our_max))
print('t_final = %.1f' % (T))
'''
lb = np.array([0.005, 0.3, 0.02, 0.2, 0.01,  5.,  1.04])
ub = np.array([ 0.07, 0.7, 0.25, 0.7, 0.10, 60., 10.40 ])
X = np.array([[fis0, xG0, xX0, yF0, lmbdE, rhox0, omega]])
X = 2.*(X - lb)/(ub - lb) - 1.

W1 = np.array([[[ 0.67813919,  0.48452556,  0.71445575,  0.66378877],
                [ 0.67638728,  0.77337388,  0.60977278,  0.24320546]],
               [[ 0.32910500,  0.60687455,  0.06172572, -0.26018624],
                [-0.28508414, -0.56751307,  0.12671069, -0.40636609]],
               [[ 0.16574545, -0.02901924,  0.62966290,  0.20960751],
                [-0.27097742,  0.11567439, -0.76165208, -0.83617288]],
               [[ 0.12529040,  0.31934871, -0.03301496, -0.14228018],
                [-0.09845739, -0.18129418,  0.02277099,  0.03171999]],
               [[ 0.52107493,  0.45239377,  0.21135817,  0.48773467],
                [-0.55788726,  0.02396285,  0.17405732, -0.19686160]],
               [[-0.33937643, -0.29871514, -0.20532291, -0.34307791],
                [ 0.08177255, -0.05872661, -0.00641201,  0.10185202]],
               [[-0.04418668, -0.01488922,  0.03671227, -0.26818104],
                [-0.24531555, -0.17190280, -0.03394246, -0.16247506]]])

idx = 4*(int(T)-1)
s = -1 if (T%1 == 0) else T%1
with open('gp_EH.pkl', 'rb') as f_id:
    for i in range(idx):
        gp = load(f_id)

    gp = load(f_id)
    f0_convRate = gp.predict(X@W1[:, :, 0], return_std=False)[0]
    gp = load(f_id)
    f0_rhoG = gp.predict(X@W1[:, :, 1], return_std=False)[0]
    gp = load(f_id)
    f0_rhoX = gp.predict(X@W1[:, :, 2], return_std=False)[0]
    gp = load(f_id)
    f0_rhoL = gp.predict(X@W1[:, :, 3], return_std=False)[0]

    if (s > 0):
        gp = load(f_id)
        f1_convRate = gp.predict(X@W1[:, :, 0], return_std=False)[0]
        gp = load(f_id)
        f1_rhoG = gp.predict(X@W1[:, :, 1], return_std=False)[0]
        gp = load(f_id)
        f1_rhoX = gp.predict(X@W1[:, :, 2], return_std=False)[0]
        gp = load(f_id)
        f1_rhoL = gp.predict(X@W1[:, :, 3], return_std=False)[0]

        f_convRate = f0_convRate*(1-s) + f1_convRate*s
        f_rhoG = f0_rhoG*(1-s) + f1_rhoG*s
        f_rhoX = f0_rhoX*(1-s) + f1_rhoX*s
        f_rhoL = f0_rhoL*(1-s) + f1_rhoL*s
    else:
        f_convRate = f0_convRate
        f_rhoG = f0_rhoG
        f_rhoX = f0_rhoX
        f_rhoL = f0_rhoL

    f_convRate = np.power(10., f_convRate)
    f_rhoG = fis0*f_rhoG
    f_rhoX = f_rhoX + rhox0
    f_rhoL = fis0*f_rhoL

#f_rhoF = 

print('\nFINAL OUTPUTS (at t = %4.4g hours)' % (T))
print('Conversion Rate = %.4f' % f_convRate)
print('rho_g = %4.4g g/L' % f_rhoG)
print('rho_x = %4.4g g/L' % f_rhoX)
print('rho_sL = %4.4g g/L' % f_rhoL)
print('rho_f = %4.4g g/L' % rho_f)

# Save the outputs into a dictionary for use as inputs for bioreactor sims
output_dict = {'enzymatic_output': {}}
output_dict['enzymatic_output']['rho_g'] = float(f_rhoG)
output_dict['enzymatic_output']['rho_x'] = float(f_rhoX)
# Soluble lignin is not currently used in bioreaction models, but we might want to.
output_dict['enzymatic_output']['rho_sL'] = float(f_rhoL)
'''
# compute dilution of furfural (which does not react during EH)
fliq0 = 1 - fis0
fliqend = 1 - result["fis"].iloc[-1]
rhof = fliq0/fliqend*rhof0
output_dict['enzymatic_output']['rho_f'] = float(rhof)
'''
output_dict['enzymatic_output']['rho_f'] = float(rho_f)

dict_to_yaml([ve_params, output_dict], params_filename)
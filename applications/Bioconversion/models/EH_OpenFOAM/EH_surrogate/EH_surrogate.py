import sys
import os
import numpy as np
from vebio.Utilities import dict_to_yaml, yaml_to_dict
from joblib import dump, load
import matplotlib as mpl
import warnings

def run_eh(ve_params, verbose=True):

    fis0 = ve_params.eh_in['fis_0']
    xG0 = ve_params.pt_out['X_G']
    xX0 = ve_params.pt_out['X_X']
    yF0 = 0.2 + 0.6*ve_params.pt_out['conv']
    lmbdE = ve_params.eh_in['lambda_e']
    dilution_factor = ve_params.eh_in['fis_0']/ve_params.pt_out['fis_0']
    rhox0 = ve_params.pt_out['rho_x']*dilution_factor
    omega = 3.14
    rho_f = ve_params.pt_out['rho_f']*dilution_factor

    T = ve_params.eh_in['t_final']
    '''
    'insoluble_solids', 'solid_cellulose', 'solid_xylan', 'facile_glucan',
                        'enzyme_loading', 'xylose_conc', 'rot_vel'
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        with open(os.path.join(os.path.dirname(__file__), 'gp_EH.pkl'), 'rb') as f_id:
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
    if verbose:
        print('\nFINAL OUTPUTS (at t = %4.4g hours)' % (T))
        print('Conversion Rate = %.4f' % f_convRate)
        print('rho_g = %4.4g g/L' % f_rhoG)
        print('rho_x = %4.4g g/L' % f_rhoX)
        print('rho_sL = %4.4g g/L' % f_rhoL)
        print('rho_f = %4.4g g/L' % rho_f)

    # Save the outputs into a dictionary for use as inputs for bioreactor sims
    output_dict = {}
    output_dict['rho_g'] = float(f_rhoG)
    output_dict['rho_x'] = float(f_rhoX)
    # Soluble lignin is not currently used in bioreaction models, but we might want to.
    output_dict['rho_sL'] = float(f_rhoL)
    '''
    # compute dilution of furfural (which does not react during EH)
    fliq0 = 1 - fis0
    fliqend = 1 - result["fis"].iloc[-1]
    rhof = fliq0/fliqend*rhof0
    output_dict['enzymatic_output']['rho_f'] = float(rhof)
    '''
    output_dict['rho_f'] = float(rho_f)

    return output_dict


if __name__ == '__main__':
    if len(sys.argv) > 1:
        # TODO: fix it
        # params_filename = sys.argv[1]
        # main(params_filename)
        pass
    else:
        raise Exception("VE parameters filename not provided. This model must be called with inputs specified via filename")
        
        
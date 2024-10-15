"""
Define initial conditions and run lignocellulose EH simulation. Uses modules in
`CEH_EmpiricalModel` (a git repo submodule). This script is based of
`driver_batch_lignocellulose.py` in that repo but has been modified for use in
Virtual Engineering.
"""

# Jim Lischeske and Jonathan Stickel, 2017-2021


# TODO:
# - presume soluble lignin was produced during pretreatment in proportion to fufural?

import numpy as np
import sys
import matplotlib as mpl
from matplotlib import pyplot as plt

import ehk_batch as ehk


def run_eh_lingocell(ve_params, show_plots=False):

    ### Set kinetics parameters - these are from Lischeske and Stickel, 2019
    fitOuts = {'KdR': 0.05000000000000001,
               'kL': 729.4513412205487,
               'kR': 14712.525849904296,
               'kF': 14712.525849904296,
               'kX': 9999.999988133452,
               'kappaRF': 9.33804072835234,
               'kappaRL': 49.99999999999999,
               'kappaRX': 11.281636078910811,
               'kappaRs': 49.99999999971725}
    ### Initialize batch object
    batch = ehk.Batch(**fitOuts) # initial creation

    ### Set conditions for modeled hydrolysate systemfrom user input
    lmbdE = ve_params.eh_in['lambda_e'] # [g/g], enzyme loading
    fis0 = ve_params.eh_in['fis_0'] # this is a target, not the output
                                                # from pretreatment
    tfin = ve_params.eh_in['t_final']
    # from pretreatment
    XG0 = ve_params.pt_out['X_G']
    XX0 = ve_params.pt_out['X_X']
    # Compute the amount of dilution required to reach the fis_0_target based on
    # the output from the pretreatment step
    dilution_factor = fis0/ve_params.pt_out['fis_0']
    rhog0 = 0.0 # g/L; no conversion of glucan in current PT -- this should be specified in
                # params file! JJS 3/14/21
    rhox0 = ve_params.pt_out['rho_x']*dilution_factor
    rhof0 = ve_params.pt_out['rho_f']*dilution_factor # furfural
    conversion_xylan = ve_params.pt_out['conv'] 
    yF0 = 0.2 + 0.6*conversion_xylan # our simple "empirical guess" model

    print('\nINPUTS')
    print('Lambda_e = %4.4g' % (lmbdE))
    print('FIS_0 = %4.4g' % (fis0))
    print('yF0 = %4.4g' % (yF0))
    print('t_final = %4.4g' % (tfin))

    conditions = {'fis0': fis0,
                  'XG0': XG0,
                  'XX0': XX0,
                  'XL0': 1 - XX0 - XG0,
                  'yF0': yF0,
                  'lmbdE': lmbdE,
                  'rhog0': rhog0,
                  'rhox0': rhox0,
                  'rhosL0': 0}
    batch.SetParams(conditions)


    # run simulation
    N = 200
    t = np.linspace(0, tfin, N)

    batch.BatchSimulate(t)
    result = batch.SimOut
    # use `result.keys()` to see the keys of the data frame

    # Save the outputs into a dictionary for use as inputs for bioreactor sims
    output_dict = {}
    output_dict['rho_g'] = float(result["rhog"].iloc[-1])
    output_dict['rho_x'] = float(result["rhox"].iloc[-1])
    # Soluble lignin is not currently used in bioreaction models, but we might want to.
    output_dict['rho_sL'] = float(result["rhosL"].iloc[-1])
    # compute dilution of furfural (which does not react during EH)
    fliq0 = 1 - fis0
    fliqend = 1 - result["fis"].iloc[-1]
    rhof = fliq0/fliqend*rhof0
    output_dict['rho_f'] = float(rhof)

    print('\nFINAL OUTPUTS (at t = %4.4g hours)' % (tfin))
    print('rho_g = %4.4g g/L' % (result["rhog"].iloc[-1]))
    print('rho_x = %4.4g g/L' % (result["rhox"].iloc[-1]))
    print('rho_f = %4.4g g/L' % rhof)
    # print('Facile Conversion = ', convF[-1])
    # print('Recalcitrant Conversion = ', convR[-1])
    print("Final FIS = %4.4g" % (result['fis'].iloc[-1]))
    print("Glucan conversion = %4.4g" % result["Gconv"].iloc[-1])
    print('Total carbohydrate conversion = %4.4g' % (result["Tconv"].iloc[-1]))
    mbmax = max(result['mbG'].max(), result['mbX'].max(), result['mbL'].max(),
                result['mbE'].max())
    print("Max mass balance error:  %4.4g" % mbmax)
    print() 

    if show_plots:

        # adjust fonts for nice display in Jupyter notebook
        font={'family':'Sans', 'size':'15'}
        mpl.rc('font',**font)
        mpl.rc('xtick',labelsize=14)
        mpl.rc('ytick',labelsize=14)

        plt.figure(1)
        plt.clf()
        plt.plot(t, result['Gconv'], lw=2, label='total glucan')
        plt.plot(t, result['GconvF'], label='facile')
        plt.plot(t, result['GconvR'], label='recalc')
        plt.plot(t, result['Tconv'], lw=2, label='total carbohydrate')
        plt.xlabel('t [h]')
        plt.ylabel('conversion')
        plt.legend(loc='best')
        
        plt.figure(2)
        plt.clf()
        plt.plot(t, result['rhog'], lw=2, label='glucose')
        plt.plot(t, result['rhox'], lw=2, label='xylose')
        plt.plot(t, result['rhosL'], lw=2, label='sol lignin')
        plt.xlabel('t [h]')
        plt.ylabel(r'$\rho_i$ [g/L]')
        plt.legend(loc='best')

        # plt.figure(3)
        # plt.clf()
        # plt.plot(t, result['mbG'], label='glucan')
        # plt.plot(t, result['mbX'], label='xylan')
        # plt.plot(t, result['mbL'], label='lignin')
        # plt.plot(t, result['mbE'], label='enzyme')
        # plt.xlabel('t [h]')
        # plt.ylabel('mass balance')
        # plt.legend(loc='best')

        plt.show()
        
        # plt.figure(4)
        # plt.clf()
        # plt.plot(result['Tconv'], result['fis'])
        # plt.xlabel('carbohydrate conversion')
        # plt.ylabel(r'$f_{is}$')

        # plt.figure(5)
        # plt.clf()
        # plt.plot(t, result['Tconv'], label='total carbohydrate')
        # plt.plot(t, result['Lconv'], label='lignin')
        # plt.xlabel('t [h]')
        # plt.ylabel('conversion')
        # plt.legend(loc='best')

        # plt.figure(6)
        # plt.clf()
        # plt.plot(t, fEA/f[:,7], label='adsorbed enzyme')
        # plt.xlabel('t [h]')
        # plt.ylabel('fraction adsorbed')
        # plt.legend(loc='best')

    return output_dict

if __name__ == '__main__':
    if len(sys.argv) > 1:
        # TODO: fix it
        # params_filename = sys.argv[1]
        # show_plots=False
        # if len(sys.argv) > 2:
        #     show_plots = sys.argv[2]
        # main(params_filename, show_plots)
        pass
    else:
        raise Exception("VE parameters filename not provided. This model must be called with inputs specified via filename")
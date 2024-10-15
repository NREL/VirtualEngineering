"""
Define initial conditions and run lignocellulose EH simulation


    #### work for CEHD FY21 Q4 milestone simulating CEH for TEA ####
    # from 2021-03 CEHD milestone report table 2

    (copy of driver_batch_lignocellulose.py)

"""

# Jim Lischeske and Jonathan Stickel, 2017-2021

# updated to use new reaction-classes approach, JJS 11/2/17

# updated to work with Jim's changes to `ehk_batch`; should be able to test output against the version in `two_phase_approach`

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

from pprint import pprint # pretty print

import ehk_batch as ehk

from vebio.Utilities import dict_to_yaml, yaml_to_dict

if len(sys.argv) > 1:
    params_filename = sys.argv[1]
    ve_params = yaml_to_dict(params_filename)
else:
    raise Exception("VE parameters filename not provided")


### Initialize batch object
#batch = ehk.Batch()

#batch = ehk.Batch(kL=100, kappaRF=5.) # this works to set kinetics parameters
#batch = ehk.Batch(XG0=0.7) # this doesn't work

# # or create parameters object and provide that at initiation of the simulation object
# P = ehk.KinParams(kL=100)
# pprint(P.__dict__) # inspect P
# batch = ehk.Batch(P)

### Set kinetics parameters
fitOuts = { 'KdR': 0.5,  # typo? should this be 10x smaller?
            'kL': 500, #1467, # Jim's value seems too large
            'kR': 3843,
            'kF': 3843,  # same as recalc?
            'kX': 1e4,
            'kappaRF': 2,
            'kappaRL': 50,
            'kappaRX': 0.84,
            'kappaRs': 1.52,
            'mGR1': 1.,
            'mX1': 1.}
#batch.SetParams(fitOuts) # if already created batch object, kinda redundant
batch = ehk.Batch(**fitOuts) # initial creation


# ### Set conditions for modeled hydrolysate system
conditions = {'fis0': ve_params['CEH_input']['fis_0'],  # initial FIS, g IS / g slurry
              'XG0': ve_params['CEH_input']['glucan_solid_fraction'],  # initial solids glucan fraction, g glucan / g IS (includes all hexoses)
              'XX0': ve_params['CEH_input']['xylan_solid_fraction'],  # initial xylan fraction (includes all pentoses)
              'XL0': ve_params['CEH_input']['lignin_solid_fraction'],  # initial lignin fraction (note that currently it is overwritten, something I want to fix)
              'yF0': ve_params['CEH_input']['facile_fraction_glucan'],  # initial fraction of glucan that is "facile", g facile / g glucan
              'lmbdE': ve_params['CEH_input']['lambda_e'], # [g/g], enzyme loading,  # enzyme loading, g enzyme / g glucan (might be better to change this definition to be g enzyme / g carbohydrates, but currently it is a glucan-only basis)
              'rhog0': 0.0,  # these are initial soluble species concentrations ("densities"), g / L liquid
              'rhox0': 0.0,
              'rhosL0': 0.0}
batch.SetParams(conditions)

# run simulation
# tfin = 5*24 # h
tfin = ve_params['CEH_input']['t_final'] # h
N = 100
t = np.linspace(0, tfin, N)

batch.BatchSimulate(t)


result = batch.SimOut


if True:
    PFR_result = batch.SimOut.iloc[-1,:]  # last line is result
    PFR_result.to_pickle('PFR_result.pkl')  # save for CEH sim inputs

    import pickle
    pickle.dump(fitOuts, open('FY21_kin_params.pkl', 'wb'))
    pickle.dump(conditions, open('FY21_init_conditions.pkl', 'wb'))


## I don't remember what this block is for -- JJL 3/4/21
# fEA = np.zeros(N)
# for i in range(N):
#     CEA = batch.Adsorption(f[i,:])
#     fEA[i] = ehk.MwE/ehk.rhoT*np.sum(CEA)
# fEf = f[:,7]-fEA


if ve_params['CEH_input']['show_plots']:
    plt.figure(1)
    plt.clf()

    plt.plot(t, result['Gconv'], lw=2, label='total')
    plt.plot(t, result['GconvF'], label='facile')
    plt.plot(t, result['GconvR'], label='recalc')
    plt.xlabel('t [h]')
    plt.ylabel('glucan conversion')
    plt.legend(loc='best')
    plt.draw()
    plt.show()

    plt.figure(2)
    plt.clf()
    plt.plot(t, result['mbG'], label='glucan')
    plt.plot(t, result['mbX'], label='xylan')
    plt.plot(t, result['mbL'], label='lignin')
    plt.plot(t, result['mbE'], label='enzyme')
    plt.xlabel('t [h]')
    plt.ylabel('mass balance')
    plt.legend(loc='best')
    plt.draw()
    plt.show()

    plt.figure(3)
    plt.clf()
    plt.plot(t, result['rhog'], lw=2, label='glucose')
    plt.plot(t, result['rhox'], lw=2, label='xylose')
    plt.plot(t, result['rhosL'], lw=2, label='sol lignin')
    plt.xlabel('t [h]')
    plt.ylabel(r'$\rho_i$ [g/L]')
    plt.legend(loc='best')
    plt.draw()
    plt.show()
    #plt.savefig('fig.pdf', bbox_inches='tight')

    plt.figure(4)
    plt.clf()
    # plt.plot(result['Tconv'], result['fis'])
    # plt.xlabel('carbohydrate conversion')
    plt.plot(t, result['fis'])
    plt.xlabel('time (h)')
    plt.ylabel(r'$f_{is}$')
    plt.draw()
    plt.show()

    plt.figure(5)
    plt.clf()
    plt.plot(t, result['Tconv'], label='total carbohydrate')
    plt.plot(t, result['Lconv'], label='lignin')
    plt.xlabel('t [h]')
    plt.ylabel('conversion')
    plt.legend(loc='best')
    plt.ylim(0,1)
    plt.draw()
    plt.show()

    # plt.figure(6)
    # plt.clf()
    # plt.plot(t, fEA/f[:,7], label='adsorbed enzyme')
    # plt.xlabel('t [h]')
    # plt.ylabel('fraction adsorbed')
    # plt.legend(loc='best')
    # plt.draw()
    # plt.show()



if False: # save results?
    np.savez('simresults.npz', **result)

"""
Define initial conditions and run lignocellulose EH simulation
"""

# Jim Lischeske and Jonathan Stickel, 2017

# updated to use new reaction-classes approach, JJS 11/2/17

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

from pprint import pprint # pretty print

import ehk_batch as ehk




# create simulation object with default or keyword parameters:

# Set constants equal to recent fitting attempt (2019--JL)
fitOuts = {'KdR': 0.05000000000000001,
           'kL': 729.4513412205487,
           'kR': 14712.525849904296,
           'kX': 9999.999988133452,
           'kappaRF': 9.33804072835234,
           'kappaRL': 49.99999999999999,
           'kappaRX': 11.281636078910811,
           'kappaRs': 49.99999999971725}
fitOuts['kF'] = fitOuts['kR']
fitOuts['mGR1'] = 1.
fitOuts['mX1'] = 1.

# Initial values for system

bchConditions = {'fis': 0.1,
                 'XG0': 0.65,
                 'XX0': 0.07,
                 'XL0': 1-0.72,
                 'rhog0': 5.0,
                 'rhox0': 30.0,
                 'rhosL0': 10.,
                 'lmbdE': 0.02,
                 'yF0': 0.5999999999999999}


batch = ehk.Batch()

# assign fitOuts to batch-model object
for key, value in fitOuts.items():
    if hasattr(batch, key):
        setattr(batch, key, value)
    else:
        print('not setting a value: '+key)
# assign bchConditions to conditions object within batch-model object  
for key, value in bchConditions.items():
    if hasattr(batch.Conditions, key):
        setattr(batch.Conditions, key, value)
    else:
        print('not setting a value: '+key)

#f0 = batch.GetF0(fis0, XG0, XX0, XL0, yF0, lmbdE, rhog0, rhox0, rhosL0)

tfin = 200. # h
N = 200
t = np.linspace(0, tfin, N)

batch.BatchSimulate(t)

batch.CalcResultValues()
result = batch.SimOut

# fEA = np.zeros(N)
# for i in range(N):
#     CEA = batch.Adsorption(f[i,:])
#     fEA[i] = ehk.MwE/ehk.rhoT*np.sum(CEA)
# fEf = f[:,7]-fEA

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
plt.plot(result['Tconv'], result['fis'])
plt.xlabel('carbohydrate conversion')
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



# if False: # save results?
#     np.savez('simresults.npz', **result)

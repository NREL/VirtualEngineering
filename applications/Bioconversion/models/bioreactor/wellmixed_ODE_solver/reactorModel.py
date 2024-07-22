##############################################################
#
# Here, we assume a constant oxygen distribution for all times
# within the reactor. We then use three different paradigms to
# simulate the reactor: (1) simple mean O2 concentration; (2)
# all species (except oxygen) are held constant for all times
# across the reactor, and the instantanesou reactor rate is
# the weighted sum of rates; and (3) each point on the
# distribution is considered it's own rate bucket, and is
# simulated from t_0 to t_f. This final scenario is analogous
# to the "worst-case" traditional leapfrogging scenario.
#
# James Lischeske
# 3/3/2020
#
#
##############################################################


import numpy as np
from scipy.stats import beta
import batchModel as bm
import matplotlib.pyplot as plt

import time

# ok, the beta distribution has domain (0,1). So, to map
# from the distribution to our possible concentrations (0,mx),
# we do C = C_mx * X, where X~beta(a,b)
#
# we also have to center our concentrations.
X_mean = 0.2
betaN = 10
alphaN = X_mean * betaN/ (1-X_mean)
xspace = np.linspace(0,1, num=20)
x_mp = (xspace[1:]+xspace[:-1])/2
X_ppf = beta.ppf(x_mp, alphaN, betaN)

# Equal-weighted concentrations
c_mx = 0.224
conc = X_ppf*c_mx
n = len(conc)

plt.figure(1)
plt.clf()

plt.plot(conc, beta.pdf(conc/c_mx, alphaN, betaN), ls='', marker='o')
plt.plot(conc, beta.cdf(conc/c_mx, alphaN, betaN))

plt.draw()
plt.show()
plt.savefig('fig1-O2_Dist.pdf', bbox_inches='tight')


#### simulation header ####
f0 = np.array([0.5, 500, 250, 0,0]) # X, G, Xy, A, B
mc = bm.modelConstants


dt = 0.01
t0 = 0.
tf = 36


#### Do simulation, o2=conc_avg ####
cpu_t0 = time.time()
t = [t0]
fmat = f0
fold = f0
c_avg = X_mean*c_mx
while t[-1]<tf:
    fnew = fold + dt*bm.dfdt(np.hstack((c_avg,fold)),0,mc)[1:]

    fmat = np.append(fmat, fnew)
    fold = fnew
    t.append(t[-1]+dt)

fmat_avg = np.reshape(fmat, (-1,5))
t_avg = np.asarray(t)

cpu_avg_tblock = time.time()-cpu_t0


#### Do simulation, summed ####
cpu_t0 = time.time()

t = [t0]
fmat = f0
fold = f0
while t[-1]<tf:
    dfdt = np.zeros(5)
    for c in conc:
        dfdt += bm.dfdt(np.hstack((c,fold)), 0, mc)[1:]
    dfdt /= n 
    fnew = fold + dt*dfdt

    fmat = np.append(fmat, fnew)
    fold = fnew
    t.append(t[-1]+dt)

fmat_sum = np.reshape(fmat, (-1,5))
t_sum = np.asarray(t)

cpu_sum_tblock = time.time()-cpu_t0


#### Do simulation, in "buckets", summed at the end ####
cpu_t0 = time.time()
fmat = fmat*0

for c in conc:
    t = [t0]
    fmat_tmp = f0
    fold = f0
    while t[-1]<tf:
        fnew = fold + dt*bm.dfdt(np.hstack((c,fold)),0,mc)[1:]
        
        fmat_tmp = np.append(fmat_tmp, fnew)
        fold = fnew
        t.append(t[-1]+dt)
    fmat += fmat_tmp
fmat /= n

fmat_buck = np.reshape(fmat, (-1,5))
t = np.asarray(t)

cpu_buck_tblock = time.time()-cpu_t0

print('number of discretizations: ')
print(n)

print('t_avg: %.2f\n' %cpu_avg_tblock)
print('t_sum: %.2f\n' %cpu_sum_tblock)
print('t_buck: %.2f\n' %cpu_buck_tblock)



plt.figure(2)
plt.clf()
plt.plot(t, fmat_avg[:,0], color='C0', ls='-', label='avg')
plt.plot(t, fmat_sum[:,0], color='C0', ls='--', label='sum')
plt.plot(t, fmat_buck[:,0], color='C0', ls=':', label='bucket')
plt.legend(loc='best')

plt.xlabel('time (h)')
plt.ylabel('biomass (kg/$m^3$)')
plt.show()
plt.savefig('fig2-Biomass.pdf', bbox_inches='tight')


plt.figure(3)
plt.clf()
plt.plot(t, fmat_avg[:,1], color='C0', label='glu')
plt.plot(t, fmat_avg[:,2], color='C1', label='xyl')
plt.plot(t, fmat_avg[:,3], color='C2', label='ace')
plt.plot(t, fmat_avg[:,4], color='C3', label='bdo')

plt.plot(t, fmat_sum[:,1], color='C0', ls='--')
plt.plot(t, fmat_sum[:,2], color='C1', ls='--')
plt.plot(t, fmat_sum[:,3], color='C2', ls='--')
plt.plot(t, fmat_sum[:,4], color='C3', ls='--')

plt.plot(t, fmat_buck[:,1], color='C0', ls=':')
plt.plot(t, fmat_buck[:,2], color='C1', ls=':')
plt.plot(t, fmat_buck[:,3], color='C2', ls=':')
plt.plot(t, fmat_buck[:,4], color='C3', ls=':')


plt.legend(loc='best')
plt.xlabel('time (h)')
plt.ylabel('substrates (mM)')
plt.show()
plt.savefig('fig3-soluble_substrates.pdf', bbox_inches='tight')


"""
Simple stepping model for fed-batch bioreaction to make lipid
"""

import numpy as np
from matplotlib.pyplot import *


# approximate experimental data from Fei et al. (2016) 
texp = np.array([0, 24, 48, 72, 96])
DCWexp = np.array([0, 20, 33, 51, 53])
lipexp = np.array([0, 7, 14, 28, 30])

# hack to manually bring in CFD results -- need to setup options to run CFD
# periodically or at specified times to get back reactor-average OUR
OURcfd = np.array([4.54, 8.21]) # mol/m^3/h, per Hari's runs
fOUR = 0.93 # OUR_actual/OUR_max; average of [4.54/4.8, 8.21/8.96]

# stoichiometry
ng = 1 # mol glucose
nl = 0.2 # mol lipid
no = 1.75 # mol oxygen
nr_l_o = nl/no # mole ratio lipid/oxygen
nr_g_o = ng/no # mole ratio glucose/oxygen

# molecular weight
mw_g = 180 # g/mol glucose
mw_l = 282 # g/mol lipid

# microorganism OUR and yield properties
lmdour = 0.48 # (mol.L)/(m^3.h.g), ratio OUR/DCW
f_l = 0.5 # g lipid / g cell, fration DCW that is lipid

# starting conditions -- these should be provided as input by user
t0 = 12 # h
DCW0 = 10 # g/L

# end time should be determined by depletion of substrate in the for loop --
# not implemented yet, JJS 5/5/21
t = np.linspace(t0, 72, 10)
#t = np.array([t0, 40, 84])
#t = np.array([t0, 36, 60, 84])
#t = np.array([t0, 48, 84])
N = t.size
OUR = np.zeros(N) # molar oxygen uptake rate
DCW = np.zeros(N) # dry cell weight (density)
rho_l = np.zeros(N) # density of lipid

#time zero
DCW[0] = DCW0
rho_l[0] = f_l*DCW0
OUR[0] = fOUR*lmdour*DCW0

# later times
for i in range(1,N):
    r_l = OUR[i-1]*nr_l_o*mw_l/1000 # mass rate of lipid production
    rho_l[i] = rho_l[i-1] + r_l*(t[i] - t[i-1])
    DCW[i] = rho_l[i]/f_l
    OUR[i] = fOUR*lmdour*DCW[i]


figure(1)
clf()
plot(t, DCW, 'o-C0', label="DCW")
plot(texp, DCWexp, 'o--C0', mfc="none")
plot(t, rho_l, 's-C1', label="lipid")
plot(texp, lipexp, 's--C1', mfc="none")
plot(t, OUR, "d-C2", label="OUR")
xlabel("t")
ylabel("DCW, lipid [g/L]; OUR [mol/($m^3$.h)]")
legend(loc="best")
#savefig("fed-batch_compare.pdf", bbox_inches="tight")

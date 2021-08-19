"""
post-process results from two-phase lignocellulose EH model (presumed
generated from c++ implementation)

"""

# Jonathan Stickel, 2021


import numpy as np

#from matplotlib.pyplot import *
import matplotlib.pyplot as plt
plt.ion()

# hard-coded parameters
MwE = 65000.0
MwG = 162.0
Mwg = 180.0
MwX = 166.129
Mwx = 150.13
rhol = 1000.0
rhoT = rhol
rGg = MwG/Mwg
rXx = MwX/Mwx

# hard-coded intial values in the c++ file, and derived values
fis0 = 0.1
XG0 = 0.62
XX0 = 0.06
XL0 = 0.32
yF0 = 0.6 
lmbdE = 0.02
rhog0 = 4.3
rhox0 = 29.3
rhosL0 = 0

fG0 = XG0*fis0
fGF0 = yF0*fG0
fGR0 = (1-yF0)*fG0
fX0 = XX0*fis0
fL0 = XL0*fis0
fliq0 = 1 - fis0
fg0 = rhog0*fliq0/rhol
fx0 = rhox0*fliq0/rhol
fsL0 = rhosL0*fliq0/rhol
fET = lmbdE*fG0

### read in the results ###
data = np.loadtxt("eh_ligno_results.csv", skiprows=1, delimiter=",")
t = data[:,0]
fGF = data[:,1]
fGR = data[:,2]
fX = data[:,3]
fL = data[:,4]
fg = data[:,5]
fx = data[:,6]
fsL = data[:,7]


# derived results

fG0 = fGF0 + fGR0
fG = fGF + fGR
fis = fG + fX + fL
fliq = 1 - fis

GconvF = (fGF[0] - fGF)/fG[0]
GconvR = (fGR[0] - fGR)/fG[0]
Gconv = 1 - fG/fG[0]

Xconv = 1 - fX/fX0
Tconv = 1 - (fG + fX)/(fG0 + fX0)

rhog = rhol/fliq*fg
rhox = rhol/fliq*fx
rhosL = rhol/fliq*fsL
rGg = 162./180. # ratio of glucan/glucose molecular weights
rXx = 166.129/150.13 # ratio of xylan/xylose molecular weights

# lignin conversion -- just to check contribution to FIS reduction
Lconv = 1 - fL/fL0

# mass balance calculations
#mbE = 1 - fET/fET0

mG0 = fG[0] + fg[0]*rGg
mG = fG + fg*rGg
mbG = 1 - mG/mG0

mX0 = fX[0] + fx[0]*rXx
mX = fX + fx*rXx
mbX = 1 - mX/mX0

mL0 = fL[0] + fsL[0]
mL = fL + fsL
mbL = 1 - mL/mL0



plt.figure(1)
plt.clf()

plt.plot(t, Gconv, lw=2, label='total')
plt.plot(t, GconvF, label='facile')
plt.plot(t, GconvR, label='recalc')
plt.xlabel('t [h]')
plt.ylabel('glucan conversion')
plt.legend(loc='best')
plt.draw()
plt.show()

plt.figure(2)
plt.clf()
plt.plot(t, mbG, label='glucan')
plt.plot(t, mbX, label='xylan')
plt.plot(t, mbL, label='lignin')
plt.xlabel('t [h]')
plt.ylabel('mass balance')
plt.legend(loc='best')
plt.draw()
plt.show()

plt.figure(3)
plt.clf()
plt.plot(t, rhog, lw=2, label='glucose')
plt.plot(t, rhox, lw=2, label='xylose')
plt.plot(t, rhosL, lw=2, label='sol lignin')
plt.xlabel('t [h]')
plt.ylabel(r'$\rho_i$ [g/L]')
plt.legend(loc='best')
plt.draw()
plt.show()
#plt.savefig('fig.pdf', bbox_inches='tight')

plt.figure(4)
plt.clf()
plt.plot(Tconv, fis)
plt.xlabel('carbohydrate conversion')
plt.ylabel(r'$f_{is}$')
plt.draw()
plt.show()

plt.figure(5)
plt.clf()
plt.plot(t, Tconv, label='total carbohydrate')
plt.plot(t, Lconv, label='lignin')
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

"""
post-process results from two-phase cellulose-only EH model (presumed
generated from c++ implementation)

"""

# Jonathan Stickel, 2021


import numpy as np

from matplotlib.pyplot import *
ion()

# hard-coded parameters
MwE = 65000.0
MwG = 162.0
Mwg = 180.0
rhol = 1000.0
rhoT = rhol
rGg = MwG/Mwg

# hard-coded intial values in the c++ file, and derived values
lmbde = 0.03
fis0 = 0.05
dilution_factor = 0.5
xG0 = 1.0
xX0 = 0.0
rhog0 = 0.0
rhox0 = 1.0*dilution_factor
rhof0 = 0.0*dilution_factor

conversion_xylan = 0.5
yF0 = 0.2 + 0.6*conversion_xylan

fG0 = xG0*fis0
fGF0 = yF0*fG0
fGR0 = (1-yF0)*fG0
fX0 = xX0*fis0
fL0 = (1 - xG0 - xX0)*fis0
fliq0 = 1 - fis0
fg0 = rhog0*fliq0/rhol
fET = lmbde*fG0

### read in the results ###
data = np.loadtxt("eh_results.csv", skiprows=1, delimiter=",")
t = data[:,0]
fGF = data[:,1]
fGR = data[:,2]
fg = data[:,3]  # labeled `fg0`, but I'm presuming it is fg
rhog = data[:,4]

# derived results
fG = fGF + fGR
fis = fG + fX0 + fL0
fliq = 1 - fis
#epsl = rhoT/rhol*fliq
rhog2 = rhol/fliq*fg
print(np.allclose(rhog, rhog2))
      
convF = (fGF0 - fGF)/fG0
convR = (fGR0 - fGR)/fG0
conv = 1 - fG/fG0
mG0 = fG0 + fg0*rGg
mG = fG + fg*rGg
mb = 1 - mG/mG0


figure(1)
clf()
xlim((-5,50))
ylim((-0.02,1))
plot(t, conv, lw=2,label='total')
#plot(t, convF, lw=2,marker='*',markevery=10,markersize=10,label='facile')
plot(t, convF, '--', lw=2, label='facile')
#plot(t, convR, lw=2,marker='^',markevery=10,markersize=10,label='recalcitrant')
plot(t, convR, '-.', lw=2, label='recalcitrant')
#    plot(tdat, convdat, 'o', markersize=10,label='exp. data')
xlabel('time (hours)')
ylabel('conversion')
legend(loc='best')
#savefig("two-phase_conversion.pdf", bbox_inches='tight')


figure(2)
clf()
plot(t, mb)
xlabel('t [h]')
ylabel('mass balance')

figure(3)
clf()
plot(t, rhog, lw=2, label='model')
#    plot(tdat, rhosdat, 'o', label='data')
xlabel('t [h]')
ylabel(r'$\rho_g$ [g/L]')
legend(loc='best')

figure(4)
clf()
xlim((-5,50))
ylim((-0.005,0.07))
#plot(t, fGF+fGR, lw=2, marker='*',markevery=10,markersize=10,label='total glucan')
#plot(t, fGF, lw=2,marker='s',markevery=10,markersize=10,label='facile')
#plot(t, fGR, lw=2,marker='^',markevery=10,markersize=10,label='recalcitrant')
plot(t, fGF+fGR, '--', lw=2, label='total glucan')
plot(t, fGF, '-.', lw=2, label='facile')
plot(t, fGR, lw=2, ls=(0,(6,2)), label='recalcitrant')
#plot(t, fGR, lw=1, ls=(5,(5,2)), label='recalcitrant')
plot(t, fg, lw=2, label='glucose')
xlabel('time (hours)')
ylabel('mass fraction')
legend(loc='best')
#savefig("two-phase_species.pdf", bbox_inches='tight')

# figure(5)
# clf()
# plot(t, fGF/fGR, lw=2, label='model')
# xlabel('t [h]')
# #ylabel('?')
# legend(loc='best')

#show()

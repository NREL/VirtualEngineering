"""
Test two-phase reaction rate model for EH. This script applies it to a batch
reaction

"""

# Jonathan Stickel, 2017


#TODO:  put the model in a module!


# variable convention:
# Ci = c^tilde = mol i per total slurry volume (mol/L)
# ci = mol i per liquid volume (mol/L)
# rhoi = mass i per liquid volume (g/L)
# fi = mass i per total slurry mass (g/g)


import numpy as np
from scipy import integrate as igt
from scipy import optimize as opt

from matplotlib.pyplot import *
import matplotlib as mpl

import sys
from vebio.Utilities import dict_to_yaml, yaml_to_dict


if len(sys.argv) > 1:
    params_filename = sys.argv[1]
    ve_params = yaml_to_dict(params_filename)

else:
    ve_params = {} # won't this just cause errors? JJS 3/14/21


font={'family':'Helvetica', 'size':'15'}
mpl.rc('font',**font)
mpl.rc('xtick',labelsize=14)
mpl.rc('ytick',labelsize=14)


# parameters that are "fixed"
MwE = 65e3 # mol weight enzymes (CBH, but close to EG and BG)
# Jeoh, T., Michener, W., Himmel, M. E., Decker, S. R., & Adney, W. S. (2008). Implications of cellobiohydrolase glycosylation for use in biomass conversion. Biotechnology for Biofuels, 1, 10. https://doi.org/10.1186/1754-6834-1-10
MwG = 162. # mol weight glucan
Mwg = 180. # mol weight glucose
rhol = 1000. # kg/m^3 = g/L
rhoT = rhol # reasonable approximation -- most terms that need rhoT cancel
            # errors, so that this approximation provides near machine
            # precision accuracy
rGg = MwG/Mwg

# kinetics model parameters
#kF = 1000. # 1/h, trial-and-error
kF = 5e3 # best fit (constrained)
kR = 1e0*kF # 1/h
## Kdi and KI are "desorption" equilibrium constants -- lower values means more
## adsorbed enzyme
#KdF = 0.05 # mol/L = kmol/m^3 , trial-and-error
KdF = 0.03773 # kmol/m^3, best fit
#KdR = 10*KdF 
KdR = 0.5079 # kmol/m^3, best fit
#KI = 1e6 # kmol/m^3; high value means less inhibition -- 1e6 = no inhibition
#KI = 0.05 # kmol/m^3; trial-and-error
KI = 1e-2 # best fit (constrained)
# lignin solubilization parameters -- not used yet
#kL = 0.6

# reactor starting conditions
#mT0 = 10. # kg -- not used


#lmbde = 0.03 # kg/kg; *1000 to get mg/g
# Direct User Input
lmbde = ve_params['enzymatic_input']['lambda_e']

#fis0 = 0.05 # kg/kg; initial fraction insoluble solids
# Direct User Input (note, this is a target, not the output from pretreatment)
fis0 = ve_params['enzymatic_input']['fis_0']

# Compute the amount of dilution required to reach the fis_0_target
# based on the output from the pretreatment step
dilution_factor = fis0/ve_params['pretreatment_output']['fis_0']

#xG0 = 1.0 # initial glucan fraction of insoluble solids -- 100% here
xG0 = ve_params['pretreatment_output']['X_G']
xX0 = ve_params['pretreatment_output']['X_X']
rhog0 = 0.0 # g/L; initial glucose concentration in the liquid
#rhox0 = 0.0 # g/L
rhox0 = ve_params['pretreatment_output']['rho_x']*dilution_factor
rhof0 = ve_params['pretreatment_output']['rho_f']*dilution_factor # furfural

#yF0 = 0.4 # fraction of glucan that is "facile" -- best fit (constrained)
conversion_xylan = ve_params['pretreatment_output']['conv']
yF0 = 0.2 + 0.6*conversion_xylan
# yF0 = 1.0*conversion_xylan

print('\nINPUTS')
print('Lambda_e = %.4f' % (lmbde))
print('FIS_0 = %.4f' % (fis0))
print('yF0 = %.4f' % (yF0))


# initial conditions in model variables
fG0 = xG0*fis0
fGF0 = yF0*fG0
fGR0 = (1-yF0)*fG0
fX0 = xX0*fis0 # used only for fis calculations
fL0 = (1 - xG0 - xX0)*fis0 # ditto
fliq0 = 1 - fis0
#epsl0 = rhoT/rhol*fliq0 # not needed
fg0 = rhog0*fliq0/rhol
#fx0 = rhox0*fliq0/rhoT
#fsl0 = 0
fET = lmbde*fG0
f0 = np.array([fGF0, fGR0, fg0])


# enzymatic hydrolysis rate equation
def adsorption(CGF, CGRA, CET, cg, epsl, KdF, KdR, KI):
    """
    Adsorption equilibrium of enzyme between facile and recalcitrant glucan,
    and accounting for inhibition by glucose (and other sugars if present)

    """
    CEGF = CET/(1 + KdF/KdR*CGRA/CGF + epsl*KdF/CGF*(1 + cg/KI))
    CEGR = CET/(1 + KdR/KdF*CGF/CGRA + epsl*KdR/CGRA*(1 + cg/KI))
    return CEGF, CEGR


def rhs(f, t, rhoT, fET, kF, kR, KdF, KdR, KI):
    """ RHS of system of ODES"""
    # unpack
    fGF, fGR, fg = f

    # convenience term
    mwGrho = MwG/rhoT
    
    # molar concentrations
    CGF = fGF/mwGrho # rhoT/MwG*fGF
    CGR = fGR/mwGrho # rhoT/MwG*fGR
    CET = rhoT/MwE*fET
    
    # calculate liquid molar sugar concentration
    fG = fGF + fGR
    fis = fG + fX0 + fL0
    fliq = 1 - fis
    epsl = rhoT/rhol*fliq 
    cg = rhol/fliq/Mwg*fg

    # enzyme adsorption -- with CGRA = CGR (no limited exposure)
    CEGF, CEGR = adsorption(CGF, CGR, CET, cg, epsl, KdF, KdR, KI)
    rF = kF*CEGF
    rR = kR*CEGR

    # ODEs
    dfGF = -mwGrho*rF
    dfGR = -mwGrho*rR
    dfg = Mwg/rhoT*(rF + rR)
    
    return np.array([dfGF, dfGR, dfg])

# deleted model fitting -- would not do that with virtual engineering sims, JJS 3/23/20


# run simulation via igt.odeint
#tfin = 200. # h
# tfin = 100. # h
tfin = ve_params['enzymatic_input']['t_final']
print('t_final = %.4f' % (tfin))

N = 200
t = np.linspace(0, tfin, N)
f = igt.odeint(rhs, f0, t, args=(rhoT, fET, kF, kR, KdF, KdR, KI))

fGF = f[:,0]
fGR = f[:,1]
fg = f[:,2]

fG = fGF + fGR
fis = fG + fX0 + fL0
fliq = 1 - fis
#epsl = rhoT/rhol*fliq 
rhog = rhol/fliq*fg 

convF = (fGF0 - fGF)/fG0
convR = (fGR0 - fGR)/fG0
conv = 1 - fG/fG0
mG0 = fG0 + fg0*rGg
mG = fG + fg*rGg
mb = 1 - mG/mG0


# Save the outputs into a dictionary for use as inputs for bioreactor sims
output_dict = {'enzymatic_output': {}}
output_dict['enzymatic_output']['rho_g'] = float(rhog[-1])

# fixed dilution_EH factor, JJS 3/14/21
dilution_EH = fliq[0]/fliq[-1]
output_dict['enzymatic_output']['rho_x'] = float(rhox0*dilution_EH)
output_dict['enzymatic_output']['rho_f'] = float(rhof0*dilution_EH)

if len(sys.argv) > 1:
    dict_to_yaml([ve_params, output_dict], params_filename)

print('\nFINAL OUTPUTS (at t = %.1f hours)' % (tfin))
print('rho_g = %.4f' % (rhog[-1]))
# print('rho_x = ', rhox0*dilution_EH)
# print('rho_f = ', rhof0*dilution_EH)
# print('Facile Conversion = ', convF[-1])
# print('Recalcitrant Conversion = ', convR[-1])
print('Total Conversion = %.4f' % (conv[-1]))
print()


# save data to compare with other modeling approaches -- probably not needed for VE use
if False:
    #np.savez('two-phase_no_struc.npz', t=t, conv=conv, rhog=rhog)
    filename='two-phase_no_struc_'+str(fis0)+'.npz'
    np.savez(filename, t=t, conv=conv, rhog=rhog, fis0=fis0)
    
    outfile=open("wellmix_sim.dat","w")
    for i in range(len(conv)):
        outfile.write("%e\t%e\n"%(t[i],conv[i]))
    outfile.close()

if ve_params['enzymatic_input']['show_plots']:
    figure(1)
    clf()
    xlim((-5,50))
    ylim((-0.02,1))
    plot(t, conv, lw=3,label='total')
    #plot(t, convF, lw=3,marker='*',markevery=10,markersize=10,label='facile')
    plot(t, convF, '--', lw=3, label='facile')
    #plot(t, convR, lw=3,marker='^',markevery=10,markersize=10,label='recalcitrant')
    plot(t, convR, '-.', lw=3, label='recalcitrant')
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
    #plot(t, fGF+fGR, lw=3, marker='*',markevery=10,markersize=10,label='total glucan')
    #plot(t, fGF, lw=3,marker='s',markevery=10,markersize=10,label='facile')
    #plot(t, fGR, lw=3,marker='^',markevery=10,markersize=10,label='recalcitrant')
    plot(t, fGF+fGR, '--', lw=3, label='total glucan')
    plot(t, fGF, '-.', lw=3, label='facile')
    plot(t, fGR, lw=3, ls=(0,(6,2)), label='recalcitrant')
    #plot(t, fGR, lw=1, ls=(5,(5,2)), label='recalcitrant')
    plot(t, fg, lw=3, label='glucose')
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

    show()

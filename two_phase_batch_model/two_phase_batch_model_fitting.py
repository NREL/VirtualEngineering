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
import yaml


if len(sys.argv) > 1:
    input_filename = sys.argv[1]
    with open(input_filename) as fp:
        input_dict = yaml.load(fp, Loader = yaml.FullLoader)
    # print(input_dict)
else:
    input_dict = {}


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

lmbde = 0.03 # kg/kg; *1000 to get mg/g
lmbde = input_dict.get('lambda_e', lmbde)

fis0 = 0.05 # kg/kg; initial fraction insoluble solids
fis0 = input_dict.get('fis_0_target', fis0)
if input_dict.get('fis_0'):
    print('\nWARNING: Found a pretreatment output for initial fraction insoluble solids')
    print('Overriding widget value: %f' % (fis0))
    print('With pretreatment value: %f' % (input_dict['fis_0']))
    dilution = input_dict['dilution_strength']*fis0 + (1.0 - input_dict['dilution_strength'])*input_dict['fis_0']
    print('After dilutions: %f' % (dilution))
    fis0 = dilution
    print('\n')

xG0 = 1.0 # initial glucan fraction of insoluble solids -- 100% here
rhog0 = 0.0 # g/L; initial glucose concentration in the liquid
rhox0 = 0.0 # g/L
yF0 = 0.4 # fraction of glucan that is "facile" -- best fit (constrained)

conversion_xylan = (1.0 - input_dict['xf'])/input_dict['xi']
new_yF0 = 0.2 + 0.6*conversion_xylan
yF0 = new_yF0


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
    fis = fG
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


# read in experimental data
import xlrd
wb = xlrd.open_workbook('JimFermEnz_modelvalidate120208_JJS.xlsx')
ws = wb.sheet_by_name('Conversion Calcs')
strow = 30
enrow = 36
tdat = np.array(ws.col_values(1, strow, enrow))
rhogdat = np.array(ws.col_values(3, strow, enrow))
rhoxdat = np.array(ws.col_values(4, strow, enrow))
convgdat = np.array(ws.col_values(11, strow, enrow))
convdat = np.array(ws.col_values(13, strow, enrow))
rhosdat = rhogdat + rhoxdat


if False:
    # fit the model to data
    def residual(p, tdat, convdat, fis0, xG0, rhog0, lmbde):
        """
        residual function computing difference between predicted and data values
        """
        yF0 = p[0]
        kF = p[1] 
        KdF = p[2]
        KdR = p[3]
        KI = p[4]
        kR = kF # use same rate for facile and recalcitrant
        fG0 = xG0*fis0
        fGF0 = yF0*fG0
        fGR0 = (1-yF0)*fG0
        fliq0 = 1 - fis0
        fg0 = rhog0*fliq0/rhol
        fET = lmbde*fG0
        f0 = np.array([fGF0, fGR0, fg0])
        t = np.linspace(0, tdat[-1], 200)
        f = igt.odeint(rhs, f0, t, args=(rhoT, fET, kF, kR, KdF, KdR, KI))
        fGF = f[:,0]
        fGR = f[:,1]
        fg = f[:,2]
        fG = fGF + fGR
        fis = fG
        fliq = 1 - fis
        rhog = rhol/fliq*fg 
        convF = (fGF0 - fGF)/fG0
        convR = (fGR0 - fGR)/fG0
        conv = 1 - fG/fG0
        convi = np.interp(tdat, t, conv)
        return convi - convdat

    p0 = np.array([yF0, kF, KdF, KdR, KI])
    optresult = opt.least_squares(residual, p0,\
                                  bounds=([.4,0,1e-4,1e-4,1e-2],\
                                          [.6, 5e3, np.inf, np.inf, np.inf]),\
                                  args=(tdat, convdat, fis0, xG0, rhog0, lmbde))
    popt = optresult.x
    print(popt)
    yF0 = popt[0]
    kF = popt[1]
    KdF = popt[2]
    KdR = popt[3]
    KI = popt[4]
    kR = kF


# initial conditions in model variables
fG0 = xG0*fis0
fGF0 = yF0*fG0
fGR0 = (1-yF0)*fG0
#fL0 = (1 - xG0)*fis0
fliq0 = 1 - fis0
#epsl0 = rhoT/rhol*fliq0 # not needed
fg0 = rhog0*fliq0/rhol
#fx0 = rhox0*fliq0/rhoT
#fsl0 = 0
fET = lmbde*fG0
f0 = np.array([fGF0, fGR0, fg0])

# run simulation via igt.odeint
#tfin = 200. # h
tfin = 100. # h
tfin = input_dict.get('t_final', tfin)

N = 200
t = np.linspace(0, tfin, N)
f = igt.odeint(rhs, f0, t, args=(rhoT, fET, kF, kR, KdF, KdR, KI))

fGF = f[:,0]
fGR = f[:,1]
fg = f[:,2]

fG = fGF + fGR
fis = fG
fliq = 1 - fis
#epsl = rhoT/rhol*fliq 
rhog = rhol/fliq*fg 

convF = (fGF0 - fGF)/fG0
convR = (fGR0 - fGR)/fG0
conv = 1 - fG/fG0
mG0 = fG0 + fg0*rGg
mG = fG + fg*rGg
mb = 1 - mG/mG0

outfile=open("wellmix_expt.dat","w")
for i in range(len(convdat)):
    outfile.write("%e\t%e\n"%(tdat[i],convdat[i]))
outfile.close()

outfile=open("wellmix_sim.dat","w")
for i in range(len(conv)):
    outfile.write("%e\t%e\n"%(t[i],conv[i]))
outfile.close()


# save data to compare with other modeling approaches
if False:
    #np.savez('two-phase_no_struc.npz', t=t, conv=conv, rhog=rhog)
    filename='two-phase_no_struc_'+str(fis0)+'.npz'
    np.savez(filename, t=t, conv=conv, rhog=rhog, fis0=fis0)

if input_dict['show_plots']:
    figure(1)
    clf()
    xlim((-5,50))
    ylim((-0.02,1))
    plot(t, conv, lw=3,label='total')
    #plot(t, convF, lw=3,marker='*',markevery=10,markersize=10,label='facile')
    plot(t, convF, '--', lw=3, label='facile')
    #plot(t, convR, lw=3,marker='^',markevery=10,markersize=10,label='recalcitrant')
    plot(t, convR, '-.', lw=3, label='recalcitrant')
    plot(tdat, convdat, 'o', markersize=10,label='exp. data')
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
    plot(tdat, rhosdat, 'o', label='data')
    xlabel('t [h]')
    ylabel(r'$\rho_s$ [g/L]')
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
    
    figure(5)
    clf()
    plot(t, fGF/fGR, lw=2, label='model')
    xlabel('t [h]')
    ylabel(r'$\rho_s$ [g/L]')
    legend(loc='best')

    show()

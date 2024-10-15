import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammainc


# all in mol/m^3, except biomass, which is Kg/m^3=g/L
modelConstants = {}
modelConstants['qsmx'] = 17 #molS/m^3 / kgBio/m^3
modelConstants['bio_mx'] = 11 #kg/m^3

modelConstants['o2sat'] = 0.214 #mol/m^3
modelConstants['K_e'] = 0.0214 #mol/m^3
modelConstants['K_s'] = 31 #mol/m^3

modelConstants['Y_xs'] = 0.009 # g/molS
modelConstants['Y_as'] = 1.02 # molA/molS
modelConstants['Y_bs'] = 0.88 # molB/molS
#modelConstants['Y_os'] = 0.00467 # molO/molS
modelConstants['Y_os'] = 0.00467*10 # molO/molS

modelConstants['alpha_s'] = 3 # --
modelConstants['beta_s'] = 12 # --
modelConstants['alpha_e'] = 1 # --
modelConstants['beta_e'] = 1e3 # --



def dfdt(f, kLa, mc):
    """
    Definition of simplified metabolic model, as defined
    in the document. 
    
    Inputs:
    --f is the concentrations of the various species:
        o2, bio, glu, xyl, ace, bdo = f
        all in mol/m^3, except biomass, which is kg/m^3
    --otr is the oxygen transfer rate, in mM/h
    --mc is the set of model constants, in a dict
    Outputs:
    --dfdt is the model output
    """
    o2, bio, glu, xyl, ace, bdo = f

    F_s = (glu+xyl)/(glu+xyl+mc['K_s'])
    F_e = (o2+ace/mc['beta_e'])/(o2+ace/mc['beta_e']+mc['K_e'])
    
    qs = mc['qsmx']*F_s*F_e

    rbio = mc['Y_xs']*qs*bio*(1-bio/mc['bio_mx'])

    sRatio = np.max(np.array([glu/(xyl+1e-8), 0]))
    chi_s = gammainc(mc['alpha_s'], mc['beta_s']*sRatio)
    eRatio = np.max(np.array([o2/(ace+1e-8), 0]))
    chi_e = gammainc(mc['alpha_e'], mc['beta_e']*eRatio)
    #chi_e = 1
    chi_p = 0.3

    rg =     -chi_s *qs * bio
    rxy = -(1-chi_s)*qs * bio

    rar =    chi_p *mc['Y_as']*qs*bio
    rbr = (1-chi_p)*mc['Y_bs']*qs*bio

    otr = kLa*(mc['o2sat'] - o2)
    ro =     -chi_e *mc['Y_os']*qs*bio + otr
    rae = -(1-chi_e)*(mc['Y_as'])*qs*bio
    rbe = -rae

    ra = rar+rae
    rb = rbr+rbe

    return np.array([ro, rbio, rg, rxy, ra, rb])



if __name__=="__main__":
    f0 = np.array([0.214, 0.5, 500, 250, 0,0])
    kLa = 5

    #print(dfdt(f0, 50, modelConstants))

    dt = 0.01
    t0 = 0.
    tf = 24

    t = [t0]
    fmat = f0
    fold = f0
    while t[-1]<tf:
        fnew = fold + dt*dfdt(fold, kLa, modelConstants)

        fmat = np.append(fmat, fnew)
        fold = fnew
        t.append(t[-1]+dt)

    fmat = np.reshape(fmat, (-1,6))
    t = np.asarray(t)

    plt.figure(1)
    plt.clf()
    plt.plot(t, fmat[:,0])
    plt.xlabel('time (h)')
    plt.ylabel('oxygen (mM)')
    plt.show()

    plt.figure(2)
    plt.clf()
    plt.plot(t, fmat[:,1])
    plt.xlabel('time (h)')
    plt.ylabel('biomass (kg/$m^3$)')
    plt.show()


    plt.figure(3)
    plt.clf()
    plt.plot(t, fmat[:,2], label='glu')
    plt.plot(t, fmat[:,3], label='xyl')
    plt.plot(t, fmat[:,4], label='ace')
    plt.plot(t, fmat[:,5], label='bdo')
    plt.legend(loc='best')
    plt.xlabel('time (h)')
    plt.ylabel('substrates (mM)')
    plt.show()

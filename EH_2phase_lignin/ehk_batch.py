#############################################################################
#
# Module to define and run simulation of batch enzymatic hydrolysis, using a
# reaction kinetics model defined in a separate file (currently,
# enz_hyd_lignocell_two_phase.py).
#
# Jonathan Stickel and James Lischeske 2017
#
##############################################################################

import numpy as np
import pandas as pd
import scipy.integrate as igt

from enz_hyd_lignocell_two_phase import *

class BatchConditions():
    def __init__(self):
        #define default parameters
        self.fis = 0.1 #insoluble solids fraction
        self.XG0 = 0.6 #solids fraction of glucan
        self.XX0 = 0.05 #solids fraction of xylan
        self.XL0 = 0.2 #solids fraction of lignin
        self.yF0 = 0.3 #fraction of initial stuff in fascile portion
        #JL: should yF0 be in the conditions or parameters???
        #JL: I think this is fine. Child script for fitting now has access
        #to both conditions and reaction parameters
        self.lmbdE = 20./1000. #mg/g / 1000 = g/g
        self.rhog0 = 0. #g/L liquid-phase glucose
        self.rhox0 = 0. #g/L liquid-phase xylose
        self.rhosL0 = 0. #g/L soluble lignin


class Batch(Reaction):
    """
    Class containing the functions to perform batch simulation of enzymatic
    hydrolysis. This class inherits the reaction kinetics model from Reaction
    contained in enz_hyd_lignocell_two_phase.py

    """

    def __init__(self, **kwargs): #function inherited from Reaction
        super(Batch, self).__init__(**kwargs)
        self.Conditions = BatchConditions()
    
    
    def BatchRHS(self,f,t):
        return self.KineticRates(f)

    
    def BatchSimulate(self, t):
        """
        runs a simulation according to input parameters
        input:
           f0: type=np.array, containing fGF, fGR, fX, fL, fg, fx, fsL,
               fET at t=0
           t: type=np.array, containing times at which f will be returned
           p: type=array, containing [kR, kF, kX, kL, KdR, kappaRF, kappaRX,
              kappaRL, kappaRs, mGR1, mX1]; set to 'default' to reset to
              defaults; do not provide to use existing values
        output:
           f
           t
        """
        #if not hasattr(self, 'f0'):
        self.GetF0() #probably shouldn't use if statement; best to reset f0 based on current state
        
        #run ODE solver
        self.f = igt.odeint(self.BatchRHS, self.f0, t)
        self.t = t
        self.CalcResultValues()
    

    def GetF0(self):
        """
        Given fis, XG0, XX0, XL0, yF0, lmbdE, rhog0, rhox0, rhosL0, return the
        system in terms of mass fractions of the total slurry

        """

        # FIXME:  put in a check that XG0 + XX0 + XL0 =< 1
        #     (remember, ash/protien also exists in this context)

        if ((self.Conditions.XG0 + self.Conditions.XX0 + self.Conditions.XL0)>1.):
            raise RuntimeWarning('Sum of solids fractions is greater than 1.')
        
        fG = self.Conditions.XG0*self.Conditions.fis
        fGF = self.Conditions.yF0*fG
        fGR = (1-self.Conditions.yF0)*fG

        fX = self.Conditions.XX0*self.Conditions.fis
        fL = self.Conditions.XL0*self.Conditions.fis

        #### FIXME: hack introduced to force sum(X)=1 by increasing XL ####
        # This is important because fis is calculated in multiple points as the
        # sum of fX, fG and fL. Doing  this properly probably means introducing
        # a new fraction that is not available for solubilization. Probably
        # unnecessary
        self.Conditions.XL0 = 1 - self.Conditions.XX0 - self.Conditions.XG0
        fL = self.Conditions.XL0*self.Conditions.fis
        
        fliq = 1 - self.Conditions.fis
        fliqorhol = fliq/rhol #rhol is a global from reaction class
        fg = fliqorhol*self.Conditions.rhog0
        fx = fliqorhol*self.Conditions.rhox0
        fsL = fliqorhol*self.Conditions.rhosL0

        fET = self.Conditions.lmbdE*fG



        self.f0 = [fGF, fGR, fX, fL, fg, fx, fsL, fET]


    def CalcResultValues(self):
        """
        doc string TBD
        """

        # FIXME:  get rid of f0 -- just use f[0,j] instead, JJS 9/25/17
        # FIXME: all these conversions are in terms of mass rather than mole

        # unpack
        fGF0, fGR0, fX0, fL0, fg0, fx0, fsL0, fET0 = self.f0

        fGF = self.f[:,0]
        fGR = self.f[:,1]
        fX = self.f[:,2]
        fL = self.f[:,3]
        fg = self.f[:,4]
        fx = self.f[:,5]
        fsL = self.f[:,6]
        fET = self.f[:,7]


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
        rGg = 162./180. # ratio of glucan/glucose molecular weights
        rXx = 166.129/150.13 # ratio of xylan/xylose molecular weights

        GconvFIX = (rGg) * (1-fis) * (rhog - rhog[0])/(fis[0]*self.Conditions.XG0*rhol)
        XconvFIX = (rXx) * (1-fis) * (rhox - rhox[0])/(fis[0]*self.Conditions.XX0*rhol)
        TconvFIX = (GconvFIX*self.Conditions.XG0 + XconvFIX*self.Conditions.XX0) / \
                   (self.Conditions.XG0 + self.Conditions.XX0)

        # lignin conversion -- just to check contribution to FIS reduction
        Lconv = 1 - fL/fL0
        
        # mass balance calculations
        mbE = 1 - fET/fET0

        mG0 = fG[0] + fg[0]*rGg
        mG = fG + fg*rGg
        mbG = 1 - mG/mG0

        mX0 = fX[0] + fx[0]*rXx
        mX = fX + fx*rXx
        mbX = 1 - mX/mX0

        mL0 = fL[0] + fsL[0]
        mL = fL + fsL
        mbL = 1 - mL/mL0

        retDict = pd.DataFrame()
        retDict['t'] = self.t

        retDict['fGF'] = fGF
        retDict['fGR'] = fGR
        retDict['fX'] = fX
        retDict['fL'] = fL
        retDict['fg'] = fg
        retDict['fx'] = fx
        retDict['fsL'] = fsL
        retDict['fET'] = fET
        
        retDict['fis'] = fis
        retDict['Gconv'] = Gconv
        retDict['GconvF'] = GconvF
        retDict['GconvR'] = GconvR
        retDict['Xconv'] = Xconv
        retDict['Tconv'] = Tconv
        retDict['Lconv'] = Lconv

        retDict['GconvFIX'] = GconvFIX
        retDict['XconvFIX'] = XconvFIX
        retDict['TconvFIX'] = TconvFIX
        
        retDict['rhog'] = rhol/fliq*fg
        retDict['rhox'] = rhol/fliq*fx
        retDict['rhosL'] = rhol/fliq*fsL

        retDict['mbE'] = mbE
        retDict['mbG'] = mbG
        retDict['mbX'] = mbX
        retDict['mbL'] = mbL

        self.SimOut = retDict

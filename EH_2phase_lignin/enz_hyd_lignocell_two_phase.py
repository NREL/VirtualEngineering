#############################################################################
#
# Module to define two-phase EH model, containing xylose and lignin. This
# follows development and documentation by JJS.
#
# James Lischeske and Jonathan Stickel 2017
#
##############################################################################

import numpy as np

## variable convention:
# Ci = c^tilde = mol i per total slurry volume (mol/L)
# ci = mol i per liquid volume (mol/L)
# rhoi = mass i per liquid volume (g/L)
# fi = mass i per total slurry mass (g/g)

## parameters that are "fixed" ##
# mol weight enzymes (CBH, but close to EG and BG) Jeoh, T., Michener, W.,
# Himmel, M. E., Decker, S. R., & Adney, W. S. (2008). Implications of
# cellobiohydrolase glycosylation for use in biomass conversion.
# Biotechnology for Biofuels, 1, 10. https://doi.org/10.1186/1754-6834-1-10
MwE = 65e3 
MwG = 162. # mol weight glucan
Mwg = 180. # mol weight glucose
MwX = 166.129 # mol weight xylan
Mwx = 150.13 # mol weight xylose
MwL = 200. # approx mol weight of lignin
# This turns out to be pretty complicated. "Lignin" is more a type than a
# specific chemical, so to define a single molecular weight is simply not
# valid. Yet, for the sake of simplicity, we treat lignin as a single
# species with a single K_I, and thus we need a single MW.  200 Da
# represents a MW within the range reported by Chua and Wayman 1979
# (reported 188-211 Da)

rhol = 1000. # kg/m^3 = g/L
rhoT = rhol # reasonable approximation -- most terms that need rhoT cancel
            # errors, so that this approximation provides near machine
            # precision accuracy
rGg = MwG/Mwg
rXx = MwX/Mwx


# number of species in this model -- can be used as a check for use in programs
# that import this module
nspecies = 8


class KinParams(object):
    """
    Class for EH kinetic model parameters (coefficients). Parameters
    are initialized to default values. Set individual parameters different from
    the defaults by providing them as keyword arguments when intializing.  The
    parameters are:

    kF
    kR
    kX
    kL
    KdR
    kappaRF
    kappaRX
    kappaRL
    kappaRs
    mGR1
    mX1

    """
    def __init__(self, **kwargs):
        # model parameters/coefficients initialized to default parameters
        # kinetics model parameters
        self.kR = 5e3 # 1/h; rate for recalcitrant glucan
        self.kF = 5e3 # 1/h; rate for facile glucan
        self.kX = 5e3 # xylan
        self.kL = 5e2 # L/mol = m^3/kmol; lignin solubilization rate coefficient
        
        self.KdR = 0.5 # kmol/m^3
        self.kappaRF = 10.
        self.kappaRX = 5.
        self.kappaRL = 10.
        self.kappaRs = 50.
        
        # accessibility model parameters; mj1 = 1 means all substrate is
        # accessible
        self.mGR1 = 2.
        self.mX1 = 2.

        #print(self.__dict__.keys())
        allowed_keys = self.__dict__.keys()
        
        for key, value in kwargs.items():
            if key in allowed_keys:
                setattr(self, key, value)
            else:
                print("warning:  the provided parameter '%s' is not recognized "
                      "and has been ignored" % key)
                

class Reaction(KinParams):
    """
    Class containing the functions used to evaluate instantaneous enzyme
    adsorption and reaction rates. It inherits the default model parameters
    from the Params class. Set individual parameters different from the
    defaults by providing them as keyword arguments when intializing OR by
    providing a Params object. If a Params object is provided, any provided
    keyword arguments will be ignored. See the docstring of Params for a list
    of the parameter keywords.

    """

    def __init__(self, ehkp=None, **kwargs):
        if ehkp is not None: # use the provided parameters object
            if len(kwargs) is not 0:
                print("warning:  the provided keyword arguments for parameters "
                      "are ignored in favor of using the provided Params "
                      "object")
            for key, value in ehkp.__dict__.items():
                setattr(self, key, value)
        else: # set defaults and process any keyword arguments
            KinParams.__init__(self, **kwargs)

            
    def Accessibility(self, XRj, mj1):
        """
        Calculate the accessibility of recalcitrant structural carbohydrate.
        """
        return 0.5*(3 - mj1)*XRj + 0.5*(mj1 - 1)*XRj**3

    
    def DynCalcValues(self, f):
        """
        Gets dynamically calculated values of lambda and D for a given f-vec
        """
        
        #lamFR = CGF/CGRA, lamXR = CXA/CGRA
        fGF, fGR, fX, fL, fg, fx, fsL, fET = f

        #Calculate Accessibility
        fRis = fGR + fX + fL
        X_RGR = fGR/fRis
        X_RX = fX/fRis
        X_RGRA = self.Accessibility(X_RGR, self.mGR1)
        X_RXA = self.Accessibility(X_RX, self.mX1)

        CGF = (rhoT/MwG) * fGF
        CGRA = (rhoT/MwG) * fRis * X_RGRA
        CXA = (rhoT/MwX) * fRis * X_RXA

        lamFR = CGF/CGRA
        lamXR = CXA/CGRA


        #calculate denominator
        fis = fRis + fGF
        eps_l = rhoT/rhol * (1 - fis)
        csL = rhoT/MwL * fsL
        css = rhoT/Mwg*fg + rhoT/Mwx*fx
        
        denom = 1 + self.kappaRF *lamFR + self.kappaRX *lamXR + eps_l/CGRA * \
                (self.KdR + self.kappaRL*csL + self.kappaRs*css)
        
        return lamFR, lamXR, denom


    def Adsorption(self, f):
        """
        Determine the equilibrium adsorption of the enzymes
        """
        # unpack f
        fGF, fGR, fX, fL, fg, fx, fsL, fET = f

        # Calculate sets of common constants Lambda and Denom
        lamFR, lamXR, denom = self.DynCalcValues(f)

        CET = (rhoT/MwE) * fET
        
        # now calculate adsorbed enzyme concentrations
        CEGR = CET/denom
        CEGF = self.kappaRF * lamFR * CET/denom
        CEX = self.kappaRX * lamXR * CET/denom
        return CEGR, CEGF, CEX
    

    def KineticRates(self, f, CEA=None):
        """
        takes in parameters and current mass-fraction values, spits out rates;
        values for adsorbed enzyme concentrations may be optionally passed in
        as a 3-element array

        """
        
        # unpack f
        fGF, fGR, fX, fL, fg, fx, fsL, fET = f
        
        if CEA is not None:
            CEGR, CEGF, CEX = CEA
        else:
            CEGR, CEGF, CEX = self.Adsorption(f)

        #CEGR, CEGF, CEX = self.Adsorption(f)
        
        # moles of lignin
        CL = (rhoT/MwE) * fL

        # calculate molar rates (rt = "r-tilda")
        rtGF = self.kF*CEGF
        rtGR = self.kR*CEGR
        rtX = self.kX*CEX
        rtL = self.kL*CL*(rtGR + rtX)
        
        # calculate rates
        R_GF = - MwG/rhoT * rtGF
        R_GR = - MwG/rhoT * rtGR
        R_X  = - MwX/rhoT * rtX
        R_g = Mwg/rhoT * (rtGF + rtGR)
        R_x = Mwx/rhoT * rtX
        R_L = - MwL/rhoT * rtL
        R_sL = - R_L
        R_ET = 0. # reaction term for enzymes is zero
        
        # pack and return rates
        return np.array([R_GF, R_GR, R_X, R_L, R_g, R_x, R_sL, R_ET])

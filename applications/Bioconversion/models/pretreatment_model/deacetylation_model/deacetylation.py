# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 14:33:24 2025

@author: spittman
"""

# Set up Excel Relation and Data Frame
# knowns should be in excel file at this point
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
import math
import matplotlib.pyplot as plt

###########################################
# PHYSICAL AND CHEMICAL CONSTANTS
###########################################
###
# Equation Constants
F = 96485 # Faraday Constant -Amps*Sec/mol- -C/mol-
R = 0.008314 # Ideal gas constant -8.314 J/mol*K- -kJ/mol-
R_gas = 0.08206 # -L*atm/mol*K-
rho_h2o = 1000 # density of water -g/L-
mrho_h2o = 55.56 # density of water -mol/m3-
Ea = 38 # Arrhenius activation energy for acetyl groups -J/mol-
Arc = 1.1e4 # Arrhenius pre-exponent factor  for acetyl groups -L/(mol-s)-
Ea_L = 55 # Arrhenius activation energy for lignin -J/mol-
Arc_L = 1.5e6 # Arrhenius pre-exponent factor for lignin -L/(mol-s)-
Ea_X = 45 # Arrhenius activation energy for xylan -J/mol-
Arc_X = 7e3 # Arrhenius pre-exponent factor for xylan -L/(mol-s)-
###
# Stoichiometric Equilibrium Constants
pKa_ace = 4.76
pK_Naace = 4.76

###
# Molar Mass - g/mol - 
M_NaOH = 40  
M_OH = 17.01
M_Na = 22.99
M_ace = 59.04 # acetate - C2H3O2
M_hace = 60.05 # acetic acid - C2H4O2
M_Naace = 82.03 # sodium acetate
M_Cel = 162.14 # cellulose - C6H10O5
M_Hem = 174.15 # hemicellulose - modeled as acetylated xylan
M_acy = 43.045 # acetyl group - C2H3O
M_X = 132.11 # deacetylated xylan
M_Lig = 188 # lignin - common approximation
M_H = 1

# Create matrix -g/mol-
loops = (3,12)
const = np.zeros(loops)
const[0,0] = M_ace
const[0,1] = M_hace
const[0,2] = M_Naace
const[0,3] = M_Cel
const[0,4] = M_Hem
const[0,5] = M_Lig
const[0,6] = M_acy
const[0,7] = M_X
const[0,8] = M_NaOH
const[0,9] = M_Na
const[0,10] = M_OH
const[0,11] = M_H


#########################################
# INITIAL CONDITIONS
#########################################
###
# VARIABLE SET INITIAL CONDITIONS 
t_final = 7200 # Total reaction time -Sec- -120 min- #7200 #5400
Tf = 363.15 # Operating Fluid Temperature -K- -363.15|90 C- -323.15|50 deg C-
V_DA = 20 # Volume of Deacetylation Reactor -L- #30 #20
pHi_DA = 7 # Inital pH
pOHi_DA = 14-pHi_DA # Initial pOH
Part = 0.75 # particle size -inch-
Decomp = 2 # storage time -unitless- -(months)-

###
# Initial Influent Content -g/L-
C_Husk = 10 # corn husk
C_Cob = 10 # corn cob
C_Stalk = 10 # corn stalk
C_NaOH = 2.01 # sodium hydoxide # 4.55 6.36 9.09 # 9.47 5.35 # 2.61 2.01 1.41
# Biomass
W_BM = 600 # biomass -g- #600 #5000

# Acetyl Groups
Acy_P = 0.023 #0.023 #0.0223 # weight percentage of acetate in biomass
# Lignin Groups
Lig_P = 0.153 # weight percentage of lignin in biomass
# Xylan Groups
Xy_P = 0.213 # corn stover


###############################################
### DEFINE FUNCTION
###############################################
###############################################
### INTEGRATE
def deacetylate(ve_params,verbose=True,show_plots=True):
    print("doing deacetylation..")

    # Read VE Inputs
    Tf=ve_params.pt_in['deacetylation_temperature']
    Acy_P = ve_params.pt_in['acetylfrac']
    fis_in = float(ve_params.pt_in['initial_solid_fraction'])
    fis_dac = float(ve_params.pt_in.get('deacetylation_solid_fraction', fis_in))
    t_final = ve_params.pt_in['deacetylation_retention_time']
    V_DA = ve_params.pt_in['deacetylation_volume']
    pHi_DA = ve_params.pt_in['deacetylation inital pH']
    W_BM = ve_params.pt_in['deacetylation_biomass_loading']
    C_NaOH = float(ve_params.pt_in['deacetylation_NaOH_loading'])

    # Acetyl Groups
    C_acy = W_BM*Acy_P/V_DA # acetyl groups -g/L-2
    Molar_acy = C_acy/M_acy
    C_ace_max = C_acy*(M_ace/M_acy)
    # Lignin Groups
    C_lig = W_BM*Lig_P/V_DA # -g/L-
    C_lig_max = C_lig
    Molar_lig = C_lig/M_Lig
    # Xylan Groups
    C_xy = W_BM*Xy_P/V_DA # acetyl groups -g/L-2
    C_xyo_max = C_xy*(M_X/M_Hem)
    Molar_xy = C_xy/M_X

    #######################################
    # Set UP Matrices
    #######################################
    # Calculation Matrix -row 0 g/L- -row 1 mol/L-
    loops=(1,12)
    DAC = np.zeros(loops)
    M_DAC = np.zeros(loops)
    dM_DAC = np.zeros(loops)
    DAC[0,0] = 0                                               # Acetate
    DAC[0,1] = 0                                               # Acetic Acid
    DAC[0,2] = 0                                               # Sodium Acetate
    DAC[0,3] = 0                                               # Cellulose
    DAC[0,4] = C_xy                                            # Hemicellulose
    DAC[0,5] = C_lig                                           # Lignin
    DAC[0,6] = C_acy                                           # Acetyl Groups
    DAC[0,7] = 0                                               # Deacetylated Xylan
    DAC[0,8] = C_NaOH                                          # Sodium Hydroxide
    DAC[0,9] = 0                                               # Sodium
    DAC[0,10] = 10**(-pOHi_DA)*M_OH                            # Hydroxide
    DAC[0,11] = 10**(-pHi_DA)*M_H                              # Hydrogen
    # Convert
    M_DAC[0,:] = DAC[0,:]/const[0,:] # -mol/L-
    # NaOH Dissociates Completely
    M_DAC[0,9] = M_DAC[0,9] + (1*M_DAC[0,8]) # Sodium
    M_DAC[0,10] = M_DAC[0,10] + (1*M_DAC[0,8]) # Hydroxide
    M_DAC[0,8] =M_DAC[0,8] - M_DAC[0,8] # Sodium Hydroxide
    # Convert
    DAC[0,:] = M_DAC[0,:]*const[0,:] # -g/L-

    # Free Acids
    M_DAC[0,10] = M_DAC[0,10] - 0.02

    # Arrehnius Rate Kinetics
    # Lignin
    kT_arh_L = Arc_L*np.exp(-Ea_L/(R*Tf)) # -L/(mol-s)-
    # Acetate
    kT_arh = Arc*np.exp(-Ea/(R*Tf)) # -L/(mol-s)-
    # Xylan
    kT_arh_X = Arc_X*np.exp(-Ea_X/(R*Tf)) # -L/(mol-s)-

    # Initial Conditions
    Cacetyl = M_DAC[0,6]
    COH = M_DAC[0,10]
    Cace = M_DAC[0,0]
    Cxy = M_DAC[0,4]
    Cxyo = 0
    Cxyoh = 0
    Clig = M_DAC[0,5]
    Cslig = 0
    t_span = (0,t_final)
    y0 = [Cacetyl, COH, Cace, Clig, Cslig, Cxy, Cxyo]

    def DACy(t, y):
        # Arrehnius Rate Kinetics
        ################################
        # Lignin
        kT_arh_L = Arc_L*np.exp(-Ea_L/(R*Tf)) # -L/(mol-s)-
        # Acetate
        kT_arh = Arc*np.exp(-Ea/(R*Tf)) # -L/(mol-s)-
        # Xylan
        kT_arh_X = Arc_X*np.exp(-Ea_X/(R*Tf)) # -L/(mol-s)-

        # Define Rate Equation
        ################################
        # Lignin
        dligdt = kT_arh_L*y[3]*y[1] # -mol/L-s-
        # Change value
        dClig_dt = -dligdt
        dCslig_dt = +dligdt
        # Acetate
        dacedt = kT_arh*y[0]*y[1] # -mol/L-s-
        # Change value
        dCacy_dt = -dacedt
        dCace_dt = +dacedt
        # Xylan
        dxydt = kT_arh_X*y[5]*y[1] # -mol/L-s-
        dCxyo_dt = +dxydt
        dCxy_dt = -dxydt
        # Final Hydroxide Stoichiometry
        dCOH_dt  =-(dacedt+(1.2*dligdt)+(dxydt))
        return [dCacy_dt, dCOH_dt, dCace_dt, dClig_dt, dCslig_dt, dCxy_dt, dCxyo_dt]

    # Integrate
    sol_DA = solve_ivp(DACy,t_span, y0, method='BDF', rtol=1e-6, atol=1e-9)

    # Retreive Final Value
    Cacetyl_f, COH_f, Cace_f, Clig_f, Cslig_f, Cxy_f, Cxyo_f = [float(sol_DA.y[i, -1]) for i in range(7)]

    # pH Calculation
    #############################################
    # Calculate pKw
    Hnot0 = 55.83
    Kw2 = 10**(-14)*((Hnot0/R)*((1/289.15)-(1/Tf)))
    pKw = -np.log10(Kw2)
    OH_final = sol_DA.y[1, -1] * M_OH
    pH = pKw+np.log10(OH_final)

    # Yield Calculations
    #################################################
    Cace_f_gL = Cace_f * M_ace
    ace_yield = Cace_f_gL / C_ace_max # acetate yield
    Clig_f_gL = Clig_f * M_Lig
    lig_yield = Clig_f_gL / C_lig_max # total lignin yield
    Cxy_f_gL = Cxy_f * M_X
    xyl_yield = Cxy_f_gL / C_xyo_max # total xylose yield

    # Dewatering Calculations
    #################################################
    # Avoid divide-by-zero / nonsense
    if fis_dac <= 0 or fis_dac >= 1 or fis_in <= 0 or fis_in >= 1:
        s = 1.0
    else:
        s = ((1.0 - fis_dac) / (1.0 - fis_in)) * (fis_in / fis_dac)
    # Liquid concentrations increase when liquid decreases
    Cxy_f_dw  = Cxy_f  / s
    Cxyo_f_dw = Cxyo_f / s
    Cace_f_dw = Cace_f / s


    if not hasattr(ve_params, "pt_out") or ve_params.pt_out is None:
        ve_params.pt_out = {}

        # standard VE pretreatment outputs expected downstream
        ve_params.pt_out["fis_0"] = fis_dac
        ve_params.pt_out["X_X"]   = float(ve_params.feedstock["xylan_solid_fraction"])
        ve_params.pt_out["X_G"]   = float(ve_params.feedstock["glucan_solid_fraction"])
        ve_params.pt_out["conv"]  = 0.0
        # Furfural
        ve_params.pt_out["rho_f"] = 0.0
        # pH
        #ve_params.pt_out["pH_final"] = float(pH_final)
        # Theoretical Maximum Yield
        ve_params.pt_out["C_ace_max_gL"] = float(C_ace_max)
        ve_params.pt_out["C_lig_max_gL"] = float(C_lig_max)
        ve_params.pt_out["C_xyo_max_gL"] = float(C_xyo_max)
        # Actual Yield 
        ve_params.pt_out["ace_yield"] = float(ace_yield)
        ve_params.pt_out["lig_yield"] = float(lig_yield)
        ve_params.pt_out["xyl_yield"] = float(xyl_yield)
        # Final Concentration - before dewatering
        ve_params.pt_out["Cace_final_gL_pre_dewater"] = float(Cace_f_gL)
        ve_params.pt_out["Clig_final_gL_pre_dewater"] = float(Clig_f_gL)
        ve_params.pt_out["Cxy_final_gL_pre_dewater"]  = float(Cxy_f_gL)
        # Final Concentration - after dewatering
        ve_params.pt_out["Cxy_final_dewatered"]  = float(Cxy_f_dw)
        ve_params.pt_out["Cxyo_final_dewatered"] = float(Cxyo_f_dw)
        ve_params.pt_out["Cace_final_dewatered"] = float(Cace_f_dw)
        ve_params.pt_out['rho_x']=float(Cxy_f_dw)
    
    return ve_params.pt_out

if __name__ == "__main__":
    class DummyVE:
        def __init__(self):
            self.pt_in = {
                "deacetylation_temperature": 350.0,
                "acetylfrac": 0.5,
                "initial_solid_fraction": 0.2,
                "deacetylation_solid_fraction": 0.3,
            }
            self.feedstock = {
                "xylan_solid_fraction": 0.25,
                "glucan_solid_fraction": 0.35,
            }
            self.pt_out = {}

    ve = DummyVE()
    deacetylate(ve, verbose=True, show_plots=False)
    print("pt_out keys:", ve.pt_out.keys())
    for k, v in ve.pt_out.items():
        print(k, v)

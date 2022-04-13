"""
This script will be used to solve for CEH steady-state outcomes given
specified reactor *design* parameters.

"""

# Jonathan Stickel, 2017, 2022

"""
Based on ceh_stead_multi_design_solve.py (03/16/2022).
Only cstr reactors (no batch reactors prior to ceh cstrs).
Add VE parameters to allow user input from notebook interface.
---- YL 03/16/2022
"""

import numpy as np

import sys
sys.path.append("../src/core/")

import ceh_cstr_model_multi as ceh

from vebio.Utilities import dict_to_yaml, yaml_to_dict

if len(sys.argv) > 1:
    params_filename = sys.argv[1]
    ve_params = yaml_to_dict(params_filename)
else:
    raise Exception("VE parameters filename not provided")

# 8 species, in this order:
# GF, GR, X, L, g, x, sL, ET
#  0,  1, 2, 3, 4, 5,  6,  7

# 4 flow streams per reactor, in this order (UF recycle stream equal to 0 if not using)
# feed, makeup buffer, permeate, purge, UF recycle
#    1,             2,        3,     4,          5  (human accounting, e.g., `F1`)
#    0,             1,        2,     3,          4  (numpy array accounting, e.g., `F[0]`)

## TODO: 03162022, hardcode this to DMR for now. Can be changed to user input.
if False: # biomass feed, *** DDA ***
        # from Lischeske and Stickel (2017) Biotechnol Biofuels 12:299
    kinParams = {'KdR': 0.05,
                 'kL': 729.5,
                 'kR': 14713,
                 'kX': 1e4,
                 'kappaRF': 9.34,
                 'kappaRL': 50,
                 'kappaRX': 11.3,
                 'kappaRs': 50}
    kinParams['kF'] = kinParams['kR']
    f1is = 0.075
    F1 = 100.0  # kg/h -- simple basis, scale later
    #F1 = 326785*.129/f1is  # based on aspen flowchart out of pretreatment:
                           # insol solids flow in kg/h
    x1G = 0.621 # fraction of solids that is glucan, from aspen flowchart out
                # of pretreatment
    y1F = 0.6 # fraction of glucan that is facile
    x1X = 0.054 # fraction of solids that is xylan, from aspen flowchart out of
                # pretreatment
    rho1g = 11.9 * f1is/0.129  # based on aspen flowchart out of pretreatment:
                               # glucose conc in g/L
    rho1x = 49.9 * f1is/0.129  # based on aspen flowchart out of pretreatment:
                               # xylose conc in g/L
    rho1sL = 10. # g/L; initial soluble lignin -- is this measured by analytical?

else: # biomass feed, *** DMR ***
    # from 2021-03 CEHD milestone report table 2, for DMR material
    kinParams = {'KdR': 0.5,
                 'kL': 1467,
                 'kR': 3843,
                 'kX': 1e4,
                 'kappaRF': 2,
                 'kappaRL': 50,
                 'kappaRX': 0.84,
                 'kappaRs': 1.52}
    # adjust params as needed
    kinParams['kF'] = kinParams['kR']  # add to dict (kF=kR for now)
    # kinParams['kX'] = kinParams['kR']  # the fit kX is way too high
    # kinParams['kappaRL'] = 10  # the fit kappaRL is too high

    #f1is = 0.100
    f1is = ve_params['CEH_input']['f1_is'] # inflow stream insoluble solids
    F1 = ve_params['CEH_input']['inflow_mass_flowrate'] # kg/h -- simple basis, scale later
    #F1 = 282628*.198/f1is  # based on aspen flowchart out of deacetylation:
                           # insol solids flow in kg/h
    x1G = ve_params['CEH_input']['glucan_solid_fraction'] # fraction of solids that is glucan, from aspen flowchart out
                # of deacetylation
    y1F = ve_params['CEH_input']['facile_fraction_glucan'] # fraction of glucan that is facile
    x1X = ve_params['CEH_input']['xylan_solid_fraction'] # fraction of solids that is xylan, from aspen flowchart out of deacetylation
    rho1g = 0 * f1is/0.198  # based on aspen flowchart out of deacetylation:
                            # glucose conc in g/L
    rho1x = 0 * f1is/0.198  # based on aspen flowchart out of deacetylation:
                            # xylose conc in g/L
    rho1sL = ve_params['CEH_input']['residue_soluble_lignin'] # g/L; initial soluble lignin -- is this measured by analytical?

## TODO: This is being hardcoded to 3 for now. 03162022
# number reactors
nr = 3  # (this value is updated by the solver later)

# target enzyme loading -- enzyme is only fed to the first reactor
lmbdE = ve_params['CEH_input']['lambda_e'] # kg/kg; *1000 to get mg/g

# target reactor fis, each reactor
#fis = np.array([0.05, 0.05, 0.05])
# length must match nr
fis = np.array([
                ve_params['CEH_input']['fis_1'],
                ve_params['CEH_input']['fis_2'],
                ve_params['CEH_input']['fis_3']
                ])

# target *carbohydrate* conversion yield *for each reactor* -- future work:
# work out a way to specify an overall conversion yield to be balanced between
# the reactors
# xi = 0.75
# xi = np.array(nr*[xi])
# different target conversions for each reactor; resulting relative sizes of
# reactors were a bit unintuitive (to me) -- it is a balance of residence times
# and amounts of solids entering each reactor
xi = np.array([
                ve_params['CEH_input']['target_conv_1'],
                ve_params['CEH_input']['target_conv_2'],
                ve_params['CEH_input']['target_conv_3'],
                ]) # length must match nr

# outlet flow criteria of each reactor
# use "theta2" now!
theta2set = 0.50
theta2 = np.array([
                ve_params['CEH_input']['theta2_1'],
                ve_params['CEH_input']['theta2_2'],
                ve_params['CEH_input']['theta2_3']
                ])# length must match nr

# enzyme rejection by membrane -- assume common for all CSTRs
# rejE = 0.5
rejE = 1.  # all enzymes stay with solids for practical purposes (FY19 Q2 CEHD experiments)

# initiate solver object, setting some kinetics parameters
#solver = ceh.CEHSolver() # using builtin default parameter values
#solver = ceh.CEHSolver(kappaRs=1., kR=1e4, kF=1e4, kL=2.5e3, kappaRF=2.) # manually specify
solver = ceh.CEHSolver(**kinParams) # use set provided above
# solve the system
cstr = solver.designsolve(F1, f1is, x1G, y1F, x1X, rho1g, rho1x, rho1sL,
                          lmbdE, fis, xi, theta2, rejE)
ceh.post(cstr)

nr = cstr.nr

print(f'++++ Pump Power Requirement for membrane loops of {nr} reactors ++++\n')
from membrane_loop_system import MembraneLoopSystem
p_m_dot = -cstr.F[:,2] ## mass flow rate of permeate kg/h
fis = cstr.fis ####
# assume density of sugar stream is equal to water?
p_rou = 1 #kg/L

membrane_loop_systems = [[]]*nr
power_consumptions = np.zeros(nr)
membrane_units = np.zeros(nr).astype(int)
for i in range(nr):
    membrane_loop_systems[i] = MembraneLoopSystem(p_m_dot[i], p_rou, fis[i])
    power_consumptions[i] = membrane_loop_systems[i].pump_power
    membrane_units[i] = membrane_loop_systems[i].membrane_units
    print('Reactor {} needs {} membrane units with pump power {:0.2f} kW.'.format(
        i+1,membrane_units[i],power_consumptions[i]))

print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

strg = "\nreactor sizes are" + (nr-1)*" %g," + " %g kg"
print(strg % tuple(cstr.mT))
strg = "Feed streams are" + (nr-1)*" %g," + " %g kg/h"
print(strg % tuple(cstr.F[:,0]))
strg = "makeup buffer streams are" + (nr-1)*" %g," + " %g kg/h"
print(strg % tuple(cstr.F[:,1]))
strg = "permeate streams are" + (nr-1)*" %g," + " %g kg/h"
print(strg % tuple(cstr.F[:,2]))
strg = "exit/purge streams are" + (nr-1)*" %g," + " %g kg/h"
print(strg % tuple(cstr.F[:,3]))
strg = "FIS in each reactor is" + (nr-1)*" %g," + " %g"
print(strg % tuple(cstr.fis))
strg = "Fraction glucan in solids:" + (nr-1)*" %g," + " %g"
print(strg % tuple(cstr.xG))
strg = "permeate glucose concentrations are" + (nr-1)*" %g," + " %g g/L"
print(strg % tuple(cstr.rho3g))
strg = "residence times in each reactor are" + (nr-1)*" %g," + " %g h"
print(strg % tuple(cstr.tau))
strg = "glucan conversions each reactor are" + (nr-1)*" %g," + " %g"
print(strg % tuple(cstr.convGR))
strg = "xylan conversions each reactor are" + (nr-1)*" %g," + " %g"
print(strg % tuple(cstr.convXR))
strg = "facile fractions of glucan are" + (nr-1)*" %g," + " %g"
print(strg % tuple(cstr.yF))

# this is not glucose recovery yeild!!!!
print("\nglucan conversion yield is %g" % cstr.convG)
print("xylan conversion yield is %g" % cstr.convX)
print("total carbohydrate conversion yield is %g" % cstr.convCarb)
print( "Glucan/glucose mass balance is %g" % cstr.mbg )
print( "Xylan/xylose mass balance is %g" % cstr.mbx )
print("Enzyme mass balance is %g" % cstr.mbE)

# Save the outputs into a dictionary for Aspen
# system level result saved in a dictionary
system_level_result_dict = {
                'Total membrane units'                          : sum(membrane_units).tolist(),
                'Total membrane loop power consumption (kW)'    : sum(power_consumptions).tolist(),
                'Glucan conversion yield'                       : cstr.convG.tolist(),
                'Xylan conversion yield'                        : cstr.convX.tolist(),
                'Total carbohydrate conversion yield'           : cstr.convCarb.tolist(),
                }
# reactor level result
reactor_level_result_dict = {
                'Reactor Size (kg)'                     : cstr.mT,
                'Membrane units'                        : membrane_units,
                'Membrane loop pump power (kW)'         : power_consumptions,
                'Feed stream (kg/h)'                    : cstr.F[:,0],
                'Makeup buffer stream (kg/h)'           : cstr.F[:,1],
                'Permeate stream kg/h'                  : cstr.F[:,2],
                'Exit stream kg/h'                      : cstr.F[:,3],
                'FIS'                                   : cstr.fis,
                'Permeate glucose concentration (g/L)'  : cstr.rho3g,
                'Residence time (h)'                    : cstr.tau,
                'Glucan conversion'                     : cstr.convGR,
                'Xylan conversion'                      : cstr.convXR
                }

output_dict = {'CEH_output': {}}

# output system level info
output_dict['CEH_output']['System level'] = system_level_result_dict

# output reactor level info
for key in reactor_level_result_dict:
    for count,value in enumerate(reactor_level_result_dict[key]):
        if 'CEH Reactor {}'.format(count+1) in output_dict['CEH_output']:
            pass
        else:
            output_dict['CEH_output']['CEH Reactor {}'.format(count+1)] = {}

        output_dict['CEH_output']['CEH Reactor {}'.format(count+1)][key] = value.tolist()

dict_to_yaml([ve_params, output_dict], params_filename)

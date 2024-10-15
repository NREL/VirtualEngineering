"""
This script will be used to solve for CEH steady-state outcomes given
specified reactor *design* parameters.

"""

# Jonathan Stickel, 2017

import numpy as np

import ceh_cstr_model_multi as ceh

import sys
#####  updated by DAS 2021-09-29

##### YL 2022-03-02: include membrane loop system calculations
from membrane_loop_system import MembraneLoopSystem

from vebio.Utilities import dict_to_yaml, yaml_to_dict

if len(sys.argv) > 1:
    params_filename = sys.argv[1]
    ve_params = yaml_to_dict(params_filename)
else:
    raise Exception("VE parameters filename not provided")

## biomass feed, *** DMR *** inputs from 24h PFR reactor simulated by ehk_batch.py
if True:
    name = 'FY21Q4'
    # load kin params used for PFR sim
    import pickle
    fitOuts = pickle.load(open('FY21_kin_params.pkl', 'rb'))
    init_conds = pickle.load(open('FY21_init_conditions.pkl', 'rb'))

    # load PFR results from sim
    import pandas as pd
    PFR_result = pd.read_pickle('PFR_result.pkl')
    f1is = PFR_result['fis']
    f1G = (PFR_result['fGF'] + PFR_result['fGR'])
    x1G = f1G / f1is  # mass frac of IS that is ALL glucan
    y1F = PFR_result['fGF'] / f1G  # frac of all glucan that is facile
    x1X = PFR_result['fX'] / f1is  # mass frac of IS that is xylan
    rho1g = PFR_result['rhog']
    rho1x = PFR_result['rhox']
    rho1sL = PFR_result['rhosL']

    # F1 = 1 # kg/h -- simple basis, scale later
    # F1 = 307633 # kg/h -- 2019 design spreadsheet
    F1 = 328996 + 6918 # kg/h -- from Chris Kinchin

    # overall EH yield for 2019 SOT was 83% (xylan+glucan)

# number reactors
nr = 2  # (this value is updated by the solver later)
MFfluxes = [50, 30]  # LMH from 2020-2021 experiments, estimates for different levels of conversion
MFrps = [21.9, 36.5]  # calculated in 'flux flow.xlsx'

# nr = 5
# target enzyme loading -- enzyme is only fed to the first reactor
# lmbdE = 0.01 # kg/kg; *1000 to get mg/g
lmbdE = 0.012 # kg/kg; *1000 to get mg/g
# target reactor fis, each reactor
#fis = np.array([0.05, 0.05, 0.05])
"""
Yudong 03/11/2022:
Change from:
the fis of the target 2 CEH reactors match upstream batch purge stream fis
to:
the fis of the target 2 CEH reactors are now user input.
(this means that, some extra care may be needed to dilute/concentrate the Batch
purge stream fis to the fis target of the 1st reactor?)
"""
#fis = np.array(nr*[f1is])
fis = np.array([
                ve_params['CEH_input']['fis_1'],
                ve_params['CEH_input']['fis_2'],
                ])

# target *carbohydrate* conversion yield *for each reactor* -- is there an
# approach to specify overall conversion yield, allowing optimization of each
# reactor?
#xi = 0.54
#xi = 0.0001
#xi = 0.8 # this seems good if set for both reactors
#xi = np.array(nr*[xi])
# different target conversions for each reactor; resulting relative sizes of
# reactors were a bit unintuitive (to me) -- it is a balance of residence times
# and amounts of solids entering each reactor

# make target conversion yield user input
#xi = np.array([0.3, 0.8])
xi = np.array([
                ve_params['CEH_input']['target_conv_1'],
                ve_params['CEH_input']['target_conv_2'],
                ])

# outlet flow criteria of each reactor
#theta3 = np.array([0.75, 0.75, 0.75])  # permeate/(sum inlet streams)
# use "theta2" now!

#make theta2 a user input
#theta2set = 0.50
#theta2 = np.array(nr*[theta2set])
theta2 = np.array([
                ve_params['CEH_input']['theta2_1'],
                ve_params['CEH_input']['theta2_2'],
                ])

# enzyme rejection by membrane -- assume common for all CSTRs
# rejE = 0.
# rejE = 0.5
rejE = 1.  # all enzymes stay with solids for practical purposes (FY19 Q2 CEHD experiments)


# initiate solver object
reactor = ceh.CEHSolver(**fitOuts)

cstr = reactor.designsolve(F1, f1is, x1G, y1F, x1X, rho1g, rho1x, rho1sL,
                           lmbdE, fis, xi, theta2, rejE)
ceh.post(cstr)

nr = cstr.nr
#%%
##########
# TODO
# YL 03/02/2022
# calculate power consumption based on FIS, permeate flow rate, with
# assumed membrane area per membrane unit
#
# equations were obtained from experiments with fis = 5%, 7%
# let's say these equations are still valid for the fis conditions in this calc.
##########
p_m_dot = -cstr.F[:,2] ## mass flow rate of permeate kg/h
fis = cstr.fis ####
# assume density of sugar stream is equal to water?
p_rou = 1 #kg/L


print(f'++++ Pump Power Requirement for membrane loops of {nr} reactors ++++\n')
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
#%%
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

# system level result saved in a dictionary
system_level_result_dict = {
                'Total membrane units'                          : sum(membrane_units).tolist(),
                'Total membrane loop power consumption (kW)'    : sum(power_consumptions).tolist(),
                'Glucan conversion yield'                       : cstr.convG.tolist(),
                'Xylan conversion yield'                        : cstr.convX.tolist(),
                'Total carbohydrate conversion yield'           : cstr.convCarb.tolist(),
                }

# extra calcs for FY21 Q4
if 'name' in locals():
    if name=='FY21Q4':
        # 5 flow streams per CEH unit:
        # inlet, water, MF perm, purge, UF perm (if applicable)
        #     0,     1,       2,     3,       4

        # 8 species, in this order:
        # GF, GR, X, L, g, x, sL, ET
        #  0,  1, 2, 3, 4, 5,  6,  7

        """
        FIXME: something is messed up here on the sol-sugar based calcs.

        1.  Calcs for the cstr object vs incoming streams for glucose check out.

        2.  g_in vs g_out are almost identical after cstr. This means there
            is almost no reaction going on! Furthermore, the carb content
            coming in from the pfr seems very low compared with the actual
            yield reported there. >>>realized that x1G/y1F/x1X were being
            taken in as frac of IS, but were in fact saved and imported as
            fraction of total slurry. Now divided by incoming FIS.

        3.  Now carbs are looking good in/out PFR. Still
            can't get to a reasonable overall sugar yield though. Seems
            like carb calculated yield is through the roof but sugars
            are way behind even with no losses in downstream SLS.
            Setting reactor to essentially 0% conversion still
            is calculating positive yield of carbs! Something's wrong.

        4.  Calculated flow of carbs from PFR and into CSTR
            don't agree.

        5.  Next issue is overall_conv != overall_conv_carbs, even
            when rec=1 (no loss of sugars in purge SLS process).
            These two yields go about it different ways but should
            be equal. Are we losing sol sugars somewhere?

        """
        g_out_pfr = F1 * (1-f1is) * rho1g / 1000 / ceh.Mwg
        g_in_cstr = cstr.F[0, 0] * cstr.f[0, 0, 4] / ceh.Mwg
        # ^^^ equal
        G_out_pfr = F1 * PFR_result['fis'] * x1G / ceh.MwG
        G_in_cstr = cstr.F[0, 0] * (cstr.f[0, 0, 0] + cstr.f[0, 0, 1]) / ceh.MwG
        # ^^^ not equal!  # now equal, JJS 12/1/21 (I think you missused y1F)

        x_out_pfr = F1 * (1-f1is) * rho1x / 1000 / ceh.Mwx
        x_in_cstr = cstr.F[0, 0] * cstr.f[0, 0, 5] / ceh.Mwx
        X_out_pfr = F1 * PFR_result['fis'] * x1X / ceh.MwX
        X_in_cstr = cstr.F[0, 0] * cstr.f[0, 0, 2] / ceh.MwX

        # feedstock going into the PFR (separate script, values saved then loaded here)
        G_in = F1 * init_conds['fis0'] * init_conds['XG0'] / ceh.MwG
        X_in = F1 * init_conds['fis0'] * init_conds['XX0'] / ceh.MwX
        L_in = F1 * init_conds['fis0'] * init_conds['XL0'] / ceh.MwL

        # Get total sugars out of all MF permeates
        g_out_perm = -np.sum(cstr.F[:, 2] * cstr.f[:, 2, 4]) / ceh.Mwg # agrees 'glucoseperm'
        x_out_perm = -np.sum(cstr.F[:, 2] * cstr.f[:, 2, 5]) / ceh.Mwx
        # seprate SLS recovery of solutes from last purge stream
        rec = 0.95  # recovery rate of solutes from SLS operation (lignin press, etc)
        g_out_purge = rec * -cstr.F[-1, 3] * cstr.f[-1, 3, 4] / ceh.Mwg
        x_out_purge = rec * -cstr.F[-1, 3] * cstr.f[-1, 3, 5] / ceh.Mwx

        g_out = g_out_perm + g_out_purge
        x_out = x_out_perm + x_out_purge
        G_out = -cstr.F[-1, 3] * (cstr.f[-1, 3, 0] + cstr.f[-1, 3, 1]) / ceh.MwG
        X_out = -cstr.F[-1, 3] * cstr.f[-1, 3, 2] / ceh.MwX
        L_out = -cstr.F[-1, 3] * cstr.f[-1, 3, 3] / ceh.MwL

        # carbs to sugars yield+recovery calc
        overall_conv = (g_out + x_out) / (G_in + X_in)
        glucan_conv = g_out/G_in
        print('\nOverall conversion for FY21Q4 (sugar method w/ recovery considered) is {:.3f}'.format(overall_conv))
        cstr.convOverall = overall_conv


        # glucan + glucose mass balance, seems to be good; JJS 12/20/21
        mbGCEH = 1 - (g_out + G_out)/(G_in_cstr + g_in_cstr)
        mbGoverall = 1 - (g_out + G_out)/G_in
        # xylan + xylose mass balance -- problem here
        mbXCEH = 1 - (x_out + X_out)/(X_in_cstr + x_in_cstr)
        mbXoverall = 1 - (x_out + X_out)/X_in

        # carbs to leftover carbs yield calc (doesn't account for purge SLS efficiency)
        overall_conv_carbs = 1 - (G_out + X_out) / (G_in + X_in)
        print('\nOverall conversion for FY21Q4 (carbs method) is {:.3f}'.format(overall_conv_carbs))

        # cstr.convOverall = overall_conv

        overall_G_conv_carbs = 1 - (G_out) / (G_in)
        overall_X_conv_carbs = 1 - (X_out) / (X_in)
        overall_L_conv = 1 - (L_out) / (L_in)
        print('\nOverall conversion for FY21Q4 (carbs method) is {:.3f} glucan and {:.3f} xylan'.format(overall_G_conv_carbs, overall_X_conv_carbs))
        print('\nOverall conversion for FY21Q4 lignin is {:.3f}'.format(overall_L_conv))

        # add entries to output dict
        system_level_result_dict['Conversion for FY21Q4 (sugar method w/ recovery considered)'] = overall_conv.tolist()
        system_level_result_dict['Conversion for FY21Q4 (carbs method)'] = {}
        system_level_result_dict['Conversion for FY21Q4 (carbs method)']['Total'] = overall_conv_carbs.tolist()
        system_level_result_dict['Conversion for FY21Q4 (carbs method)']['Glucan'] = overall_G_conv_carbs.tolist()
        system_level_result_dict['Conversion for FY21Q4 (carbs method)']['Xylan'] = overall_X_conv_carbs.tolist()
        system_level_result_dict['Conversion for FY21Q4 (carbs method)']['Lignin'] = overall_L_conv.tolist()

if True:
    """
    Created on Thu Dec  7 11:51:43 2017
    @author: dsievers

    This script is for exporting previously computed CEH results to a stream flow worksheet.
    """

    # =============================================================================
    # membrane settings - these values were calculated in "flux flow.xlsx"
    # =============================================================================

    #MFflux = 20  # LMH, 100% enzyme retaining material (value estimate from Jonathan)
    #MFrp = 16.7

    # MFflux = 50  # LMH, Koch membrane, design value from 2016 milestone
    # MFrp = 7.5

    # MFflux = 100  # LMH, pure MF, case 3
    # MFrp = 7.49

    # UFflux = 10 * 2*50/15.  # LMH, see permeance figure from Jim's data (50 psi, x2 for temperature)



    # =============================================================================
    #
    # =============================================================================

    textdict = {
                'header' : ['1.1 (feed)'],
                'm' : [cstr.F[0,0]],
                'FIS' : [f1is],
                'xGt' : [cstr.f[0,0,0]+cstr.f[0,0,1]],
                'xXt' : [cstr.f[0,0,2]],
                'xLt' : [cstr.f[0,0,3]],
                'xsLt' : [cstr.f[0,0,6]],
                'xgt' : [cstr.f[0,0,4]],
                'xxt' : [cstr.f[0,0,5]],
                'xE' : [cstr.f[0,0,7]],
                'xG' : [x1G],
                'xX' : [x1X],
                'xL' : [1-(x1G+x1X)],
                'cg' : [rho1g],
                'cx' : [rho1x],
                'csL' : [cstr.f[0,0,6]],
                'ce' : [1e3*cstr.f[0,0,7]],
                'RT' : [''],
                'eG' : [''],
                'eX' : [''],
                }

    # row descriptors with dictionary keys
    lefter = [
            ['stream/reactor ID,','header'],
            ['size or total flow,(kg or kg/h)','m'],
            ['FIS,','FIS'],
            ['mass frac glucan overall,','xGt'],
            ['mass frac xylan overall,','xXt'],
            ['mass frac lignin overall,','xLt'],
            ['mass frac sol. lignin overall,','xsLt'],
            ['mass frac glucose overall,','xgt'],
            ['mass frac xylose overall,','xxt'],
            ['mass frac enzyme overall,','xE'],
            ['mass frac glucan in solids,','xG'],
            ['mass frac xylan in solids,','xX'],
            ['mass frac lignin in solids,','xL'],
            ['glucose concentration,(g/L)','cg'],
            ['xylose concentration,(g/L)','cx'],
            ['sol. lignin conc.,(g/L)','csL'],
            ['enzyme concentration,(g/L)','ce'],
            ['residence time,(h)','RT'],
            ['glucan conversion,','eG'],
            ['xylan conversion,','eX'],
            ]


    for i in range (cstr.nr):
        textdict['header'].append('CSTR '+str(i+1))
        textdict['header'].append(str(i+1)+'.2 (makeup)')
        textdict['header'].append(str(i+1)+'.3 (MF perm.)')
        textdict['header'].append(str(i+1)+'.5 (UF recycle)')
        textdict['header'].append(str(i+1)+'.4/1 (outlet/feed)')

        textdict['m'].append(cstr.mT[i])
        textdict['m'].append(cstr.F[i,1])
        textdict['m'].append(-cstr.F[i,2])
        textdict['m'].append(cstr.F[i,4])
        textdict['m'].append(-cstr.F[i,3])

        textdict['FIS'].append(cstr.fis[i])
        textdict['FIS'].append(0)
        textdict['FIS'].append(0)
        textdict['FIS'].append(0)
        textdict['FIS'].append(cstr.fis[i])

        textdict['xGt'].append(cstr.fk[i,0]+cstr.fk[i,1])
        textdict['xGt'].append(cstr.f[i,1,0]+cstr.f[i,1,1])
        textdict['xGt'].append(cstr.f[i,2,0]+cstr.f[i,2,1])
        textdict['xGt'].append(cstr.f[i,4,0]+cstr.f[i,4,1])
        textdict['xGt'].append(cstr.f[i,3,0]+cstr.f[i,3,1])

        textdict['xXt'].append(cstr.fk[i,2])
        textdict['xXt'].append(cstr.f[i,1,2])
        textdict['xXt'].append(cstr.f[i,2,2])
        textdict['xXt'].append(cstr.f[i,4,2])
        textdict['xXt'].append(cstr.f[i,3,2])

        textdict['xLt'].append(cstr.fk[i,3])
        textdict['xLt'].append(cstr.f[i,1,3])
        textdict['xLt'].append(cstr.f[i,2,3])
        textdict['xLt'].append(cstr.f[i,4,3])
        textdict['xLt'].append(cstr.f[i,3,3])

        textdict['xsLt'].append(cstr.fk[i,6])
        textdict['xsLt'].append(cstr.f[i,1,6])
        textdict['xsLt'].append(cstr.f[i,2,6])
        textdict['xsLt'].append(cstr.f[i,4,6])
        textdict['xsLt'].append(cstr.f[i,3,6])

        textdict['xgt'].append(cstr.fk[i,4])
        textdict['xgt'].append(cstr.f[i,1,4])
        textdict['xgt'].append(cstr.f[i,2,4])
        textdict['xgt'].append(cstr.f[i,4,4])
        textdict['xgt'].append(cstr.f[i,3,4])

        textdict['xxt'].append(cstr.fk[i,5])
        textdict['xxt'].append(cstr.f[i,1,5])
        textdict['xxt'].append(cstr.f[i,2,5])
        textdict['xxt'].append(cstr.f[i,4,5])
        textdict['xxt'].append(cstr.f[i,3,5])

        textdict['xE'].append(cstr.fk[i,7])
        textdict['xE'].append(cstr.f[i,1,7])
        textdict['xE'].append(cstr.f[i,2,7])
        textdict['xE'].append(cstr.f[i,4,7])
        textdict['xE'].append(cstr.f[i,3,7])

        textdict['xG'].append(cstr.xG[i])
        textdict['xG'].append(0)
        textdict['xG'].append(0)
        textdict['xG'].append(0)
        textdict['xG'].append(cstr.xG[i])

        textdict['xX'].append(cstr.xX[i])
        textdict['xX'].append(0)
        textdict['xX'].append(0)
        textdict['xX'].append(0)
        textdict['xX'].append(cstr.xX[i])

        textdict['xL'].append(1-(cstr.xG[i]+cstr.xX[i]))
        textdict['xL'].append(0)
        textdict['xL'].append(0)
        textdict['xL'].append(0)
        textdict['xL'].append(1-(cstr.xG[i]+cstr.xX[i]))

        textdict['cg'].append(cstr.rhog[i])
        textdict['cg'].append(0)
        textdict['cg'].append(cstr.rho3g[i])
        textdict['cg'].append(cstr.rho5g[i])
        textdict['cg'].append(cstr.rhog[i])

        textdict['cx'].append(cstr.rhox[i])
        textdict['cx'].append(0)
        textdict['cx'].append(cstr.rho3x[i])
        textdict['cx'].append(cstr.rho5x[i])
        textdict['cx'].append(cstr.rhox[i])

        textdict['csL'].append(1e3*cstr.fk[i,6])
        textdict['csL'].append(1e3*cstr.f[i,1,6])
        textdict['csL'].append(1e3*cstr.f[i,2,6])
        textdict['csL'].append(1e3*cstr.f[i,4,6])
        textdict['csL'].append(1e3*cstr.f[i,3,6])

        textdict['ce'].append(cstr.rhoE[i])
        textdict['ce'].append(1e3*cstr.f[i,1,7])
        textdict['ce'].append(1e3*cstr.f[i,2,7])
        textdict['ce'].append(1e3*cstr.f[i,4,7])
        textdict['ce'].append(1e3*cstr.f[i,3,7])

        textdict['RT'].append(cstr.tau[i])
        textdict['RT'].append('')
        textdict['RT'].append('')
        textdict['RT'].append('')
        textdict['RT'].append('')

        textdict['eG'].append(cstr.convGR[i])
        textdict['eG'].append('')
        textdict['eG'].append('')
        textdict['eG'].append('')
        textdict['eG'].append('')

        textdict['eX'].append(cstr.convXR[i])
        textdict['eX'].append('')
        textdict['eX'].append('')
        textdict['eX'].append('')
        textdict['eX'].append('')


    text = ''
    text += 'Nr. CSTRs,{:d}\n'.format(nr)
    text += 'FIS,{:.1f},%\n'.format(f1is*100)
    text += 'enzyme loading,{:.0f},mg/g,(independent calculation)\n'.format((cstr.F[0,1]*cstr.f[0,1,7])/(cstr.F[0,0]*np.sum(cstr.f[0,0,0:2]))*1e3)  # lmbdE*1e3
    text += 'theta2 (makeup/feed) for reactor1 :{:.2f}; for reactor2: {:.2f},\n'.format(theta2[0],theta2[1]) # hardcode reactor number=2 for now
    text += 'glucan conversion,{:.1f},%\n'.format(cstr.convG*100)
    text += 'xylan conversion,{:.1f},%\n'.format(cstr.convX*100)
    text += 'total carbohydrate conversion,{:.1f},%\n'.format(cstr.convCarb*100)
    text += 'glucan/glucose mass balance,' + str(cstr.mbg) + '\n'
    text += 'enzyme mass balance,' + str(cstr.mbE) + '\n'
    text += '\n\n'

    for i in range(len(lefter)):
        text += lefter[i][0] + ','
        for t in textdict[lefter[i][1]]:
            text += str(t) + ','
        text = text[:-1] + '\n'

    #text += '\n\n'
    #text += 'enzyme feed concentration,{:.1f},g/L\n'.format(cstr.rho2ET)
    #text += 'enzyme recycle concentration,{:.1f},g/L\n'.format(cstr.rhoretE)
    text += '\n\n'
    text += 'total MF permeate,{:.1f},kg/h\n'.format(cstr.Fperm)
    if cstr.F[:,4].sum()/cstr.F[0,0] > 1e-6:
        text += 'UF retentate,{:.1f},kg/h\n'.format(cstr.F[:,4].sum())
        text += 'UF permeate,{:.1f},kg/h\n'.format(cstr.Fufperm)
        text += 'UF permeate glucose concentration,{:.1f},g/L\n'.format(cstr.rhogufperm)
        text += 'UF permeate xylose concentration,{:.1f},g/L\n'.format(cstr.rhoxufperm)
    else:
        text += 'permeate glucose concentration,{:.1f},g/L\n'.format(cstr.rhogufperm)
        text += 'permeate xylose concentration,{:.1f},g/L\n'.format(cstr.rhoxufperm)
    text += '\n\n'
    text += 'MF area reactor 1,{:.0f},m2\n'.format(-cstr.F[0][2]/MFfluxes[0])
    text += 'MF area reactor 2,{:.0f},m2\n'.format(-cstr.F[1][2]/MFfluxes[1])
    text += 'MF recirc rate reactor 1,{:.0f},kg/h\n'.format(-cstr.F[0][2]*MFrps[0])
    text += 'MF recirc rate reactor 2,{:.0f},kg/h\n'.format(-cstr.F[1][2]*MFrps[1])

    with open("flowsheet.csv", "w") as f:
        f.write(text)


# Save the outputs into a dictionary for Aspen

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

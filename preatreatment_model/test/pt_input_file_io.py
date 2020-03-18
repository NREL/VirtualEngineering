"""
Functions to read in and write out the 'pretreat_defs.inp' file
"""
# Jonathan Stickel, 2014


import numpy as np


def readinpfile(filename):
    inpfile = open(filename, 'r')

    # read the entire file into a list of strings
    lines = inpfile.readlines()
    inpfile.close()
    
    # mesh parameters
    meshp = dict()
    i=1 # skip first two lines
    i+=1; meshp['coord_system'] = int(lines[i].split()[-1])
    i+=1; meshp['minx'] = float(lines[i].split()[-1])
    i+=1; meshp['maxx'] = float(lines[i].split()[-1])
    i+=1; meshp['enum'] = int(lines[i].split()[-1])
    i+=1; meshp['pord'] = int(lines[i].split()[-1])
    i+=1; meshp['time_order'] = int(lines[i].split()[-1])
    i+=1; meshp['tstep'] = float(lines[i].split()[-1])
    i+=1; meshp['ftime'] = float(lines[i].split()[-1])
    i+=1; meshp['ptime'] = float(lines[i].split()[-1])
    i+=1; meshp['ppnts'] = int(lines[i].split()[-1])

    #scales
    scales = dict()
    i+=1; scales['xscl'] = float(lines[i].split()[-1])
    i+=1; scales['tscl'] = float(lines[i].split()[-1])
    i+=1; scales['stsc'] = float(lines[i].split()[-1])
    i+=1; scales['lisc'] = float(lines[i].split()[-1])
    i+=1; scales['Tesc'] = float(lines[i].split()[-1])
    i+=1; scales['xysc'] = float(lines[i].split()[-1])
    i+=1; scales['xgsc'] = float(lines[i].split()[-1])
    i+=1; scales['xssc'] = float(lines[i].split()[-1])
    i+=1; scales['fusc'] = float(lines[i].split()[-1])

    # initial and boundary conditions
    IBCs = dict()
    i+=2 # skip two lines 
    i+=1; IBCs['xyfr'] = float(lines[i].split()[-1])
    i+=1; IBCs['poro'] = float(lines[i].split()[-1])
    i+=1; IBCs['bkst'] = float(lines[i].split()[-1])
    i+=1; IBCs['acid'] = float(lines[i].split()[-1])
    i+=1; IBCs['lifr'] = float(lines[i].split()[-1])
    i+=1; IBCs['Temp'] = float(lines[i].split()[-1])
    i+=1; IBCs['stmT'] = float(lines[i].split()[-1])

    # Energetics
    rrates = {'stdf':np.zeros(2), 'xgdf':np.zeros(2), 'xsdf':np.zeros(2),\
              'fudf':np.zeros(2), 'kcond':np.zeros(3), 'kxylog':np.zeros(3),\
              'kxyl1':np.zeros(3), 'kxyl2':np.zeros(3), 'kfurf':np.zeros(3)}
    Egtcs = dict()
    i+=2 #skip two lines 
    i+=1; rrates['stdf'][1] = float(lines[i].split()[-1])
    i+=1; rrates['xgdf'][1] = float(lines[i].split()[-1])
    i+=1; rrates['xsdf'][1] = float(lines[i].split()[-1])
    i+=1; rrates['fudf'][1] = float(lines[i].split()[-1])
    i+=1; Egtcs['lhst'] = float(lines[i].split()[-1])
    i+=1; Egtcs['ksol'] = float(lines[i].split()[-1])
    i+=1; Egtcs['kliq'] = float(lines[i].split()[-1])
    i+=1; Egtcs['kgas'] = float(lines[i].split()[-1])
    i+=1; Egtcs['shso'] = float(lines[i].split()[-1])
    i+=1; Egtcs['shli'] = float(lines[i].split()[-1])
    i+=1; Egtcs['shgs'] = float(lines[i].split()[-1])
    i+=1; Egtcs['chtr'] = float(lines[i].split()[-1])

    # Densities and tortuosities
    deto = dict()
    i+=2 #skip two lines 
    i+=1; deto['sode'] = float(lines[i].split()[-1])
    i+=1; deto['lide'] = float(lines[i].split()[-1])
    i+=1; deto['gsde'] = float(lines[i].split()[-1])
    i+=1; deto['xymw'] = float(lines[i].split()[-1])
    i+=1; deto['xomw'] = float(lines[i].split()[-1])
    i+=1; deto['xsmw'] = float(lines[i].split()[-1])
    i+=1; deto['fumw'] = float(lines[i].split()[-1])
    i+=1; deto['xamw'] = float(lines[i].split()[-1])
    i+=1; deto['togs'] = float(lines[i].split()[-1])
    i+=1; deto['toli'] = float(lines[i].split()[-1])
    i+=1; rrates['stdf'][0] = float(lines[i].split()[-1])
    i+=1; rrates['xgdf'][0] = float(lines[i].split()[-1])
    i+=1; rrates['xsdf'][0] = float(lines[i].split()[-1])
    i+=1; rrates['fudf'][0] = float(lines[i].split()[-1])

    # Reaction rates
    i+=2 #skip two lines 
    i+=1; rrates['kcond'] = np.array(lines[i].split()[-3:]).astype(float)
    i+=1; rrates['kevap'] = np.array(lines[i].split()[-3:]).astype(float)
    i+=1; rrates['kxylog'] = np.array(lines[i].split()[-3:]).astype(float)
    i+=1; rrates['kxyl1'] = np.array(lines[i].split()[-3:]).astype(float)
    i+=1; rrates['kxyl2'] = np.array(lines[i].split()[-3:]).astype(float)
    i+=1; rrates['kfurf'] = np.array(lines[i].split()[-3:]).astype(float)

    return meshp, scales, IBCs, rrates, Egtcs, deto

    
def writeinpfile(filename, meshp, scales, IBCs, rrates, Egtcs, deto):

    # convert all the elements to strings for writing out
    # meshp = {k: str(v) for k, v in meshp.iteritems()}
    # scales = {k: str(v) for k, v in scales.iteritems()}
    # IBCs = {k: str(v) for k, v in IBCs.iteritems()}
    # Egtcs = {k: str(v) for k, v in Egtcs.iteritems()}
    # deto = {k: str(v) for k, v in deto.iteritems()}
    # rrates = {k: v.astype(str) for k, v in rrates.iteritems()}
    # Change for Python 3
    meshp = {k: str(v) for k, v in meshp.items()}
    scales = {k: str(v) for k, v in scales.items()}
    IBCs = {k: str(v) for k, v in IBCs.items()}
    Egtcs = {k: str(v) for k, v in Egtcs.items()}
    deto = {k: str(v) for k, v in deto.items()}
    rrates = {k: v.astype(str) for k, v in rrates.items()}
    
    inpfile = open(filename,'w')
    
    #mesh parameters
    inpfile.write("Mesh parameters\n\n")
    inpfile.write("Coord_system\t" + meshp['coord_system'] + "\n")
    inpfile.write("Min_x\t"+ meshp['minx'] +"\n")
    inpfile.write("Max_x\t"+ meshp['maxx'] +"\n")
    inpfile.write("No_of_elements\t"+ meshp['enum'] +"\n")
    inpfile.write("Poly_order\t"+ meshp['pord'] +"\n")
    inpfile.write("Time_order\t" + meshp['time_order'] + "\n")
    inpfile.write("Timestep\t"+ meshp['tstep'] +"\n")
    inpfile.write("Finaltime\t"+ meshp['ftime'] +"\n")
    inpfile.write("printtime\t"+ meshp['ptime'] +"\n")
    inpfile.write("printpoints\t"+ meshp['ppnts'] +"\n")

    #scales
    inpfile.write("Length_scale\t"  +  scales['xscl'] +"\n")
    inpfile.write("Time_scale\t"    +  scales['tscl'] +"\n")
    inpfile.write("Steam_scale\t"   +  scales['stsc'] +"\n")
    inpfile.write("Liquid_scale\t"  +  scales['lisc'] +"\n")
    inpfile.write("Temp_scale\t"    +  scales['Tesc'] +"\n")
    inpfile.write("Xylan_scale\t"   +  scales['xysc'] +"\n")
    inpfile.write("Xylog_scale\t"   +  scales['xgsc'] +"\n")
    inpfile.write("Xylose_scale\t"  +  scales['xssc'] +"\n")
    inpfile.write("Furfural_scale\t"+  scales['fusc'] +"\n")

    #Initial and BCs
    inpfile.write("\nInitial_and_BCs\n")

    inpfile.write("Xylan_solid_fraction\t" +  IBCs['xyfr'] +"\n")
    inpfile.write("Initial_porosity\t" +  IBCs['poro'] +"\n")
    inpfile.write("Bulk_steam_conc\t" +  IBCs['bkst'] +"\n")
    inpfile.write("Initial_acid_conc\t" +  IBCs['acid'] +"\n")
    inpfile.write("Initial_liq_fraction\t" +  IBCs['lifr'] +"\n")
    inpfile.write("Initial_Temperature\t" +  IBCs['Temp'] +"\n")
    inpfile.write("Steam_Temperature\t" +  IBCs['stmT'] +"\n")

    # Energetics
    inpfile.write("\nEnergetics\n")
    inpfile.write("Diffusion_activ_energy_steam\t"+  rrates['stdf'][1] +"\n")
    inpfile.write("Diffusion_activ_energy_Xylog\t"+  rrates['xgdf'][1] +"\n")
    inpfile.write("Diffusion_activ_energy_Xylos\t"+  rrates['xsdf'][1] +"\n")
    inpfile.write("Diffusion_activ_energy_Furf\t" +  rrates['fudf'][1] +"\n")
    inpfile.write("Latent_heat_of_steam\t" +  Egtcs['lhst'] +"\n")
    inpfile.write("Thermal_Conductivity_solid\t" +  Egtcs['ksol'] +"\n")
    inpfile.write("Thermal_Conductivity_liquid\t" +  Egtcs['kliq'] +"\n")
    inpfile.write("Thermal_Conductivity_gas\t" +  Egtcs['kgas'] +"\n")
    inpfile.write("Specific_heat_solid\t" +  Egtcs['shso'] +"\n")
    inpfile.write("Specific_heat_liquid\t" +  Egtcs['shli'] +"\n")
    inpfile.write("Specific_heat_gas\t" +  Egtcs['shgs'] +"\n")
    inpfile.write("Convective_heat_transfer\t" +  Egtcs['chtr'] +"\n")

    # Densities and tortuosities
    inpfile.write("\nDensities_and_tortuosities\n")
    inpfile.write("Solid_density\t"     +  deto['sode'] +"\n")
    inpfile.write("Liquid_density\t"    +  deto['lide'] +"\n")
    inpfile.write("Gas_density\t"       +  deto['gsde'] +"\n")
    inpfile.write("Mol_wt_Xylan\t"      +  deto['xymw'] +"\n")
    inpfile.write("Mol_wt_Xylog\t"      +  deto['xomw'] +"\n")
    inpfile.write("Mol_wt_Xylose\t"     +  deto['xsmw'] +"\n")
    inpfile.write("Mol_wt_Furf\t"       +  deto['fumw'] +"\n")
    inpfile.write("Mol_wt_water\t"      +  deto['xamw'] +"\n")
    inpfile.write("Tortuosity_gas\t"    +  deto['togs'] +"\n")
    inpfile.write("Tortuosity_liq\t"    +  deto['toli'] +"\n")
    inpfile.write("pore_diameter_nm\t"       +  rrates['stdf'][0] +"\n")
    inpfile.write("xylog_mol_size_nm\t"      +  rrates['xgdf'][0] +"\n")
    inpfile.write("xylose_mol_size_nm\t"     +  rrates['xsdf'][0] +"\n")
    inpfile.write("furf_mol_size_nm\t"       +  rrates['fudf'][0] +"\n")

    # Reaction rates
    inpfile.write("\nReaction_rates\n")
    inpfile.write("condensation\t"+ rrates['kcond'][0] +"\t")
    inpfile.write( rrates['kcond'][1] +"\t"+ rrates['kcond'][2] +"\n")
    
    inpfile.write("evaporation\t"+ rrates['kevap'][0] +"\t")
    inpfile.write( rrates['kevap'][1] +"\t"+ rrates['kevap'][2] +"\n")

    inpfile.write("xylan_to_xylog\t"+ rrates['kxylog'][0] +"\t")
    inpfile.write( rrates['kxylog'][1] +"\t"+ rrates['kxylog'][2] +"\n")

    inpfile.write("xylan_to_xylose\t"+ rrates['kxyl1'][0] +"\t")
    inpfile.write( rrates['kxyl1'][1] +"\t"+ rrates['kxyl1'][2] +"\n")

    inpfile.write("xylog_to_xylose\t"+ rrates['kxyl2'][0] +"\t")
    inpfile.write( rrates['kxyl2'][1] +"\t"+ rrates['kxyl2'][2] +"\n")

    inpfile.write("xylose_to_furf\t"+ rrates['kfurf'][0] +"\t")
    inpfile.write( rrates['kfurf'][1] +"\t"+ rrates['kfurf'][2] )

    inpfile.write('\n')
    inpfile.close()

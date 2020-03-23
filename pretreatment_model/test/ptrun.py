import pt
import matplotlib.pyplot as plt
import numpy as np
import pt_input_file_io as pt_input
import timeit as timerlib
import yaml
import sys

# run pretreatment model
inputfilename='pretreat_defs.inp'
outputfilebase = 'ptsim_'


inputfile = open(inputfilename, 'r')

meshp, scales, IBCs, rrates, Egtcs, deto =\
    pt_input.readinpfile('pretreat_defs.inp')


# If applicable, load the input file into a dictionary
if len(sys.argv) > 1:
    input_filename = sys.argv[1]
    with open(input_filename) as fp:
        input_dict = yaml.load(fp, Loader = yaml.FullLoader)
    # print(input_dict)

    # Override pretreat_defs.inp definitions with those from the pretreatment widgets
    IBCs['acid'] = input_dict['initial_acid_conc']
    IBCs['stmT'] = input_dict['steam_temperature']
    IBCs['bkst'] = input_dict['bulk_steam_conc']
    meshp['ftime'] = input_dict['final_time']
    IBCs['xyfr'] = input_dict['xylan_solid_fraction']
    IBCs['lifr'] = 1.0 - input_dict['initial_solid_fraction']
    IBCs['poro'] = input_dict['initial_porosity']
else:
    input_dict = {}


# read in number of elements from input file
nelem = meshp['enum']

ppelem = meshp['ppnts']
gneq = 7

# calculate the shape of the output array
m = nelem*(ppelem - 1) + 1
n = gneq + 1

# read in final time
finaltime = meshp['ftime']

#establish parameters for porosity and time dependent [acid] calcs
fx0 = IBCs['xyfr'] # initial fraction of xylan in the solids, a.k.a., X_X0 
ep0 = IBCs['poro']
cacid0 = IBCs['acid']
eL0 = IBCs['lifr']
l = meshp['maxx']

new_inputfilename = 'pretreat_defs_updated.inp'
pt_input.writeinpfile(new_inputfilename, meshp, scales, IBCs, rrates, Egtcs, deto)

# run the simulation
simtime=-timerlib.default_timer()
solnvec = pt.main(m,n, new_inputfilename)
simtime=simtime+timerlib.default_timer()

#solnvec=pt.ptmain.interpsoln
n=len(solnvec)

x=solnvec[:,0]
steam=solnvec[:,1]
liquid=solnvec[:,2]
Temp=solnvec[:,3]
xylan=solnvec[:,4] # what is this? units? dimensionless concentration? what
                   # basis? JJS 3/22/20
xylog=solnvec[:,5]
xylose=solnvec[:,6]
furfural=solnvec[:,7]
porosity=ep0+fx0*(1.0-ep0)-xylan
cacid = cacid0*eL0/liquid

# integrate concentrations to determine bulk xylose, xylog, and furfural concentrations
solidvfrac = 1.0-porosity
xylanweight = np.trapz(xylan,x)
solidweight = np.trapz(solidvfrac,x)
liquid_bulk = np.trapz(liquid,x)
gas_bulk    = np.trapz(porosity-liquid,x)

xylan_bulk    = xylanweight/solidweight # this looks like X_X (fraction of solids)?
xylose_bulk   = np.trapz(xylose,x)/liquid_bulk
xylog_bulk    = np.trapz(xylog,x)/liquid_bulk
furfural_bulk = np.trapz(furfural,x)/liquid_bulk
steam_bulk    = np.trapz(steam*(porosity-liquid),x)/gas_bulk

M_xylose = 150.0
M_furf   = 100.0
M_xylog  = 450.0

print( "\n**************************************************************************")
print( "xylan weight:",        xylanweight)
print( "solid weight:",        solidweight)
print( "[Xylan] (w/w) ",       xylan_bulk)
print( "[xylose] (M/ml) =" ,   xylose_bulk)
print( "[xylose] (g/L) =" ,    xylose_bulk*1000*M_xylose)
print( "[xylog] (M/ml) =" ,    xylog_bulk)
print( "[xylog] (g/L) =" ,     xylog_bulk*1000*M_xylog)
print( "[furfural] (M/ml) =" , furfural_bulk)
print( "[furfural] (g/L) =" ,  furfural_bulk*1000*M_furf)
print( "Avg. liquid = ",       liquid_bulk)
print( "Avg. gas    = ",       gas_bulk)
print( "weight of steam added per 1g of initial liquid:",(liquid_bulk-eL0*l)/(eL0*l))
print( "Fraction of Insoluble solids:", solidweight/(solidweight+liquid_bulk))
print( "simulation time: %f seconds"%(simtime))
print( "**************************************************************************\n\n")
print( "Mass balance calculations")
print( "*************************")

xylanweight0 = fx0*(1-ep0)*l
print( "initial xylan mass   (density =  1 g/ml):",    xylanweight0)
print( "final xylan mass     (density =  1 g/ml):",    xylanweight)
print( "reacted xylan mass   (density =  1 g/ml):",    fx0*(1-ep0)*l-xylanweight)

prodmass = liquid_bulk*(xylose_bulk*M_xylose + xylog_bulk*M_xylog + furfural_bulk*M_furf)
reactmass = xylanweight0 - xylanweight
conv = reactmass/xylanweight0

print( "liquid weight (density = 1 g/ml):", liquid_bulk)
print( "liquid volume (density = 1 g/ml):", liquid_bulk)
print( "xylose mass                     :", liquid_bulk*xylose_bulk*M_xylose)
print( "xylooligomer mass               :", liquid_bulk*xylog_bulk*M_xylog)
print( "furfural mass                   :", liquid_bulk*furfural_bulk*M_furf)
print( "total mass of products          :", prodmass)
print( "total xylan conversion (%)      :", 100*conv)
print( "% mass balance                  :", 100*(1.0-(prodmass-reactmass)/reactmass))

# compute an updated glucan fraction based on xylan conversion
X_G = input_dict['glucan_solid_fraction']/(1 - input_dict['xylan_solid_fraction']*conv)

# Save the outputs into a dictionary
output_dict = {}
output_dict['fis_0'] = float(solidweight/(solidweight+liquid_bulk))
output_dict['conv'] = float(conv)
output_dict['X_X'] = float(xylan_bulk) # is this correct? JJS 3/22/20
output_dict['X_G'] = float(X_G)
output_dict['rho_x'] = float(xylose_bulk*1000*M_xylose)
output_dict['rho_f'] = float(furfural_bulk*1000*M_furf) 
if len(sys.argv) > 2:
    # Save the output dictionary to a .yaml file
    output_filename = sys.argv[2]
    with open(output_filename, 'w') as fp:
        yaml.dump(output_dict, fp)



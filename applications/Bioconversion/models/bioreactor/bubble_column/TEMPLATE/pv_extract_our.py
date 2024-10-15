"""
extract data, using Paraview-python modules, to numpy arrays

This script will focus on getting volume-average data for scalar parameters in
the liquid phase

""" 
# H Sitaraman and J Stickel, 2020

import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 
import sys
from sys import argv


ofreader = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
ofreader.CaseType = 'Reconstructed Case'
ofreader.MeshRegions = ['internalMesh']
ofreader.SkipZeroTime = 0  #dont know why this is not working

ourmax=##VARIABLE_OXYGEN_UPTAKE_RATE##
ko=0.01

t = np.array(ofreader.TimestepValues)
N = t.size
mwO2 = 32.0 # g/mol
rhol = 1000.0 # kg/m^3 -- could get an average value from data...
kg_to_g = 1000.0
calc0 = pv.Calculator(Input=ofreader, AttributeType='Cell Data',\
                         ResultArrayName='o2conc',\
                         Function='O2.liquid*%5.2e*%5.2e/%5.2e'%(rhol,kg_to_g,mwO2))
calc1 = pv.Calculator(Input=calc0, AttributeType='Cell Data',\
                         ResultArrayName='OUR',\
                         Function='%5.2e*(o2conc/(o2conc+%5.2e)*alpha.liquid)'%(ourmax,ko))

# threshold filter to get only the "aerated liquid"
liquidthreshold = pv.Threshold(Input=calc1, Scalars=['CELLS', 'alpha.gas'],\
                               ThresholdRange=[0., 0.6])


integrateliq = pv.IntegrateVariables(Input=liquidthreshold)
interval = 1
idx    = range(0,N,interval)
tint   = t[idx]
nint   = tint.size
Vl     = np.zeros(nint) # liquid volume; m^3
ouravg = np.zeros(nint) # avg our

print("processing time = %g" % tint[0])
pv.UpdatePipeline(time=tint[0], proxy=integrateliq)
idat       = dsa.WrapDataObject( pv.servermanager.Fetch(integrateliq) )
Vl[0]      = idat.CellData['alpha.liquid'].item()
ouravg[0]  = idat.CellData['OUR'].item()

for i in range(1,nint):
    print("processing time = %g" % tint[i])
    pv.UpdatePipeline(time=tint[i], proxy=integrateliq)
    idat = dsa.WrapDataObject( pv.servermanager.Fetch(integrateliq) )
    Vl[i]      = idat.CellData['alpha.liquid'].item()
    ouravg[i]  = idat.CellData['OUR'].item()

# average values
ouravg = ouravg/Vl
outfile=open("our_avg.dat","w")

for i in range(len(tint)):
	outfile.write("%e\t%e\t%e\n"%(tint[i],ouravg[i],Vl[i]))

outfile.close()

# trace generated using paraview version 5.8.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
import numpy as np
from paraview import simple as pv
import vtk.numpy_interface.dataset_adapter as dsa 
import sys
from sys import argv
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection

solfoam = pv.OpenFOAMReader(FileName = './soln.foam') # just need to provide folder
solfoam.CaseType = 'Reconstructed Case'
solfoam.MeshRegions = ['internalMesh']
#solfoam.CellArrays = ['C']
solfoam.PointArrays = ['phis','phifs']
t = np.array(solfoam.TimestepValues)
N=t.size

ofvtkdata = pv.servermanager.Fetch(solfoam)
ofdata = dsa.WrapDataObject( ofvtkdata)
ofpts = np.array(ofdata.Points.Arrays[0])
ptsmin = ofpts.min(axis=0)
ptsmax = ofpts.max(axis=0)
print(ptsmin)
print(ptsmax)

print("doing time:",t[-1])

pv.UpdatePipeline(time=t[-1], proxy=solfoam)
pltline1 = pv.PlotOverLine(Input=solfoam,
Source='High Resolution Line Source')

pltline1.Source.Point1 = [0.25*ptsmin[0]+0.75*ptsmax[0], 0.5*(ptsmin[1]+ptsmax[1]),ptsmin[2]]
pltline1.Source.Point2 = [0.25*ptsmin[0]+0.75*ptsmax[0], 0.5*(ptsmin[1]+ptsmax[1]),ptsmax[2]]

idat1    = dsa.WrapDataObject(pv.servermanager.Fetch(pltline1))

phis  = abs(idat1.PointData['phis'])
phifs = abs(idat1.PointData['phifs'])
outarr=np.array([idat1.Points[:,2],phis,phifs])
np.savetxt("solids_along_line.dat",np.transpose(outarr),delimiter=" ")

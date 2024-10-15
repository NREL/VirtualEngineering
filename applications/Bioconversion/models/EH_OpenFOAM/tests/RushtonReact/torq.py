# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
import vtk.numpy_interface.dataset_adapter as dsa
from sys import argv
#### disable automatic camera reset on 'Show'

# create a new 'OpenFOAMReader'
solfoam = OpenFOAMReader(FileName='./soln.foam')
solfoam.MeshRegions = [argv[1]]
solfoam.CellArrays = ['visc', 'grad(U)']
t = np.array(solfoam.TimestepValues)
N = t.size

# create a new 'Generate Surface Normals'
surf = GenerateSurfaceNormals(Input=solfoam)

calc1 = Calculator(Input=surf, AttributeType='Point Data',\
                         ResultArrayName='fx',\
                         Function='-visc*( 2.0*grad(U)_0*Normals_X + (grad(U)_1+grad(U)_3)*Normals_Y + (grad(U)_2+grad(U)_6)*Normals_Z )')
calc2 = Calculator(Input=calc1, AttributeType='Point Data',\
                         ResultArrayName='fy',\
                         Function='-visc*( (grad(U)_1+grad(U)_3)*Normals_X + 2.0*grad(U)_4*Normals_Y + (grad(U)_5+grad(U)_7)*Normals_Z)')
calc3 = Calculator(Input=calc2, AttributeType='Point Data',\
                         ResultArrayName='fz',\
                         Function='-visc*( (grad(U)_2+grad(U)_6)*Normals_X + (grad(U)_5+grad(U)_7)*Normals_Y + 2.0*grad(U)_8*Normals_Z)')

calc4 = Calculator(Input=calc3, AttributeType='Point Data',\
                         ResultArrayName='torqx',\
                         Function='coordsY*fz-coordsZ*fy')

calc5 = Calculator(Input=calc4, AttributeType='Point Data',\
                         ResultArrayName='torqy',\
                         Function='coordsZ*fx-coordsX*fz')

calc6 = Calculator(Input=calc5, AttributeType='Point Data',\
                         ResultArrayName='torqz',\
                         Function='coordsX*fy-coordsY*fx')

int1=IntegrateVariables(Input=calc6)

outfile=open("torq.dat","w")
for i in range(N):
    UpdatePipeline(time=t[i], proxy=int1)
    idat = dsa.WrapDataObject( servermanager.Fetch(int1) )
    area = idat.CellData['Area'].item()
    torqx=idat.PointData['torqx'].item()
    torqy=idat.PointData['torqy'].item()
    torqz=idat.PointData['torqz'].item()
    print("processing time = %e\t%e\t%e\t%e\t%e" % (t[i],area,torqx,torqy,torqz))
    outfile.write("%e\t%e\n"%(t[i],torqz))
outfile.close()

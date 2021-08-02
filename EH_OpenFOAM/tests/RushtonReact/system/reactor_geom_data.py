import numpy as np

#geometry ========
Dt = 8.8               # Tank Diameter
Da = 4.0                 # impeller tip Diameter
H  = 9.0              # height of reactor (includes D/4 with the air phase only)
nimpellers = 1
C = [Dt/4,Dt*3/4,Dt*5/4,Dt*7/4,Dt*9/4]      # height of the center of impellers
W = 1.5*Da/5               # height of impeller (Width) 
L = Da/8               # impeller blade width (beyond the hub)
Dh =Da-2*L             # Hub Diameter
Lin = L                # impeller blade width (inside the hub)
J =  0.5*Dt/10             # Baffle Width
Wh = W/6               # Hub height (Width) 
polyrad=1.5*Da/8/2           # Stem radius

Z0 = 0.0               # bottom of reactor
Dmrf = (Da+Dt-2*J)/2   # MRF region Diameter

#mesh ========
nr  = 5  	  # mesh points per unit radial length
nz  = 10          # mesh points per unit axial length
Npoly = 1         # mesh points in the polygon at the axis
Na = 3            # mesh points in the azimuthal direction

nbaffles = 6          # number of baffles and impeller fins

nsplits=2*nbaffles    #we need twice the number of splits
dangle=2.0*np.pi/float(nsplits)

circradii=np.array([Dh/2-Lin,Dh/2,Da/2,Dmrf/2,Dt/2-J,Dt/2])
ncirc = len(circradii)
hub_circ   = 1 
inhub_circ = hub_circ-1  #circle inside hub
rot_circ   = hub_circ+1
mrf_circ   = rot_circ+1 
tank_circ = ncirc-1 

reacthts = [Z0]
baff_sections = []
baff_volumes = []
hub_volumes=[]
count=1
for n_imp in range(nimpellers):
    reacthts.append(Z0 + C[n_imp] -  W/2)

    baff_sections.append(count)
    baff_volumes.append(count)
    count=count+1

    reacthts.append(Z0 + C[n_imp] - Wh/2)

    baff_sections.append(count)
    baff_volumes.append(count)
    hub_volumes.append(count)
    count=count+1

    reacthts.append(Z0 + C[n_imp] + Wh/2)

    baff_sections.append(count)
    baff_volumes.append(count)
    count=count+1

    reacthts.append(Z0 + C[n_imp] +  W/2)
    baff_sections.append(count)
    count=count+1

reacthts.append(Z0+H)


nsections = len(reacthts)
nvolumes = nsections-1
meshz = nz*np.diff(reacthts)
meshz = meshz.astype(int)+1 #avoid zero mesh elements

all_volumes=range(nvolumes)
nonbaff_volumes=[sec for sec in all_volumes if sec not in baff_volumes]
nonstem_volumes=[0,1] #this is 0,1 no matter how many impellers are there


#note: stem_volumes include hub volumes also
#these are volumes where we miss out polygon block
stem_volumes=[sec for sec in all_volumes if sec not in nonstem_volumes]

#removes hub_volumes here for declaring patches
only_stem_volumes=[sec for sec in stem_volumes if sec not in hub_volumes]

#to define mrf region
#note that [1] is not a stem volume but baffles are there
all_mrf_volumes=[1]+stem_volumes

#for avoiding the stem
#mrf_volumes=all_mrf_volumes[0:-1]
mrf_volumes=all_mrf_volumes
print(mrf_volumes)

#increase grid points in the impeller section
for i in baff_volumes:
    meshz[i] *=2 

meshr = nr*np.diff(circradii)

#adding polygon to hub mesh resolution
meshr = np.append(nr*(circradii[0]-polyrad),meshr)
meshr = meshr.astype(int)

centeroffset     = 1 #one point on the axis
polyoffset       = nsplits #number of points on polygon
npts_per_section = centeroffset + polyoffset + ncirc*nsplits #center+polygon+circles

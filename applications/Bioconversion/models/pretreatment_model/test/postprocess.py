import numpy as np
import pt_input_file_io as pt_input
import matplotlib.pyplot as plt
import timeit as timerlib
import sys
from sys import argv
import glob
#======================================================================
def readtimehist(filename):
    it=0
    infile=open(filename,'r')

    t        = np.array([])
    xylose   = np.array([])
    xylog    = np.array([])
    furfural = np.array([])
    liquid   = np.array([])
    xylan    = np.array([])
    FIS      = np.array([])

    for line in infile:

        values=line.split()
        t         =  np.append(t,float(values[0]))
        xylose    =  np.append(xylose,float(values[1]))
        xylog     =  np.append(xylog,float(values[2]))
        furfural  =  np.append(furfural,float(values[3]))
        liquid    =  np.append(liquid,float(values[4]))
        xylan     =  np.append(xylan,float(values[5]))
        FIS       =  np.append(FIS,float(values[6]))

    infile.close()
    return(t,xylose,xylog,furfural,liquid,xylan,FIS)
#======================================================================
def findporosity(fX0,ep0,fXtilde):

    #porosity=ep0+fx0*(1.0-ep0)-fXtilde
    dratio=2.0

    k = (fX0/(1-fX0) + dratio)/(1-ep0)

    a = k;
    b = -(dratio+k*fXtilde);
    c = dratio*fXtilde - fXtilde;

    D = b*b - 4.*a*c

    r1 = 0.5*(-b + np.sqrt(D))/a
    r2 = 0.5*(-b - np.sqrt(D))/a

    porosity=1-r1

    return(porosity)
#======================================================================


inputfilename='pretreat_defs_updated.inp'
meshp, scales, IBCs, rrates, Egtcs, deto =\
        pt_input.readinpfile(inputfilename)

l = meshp['maxx']
nelem = meshp['enum']

finaltime = meshp['ftime']
ptime	  = meshp['ptime']
#establish parameters for porosity and time dependent [acid] calcs
fx0 = IBCs['xyfr']
ep0 = IBCs['poro']
cacid0 = IBCs['acid']
eL0 = IBCs['lifr']
rho_s = deto['sode']
rho_l = deto['lide']

M_xylose = deto['xsmw'] 
M_furf   = deto['fumw']
M_xylog  = deto['xomw']
#print(M_xylose,M_furf,M_xylog) % molecular weights? JJS 3/21/21

filenames = sorted(glob.glob(argv[1]), key=lambda f: int(f.split(".")[0].split("_")[1]))
numfiles=len(filenames)

expfilename=argv[2]
(texp,xyexp,xoexp,fexp,liqexp,xylexp,FISexp)=readtimehist(expfilename)

tsim=np.zeros(numfiles)
xysim=np.zeros(numfiles)
xosim=np.zeros(numfiles)
fsim=np.zeros(numfiles)
liqsim=np.zeros(numfiles)
xylsim=np.zeros(numfiles)
FISsim=np.zeros(numfiles)

for i in range(numfiles):

    infile=open(filenames[i],'r')
    it=0
    #time in minutes
    time = float(filenames[i].split('.')[0].split('_')[-1])*ptime/60.0

    x        = np.array([])
    steam    = np.array([])
    liquid   = np.array([])
    Temp     = np.array([])
    xylan    = np.array([])
    xylog    = np.array([])
    xylose   = np.array([])
    furfural = np.array([])

    for line in infile:

        if(it>0):
            values=line.split()
            x         =  np.append(x,float(values[0]))
            steam     =  np.append(steam,float(values[1]))
            liquid    =  np.append(liquid,float(values[2]))
            Temp      =  np.append(Temp,float(values[3]))
            xylan     =  np.append(xylan,float(values[4]))
            xylog     =  np.append(xylog,float(values[5]))
            xylose    =  np.append(xylose,float(values[6]))
            furfural  =  np.append(furfural,float(values[7]))

        it=it+1

    n=it-1
    porosity = np.array([])

    for j in range(n):
        porosity = np.append(porosity,findporosity(fx0,ep0,xylan[j]))


    # integrate concentrations to determine bulk xylose, xylog, and furfural concentrations
    solidvfrac = 1.0-porosity
    xylanweight = np.trapz(xylan,x)
    solidweight = np.trapz(solidvfrac,x)
    liquid_bulk = np.trapz(liquid,x)
    gas_bulk    = np.trapz(porosity-liquid,x)

    xylan_bulk    = xylanweight/solidweight
    xylose_bulk   = np.trapz(xylose,x)/liquid_bulk
    xylog_bulk    = np.trapz(xylog,x)/liquid_bulk
    furfural_bulk = np.trapz(furfural,x)/liquid_bulk
    steam_bulk    = np.trapz(steam*(porosity-liquid),x)/gas_bulk

    tsim[i]=time
    xysim[i]=xylose_bulk*1000*M_xylose
    xosim[i]=xylog_bulk*1000*M_xylog
    fsim[i]=furfural_bulk*1000*M_furf
    liqsim[i]=(liquid_bulk-eL0*l)/(eL0*l)
    xylsim[i]=xylan_bulk*100
    FISsim[i]=solidweight*rho_s/(solidweight*rho_s+liquid_bulk*rho_l)

    infile.close()

plt.close('all')
fig, ax=plt.subplots(2,3)
fig.tight_layout()

plt.ticklabel_format(style='sci', axis='y')

ax[0,0].plot(tsim,xysim,label='sim',color='red',linewidth=2)
#ax[0,0].plot(texp,xyexp,label='exp',color='blue',linewidth=2)
#ax[0,0].axis([0,1.2*max(texp),0,2.5*max(xyexp)])
ax[0,0].set_title("Xylose (g/L)")

ax[0,1].plot(tsim,xosim,label='sim',color='red',linewidth=2)
#ax[0,1].plot(texp,xoexp,label='exp',color='blue',linewidth=2)
#ax[0,1].axis([0,1.2*max(texp),0,2.5*max(xosim)])
ax[0,1].set_title("Xylog (g/L)")

ax[0,2].plot(tsim,fsim,label='sim',color='red',linewidth=2)
#ax[0,2].plot(texp,fexp,label='exp',color='blue',linewidth=2)
#ax[0,2].axis([0,1.2*max(texp),0,2.5*max(fexp)])
ax[0,2].set_title("Furfural (g/L)")

ax[1,0].plot(tsim,liqsim,label='sim',color='red',linewidth=2)
#ax[1,0].plot(texp,liqexp,label='exp',color='blue',linewidth=2)
#ax[1,0].axis([0,1.2*max(texp),0,2.5*max(liqexp)])
ax[1,0].set_title("Dilution (g/g)")

ax[1,1].plot(tsim,xylsim,label='sim',color='red',linewidth=2)
#ax[1,1].plot(texp,xylexp,label='exp',color='blue',linewidth=2)
#ax[1,1].axis([0,1.2*max(texp),0,2.5*max(xylexp)])
ax[1,1].set_title("Xylan fraction (%)")

ax[1,2].plot(tsim,FISsim,label='sim',color='red',linewidth=2)
#ax[1,2].plot(texp,FISexp,label='exp',color='blue',linewidth=2)
#ax[1,2].axis([0,1.2*max(texp),0,2.5*max(FISexp)])
ax[1,2].set_title("FIS")

#plt.legend(loc=1)
plt.show()

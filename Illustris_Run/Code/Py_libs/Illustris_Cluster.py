# This script reads Illustris1,2 &3 to produce files containing phaseSpace as well as the masses info of subhalos and 
# Phase space info is kept in ..../Illustris_Run/data/phaseSpace/Illustris{1..3}.csv
# Masses are kept in ..../Illustris_Run/data/Galaxy{1..3}.csv
# Phase space format is:     x,y,z,vx,vy,vz,mG,mDM
# Masses format is:    mG,mDM,m_,m_,mStellar,mBH

import matplotlib
matplotlib.use('Agg')


import numpy as np
from arepo import *
import matplotlib.pyplot as plt

# Path of the simulations
path1 = "/hits/universe/Illustris/L75n1820C/" 
path2 = "/hits/universe/Illustris/L75n910FP/" 
path3 = "/hits/universe/Illustris/L75n455FP/"
# Filenames for phaseSpace
filenamesALFALFA = ['Illustris1G','Illustris2G','Illustris3G']
massFilenamesALFALFA = ['Galaxy1','Galaxy2','Galaxy3']
filenamesSDSS = ['Illustris1H','Illustris2H','Illustris3H']
massFilenamesSDSS = ['Halo1','Halo2','Halo3']


# Snapshot number
snapnum = 135
for path,filenameG,massfnG,filenameH,massfnH 
in zip([path1,path2,path3],filenamesALFALFA,massFilenamesALFALFA, filenamesSDSS,massFilenamesSDSS): 
 
    # Loads groups to select the principal halo
    subSn = Subfind(path,snapnum,combineFiles=True, verbose = True)
    
    # Loads important data
    pos  = subSn.SubhaloPos
    vel  = subSn.SubhaloVel
    mass = subSn.SubhaloMassType[:,(0,1,4,5)]
    sph  = subSn.SubhaloStellarPhotometrics[:,5]
    # Have into account only subhalos with mass
    #for mag in sph[np.where(sph.T>0)[0]]:
    #    print(mag)
    indSDSS, = np.where((mass.T[0]>0) & (sph < -19)) # Index for SDSS galaxies
    indALFALFA, = np.where((mass.T[0]>0) )) # Index for ALFALFA galaxies
    #ind, = np.where((mass.T[0]>0) )
    print(len(ind))
    posALFALFA  = pos[indALFALFA]
    velALFALFA  = vel[indALFALFA]
    massALFALFA = mass[indALFALFA]
    sphALFALFA  = sph[indALFALFA]
    posSDSS  = pos[indSDSS]
    velSDSS  = vel[indSDSS]
    massSDSS = mass[indSDSS]
    sphSDSS  = sph[indSDSS]
	 
    # Histogram of luminosity band r
    sph[np.where(abs(sph) > 1e10)[0]] = 10
    hist,bins = np.histogram(sph,bins = 20)
    x = 0.5*(bins[:-1]+bins[1:])
    y = np.log10(hist+1)
    plt.plot(x,y,marker = '.', linewidth = 0, label = filename)
    plt.plot([-19,-19],[0,max(y)+1])
    plt.plot([0,0],[0,max(y)+1])
    plt.xlabel("R-Band Luminosity")
    plt.ylabel("Log10(n)")
    plt.savefig("../data/pics/Rband_lumFunc_19cut"+filename+".png")
    	

    # Creates and fills the numpy array to save
    arrALFALFA = np.zeros((len(posALFALFA),10),dtype = np.float)
    arrSDSS = np.zeros((len(posSDSS),10),dtype = np.float)
    arrALFALFA[:,:3]  = posALFALFA
    arrALFALFA[:,3:6] = velALFALFA
    arrALFALFA[:,6:]  = massALFALFA[:,:4]
    arrSDSS[:,:3]  = posSDSS
    arrSDSS[:,3:6] = velSDSS
    arrSDSS[:,6:]  = massSDSS[:,:4]
    np.savetxt("../data/phaseSpace/"+filenameG+".csv",arrALFALFA,delimiter = ',')
    np.savetxt("../data/"+massfnG+".csv",massALFALFA,delimiter = ',')
    np.savetxt("../data/phaseSpace/"+filenameH+".csv",arrSDSS,delimiter = ',')
    np.savetxt("../data/"+massfnH+".csv",massSDSS,delimiter = ',')
    
plt.legend()
plt.savefig("../data/pics/Rband_lumFunc.png")





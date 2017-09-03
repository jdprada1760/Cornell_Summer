# This script reads Illustris1,2 &3 to produce files containing phaseSpace as well as the masses info of subhalos and 
# Phase space info is kept in ..../Illustris_Run/data/phaseSpace/Illustris{1..3}.csv
# Masses are kept in ..../Illustris_Run/data/Galaxy{1..3}.csv
# Phase space format is:     x,y,z,vx,vy,vz,mG,mDM
# Masses format is:    mG,mDM,m_,m_,mStellar,mBH

import numpy as np
from arepo import *

# Path of the simulations
path1 = "/hits/universe/Illustris/L75n1820C/" 
path2 = "/hits/universe/Illustris/L75n910FP/" 
path3 = "/hits/universe/Illustris/L75n455FP/"
# Filenames for phaseSpace
filenames = ['Illustris1','Illustris2','Illustris3']
massFilenames = ['Galaxy1','Galaxy2','Galaxy3']
# Snapshot number
snapnum = 135
for path,filename,massfn in zip([path1,path2,path3],filenames,massFilenames): 
 
    # Loads groups to select the principal halo
    subSn = Subfind(path,snapnum,combineFiles=True, verbose = True)
    
    # Loads important data
    pos  = subSn.SubhaloPos
    vel  = subSn.SubhaloVel
    mass = subSn.SubhaloMassType
    sph  = subSn.SubhaloStellarPhotometrics[:,5]
    # Have into account only subhalos with mass
    #for mag in sph[np.where(sph.T>0)[0]]:
    #    print(mag)
    ind, = np.where((mass.T[0]>0) & (sph < 0))
    print(len(pos))
    pos  = pos[ind]
    vel  = vel[ind]
    mass = mass[ind]

    # Creates and fills the numpy array to save
    arr = np.zeros((len(pos),8),dtype = np.float)
    arr[:,:3]  = pos
    arr[:,3:6] = vel
    arr[:,6:]  = mass[:,:2]
    np.savetxt("../data/phaseSpace/"+filename+".csv",arr,delimiter = ',')
    np.savetxt("../data/"+massfn+".csv",mass,delimiter = ',')
    

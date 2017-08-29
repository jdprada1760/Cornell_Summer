import numpy as np
from arepo import *

# Path of the simulations
path1 = "/hits/universe/Illustris/L75n1820C/" 
path2 = "/hits/universe/Illustris/L75n910FP/" 
path3 = "/hits/universe/Illustris/L75n455FP/"
# Filenames
filenames = ['Illustris1','Illustris2','Illustris3']
# Snapshot number
snapnum = 135
for path,filename in zip([path1,path2,path3],filenames): 
 
    # Loads groups to select the principal halo
    subSn = Subfind(path,snapnum)
    
    # Loads important data
    pos  = subSn.subhalo.SubhaloPos
    vel  = subSn.subhalo.SubhaloVel
    mass = subSn.subhalo.SubhaloMassType
    
    # Creates and fills the numpy array to save
    arr = np.zeros((len(pos),12))
    arr[:,:3]  = pos
    arr[:,3:6] = vel
    arr[:,6:]  = mass
    np.savetxt("../data/phaseSpace/"+filename+".txt",arr)

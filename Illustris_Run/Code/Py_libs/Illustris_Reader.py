import numpy as np
import illustris_python as il
import sys
import time
# Generates the basepath depending on the simulation by parameter
basePath = "./Illustris-1"

######################################################
#                   Load Subhalos                    #
######################################################
start_time = time.time()

# Selects the fields to load for subhalos
fields = ['SubhaloCM','SubhaloVel','SubhaloStellarPhotometrics','SubhaloMassType']
# Loads subhalos in g
g = il.groupcat.loadSubhalos(basePath,135,fields=fields)

# Filters galaxies by magnitude
print "Filtering galaxies by magnitude..."
print "Number of galaxies before: " + str(g['count'])
# Creates lists for galaxies that fulfill the condition
gr = []
gv = []
# Filter the galaxies for its photometrics (only rband < -19)
i = 0
for i in range(g['count']):
    if g['SubhaloStellarPhotometrics'][i][5] > -19.0:
        gr.append([g['SubhaloCM'][i][0],g['SubhaloCM'][i][1],g['SubhaloCM'][i][2]])
        gv.append([g['SubhaloVel'][i][0],g['SubhaloVel'][i][1],g['SubhaloVel'][i][2]])

gr = np.array(gr)
gv = np.array(gv)
g['count'] = len(gr)

print "Number of galaxies after: " + str(len(gr))
print "Time elapsed: " + str(time.time() - start_time)



######################################################
#                     Load Halos                     #
######################################################

# Selects the fields to load for haloes
fields = ['GroupPos','GroupVel','GroupFirstSub','GroupMassType']
h = il.groupcat.loadHalos(basePath,135,fields=fields)
#print "LOL ya.\n"

# Creates the arrays to convert to files
print "Formating info..."
Halos = np.array([[0,0,0,0,0,0] for i in range(h['count']) ])
Galaxies = np.array([[0,0,0,0,0,0] for i in range(g['count']) ])

# Fills the arrays
for i in range(h['count']):
    Halos[i][0] = h['GroupPos'][i][0]
    Halos[i][1] = h['GroupPos'][i][1]
    Halos[i][2] = h['GroupPos'][i][2]
    Halos[i][3] = h['GroupVel'][i][0]
    Halos[i][4] = h['GroupVel'][i][1]
    Halos[i][5] = h['GroupVel'][i][2]
for i in range(g['count']):
    Galaxies[i][0] = gr[i][0]
    Galaxies[i][1] = gr[i][1]
    Galaxies[i][2] = gr[i][2]
    Galaxies[i][3] = gv[i][0]
    Galaxies[i][4] = gv[i][1]
    Galaxies[i][5] = gv[i][2]

# Saves the arrays to files
print "Writing files"
np.savetxt('Halo.csv', Halos, delimiter = ',')
np.savetxt('Galaxy.csv', Galaxies, delimiter = ',')
print "Total time elapsed: " + str(time.time() - start_time)

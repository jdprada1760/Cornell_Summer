import struct
import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il
import time

def read_CIC_scalar(filename):
    f = open(filename, "rb")
    dumb = f.read(38)

    dumb = f.read(4)
    n_x = f.read(4)
    n_y = f.read(4)
    n_z = f.read(4)
    nodes = f.read(8)
    x0 = f.read(4)
    y0 = f.read(4)
    z0 = f.read(4)
    dx = f.read(4)
    dy = f.read(4)
    dz = f.read(4)
    dumb = f.read(4)

    n_x = (struct.unpack('i', n_x))[0]
    n_y = (struct.unpack('i', n_y))[0]
    n_z = (struct.unpack('i', n_z))[0]
    nodes = (struct.unpack('q', nodes))[0]
    dx = (struct.unpack('f', dx))[0]
    dy = (struct.unpack('f', dy))[0]
    dz = (struct.unpack('f', dz))[0]
    x0 = (struct.unpack('f', x0))[0]
    y0 = (struct.unpack('f', y0))[0]
    z0 = (struct.unpack('f', z0))[0]
    print n_x, n_y, n_z, nodes, dx, dy, dz

    total_nodes = n_x * n_y *n_z
    dumb = f.read(4)
    array_data = f.read(total_nodes*4)
    dumb = f.read(4)
    format_s = str(total_nodes)+'f'
    array_data = struct.unpack(format_s, array_data)
    f.close()
    array_data  = np.array(array_data)
    array_data.resize(n_z,n_y,n_x)
    array_data = array_data.transpose()
    return {'eigenval':array_data, 'delta_x':dx}

data = read_CIC_scalar("../data/tweb/snap_135.s1.00.eigen_1")
eigen1 = data['eigenval']
data = read_CIC_scalar("../data/tweb/snap_135.s1.00.eigen_2")
eigen2 = data['eigenval']
data = read_CIC_scalar("../data/tweb/snap_135.s1.00.eigen_3")
eigen3 = data['eigenval']
#print("delta_x", data['delta_x'])
#print("shape", np.shape(data['eigenval']))


# Size of grid cells
delx = data['delta_x']

######################################################
#             Environment from eigen vals            #
######################################################
def getEnv(x,y,z):
    # Defines the eigenvals cut
    cut = 0.2
    mx = x >= cut
    my = y >= cut
    mz = z >= cut

    return int(mx)+int(my)+int(mz)

######################################################
#                   Load Subhalos                    #
######################################################
# cut: the cut in photometrics r-band
def getGalaxies(basePath,cut):
    start_time = time.time()
    # Selects the fields to load for subhalos
    fields = ['SubhaloCM','SubhaloStellarPhotometrics','SubhaloMassType']
    # Loads subhalos to variable g
    g = il.groupcat.loadSubhalos(basePath,135,fields=fields)

    # Counts the galaxies by magnitude and not null gas mass
    # filter the galaxies for its photometrics (only rband < -19)
    print "Number of galaxies before: " + str(g['count'])
    ngalaxies = len(np.where((g[fields[2]][:,0]!=0)*(g[fields[1]][:,5]<cut))[0])
    print "Number of galaxies after: " + str(ngalaxies)

    # Creates arrays for the candidate galaxies and fills it
    Galaxies = np.array([[0.0,0.0,0.0] for i in range(ngalaxies)])
    i = 0
    imax = 0
    for j in range(g['count']):
        if( (g['SubhaloMassType'][j][0] != 0) and ((g[fields[1]][j][5] < cut)) ):
        #if(True):
            #print g['SubhaloMassType'][j][0]
            ix = int(g['SubhaloCM'][j][0]/delx)
            iy = int(g['SubhaloCM'][j][1]/delx)
            iz = int(g['SubhaloCM'][j][2]/delx)
            ex = eigen1[ix,iy,iz]
            ey = eigen2[ix,iy,iz]
            ez = eigen3[ix,iy,iz]
            #print("ix+iy+iz", ex, ey, ez, int( int(ex > 0.2) + int(ey > 0.2) + int(ez > 0.2) ) )
            #Galaxies[i][0] = getEnv(ix,iy,iz)
            Galaxies[i][0] = getEnv(ex,ey,ez);
            Galaxies[i][1] = g['SubhaloMassType'][j][1]
            Galaxies[i][2] = g['SubhaloMassType'][j][0]
            i+=1
    #print "Time elapsed: " + str(time.time() - start_time)
    #print imax
    return Galaxies

basePath = "../Illustris-1"
Galaxies = getGalaxies(basePath,0)
np.savetxt('../data/twEnv/Tweb1.csv', Galaxies, delimiter = ',')
basePath = "../Illustris-2"
Galaxies = getGalaxies(basePath,0)
np.savetxt('../data/twEnv/Tweb2.csv', Galaxies, delimiter = ',')
basePath = "../Illustris-3"
Galaxies = getGalaxies(basePath,0)
np.savetxt('../data/twEnv/Tweb3.csv', Galaxies, delimiter = ',')
#plt.imshow(data['eigenval'][:,:,10].T)
#plt.show()

import struct
import numpy as np
import matplotlib.pyplot as plt
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

data = read_CIC_scalar("../data/tweb/snap_135.s1.00.eigen_2")
print("delta_x", data['delta_x'])
print("shape", np.shape(data['eigenval']))

plt.imshow(data['eigenval'][:,:,10].T)
plt.show()

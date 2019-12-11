import numpy as np

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

def obj_data_to_mesh3d(odata):
    # odata is the string read from an obj file
    vertices = []
    faces = []
    lines = odata.splitlines()

    for line in lines:
        slist = line.split()
        if slist:
            if slist[0] == 'v':
                vertex = list(map(float, slist[1:]))
                vertices.append(vertex)
            elif slist[0] == 'f':
                face = []
                for k in range(1, len(slist)):
                    face.append([int(s) for s in
                                 slist[k].replace('//', '/').split('/')])

                if len(face) > 3:  # triangulate the n-polyonal face, n>3
                    faces.extend([[face[0][0] - 1, face[k][0] - 1,
                                   face[k + 1][0] - 1] for k in
                                  range(1, len(face) - 1)])
                else:
                    faces.append([face[j][0] - 1 for j in range(len(face))])
            else:
                pass

    return np.array(vertices), np.array(faces)


def obj_viewer(obj_file):

    obj_data = open(obj_file,'r').read()
    vertices, faces = obj_data_to_mesh3d(obj_data)

    # Check whether all faces are triangular:
    if sum([1 if len(face) > 3 else 0 for face in faces]) == 0:
        print('True')
    else:
        print('False')

    x,y,z = vertices[:,:3].T
    I,J,K = faces.T

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z, c='r', marker='o', s=0.01)
    ax.scatter(-0.00021, -0.00061,  0.00410, c='b', marker='o', s=1)
    ax.scatter(0.00165, 0.00051, -0.000128, c='g', marker='o', s=1)
    ax.view_init(0, 0)
    plt.show()

    ax.view_init(90, 0)
    plt.show()



obj_viewer('/Users/mcosta/SPICE/JUICE/kernels/dsk/juice_sc_bus_v02.obj')
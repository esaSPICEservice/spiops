import numpy as np
import plotly.graph_objs as go
import urllib.request
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot

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


def draw_obj_and_points(obj, points=[], time='', notebook=False, title=''):

    req = urllib.request.Request('file://'+obj)
    try:
        response=urllib.request.urlopen(req)
        obj_data = response.read().decode('utf-8')
    except urllib.error.URLError as e:
        print(e.reason)

    vertices, faces = obj_data_to_mesh3d(obj_data)

    if sum([1 if len(face) > 3 else 0 for face in faces]) == 0:
        print('True')
    else:
        print('False')

    if notebook:
        init_notebook_mode(connected=True)

    vertices[:, 3:][0]  # color codes for the vertex 0

    x, y, z = vertices[:, :3].T
    I, J, K = faces.T

    mesh = dict(type='mesh3d',
                x=x,
                y=y,
                z=z,
                vertexcolor=vertices[:, 3:],
                i=I,
                j=J,
                k=K,
                name='',
                showscale=False
                )

    mesh.update(lighting=dict(ambient=0.18,
                              diffuse=1,
                              fresnel=.1,
                              specular=1,
                              roughness=.1),

                lightposition=dict(x=100,
                                   y=200,
                                   z=150))

    layout = dict(title=title,
                  font=dict(size=14, color='black'),
                  width=900,
                  height=800,
                  scene=dict(xaxis=dict(visible=False),
                             yaxis=dict(visible=False),
                             zaxis=dict(visible=False),
                             aspectratio=dict(x=1.2,
                                              y=0.9,
                                              z=1.5
                                              ),
                             camera=dict(eye=dict(x=1., y=1., z=0.5)),
                             ),
                  paper_bgcolor='rgb(235,235,235)',
                  margin=dict(t=175)
                  )

    pl_mygray= [[0.0, 'rgb(112, 123, 143)'],
               [0.25, 'rgb(133, 153, 165)'],
               [0.5,  'rgb(156, 184, 188)'],
               [0.75, 'rgb(185, 210, 210)'],
               [1.0,  'rgb(220, 233, 233)'],
               ]

    mesh.update(intensity=z, colorscale=pl_mygray)
    mesh.pop('vertexcolor')
    layout.update(paper_bgcolor='rgb(50,50,50)',
                  font=dict(size=14, color='white'))

    fig = go.Figure(data=[mesh], layout=layout)

    for point in points:
        fig.add_trace(go.Scatter3d(
                    x=[point[2][0]], y=[point[2][1]], z=[point[2][2]],
                    mode="markers", marker=dict(size=5,
                                                color=point[1]),
                                                name=point[0]))

    iplot(fig)

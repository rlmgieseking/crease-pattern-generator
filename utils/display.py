# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 08:18:03 2023

@author: gieseking
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import numpy as np


#def write_svg_cp(filename, xdim, ydim, edgelist):
def write_svg_cp(filename, shape, verts=False, edges=True, rules=True):
    # Layout notes:
    # 1. All points are scaled by a factor of 100 (a gore width of 1 is 100 mm)
    # 2. SVG files use the top left corner as (0,0) instead of the bottom left.
    #    Y dimensions of all points are translated accordingly.
    
    svg = open(filename,'w')
    xdim = str((shape.ngores + shape.overlap)*shape.gorewidth*10)
    ydim = str(shape.currentheight2d*10)
    svg.write('''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg 
  xmlns="http://www.w3.org/2000/svg"
  id="svg8"
  version="1.1" 
  height="''' + ydim + '''mm" 
  width="'''  + xdim + '''mm">
<g>
''')
    #svg.write('<rect x="0" y="0" height="10" width="10" ' +
    #          'style="fill:none;stroke:#000000;stroke-width:1">\n')
    if edges:
        for edge in shape.edges:
            if rules or edge.direction != 'U':
                x1, y1 = str(edge.end1.pos2d[0,0] * 10), str(float(ydim) - edge.end1.pos2d[0,1] * 10)
                x2, y2 = str(edge.end2.pos2d[0,0] * 10), str(float(ydim) - edge.end2.pos2d[0,1] * 10)
                svg.write('  <line x1="' + x1 + 'mm" y1="' + y1 + 'mm" x2="' + x2 + 'mm" y2="' + y2 + 'mm" ') 
                if edge.direction == 'B':
                    svg.write('stroke="black" stroke-width="2"')
                elif edge.direction == 'M':
                    svg.write('stroke="red" stroke-width="1.5"')
                elif edge.direction == 'V':
                    svg.write('stroke="blue" stroke-width="1.5"')
                else:
                    svg.write('stroke="gray" stroke-width="0.3"')
                svg.write(' />\n')
    svg.write('</g>\n</svg>')
    svg.close()

def plot2d(shape, verts=False, edges=True):
    fig = plt.figure()
    ax = fig.add_subplot()
    
    if verts:
        xdata, ydata = [],[]
        for vertex in shape.verts:
            xdata.append(vertex.pos2d[0,0])
            ydata.append(vertex.pos2d[0,1])
            ax.scatter(vertex.pos2d[0,0], vertex.pos2d[0,1], color='k', s=8)
    if edges:
        for edge in shape.edges:
            edge2d = [edge.end1.pos2d,edge.end2.pos2d]
            if edge.direction == "B":
                ax.plot([edge2d[0][0,0],edge2d[1][0,0]],
                        [edge2d[0][0,1],edge2d[1][0,1]], color = 'k')
            elif edge.direction == "M":
                ax.plot([edge2d[0][0,0],edge2d[1][0,0]],
                        [edge2d[0][0,1],edge2d[1][0,1]], color = 'r', linewidth=1)
            elif edge.direction == "V":
                ax.plot([edge2d[0][0,0],edge2d[1][0,0]],
                        [edge2d[0][0,1],edge2d[1][0,1]], color = 'b', linewidth=1)
            else:
                ax.plot([edge2d[0][0,0],edge2d[1][0,0]],
                        [edge2d[0][0,1],edge2d[1][0,1]], color = '0.7', linewidth=0.3)
                    

    ax.set_aspect('equal')
    # Hide axes in 3D plot
    ax.set_axis_off()
    plt.show()

def plot3d(filename, shape, verts=False, edges=False, faces=True):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    x3d, y3d, z3d = [],[],[]
    for vertex in shape.verts:
        x3d.append(vertex.pos3dvis[0,0])
        y3d.append(vertex.pos3dvis[0,1])
        z3d.append(vertex.pos3dvis[0,2])
        if verts:
            ax.scatter(vertex.pos3dvis[0,0],vertex.pos3dvis[0,1],vertex.pos3dvis[0,2], color='k')
        else:
            ax.scatter(vertex.pos3dvis[0,0],vertex.pos3dvis[0,1],vertex.pos3dvis[0,2], color='k', s=0.0)
    
    if faces:
        lightvector = np.array([[-1.0], [0.5], [-1.5]])
        lightvector /= np.linalg.norm(lightvector)
        ls = mpl.colors.LightSource(azdeg=110, altdeg=45)

        for face in shape.faces:
            verts = np.zeros((len(face.verts), 3))
            #print(verts)
            for i in range(len(face.verts)):
                #print(face.verts[i].pos3dvis[0])
                verts[i, :] = face.verts[i].pos3dvis[0]
            #print(verts)
            #verts = list(zip(*verts))
            #print(verts)
            facet = a3.art3d.Poly3DCollection([verts], 
                                              facecolors=face.color, edgecolor='k', linewidth=0.5,
                                              shade=True, lightsource=ls)
            #color_scaling = np.dot(face.normal, lightvector) * color_norm_scaling
            #print(color_scaling)
            #color = face.color + color_scaling
            #facet.set_facecolor(color)
            #facet.set_facecolor(face.color)
            #facet.set_edgecolor('k')
            ax.add_collection3d(facet)
    if edges:
        for edge in shape.edges:
            edge3d = [edge.end1.pos3dvis,edge.end2.pos3dvis]
            if edge.direction == "B" or edge.direction == "M" or edge.direction == "V":
                ax.plot([edge3d[0][0,0],edge3d[1][0,0]],
                        [edge3d[0][0,1],edge3d[1][0,1]],
                        [edge3d[0][0,2],edge3d[1][0,2]], color = 'k')
            #colors.append((1,0,0,1))
    limits = np.array([
            ax.get_xlim3d(),
            ax.get_ylim3d(),
            ax.get_zlim3d(),
        ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])
    ax.set_box_aspect([1,1,1])
    # Hide axes in 3D plot
    ax.set_axis_off()
    plt.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()
    #plt.savefig(filename)
    
def fold2d(filename, shape, verts=True, edges=True, faces=True):
    fold = open(filename,'w')
    fold.write('''{
    "file_spec": 1,
    "file_creator": "Vase-CP-Generator (Rebecca Gieseking)",
    "file_classes": ["singleModel"],
    "file_title": "''' + filename + ''' 3D shape",
    "frame_classes": ["unfoldedForm"],
    "frame_attributes": ["3D"]''')
    if verts:
        fold.write(',\n    "vertices_coords": [\n')
        for i, vert in enumerate(shape.verts):
            if i > 0:
                fold.write(',\n')
            fold.write('      [' + ('%.6f'%vert.pos2d[0,0]).rjust(12) + ','
                                 + ('%.6f'%vert.pos2d[0,1]).rjust(12) + ' ]')
        fold.write('\n    ]')
    if edges:
        fold.write(',\n    "edges_vertices": [\n')
        for i, edge in enumerate(shape.edges):
            if i > 0:
                fold.write(',\n')
            end1 = shape.verts.index(edge.end1)
            end2 = shape.verts.index(edge.end2)
            fold.write('      [' + str(end1).rjust(5) + ',' + str(end2).rjust(5) + ' ]')
        fold.write('\n    ],\n    "edges_assignment": [\n')
        for i, edge in enumerate(shape.edges):
            if i > 0:
                fold.write(',\n')
            fold.write('      "' + edge.direction + '"')
        fold.write('    ]')
    if faces:
        fold.write(',\n    "faces_vertices": [\n')
        for j, face in enumerate(shape.faces):
            if j > 0:
                fold.write(',\n')
            fold.write('      [')
            for i, vert in enumerate(face.verts):
                vind = shape.verts.index(vert)
                if i > 0:
                    fold.write(',')
                fold.write(str(vind).rjust(5))
            fold.write(' ]')
        fold.write('    ]')
    fold.write('}')    
    fold.close()

def fold3d(filename, shape, verts=True, edges=True, faces=True):
    fold = open(filename,'w')
    fold.write('''{
    "file_spec": 1,
    "file_creator": "Vase-CP-Generator (Rebecca Gieseking)",
    "file_classes": ["singleModel"],
    "file_title": "''' + filename + ''' 3D shape",
    "frame_classes": ["unfoldedForm"],
    "frame_attributes": ["3D"]''')
    if verts:
        fold.write(',\n    "vertices_coords": [\n')
        for i, vert in enumerate(shape.verts):
            if i > 0:
                fold.write(',\n')
            fold.write('      [' + ('%.6f'%vert.pos3d[0,0]).rjust(12) + ','
                                 + ('%.6f'%vert.pos3d[0,1]).rjust(12) + ','
                                 + ('%.6f'%vert.pos3d[0,2]).rjust(12) + ' ]')
        fold.write('\n    ]')
    if edges:
        fold.write(',\n    "edges_vertices": [\n')
        for i, edge in enumerate(shape.edges):
            if i > 0:
                fold.write(',\n')
            end1 = shape.verts.index(edge.end1)
            end2 = shape.verts.index(edge.end2)
            fold.write('      [' + str(end1).rjust(5) + ',' + str(end2).rjust(5) + ' ]')
        fold.write('\n    ],\n    "edges_assignment": [\n')
        for i, edge in enumerate(shape.edges):
            if i > 0:
                fold.write(',\n')
            fold.write('      "' + edge.direction + '"')
        fold.write('    ]')
    if faces:
        fold.write(',\n    "faces_vertices": [\n')
        for j, face in enumerate(shape.faces):
            if j > 0:
                fold.write(',\n')
            fold.write('      [')
            for i, vert in enumerate(face.verts):
                vind = shape.verts.index(vert)
                if i > 0:
                    fold.write(',')
                fold.write(str(vind).rjust(5))
            fold.write(' ]')
        fold.write('    ]')
    fold.write('}')    
    fold.close()


def plotly3d(filename, shape, verts=False, edges=False, faces=True):
    import plotly.graph_objects as go
    import plotly.io as pio
    pio.renderers.default='browser'
    
    data=[]
    
    vert_x, vert_y, vert_z = [],[],[]
    for vert in shape.verts:
        vert_x.append(vert.pos3dvis[0,0])
        vert_y.append(vert.pos3dvis[0,1])
        vert_z.append(vert.pos3dvis[0,2])
    data.append(go.Scatter3d(x=vert_x, y=vert_y, z=vert_z, mode='markers'))
    
    for face in shape.faces:
        face_x, face_y, face_z = [],[],[]
        for i, vert in enumerate(face.verts):
            face_x.append(vert.pos3dvis[0,0] *(1 + i*.00001))
            face_y.append(vert.pos3dvis[0,1] *(1 + i*.00001))
            face_z.append(vert.pos3dvis[0,2])
        #print(face_x, face_y, face_z)
        data.append(go.Mesh3d(x=face_x, y=face_y, z=face_z, 
                               showscale=False, opacity=0.5, color='rgb(200,200,200)',
                               alphahull=0))
    fig = go.Figure(data)
    fig.show()
    
    '''
    a, b, d = 1.32, 1., 0.8
    c = a**2 - b**2
    u, v = np.mgrid[0:2*np.pi:100j, 0:2*np.pi:100j]
    x = (d * (c - a * np.cos(u) * np.cos(v)) + b**2 * np.cos(u)) / (a - c * np.cos(u) * np.cos(v))
    y = b * np.sin(u) * (a - d*np.cos(v)) / (a - c * np.cos(u) * np.cos(v))
    z = b * np.sin(v) * (c*np.cos(u) - d) / (a - c * np.cos(u) * np.cos(v))
    fig = go.Figure(go.Surface(x=x, y=y, z=z, colorbar_x=-0.07))
    fig.show()
    '''
    return

def vpython3d(filename, shape, verts=False, edges=False, faces=True, savefile=False):
    import vpython 
    caption = """Rotate: right-click or Ctrl+left-click and drag
Zoom: middle-click or Alt+left-click and drag
Pan: Shift+left-click and drag"""
    scene = vpython.canvas(userzoom=True, userspin=True, 
                           forward=vpython.vector(0,1,0), up=vpython.vector(0,0,1),
                           center=vpython.vector(0,0,shape.currentheight3d/2),
                           background=vpython.color.white,
                           caption=caption)
    #scene.lights = []
    #vpython.distant_light(direction=vpython.vector( 0.0,  0.0,  1.0), 
    #                      color=vpython.color.gray(0.8))
    lightvector = np.array([[1.0], [0.5], [1.5]])
    lightvector /= np.linalg.norm(lightvector)
    color_norm_scaling = 0.2
    color = []
    if edges:
        for edge in shape.edges:
            if edge.direction == "B" or edge.direction == "M" or edge.direction == "V":
                vpython.curve(pos=[vpython.vec(edge.end1.pos3dvis[0,0],
                                               edge.end1.pos3dvis[0,1],
                                               edge.end1.pos3dvis[0,2]),
                                   vpython.vec(edge.end2.pos3dvis[0,0],
                                               edge.end2.pos3dvis[0,1],
                                               edge.end2.pos3dvis[0,2])],
                              radius=0.01,
                              color=vpython.color.black)
    
    if faces:
        for face in shape.faces:
            color_scaling = np.dot(face.normal, lightvector) * color_norm_scaling
            color.append(face.color + color_scaling)
            if len(face.verts) == 3:
                vpython.triangle(v0=vpython.vertex(pos=vpython.vec(face.verts[0].pos3dvis[0,0],
                                                                   face.verts[0].pos3dvis[0,1],
                                                                   face.verts[0].pos3dvis[0,2]),
                                                   color=vpython.vector(color[-1][0,0],
                                                                        color[-1][0,1],
                                                                        color[-1][0,2]),
                                                   shininess=0),
                                 v1=vpython.vertex(pos=vpython.vec(face.verts[1].pos3dvis[0,0],
                                                                   face.verts[1].pos3dvis[0,1],
                                                                   face.verts[1].pos3dvis[0,2]),
                                                   color=vpython.vector(color[-1][0,0],
                                                                        color[-1][0,1],
                                                                        color[-1][0,2]),
                                                   shininess=0),
                                 v2=vpython.vertex(pos=vpython.vec(face.verts[2].pos3dvis[0,0],
                                                                   face.verts[2].pos3dvis[0,1],
                                                                   face.verts[2].pos3dvis[0,2]),
                                                   color=vpython.vector(color[-1][0,0],
                                                                        color[-1][0,1],
                                                                        color[-1][0,2]),
                                                   shininess=0))
            else:
                vpython.quad(v0=vpython.vertex(pos=vpython.vec(face.verts[0].pos3dvis[0,0],
                                                               face.verts[0].pos3dvis[0,1],
                                                               face.verts[0].pos3dvis[0,2]),
                                               color=vpython.vector(color[-1][0,0],
                                                                    color[-1][0,1],
                                                                    color[-1][0,2]),
                                               shininess=0),
                             v1=vpython.vertex(pos=vpython.vec(face.verts[1].pos3dvis[0,0],
                                                               face.verts[1].pos3dvis[0,1],
                                                               face.verts[1].pos3dvis[0,2]),
                                               color=vpython.vector(color[-1][0,0],
                                                                    color[-1][0,1],
                                                                    color[-1][0,2]),
                                               shininess=0),
                             v2=vpython.vertex(pos=vpython.vec(face.verts[2].pos3dvis[0,0],
                                                               face.verts[2].pos3dvis[0,1],
                                                               face.verts[2].pos3dvis[0,2]),
                                               color=vpython.vector(color[-1][0,0],
                                                                    color[-1][0,1],
                                                                    color[-1][0,2]),
                                               shininess=0),
                             v3=vpython.vertex(pos=vpython.vec(face.verts[3].pos3dvis[0,0],
                                                               face.verts[3].pos3dvis[0,1],
                                                               face.verts[3].pos3dvis[0,2]),
                                               color=vpython.vector(color[-1][0,0],
                                                                    color[-1][0,1],
                                                                    color[-1][0,2]),
                                               shininess=0))
    if savefile:
        scene.capture(filename)
    return
    
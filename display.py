# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 08:18:03 2023

@author: gieseking
"""

#import elements as el
import model
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import numpy as np
from array import array

def write_svg_cp(filename, xdim, ydim, edgelist):
    # Layout notes:
    # 1. All points are scaled by a factor of 100 (a gore width of 1 is 100 px)
    # 2. SVG files use the top left corner as (0,0) instead of the bottom left.
    #    Y dimensions of all points are translated accordingly.
    svg = open(filename,'w')
    xdim = str(xdim*100)
    ydim = str(ydim*100)
    svg.write('''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg 
  xmlns="http://www.w3.org/2000/svg"
  id="svg8"
  version="1.1" 
  height="''' + ydim + '''mm" 
  width="''' + xdim + '''mm">
<g>
''')
    #svg.write('<rect x="0" y="0" height="10" width="10" ' +
    #          'style="fill:none;stroke:#000000;stroke-width:1">\n')
    for edge in edgelist:
        x1, y1 = str(edge.end1.pos2d[0,0] * 100), str(float(ydim) - edge.end1.pos2d[0,1] * 100)
        x2, y2 = str(edge.end2.pos2d[0,0] * 100), str(float(ydim) - edge.end2.pos2d[0,1] * 100)
        svg.write('  <line x1="' + x1 + '" y1="' + y1 + '" x2="' + x2 + '" y2="' + y2 + '" ') 
        if edge.direction == 'B':
            svg.write('stroke="black" stroke-width="3"')
        elif edge.direction == 'M':
            svg.write('stroke="red" stroke-width="2"')
        elif edge.direction == 'V':
            svg.write('stroke="blue" stroke-width="2"')
        else:
            svg.write('stroke="gray" stroke-width="0.5"')
        svg.write(' />\n')
        pass
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
        x3d.append(vertex.pos3d[0,0])
        y3d.append(vertex.pos3d[0,1])
        z3d.append(vertex.pos3d[0,2])
        if verts:
            ax.scatter(vertex.pos3d[0,0],vertex.pos3d[0,1],vertex.pos3d[0,2], color='k')
        else:
            ax.scatter(vertex.pos3d[0,0],vertex.pos3d[0,1],vertex.pos3d[0,2], color='k', s=0.0)
    
    if faces:
        lightvector = np.array([[-1.0], [0.5], [-1.5]])
        lightvector /= np.linalg.norm(lightvector)
        ls = mpl.colors.LightSource(azdeg=110, altdeg=45)

        for face in shape.faces:
            verts = np.zeros((len(face.verts), 3))
            #print(verts)
            for i in range(len(face.verts)):
                #print(face.verts[i].pos3d[0])
                verts[i, :] = face.verts[i].pos3d[0]
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
            edge3d = [edge.end1.pos3d,edge.end2.pos3d]
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
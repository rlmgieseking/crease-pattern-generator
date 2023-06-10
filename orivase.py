# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 16:51:09 2023

@author: gieseking
"""

import utils.readinp as readinp
from datetime import datetime

t = []
t.append(datetime.now())
# Set up file names
basefile = 'examples/vase3'
if '.' in basefile:
    basefile = basefile[:basefile.index('.')]
inpfile = basefile+'.txt'

# Read the input
blocks = readinp.readfile(inpfile)
t.append(datetime.now())
print("Input read, time = ", t[-1] - t[-2])

# Generate the model
shape = readinp.generatemodel(blocks)
t.append(datetime.now())
print("Model generated, time = ", t[-1] - t[-2])

# Plot the output
readinp.displaymodel(shape, blocks, basefile)
t.append(datetime.now())
print("Display done, time = ", t[-1] - t[-2])
print("Total time = ", t[-1] - t[0])
"""
display.vpython3d(shape, edges=True)
'''
display.plotly3d(basefile+'_3D.html', shape)
'''

display.plot2d(shape, verts=False)
t.append(datetime.now())
print("2D plot done, time = ", t[-1] - t[-2])
'''
display.plot3d(basefile+'_3D.png', shape, verts=False, edges=False, faces=True)
t.append(datetime.now())
print("3D plot done, time = ", t[-1] - t[-2])
'''
display.write_svg_cp(basefile+'.svg', shape.ngores*shape.gorewidth, 
                     shape.currentheight2d, shape.edges)
t.append(datetime.now())
print("SVG done, time = ", t[-1] - t[-2])
print("Total time = ", t[-1] - t[0])
"""

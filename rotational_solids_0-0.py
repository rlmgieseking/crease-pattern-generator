# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 16:51:09 2023

@author: gieseking
"""

# Version 0.0  Replicate the functionality of OriRevo, Lang's Rotational Solids, etc.
import display 
import readinp
from datetime import datetime

start = datetime.now()
# Set up file names
basefile = 'vase3'
if '.' in basefile:
    basefile = basefile[:basefile.index('.')]
inpfile = basefile+'.txt'

# Read the input
blocks = readinp.readfile(inpfile)
t1 = datetime.now()
print("Input read, time = ", t1 - start)

# Generate the model
shape = readinp.generatemodel(blocks)
t2 = datetime.now()
print("Model generated, time = ", t2 - t1)

# Plot the output
display.plot2d(shape, verts=True)
t3 = datetime.now()
print("2D plot done, time = ", t3 - t2)

display.plot3d(basefile+'_3D.png', shape, verts=True, edges=True)
t4 = datetime.now()
print("3D plot done, time = ", t4 - t3)

display.write_svg_cp(basefile+'.svg', shape.ngores*shape.gorewidth, shape.currentheight2d, shape.edges)
t5 = datetime.now()
print("SVG done, time = ", t5 - t4)
print("Total time = ", t5 - start)


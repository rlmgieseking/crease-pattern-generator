'''
Welcome to OriVase 0.1
Rebecca Gieseking, 2023

Main script
See instructions for use in README.md
'''

import utils.readinp as readinp
from datetime import datetime
import sys

t = []
t.append(datetime.now())

# Program requires a text input file (details in README.md)
# Get input file name
if len(sys.argv) <= 1:
    #inpfile = input("Enter the OriVase input file name: ")
    inpfile = 'examples/diagonal_shift_vase.txt'
else:
    inpfile = sys.argv[1]

if '.' in inpfile:
    basefile = inpfile[:inpfile.rfind('.')]
else:
    basefile = inpfile

# Read the input from the file and convert into blocks of information
blocks = readinp.readfile(inpfile)
t.append(datetime.now())
print("Input read, time = ", t[-1] - t[-2])

# Generate the model using the blocks from the input file
# Construct the 2D and 3D location of points and lines
shape = readinp.generatemodel(blocks)
t.append(datetime.now())
print("Model generated, time = ", t[-1] - t[-2])

# Print verts (temp)
#for v in shape.verts:
#    print(v.pos2d[0,0], v.pos2d[0,1], v.pos3d[0,0], v.pos3d[0,1], v.pos3d[0,2])

# Convert the 2D/3D positions of points and lines into the requested outputs
readinp.displaymodel(shape, blocks, basefile)
t.append(datetime.now())
print("Display done, time = ", t[-1] - t[-2])
print("Total time = ", t[-1] - t[0])

# Welcome to Orivase 0.1!
Author: Rebecca Gieseking

Orivase is a Python 3 script that generates crease patterns and 3D visualizations for origami vase forms. As of Version 0.1, the script can create forms based on cones (similar to the capabilities of the ORI-REVO package) and diagonal shifts.

## Input file format
The script reads a text-based input file and uses the information in that file to generate the output. The input format is fairly flexible about spacing, capitalization, number of lines, etc. 

### Comments:
Comments start with the symbol #. Any text from the # to the end of the line will be treated as a comment and discarded when the script reads the input file. The # may appear at the beginning or middle of any line.

### Setup block:
In most cases, the input file should start with a setup block. If not, the default setup parameters will be used. An example of a setup block is:

>setup

>  ngores     4

>  radius     0.5

>end 

The setup block must start with the word "setup" and end with the word "end". If there is text following "end" on the same line, it will be ignored. 

Within the setup block, there are several options that can be set:

| Option      | Default | Description |
| ------      | ------- | ----------- |
| ngores      | 8       | Number of gores/facets/sides the model has (must be an integer) | 
| gorewidth   | 1.0     | Width of each gore, in cm |
| radius      | 0.0     | Initial inner radius of the cylinder |
| cwrot       | False   | Initial direction of flanges surrounding the model, as viewed from above (counterclockwise by default, clockwise if cwrot = True) |
| relativedim | True    | Sets whether later blocks should use relative dimensions (defined such that 1.0 is the  |




# Welcome to Orivase 0.1!
Author: Rebecca Gieseking

Orivase is a Python 3 script that generates crease patterns and 3D visualizations for origami vase forms. As of Version 0.1, the script can create forms based on cones (similar to the capabilities of the ORI-REVO package) and diagonal shifts.

# Setup and installation

You must have Python 3 installed to run this script.

First, follow GitHub's instructions for cloning a repository (https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) to create a copy of the code on your own computer.

Once you have a local copy, use the command line to navigate to the `crease-pattern-generator` base directory and use the command:

```
pip install -r requirements.txt
```

to install the packages that are needed to run the script.

# Running the script

As of version 0.1, the script requires a text-based input file (described below). Once you have an input file, you can run the script from the `crease-pattern-generator` base directory using the following command:

```
python orivase.py <filename.txt>
```

where `<filename.txt>` is replaced by the path and name of your input file.

## Input file format

The script reads a text-based input file and uses the information in that file to generate the output. The input format is flexible about spacing, number of options per line, etc. The input is not case-sensitive. Several example input files are included in the `examples` directory.

The input interpreter is set up to correct errors and give a final model as often as possible. If an invalid value is given (for example, a non-numerical value for an option that needs a number), the script will use the default value of that option. If the same option is defined twice within the same block, the script will only use the value that is defined first.

### Comments

Comments start with the symbol `#`. Any text from the `#` to the end of the line will be treated as a comment and discarded when the script reads the input file. The `#` may appear at the beginning or middle of any line.

### Setup block

In most cases, the input file should start with a setup block. If not, the default setup parameters will be used. An example of a setup block is:

```
setup
  ngores     4
  radius     0.5
end 
```

The setup block must start with the command `setup` and end with the command `end`. If there is text following `end` on the same line, it will be treated as a comment. I find it easier to read the input file if there are spaces before each option, but they are not necessary for the code to work. It is also possible to put more than one option on the same line, as long as there is at least one space between them.

Within the setup block, the possible options are:

| Option      | Default | Description |
| ------      | ------- | ----------- |
| ngores      | 8       | Number of gores/facets/sides the model has. Must be an integer. | 
| gorewidth   | 1.0     | Width of each gore, in cm. |
| radius      | 0.0     | Initial inner radius of the model, usually at the bottom edge. By default, defined as a fraction of the maximum possible radius. If relativedim = False, defined in cm. |
| cwrot       | False   | Initial direction of flanges surrounding the model, as viewed from above. Counterclockwise by default, clockwise if cwrot = True. |
| relativedim | True    | Sets whether later blocks should use relative dimensions (defined as a fraction of the maximum possible radius given the number of gores and their width) or absolute dimensions (cm).  |


### Model-building blocks

There are several types of blocks that can be used to add various components to the model. The code was structured with the intention of building models from bottom to top, so a positive height generally goes up in 3D space.

Each component block must start with the command `add <component>`, where `<component>` is one of the options below. Like the setup block, the component block must end with the command `end`. If there is text following `end` on the same line, it will be treated as a comment. 

#### Cones

Cones are the basic building block for rotationally symmetric solids. The 'cones' created by this code are actually truncated regular pyramids with `ngores` sides (where `ngores` was defined in the setup block), or prisms if the initial and final radius are the same. The radius at the initial edge (usually the bottom edge) of the cone is the radius from the setup block if this is the first component, or the radius from the previous component if this is not the first component. Because of this, only the height of the cone and its final radius are defined in this block. 

An example of a cone block is:

```
add cone
  height     1.0
  radius     0.5
end
```

Within the block, the possible options are:

| Option      | Default | Description |
| ------      | ------- | ----------- |
| height      | 0.0     | Height of the cone, from bottom edge to top edge. Positive values are up. Defined either as a fraction of the maximum possible radius or in cm, depending on `relativedim` from the input block. |
| radius      | 0.0     | Ending radius of the cone. Defined either as a fraction of the maximum possible radius or in cm, depending on `relativedim` from the input block. |
| startedge   | auto    | Fold format for the first edge of the cone. If `auto`, the code will use the relative angle between this cone and the previous component to assign the fold as a mountain, valley, or rule line (flat). If 'none', the fold will be assigned as a rule line; this can be helpful if you are constructing an approximate curve from a series of cones. |

#### Curved cones

Curved cones allow you to approximate a curved form using a series of cones. As for a single cone, the initial radius is set by default, and the final height and radius must be defined. In addition, you must define the shape of the curve. The curve is constructed by defining the point where the the paper is parallel to the line from the start point to the end point of the curved cone, called the vertex. Both the height and the radius of the vertex must be defined. If the vertex is inset from the start-to-end line, the entire curve is concave; if it is further out, the entire curve is convex. Constructing an S-like curve requires one curved cone for the concave part and a second for the convex part.

The curve is constructed based on two parabolas: one from the vertex to the start point, and one from the vertex to the end point. For both parabolas, the vertex of the parabola is the point defined as the vertex of the curved cone, and the directrix is parallel to the to the line from the start point to the end point of the curved cone. The width of each parabola is set to ensure that it passes through the start point or the end point as appropriate. 

An example of a curved cone block is:

```
add curvedcone
  height     2.0
  radius     0.5
  vertexheight  1.5
  vertexradius  1.0
  nsegments     8
end
```
Within the block, the possible options are:

| Option      | Default | Description |
| ------      | ------- | ----------- |
| height      | 0.0     | Height of the cone, from bottom edge to top edge. Positive values are up. Defined either as a fraction of the maximum possible radius or in cm, depending on `relativedim` from the input block. |
| radius      | 0.0     | Ending radius of the cone. Units are the same as above. |
| vertexheight | 0.0    | Height of the vertex of the curve. Must be between 0 and `height`. Units are the same as above. |
| vertexradius | 0.0    | Radius of the vertex of the curve. Units are the same as above. |
| nsegments   | 4       | Number of cones used to construct each half of the curve. For the default of 4, a total of 8 cones will be used. |
| startedge   | auto    | Fold format for the first edge of the cone. If `auto`, the code will use the relative angle between this cone and the previous component to assign the fold as a mountain, valley, or rule line (flat). If 'none', the fold will be assigned as a rule line; this can be helpful if you are constructing an approximate curve from a series of cones. |

#### Diagonal shifts

Diagonal shifts are a component that I designed around 2013 that create the appearance of slicing a cylinder along a diagonal and shifting the top half of the cylinder uphill along that diagonal, [see examples on my website](http://rebecca.gieseking.us/2013/08/test-models-diagonal-shift-part-2/). The diagonal shift is folded essentially as a twist, where the paper rotates (`ngores`/2 - 1)/`ngores` * 360 degrees. As `ngores` increases, the rotation angle approaches 180 degrees. In the center of the twist, `ngores` pleats come together to (approximately) one point, which I call the convergence point. This twist also causes the flanges to switch between clockwise and counterclockwise.

The angle of the diagonal and the distance of the shift are linked and cannot be changed independently. In this code, the extent of the shift is set up `offsetfract`, which is how far the convergence point is from the center of the cylinder. The height of the diagonal shift block is the minimum height needed to contain the shift set by the other options.

An example of a diagonal shift block is:

```
add diagshift
  radius       1.0
  offsetfract  0.4
  tiltrot      0
end
```

Within the block, the possible options are:

| Option      | Default | Description |
| ------      | ------- | ----------- |
| radius      | 0.0     | Ending radius of the shift. Defined either as a fraction of the maximum possible radius or in cm, depending on `relativedim` from the input block. |
| offsetfract | 0.0    | Offset of the convergence point of the twist from the center of the cylinder, as a fraction of the cylinder width. Must be between 0.0 (for a flat twist) and 0.8. An offset larger than 0.8 would give such a large tilt angle of the diagonal plane that it would be impractical to fold. |
| tiltrot     | 0.0    | Rotation of the convergence point in the xy plane, in degrees. |
| startedge   | auto    | Fold format for the first edge of the cone. If `auto`, the code will use the relative angle between this cone and the previous component to assign the fold as a mountain, valley, or rule line (flat). If 'none', the fold will be assigned as a rule line; this can be helpful if you are constructing an approximate curve from a series of cones. |


### Display and output

By default, the script creates 5 different outputs, two including the 2D crease pattern and three including the 3D folded form. The output types are:

1. `svpcp`: Crease pattern in .svg format, which is a common vector graphics format
2. `vpy`: Interactive 3D VPython visualization that can be rotated, zoomed, etc., which opens in your default browser
3. `vpy_script`: Prints a VPython script to visualize the 3D model that can be used at https://www.glowscript.org/ to create a web-based 3D interactive visualization
4. `fold2d`: Crease pattern in .fold format (see https://github.com/edemaine/fold for details)
5. `fold3d`: 3D folded form in .fold format

The specific options are:

| Option      | Default | Description |
| ------      | ------- | ----------- |
svgcp         |  True   | If True, prints a crease pattern in .svg format
svgcp_verts   |  False  | If True, .svg crease pattern includes points for each vertex
svgcp_edges   |  True   | If True, .svg crease pattern includes borders (black) and folded edges (mountain = red, valley = blue)
svgcp_rules   |  True   | If True, .svg crease pattern includes rule lines as light gray lines
vpy           |  True   | If True, generates an interactive 3D visualization of the folded model using VPython, which will open in your default browser. This visualization can be rotated, zoomed, etc.
vpy_script    |  False  | If True, prints a VPython script that can be used at https://www.glowscript.org/ to create a web-based 3D interactive visualization
vpy_verts     |  False  | If True, 3D visualization includes points for each vertex
vpy_edges     |  True   | If True, 3D visualization includes black lines for each border and folded edge
vpy_faces     |  True   | If True, 3D visualization includes colored/shaded faces
vpy_image     |  True   | If True, 3D visualization also saves a .png image file of the default view
fold2d        |  True   | If True, prints a .fold file of the 2D crease pattern
fold2d_verts  |  True   | If True, 2D .fold file includes vertex coordinates
fold2d_edges  |  True   | If True, 2D .fold file includes vertex indices and fold assignments (mountain, valley, etc.) for each edge
fold2d_faces  |  True   | If True, 2D .fold file includes vertex indices for each face
fold3d        |  True   | If True, prints a .fold file of the 3D folded form
fold3d_verts  |  True   | If True, 3D .fold file includes vertex coordinates
fold3d_edges  |  True   | If True, 3D .fold file includes vertex indices and fold assignments (mountain, valley, etc.) for each edge
fold3d_faces  |  True   | If True, 2D .fold file includes vertex indices for each face
                  
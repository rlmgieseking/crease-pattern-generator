# Welcome to Orivase 0.1!
Author: Rebecca Gieseking

Orivase is a Python 3 script that generates crease patterns and 3D visualizations for origami vase forms. As of Version 0.1, the script can create forms based on cones (similar to the capabilities of the ORI-REVO package) and diagonal shifts.

## Input file format
The script reads a text-based input file and uses the information in that file to generate the output. The input format is fairly flexible about spacing, capitalization, number of lines, etc. 

### Comments:
Comments start with the symbol #. Any text from the # to the end of the line will be treated as a comment and discarded when the script reads the input file. The # may appear at the beginning or middle of any line.

### Setup block:
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

Cones are the basic building block for rotationally symmetric solids. The 'cones' created by this code are actually truncated regular pyramids with ngores sides (where ngores was defined in the setup block), or prisms if the initial and final radius are the same. The radius at the initial edge (usually the bottom edge) of the cone is the radius from the setup block if this is the first component, or the radius from the previous component if this is not the first component. Because of this, only the height of the cone and its final radius are defined in this block. 

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

Curved cones allow you to approximate a curved form using a series of cones. As for a single cone, the initial radius is set by default, and the final height and radius must be defined. In addition, you must define the shape of the curve. The curve is constructed by defining the point where the the paper is parallel to the line from the start point to the end point of the curved cone, called the vertex. Both the height and the radius of the vertex must be defined. 

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

# Diagonal shifts

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





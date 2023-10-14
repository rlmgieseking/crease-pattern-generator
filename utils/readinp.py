# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:28:19 2023

@author: gieseking
"""

import utils.model as model
import utils.display as display
import re

# Define variables available for each block and their default values
setupvars =      [['ngores',         int,   8], 
                  ['gorewidth',      float, 1.0], 
                  ['radius',         float, 0.0], 
                  ['cwrot',          bool,  False],
                  ['overlap',        int,   0],
                  ['relativedim',    bool,  True]]
conevars =       [['height',         float, 0.0],
                  ['radius',         float, 0.0],
                  ['startedge',      str,   'auto']]
curvedconevars = [['height',         float, 0.0],
                  ['radius',         float, 0.0],
                  ['startedge',      str,   'auto'],
                  ['vertexheight',   float, 0.0],
                  ['vertexradius',   float, 0.0],
                  ['nsegments',      int,   4]]
diagshiftvars =  [['radius',         float, 0.0],
                  ['offsetfract',    float, 0.0],
                  ['tiltrot',        float, 0.0],
                  ['startedge',      str,   'auto']]
displayvars =    [['svgcp',          bool,  True],
                  ['svgcp_verts',    bool,  False],
                  ['svgcp_edges',    bool,  True],
                  ['svgcp_rules',    bool,  True],
                  ['vpy',            bool,  True],
                  ['vpy_script',     bool,  False],
                  ['vpy_verts',      bool,  False],
                  ['vpy_edges',      bool,  True],
                  ['vpy_faces',      bool,  True],
                  ['vpy_image',      bool,  True],
                  ['fold2d',         bool,  True],
                  ['fold2d_verts',   bool,  True],
                  ['fold2d_edges',   bool,  True],
                  ['fold2d_faces',   bool,  True],
                  ['fold3d',         bool,  True],
                  ['fold3d_verts',   bool,  True],
                  ['fold3d_edges',   bool,  True],
                  ['fold3d_faces',   bool,  True]]

def getvar(block, varname, vartype, vardefault):
    var = None
    if varname in block:
        if block.index(varname) == len(block) - 1:
            var = vardefault
            print('Error: Value could not be found for',varname)
        else:
            try:
                if vartype == bool:
                    var = eval(block[block.index(varname)+1].capitalize())
                else:
                    var = vartype(block[block.index(varname)+1])
                block.pop(block.index(varname)+1)
                block.pop(block.index(varname))
            except:
                print('Error:', varname, 'value of', block[block.index(varname)+1],
                      'is not a', vartype, '. Using the default value of', vardefault)
                var = vardefault
    else:
        var = vardefault
    return block, var

def readfile(inpfile):
    # Read input file
    ifile = open(inpfile, 'r')
    #iline = ifile.readline().lower()
    itext = ifile.read().lower().splitlines()
    
    blocks = []
    #print(itext)
    # Process text
    # Remove comments 
    for i in range(len(itext)-1,-1, -1):
        if itext[i].find('#') > -1:
            itext[i] = itext[i][:itext[i].find('#')]
        if len(itext[i]) == 0:
            itext.pop(i)
    #print(itext)
    
    # Split text into blocks
    while len(itext) > 0:
        block = []
        while 'end' not in block and len(itext) > 0:
            block.extend(re.split(" |,|=",itext.pop(0)))
            #line = line.split()
            #for word in line:
            #    block.append(word) 
        block= list(filter(None,block))
        #print(block)
        if block.index('end') != len(block) - 1:
            print("Error: Variables after end statement will be ignored. Ignoring ", 
                  block[block.index('end')+1:])    
        block = block[:block.index('end')]
        blocks.append(block)
        #itext.pop(0)
    #print(blocks)
    return blocks

def generatemodel(blocks):
    modelexists = False
    for block in blocks:
        if 'setup' in block:
            block.remove('setup')
            if modelexists:
                print('Error: Setup block found after model started. Ignoring setup block',
                      block)
            else:
                #print(setupvars)
                for [varname, vartype, vardefault] in setupvars:
                    #print(varname, vartype)
                    block, var = getvar(block, varname, vartype, vardefault)
                    #print(varname, var)
                    exec("%s = %r" % (varname, var), globals())
                #print(ngores)
                if len(block) > 0:
                    print('Error: Variables', block, 'could not be interpreted. Ignoring these values.')
                shape = model.Model(ngores=ngores, 
                                    gorewidth=gorewidth, 
                                    cwrot=cwrot, 
                                    currentradius=radius, 
                                    relativedim=relativedim,
                                    overlap=overlap)
                modelexists = True
        elif 'add' in block:
            block.remove('add')
            if not modelexists:
                print('No setup block before components are added. Using default setup parameters.')
                for [varname, vartype, vardefault] in setupvars:
                    #print(varname, vartype)
                    block, var = getvar([], varname, vartype, vardefault)
                    exec("%s = %r" % (varname, var), globals())
                shape = model.Model(ngores=ngores, 
                                    gorewidth=gorewidth, 
                                    cwrot=cwrot, 
                                    currentradius=radius, 
                                    relativedim=relativedim,
                                    overlap=overlap)
                modelexists = True
            if 'cone' in block:
                block.remove('cone')
                for [varname, vartype, vardefault] in conevars:
                    #print(varname, vartype)
                    block, var = getvar(block, varname, vartype, vardefault)
                    exec("%s = %r" % (varname, var), globals())
                if len(block) > 0:
                    print('Error: Variables', block, 'could not be interpreted. Ignoring these values.')
                shape.add_cone(radius, height, startedge)
            elif 'curvedcone' in block:
                block.remove('curvedcone')
                for [varname, vartype, vardefault] in curvedconevars:
                    #print(varname, vartype)
                    block, var = getvar(block, varname, vartype, vardefault)
                    exec("%s = %r" % (varname, var), globals())
                if len(block) > 0:
                    print('Error: Variables', block, 'could not be interpreted. Ignoring these values.')
                shape.add_curvedcone(radius, height, vertexradius, vertexheight,
                                     nsegments, startedge)
            elif 'diagshift' in block:
                block.remove('diagshift')
                for [varname, vartype, vardefault] in diagshiftvars:
                    #print(varname, vartype, vardefault)
                    block, var = getvar(block, varname, vartype, vardefault)
                    exec("%s = %r" % (varname, var), globals())
                if len(block) > 0:
                    print('Error: Variables', block, 'could not be interpreted. Ignoring these values.')
                shape.add_diagshift(radius, offsetfract, tiltrot, startedge)
            elif 'display' in block:
                pass
            else:
                print('Error: Block has undetermined component type. Ignoring block', block)
    return shape

def displaymodel(shape, blocks, basefile):
    # Set variables to defaults
    for [varname, vartype, vardefault] in displayvars:
        #print(varname, vartype)
        block, var = getvar([], varname, vartype, vardefault)
        exec("%s = %r" % (varname, var), globals())
    # Look for display block
    for block in blocks:
        if 'display' not in block:
            continue
        block.remove('display')
        for [varname, vartype, vardefault] in displayvars:
            #print(varname, vartype)
            block, var = getvar(block, varname, vartype, vardefault)
            exec("%s = %r" % (varname, var), globals())
        if len(block) > 0:
            print('Error: Variables', block, 'could not be interpreted. Ignoring these values.')
        
    if svgcp:
        display.write_svg_cp(basefile+'.svg', shape, verts=svgcp_verts, edges=svgcp_edges, rules=svgcp_rules)
    if vpy:
        display.vpython3d(basefile+'.png', shape, verts=vpy_verts, edges=vpy_edges, faces=vpy_faces, savefile=vpy_image)
    if vpy_script:
        display.glowscript(basefile+'_vpy.py', shape, verts=vpy_verts, edges=vpy_edges, faces=vpy_faces)
    if fold2d:
        display.fold2d(basefile+'_2D.fold', shape, verts=fold2d_verts, edges=fold2d_edges, faces=fold2d_faces)
    if fold3d:
        display.fold3d(basefile+'_3D.fold', shape, verts=fold3d_verts, edges=fold3d_edges, faces=fold3d_faces)
        
    
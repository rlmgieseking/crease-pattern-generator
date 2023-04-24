# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:28:19 2023

@author: gieseking
"""

import model

# Define variables available for each block and their default values
setupvars = [['ngores',    int,   8], 
             ['gorewidth', float, 1.0], 
             ['radius',    float, 0.0], 
             ['cwrot',     bool,  False]]
conevars =  [['height',    float, 0.0],
             ['radius',    float, 0.0],
             ['startedge', str,   'auto']]
curvedconevars = [['height',    float, 0.0],
                  ['radius',    float, 0.0],
                  ['startedge', str,   'auto'],
                  ['vertexheight',    float, 0.0],
                  ['vertexradius',    float, 0.0],
                  ['nsegments',       int,   2]]
diagshiftvars =  [['radius',      float, 0.0],
                  ['offsetfract', float, 0.0],
                  ['tiltrot',     float, 0.0],
                  ['startedge',   str,   'auto']]

def getvar(block, varname, vartype, vardefault):
    var = None
    if varname in block:
        try:
            if vartype == bool:
                var = eval(block[block.index(varname)+1].capitalize())
            else:
                var = vartype(block[block.index(varname)+1])
            block.pop(block.index(varname)+1)
            block.pop(block.index(varname))
        except:
            print('Error:', varname, block[block.index(varname)+1],
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
    # print(itext)
    
    # Split text into blocks
    while len(itext) > 0:
        block = itext[0].split()
        while 'end' not in itext[0] and len(itext) > 0:
            itext.pop(0)
            line = itext[0].split()
            for word in line:
                block.append(word) if word != '=' else None
        if block.index('end') != len(block) - 1:
            print("Error: Variables after end statement will be ignored. Ignoring ", 
                  block[block.index('end')+1:])    
        block = block[:block.index('end')]
        blocks.append(block)
        itext.pop(0)
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
                shape = model.Model(ngores, gorewidth, cwrot, radius)
                modelexists = True
        elif 'add' in block:
            block.remove('add')
            if not modelexists:
                print('No setup block before components are added. Using default setup parameters.')
                for [varname, vartype, vardefault] in setupvars:
                    #print(varname, vartype)
                    block, var = getvar([], varname, vartype, vardefault)
                    exec("%s = %r" % (varname, var), globals())
                shape = model.Model(ngores, gorewidth, cwrot, radius)
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
            else:
                print('Error: Block has undetermined component type. Ignoring block', block)
    return shape
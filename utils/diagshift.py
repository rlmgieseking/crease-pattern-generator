# -*- coding: utf-8 -*-
"""
DiagonalShift class
"""

from utils.molecule import Molecule
import math
import utils.elements as el
import numpy as np
from copy import deepcopy

# tiltangle = angle from horizontal
# tiltrot = angle in the xy plane where the tilt plane is at the start height
class DiagonalShift(Molecule):
    def __init__(self, ngores, gorewidth, cwrot, 
                 startradius, startheight2d, startheight3d, startverts, startedges,
                 incomingangle, xyrot, axis, center, startedgetype, 
                 endradius, offsetfract, tiltrot):
        Molecule.__init__(self, ngores, gorewidth, cwrot, 
                          startradius, startheight2d, startheight3d, startverts, startedges,
                          incomingangle, xyrot, axis, center, startedgetype)

        self.endradius = endradius
        self.offsetfract = offsetfract
        self.tiltrot = tiltrot * math.pi/180
        self.angle = math.pi/2
        #print('angle', self.angle)
        #self.tiltangle = math.atan(offsetfract / math.sqrt(1 - offsetfract**2))
        #print('tilt', self.tiltangle)
        
        if self.incomingangle != None:
            self.anglediff = self.angle - self.incomingangle
        else:
            self.anglediff = None
        
        self.generateverts()
        self.generateedges()
        self.generatefaces()
    
    def generateverts(self):
        ####
        # Generate or load bottom edge points
        ####
        visshift = 0.01
        if self.startverts and len(self.startverts) == self.ngores * 2 + 1:
            self.verts = self.startverts
            #print('startverts')
            #for vert in self.startverts:
            #    print(vert.pos2d, vert.pos3d)
            #print('end')
        else:
            if self.cwrot:
                startouter3d = [self.startradius, -self.gorewidth/2, self.startheight3d]
                self.verts.append(el.Vertex([0.0, self.startheight2d], startouter3d))
                self.verts[-1].rotate(-math.pi * 2 / self.ngores, 'z') 
                startinner3d = [self.startradius, 
                                self.startradius * math.tan(math.pi/self.ngores), 
                                self.startheight3d]
                startinner3dvis = [self.startradius - visshift, 
                                   self.startradius * math.tan(math.pi/self.ngores), 
                                   self.startheight3d]
                self.verts.append(el.Vertex([self.gorewidth/2 
                                             - self.startradius * math.tan(math.pi/self.ngores), 
                                             self.startheight2d], 
                                            startinner3d,
                                            startinner3dvis))
            else: 
                startouter3d = [self.startradius, self.gorewidth/2, self.startheight3d]
                self.verts.append(el.Vertex([0.0, self.startheight2d], startouter3d))
                startinner3d = [self.startradius, 
                                -self.startradius * math.tan(math.pi/self.ngores), 
                                self.startheight3d]
                startinner3dvis = [self.startradius - visshift, 
                                   -self.startradius * math.tan(math.pi/self.ngores), 
                                   self.startheight3d]
                self.verts.append(el.Vertex([self.gorewidth/2 
                                             + self.startradius * math.tan(math.pi/self.ngores), 
                                             self.startheight2d], 
                                            startinner3d,
                                            startinner3dvis))
            self.verts[-2].rotate(self.xyrot, 'z')
            self.verts[-1].rotate(self.xyrot, 'z')
            
            # Replicate verts for each gore, with the correct rotations and translations
            for i in range(1, self.ngores):
                self.verts.append(deepcopy(self.verts[-2]))
                self.verts[-1].translate2d([self.gorewidth, 0.0])
                self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
                #print('vertex', self.verts[-1].pos3d)
                
                self.verts.append(deepcopy(self.verts[-2]))
                self.verts[-1].translate2d([self.gorewidth, 0.0])
                self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
                #print('vertex', self.verts[-1].pos3d)
                
            # Add final edge vertex
            self.verts.append(deepcopy(self.verts[-2]))
            self.verts[-1].translate2d([self.gorewidth, 0.0])
            self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
            
            # Translate start edge based on center 
            for vert in self.verts:
                vert.translate3d([self.center[0],self.center[1], 0.0])
        
        ####
        # Generate top edge vertices, 3D height for now is arbitrary
        ####
        topverts = []
        if not self.cwrot:
            endouter3d = [self.endradius, -self.gorewidth/2, 0]
            topverts.append(el.Vertex([0.0, self.startheight2d], endouter3d))
            topverts[-1].rotate(-math.pi * 2 / self.ngores, 'z') 
            topverts[-1].rotate(-2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
            endinner3d = [self.endradius, 
                          self.endradius * math.tan(math.pi/self.ngores), 
                          0]
            endinner3dvis = [self.endradius - visshift, 
                             self.endradius * math.tan(math.pi/self.ngores), 
                             0]
            topverts.append(el.Vertex([self.gorewidth/2 
                                       - self.endradius * math.tan(math.pi/self.ngores), 
                                       self.startheight2d], 
                                      endinner3d,
                                      endinner3dvis))
            topverts[-1].rotate(-2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
        else: 
            endouter3d = [self.endradius, self.gorewidth/2, 0]
            topverts.append(el.Vertex([0.0, self.startheight2d], endouter3d))
            topverts[-1].rotate(2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
            endinner3d = [self.endradius, 
                          -self.endradius * math.tan(math.pi/self.ngores), 
                          0]
            endinner3dvis = [self.endradius - visshift, 
                             -self.endradius * math.tan(math.pi/self.ngores), 
                             0]
            topverts.append(el.Vertex([self.gorewidth/2 
                                       + self.endradius * math.tan(math.pi/self.ngores), 
                                      self.startheight2d], 
                                     endinner3d,
                                     endinner3dvis))
            topverts[-1].rotate(2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
        topverts[-2].rotate(self.xyrot, 'z')
        topverts[-1].rotate(self.xyrot, 'z')
        # Replicate verts for each gore, with the correct rotations and translations
        for i in range(1, self.ngores):
            topverts.append(deepcopy(topverts[-2]))
            topverts[-1].translate2d([self.gorewidth, 0.0])
            topverts[-1].rotate(math.pi * 2 / self.ngores, 'z')
            #print('vertex2', topverts[-1].pos3d)
            
            topverts.append(deepcopy(topverts[-2]))
            topverts[-1].translate2d([self.gorewidth, 0.0])
            topverts[-1].rotate(math.pi * 2 / self.ngores, 'z')
            #print('vertex2', topverts[-1].pos3d)
            
        # Add final edge vertex
        topverts.append(deepcopy(topverts[-2]))
        topverts[-1].translate2d([self.gorewidth, 0.0])
        topverts[-1].rotate(math.pi * 2 / self.ngores, 'z')
        
        # Translate top edge in the xy plane based on offset
        dx = (self.offsetfract * (self.startradius + self.endradius) * math.cos(self.tiltrot) 
              + self.center[0])
        dy = (self.offsetfract * (self.startradius + self.endradius) * math.sin(self.tiltrot) 
              + self.center[1])
        topcenter = [dx, dy]
        for vert in topverts:
            vert.translate3d([dx,dy, 0.0])
        
        ####
        # Set up convergence point xy position; z will be modified later
        ####
        if self.offsetfract > 0.8:
            print('Diagonal shift offset fraction of ', self.offsetfract, ' is too large')
            print('Resetting offset fraction to 0.8')
            self.offsetfraction = 0.8
        offset = self.offsetfract * self.startradius
        h = 2 * self.startradius * offset / math.sqrt(self.startradius**2 - offset**2)
        conv3d = np.array([[math.cos(self.tiltrot) * offset + self.center[0],
                            math.sin(self.tiltrot) * offset + self.center[1],
                            h * (self.startradius + offset)/(2 * self.startradius)]])
        #print('conv3d', conv3d, conv3d[0,0:2])
        
        ####
        # Generate remaining points for top and bottom half
        # Then, combine all points into one set of vertices
        ####
        self.verts, conv3d, conv2d = self.generatehalftwist(self.verts, conv3d, 
                                                            self.center, self.offsetfract, 
                                                            toptwist=True)
        self.cwrot = not self.cwrot
        topverts, conv3dtop, conv2dtop = self.generatehalftwist(topverts, conv3d, 
                                                                topcenter, self.offsetfract,
                                                                toptwist=False)
        
        iconv = self.ngores*4 + 3
        dz2d = self.verts[iconv].pos2d[0,1] - topverts[iconv].pos2d[0,1]
        dz3d = self.verts[iconv].pos3d[0,2] - topverts[iconv].pos3d[0,2]
        #print('conv2d', conv2d, self.startheight2d)
        #print('conv3d', conv3d, self.startheight3d)

        for i in range(self.ngores*2+1, self.ngores*4+2):
            self.verts.append(topverts[i])
            self.verts[-1].pos2d[0,1] += dz2d
            self.verts[-1].pos3d[0,2] += dz3d
            self.verts[-1].pos3dvis[0,2] += dz3d
        for i in range(0, self.ngores*2+1):
            self.verts.append(topverts[i])
            self.verts[-1].pos2d[0,1] += dz2d - self.startheight2d
            self.verts[-1].pos3d[0,2] += dz3d
            self.verts[-1].pos3dvis[0,2] += dz3d
        
        #print(self.verts[-1].pos3d)
        self.endheight3d = self.verts[-1].pos3d[0,2]
        self.endheight2d = self.verts[-1].pos2d[0,1]
        
        # Correct the final 3D positions of vertices
        #self.rotate_xy()
        self.rotate_axis()
        #self.translate3d(z = False)
        # Select end vertices
        self.endverts = self.verts[-2*self.ngores-1:]
        #print('endverts')
        #for vert in self.endverts:
        #    print(vert.pos2d, vert.pos3d)
        #print('end')
        self.center = topcenter
        
        #for vert in self.verts:
        #    print(vert.pos2d, vert.pos3d)

        
    def generateedges(self):
        if self.startedges and len(self.startedges) == self.ngores * 2:
            self.edges = self.startedges
        else:
            # Generate start edges
            for i in range(0, self.ngores*2, 2):
                self.edges.append(el.Edge(self.verts[i],                     self.verts[i+1], 'B'))
                self.edges.append(el.Edge(self.verts[i+1],                   self.verts[i+2], 'B'))
        # Correct directions of start edges
        for i in range(0, self.ngores*2, 2):
            if self.anglediff == None:
                pass
            elif abs(self.anglediff) < 0.01 or self.startedgetype == 'none':
                self.edges[i].direction   = 'U'
                self.edges[i+1].direction = 'U'
            elif (self.anglediff > 0 and self.cwrot) or (self.anglediff < 0 and not self.cwrot):
                self.edges[i].direction   = 'V'
                self.edges[i+1].direction = 'M'
            else:
                self.edges[i].direction   = 'M'
                self.edges[i+1].direction = 'V'
        # Connect start cylinder to start shift
        for i in range(0, self.ngores*2, 2):
            self.edges.append(el.Edge(self.verts[i],   self.verts[i+1  + self.ngores*2], 'V'))
            if i == 0: self.edges[-1].direction = 'B' 
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2  + self.ngores*2], 'M'))
        self.edges.append(el.Edge(self.verts[self.ngores*2], self.verts[self.ngores*4+1], 'B'))
        
        # Connect start shift
        for i in range(self.ngores*2+1, self.ngores*4+1, 2):
            if self.cwrot:
                self.edges.append(el.Edge(self.verts[i]  , self.verts[i+1], 'V'))
                self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'M'))
            else:
                self.edges.append(el.Edge(self.verts[i]  , self.verts[i+1], 'M'))
                self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'V'))
        
        # Connect start shift to middle shift
        for i in range(self.ngores*2+1, self.ngores*4+1, 2):
            self.edges.append(el.Edge(self.verts[i],   self.verts[i+1  + self.ngores*2], 'V'))
            if i == self.ngores*2+1: self.edges[-1].direction = 'B' 
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2  + self.ngores*2], 'M'))
        self.edges.append(el.Edge(self.verts[self.ngores*4+1], self.verts[self.ngores*6+2], 'B'))

        
        # Connect middle shift to end shift
        for i in range(self.ngores*4+2, self.ngores*6+2, 2):
            self.edges.append(el.Edge(self.verts[i],   self.verts[i+1  + self.ngores*2], 'V'))
            if i == self.ngores*4+2: self.edges[-1].direction = 'B' 
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2  + self.ngores*2], 'M'))
        self.edges.append(el.Edge(self.verts[self.ngores*6+2], self.verts[self.ngores*8+3], 'B'))
        
        # Connect end shift
        for i in range(self.ngores*6+3, self.ngores*8+3, 2):
            if self.cwrot:
                self.edges.append(el.Edge(self.verts[i]  , self.verts[i+1], 'M'))
                self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'V'))
            else:
                self.edges.append(el.Edge(self.verts[i]  , self.verts[i+1], 'V'))
                self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'M'))
        
        # Connect end shift to end cylinder
        for i in range(self.ngores*6+3, self.ngores*8+3, 2):
            self.edges.append(el.Edge(self.verts[i],   self.verts[i+1  + self.ngores*2], 'V'))
            if i == self.ngores*6+3: self.edges[-1].direction = 'B' 
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2  + self.ngores*2], 'M'))
        self.edges.append(el.Edge(self.verts[self.ngores*8+3], self.verts[self.ngores*10+4], 'B'))
       
        # Connect end cylinder
        for i in range(self.ngores*8+4, self.ngores*10+4, 2):
            self.edges.append(el.Edge(self.verts[i]  , self.verts[i+1], 'B'))
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'B'))
        
        self.endedges = self.edges[-self.ngores*2:]
        #print('edges', len(self.endedges))
        
        
    def generatefaces(self):
        # Add faces
        for i in range(self.ngores*2):
            self.faces.append(el.Face(verts=[self.verts[i], 
                                             self.verts[i + self.ngores*2 + 1],
                                             self.verts[i + self.ngores*2 + 2],
                                             self.verts[i+1]],
                                      edgelist = self.edges))
        for i in range(self.ngores*2):
            istart = i + self.ngores*2 + 1
            self.faces.append(el.Face(verts=[self.verts[istart],
                                             self.verts[istart + self.ngores*2 + 1],
                                             self.verts[istart + self.ngores*2 + 2],
                                             self.verts[istart + 1]],
                                      edgelist = self.edges))
        
        for i in range(self.ngores*2):
            istart = i + self.ngores*4 + 2
            self.faces.append(el.Face(verts=[self.verts[istart],
                                             self.verts[istart + self.ngores*2 + 1],
                                             self.verts[istart + self.ngores*2 + 2],
                                             self.verts[istart + 1]],
                                      edgelist = self.edges))
        
        for i in range(self.ngores*2):
            istart = i + self.ngores*6 + 3
            self.faces.append(el.Face(verts=[self.verts[istart],
                                             self.verts[istart + self.ngores*2 + 1],
                                             self.verts[istart + self.ngores*2 + 2],
                                             self.verts[istart + 1]],
                                      edgelist = self.edges))
        

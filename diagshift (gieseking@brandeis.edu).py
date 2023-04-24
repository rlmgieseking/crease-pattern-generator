# -*- coding: utf-8 -*-
"""
DiagonalShift class
"""

from molecule import Molecule
import math
import elements as el
import numpy as np
from copy import deepcopy
from scipy.optimize import fsolve

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
        self.tiltrot = tiltrot
        self.angle = math.pi/2
        self.tiltangle = math.atan(offsetfract / math.sqrt(1 - offsetfract**2))
        print('tilt', self.tiltangle)
        
        if self.incomingangle != None:
            self.anglediff = self.angle - self.incomingangle
        else:
            self.anglediff = None

        self.generateverts()
        #self.generateedges()
        #self.generatefaces()
    
    def generateverts(self):
        # Generate or load bottom edge points
        innerscale = 1.0
        if self.startverts and len(self.startverts) == self.ngores * 2 + 1:
            self.verts = self.startverts
        else:
            if self.cwrot:
                startouter3d = [self.startradius, -self.gorewidth/2, 0]
                self.verts.append(el.Vertex([0.0, self.startheight2d], startouter3d))
                self.verts[-1].rotate(-math.pi * 2 / self.ngores, 'z') 
                startinner3d = [self.startradius, self.startradius * math.tan(math.pi/self.ngores) * innerscale, 0]
                self.verts.append(el.Vertex([self.gorewidth/2 - self.startradius * math.tan(math.pi/self.ngores), self.startheight2d], startinner3d))
            else: 
                startouter3d = [self.startradius, self.gorewidth/2, 0]
                self.verts.append(el.Vertex([0.0, self.startheight2d], startouter3d))
                startinner3d = [self.startradius, -self.startradius * math.tan(math.pi/self.ngores) * innerscale, 0]
                self.verts.append(el.Vertex([self.gorewidth/2 + self.startradius * math.tan(math.pi/self.ngores), self.startheight2d], startinner3d))
            
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
        
        # Generate top edge vertices, 3D height for now is arbitrary
        topedge = []
        if not self.cwrot:
            endouter3d = [self.endradius, -self.gorewidth/2, 0]
            topedge.append(el.Vertex([0.0, self.startheight2d], endouter3d))
            topedge[-1].rotate(-math.pi * 2 / self.ngores, 'z') 
            topedge[-1].rotate(- 2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
            endinner3d = [self.endradius, 
                          self.endradius * math.tan(math.pi/self.ngores) * innerscale, 
                          0]
            topedge.append(el.Vertex([self.gorewidth/2 - self.endradius * math.tan(math.pi/self.ngores), 
                                       self.startheight2d], 
                                      endinner3d))
            topedge[-1].rotate(- 2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
        else: 
            endouter3d = [self.endradius, self.gorewidth/2, 0]
            topedge.append(el.Vertex([0.0, self.startheight2d], endouter3d))
            topedge[-1].rotate(2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
            endinner3d = [self.endradius, -self.endradius * math.tan(math.pi/self.ngores) * innerscale, 0]
            topedge.append(el.Vertex([self.gorewidth/2 + self.endradius * math.tan(math.pi/self.ngores), 
                                      self.startheight2d], 
                                     endinner3d))
            topedge[-1].rotate(2 * math.pi * (self.ngores/2 - 1) / self.ngores, 'z')
        
        # Replicate verts for each gore, with the correct rotations and translations
        for i in range(1, self.ngores):
            topedge.append(deepcopy(topedge[-2]))
            topedge[-1].translate2d([self.gorewidth, 0.0])
            topedge[-1].rotate(math.pi * 2 / self.ngores, 'z')
            #print('vertex2', topedge[-1].pos3d)
            
            topedge.append(deepcopy(topedge[-2]))
            topedge[-1].translate2d([self.gorewidth, 0.0])
            topedge[-1].rotate(math.pi * 2 / self.ngores, 'z')
            #print('vertex2', topedge[-1].pos3d)
            
        # Add final edge vertex
        topedge.append(deepcopy(topedge[-2]))
        topedge[-1].translate2d([self.gorewidth, 0.0])
        topedge[-1].rotate(math.pi * 2 / self.ngores, 'z')
        
        # Translate top edge in the xy plane based on offset
        dx = self.offsetfract * (self.startradius + self.endradius) * math.cos(self.tiltrot)
        dy = self.offsetfract * (self.startradius + self.endradius) * math.sin(self.tiltrot)
        for vert in topedge:
            vert.translate3d([dx,dy, 0.0])
        
        # Set up convergence point xyz position
        # z position will be fine-tuned later
        if self.offsetfract > 0.8:
            print('Diagonal shift offset fraction of ', self.offsetfract, ' is too large')
            print('Resetting offset fraction to 0.8')
            self.offsetfraction = 0.8
        offset = self.offsetfract * self.startradius
        h = 2 * self.startradius * offset / math.sqrt(self.startradius**2 - offset**2)
        conv3d = np.array([[math.cos(self.tiltrot) * offset,
                           math.sin(self.tiltrot) * offset,
                           h * (self.startradius + offset)/(2 * self.startradius)]])
        #print('conv3d', conv3d, conv3d[0,0:2])
        

        
        # Solve heights function to get the actual heights of the points
        init = np.array([1] + [0] * (self.ngores - 1))
        zpos = fsolve(self.heights, init, args=(self.verts, conv3d[0,0:2], True))
        minz = min(0, min(zpos[1:]))
        print('zpos',zpos)
        #print('error',self.heights(zpos, self.verts, conv3d[0,0:2], True))
        conv3d[0,2] = zpos[0] - minz
        conv2d = - minz + math.sqrt((self.verts[1].pos3d[0,0] - conv3d[0,0])**2 +
                                    (self.verts[1].pos3d[0,1] - conv3d[0,1])**2 +
                                    (0 - zpos[0])**2 - 
                                    (self.verts[1].pos2d[0,0] - self.gorewidth/2)**2)
        #print('conv',conv3d, conv2d)
        
        # Add new points for transition from cylinder to shift
        for i in range(self.ngores*2+1):
            self.verts.append(deepcopy(self.verts[i]))
            if i == 0 or i >= self.ngores*2 - 1:
                self.verts[-1].pos3d[0,2] = - minz
                self.verts[-1].pos2d[0,1] = - minz
            else:
                self.verts[-1].pos3d[0,2] = - minz + zpos[math.floor((i+1)/2)]
                self.verts[-1].pos2d[0,1] = - minz + zpos[math.floor((i+1)/2)]
        
        # Add vertices surrounding convergence point 
        # Inner vertices get exact positions
        # Outer vertices get placeholder 3D positions to be fixed in the next step
        for i in range(0, self.ngores * 2 + 1, 2):
            i1 = i + self.ngores*2     if i > 0 else self.ngores*4
            i3 = i + self.ngores*2 + 2 if i <= self.ngores*2 else self.ngores*2 + 2
            x13 = self.verts[i1].pos3d[0,0] - self.verts[i3].pos3d[0,0]
            y13 = self.verts[i1].pos3d[0,1] - self.verts[i3].pos3d[0,1]
            z13 = self.verts[i1].pos3d[0,2] - self.verts[i3].pos3d[0,2]
            d13 = math.sqrt(x13**2 + y13**2 + z13**2)
            if self.cwrot:
                xd = (self.gorewidth/2) * x13 / d13
                yd = (self.gorewidth/2) * y13 / d13
                zd = (self.gorewidth/2) * z13 / d13
            else:
                xd = -(self.gorewidth/2) * x13 / d13
                yd = -(self.gorewidth/2) * y13 / d13
                zd = -(self.gorewidth/2) * z13 / d13
            self.verts.append(el.Vertex([self.gorewidth*i/2, conv2d],
                                        [conv3d[0,0] + xd, conv3d[0,1] + yd, conv3d[0,2] + zd]))
            self.verts.append(el.Vertex([self.gorewidth*(i+1)/2 , conv2d],
                                        [conv3d[0,0], conv3d[0,1], conv3d[0,2]]))
        self.verts.pop(-1)
        
        # Fix positions of outer vertices surrounding the convergence point
        vertind = [0,0,0]
        for i in range(self.ngores+1):
            vertind[0] = i*2 + self.ngores*2 + 1
            vertind[1] = i*2 + self.ngores*4 + 3 if i < self.ngores else i*2 + self.ngores*4 + 1
            vertind[2] = vertind[0] + 2 if i < self.ngores else vertind[0] + 2 - self.ngores*2
            v = i*2 + self.ngores*4 + 2
            #print('vertind', i, v, vertind)
            dx = abs(self.verts[vertind[1]].pos2d[0,0] - self.verts[v].pos2d[0,0])
            dy = abs(self.verts[vertind[0]].pos2d[0,1] - self.verts[v].pos2d[0,1])
            #print(self.verts[v].pos3d)
            self.verts[v].pos3d = np.array([fsolve(self.locpt,self.verts[v].pos3d, 
                                         args=(dx, dy, self.verts, vertind))])
            #print(self.locpt(self.verts[v].pos3d, dx, dy, self.verts, vertind))
            
            
            
        
        '''
        # Add new points for top edge of shift
        for i in range(self.ngores*2+1):
            self.verts.append(deepcopy(topedge[i]))
            if i <= 1 or i > self.ngores*2 - 1:
                self.verts[-1].pos3d[0,2] = - minz + zpos[self.ngores - 1]
                self.verts[-1].pos2d[0,1] = - minz + zpos[self.ngores - 1] + zpos[-1]
            else:
                self.verts[-1].pos3d[0,2] = - minz + zpos[self.ngores + math.floor((i-2)/2)]
                self.verts[-1].pos2d[0,1] = - minz + zpos[self.ngores + math.floor((i-2)/2)] + zpos[-1]
            maxz = max(maxz, self.verts[-1].pos3d[0,2])
            maxy = max(maxy, self.verts[-1].pos2d[0,1])
        # Add new points for straight top edge
        for i in range(self.ngores*2+1):
            self.verts.append(deepcopy(topedge[i]))
            self.verts[-1].pos3d[0,2] = maxz
            self.verts[-1].pos2d[0,1] = maxy
        '''
        
        #print(self.verts[-1].pos3d)
        self.endheight3d = self.verts[-1].pos3d[0,2]
        self.endheight2d = self.verts[-1].pos2d[0,1]
        
        for vert in self.verts:
            print(vert.pos2d, vert.pos3d)
        #for vert in topedge:
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
            if not self.anglediff:
                pass
            elif self.anglediff == 0 or self.startedgetype == 'none':
                self.edges[i].direction   = 'U'
                self.edges[i+1].direction = 'U'
            elif (self.anglediff > 0 and not self.cwrot) or (self.anglediff < 0 and self.cwrot):
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
            self.edges.append(el.Edge(self.verts[i]  , self.verts[i+1], 'V'))
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'M'))
        
        # Connect start shift to end shift
        for i in range(self.ngores*2+1, self.ngores*4+1, 2):
            self.edges.append(el.Edge(self.verts[i],   self.verts[i+1  + self.ngores*2], 'V'))
            if i == self.ngores*2+1: self.edges[-1].direction = 'B' 
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2  + self.ngores*2], 'M'))
        self.edges.append(el.Edge(self.verts[self.ngores*4+1], self.verts[self.ngores*6+2], 'B'))

        '''
        # Connect middle shift to end shift
        for i in range(self.ngores*4+2, self.ngores*6+2, 2):
            self.edges.append(el.Edge(self.verts[i],   self.verts[i+1  + self.ngores*2], 'V'))
            if i == self.ngores*4+2: self.edges[-1].direction = 'B' 
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2  + self.ngores*2], 'M'))
        self.edges.append(el.Edge(self.verts[self.ngores*6+2], self.verts[self.ngores*8+3], 'B'))
        '''
        # Connect end shift
        for i in range(self.ngores*4+2, self.ngores*6+2, 2):
            self.edges.append(el.Edge(self.verts[i]  , self.verts[i+1], 'M'))
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'V'))
        
        # Connect end shift to end cylinder
        for i in range(self.ngores*4+2, self.ngores*6+2, 2):
            self.edges.append(el.Edge(self.verts[i],   self.verts[i+1  + self.ngores*2], 'V'))
            if i == self.ngores*6+3: self.edges[-1].direction = 'B' 
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2  + self.ngores*2], 'M'))
        self.edges.append(el.Edge(self.verts[self.ngores*6+2], self.verts[self.ngores*8+3], 'B'))
       
        # Connect end cylinder
        for i in range(self.ngores*6+3, self.ngores*8+3, 2):
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
        '''
        for i in range(self.ngores*2):
            istart = i + self.ngores*6 + 3
            self.faces.append(el.Face(verts=[self.verts[istart],
                                             self.verts[istart + self.ngores*2 + 1],
                                             self.verts[istart + self.ngores*2 + 2],
                                             self.verts[istart + 1]],
                                      edgelist = self.edges))
        '''
    
       
        
        
        '''
        # Set up convergence point xyz position
        # z position will be fine-tuned later
        if self.offsetfract > 0.8:
            print('Diagonal shift offset fraction of ', self.offsetfract, ' is too large')
            print('Resetting offset fraction to 0.8')
            self.offsetfraction = 0.8
        offset = self.offsetfract * self.startradius
        h = 2 * self.startradius * offset / math.sqrt(self.startradius**2 - offset**2)
        conv3d = np.array([[math.cos(self.tiltrot) * offset,
                           math.sin(self.tiltrot) * offset,
                           h * (self.startradius + offset)/(2 * self.startradius)]])
        hd = h * (self.startradius - offset)/(2 * self.startradius)
        
        # For each outer point, compute the xy distance to the convergence point
        xydist = []
        for vertex in self.verts:
            xydist.append(math.sqrt((vertex.pos3d[0,0] - conv3d[0,0])**2 
                                    + (vertex.pos3d[0,1] - conv3d[0,1])**2))
        #print(xydist)
        minxydist = min(xydist)
        denom = 2 * (math.sqrt(minxydist**2 + hd**2) + hd)
        # For each outer point, compute the z that equalizes lengths & distance to conv3d
        zdist = []
        totdist = []
        for i in range(len(xydist)):
            zdist.append(conv3d[0,2] + hd 
                         - (xydist[i]**2 - minxydist**2)/denom)
            totdist.append(math.sqrt(xydist[i]**2 + zdist[i]**2))
            
        # z coordinates of outer vertices are based on a straight line between the
        # adjacent inner vertices
        for i in range(0, len(zdist), 2):
            i1 = i - 1 if i > 0 else len(zdist) - 2
            i3 = i + 1 if i + 1 < len(zdist) else 1
            d13 = math.sqrt((self.verts[i1].pos3d[0,0] - self.verts[i3].pos3d[0,0])**2 
                           +(self.verts[i1].pos3d[0,1] - self.verts[i3].pos3d[0,1])**2 )
            z13 = zdist[i1] - zdist[i3]
            #print(i, i1, i3, d13, z13)
            if self.cwrot:
                d23 = math.sqrt((self.verts[i].pos3d[0,0] - self.verts[i3].pos3d[0,0])**2 
                               +(self.verts[i].pos3d[0,1] - self.verts[i3].pos3d[0,1])**2 )
                zdist[i] = zdist[i3] + d23 * z13 / d13
            else:
                d21 = math.sqrt((self.verts[i].pos3d[0,0] - self.verts[i1].pos3d[0,0])**2 
                               +(self.verts[i].pos3d[0,1] - self.verts[i1].pos3d[0,1])**2 )
                zdist[i] = zdist[i1] + d21 * z13 / d13
        
        # Correct z distances based on minimum z distance
        minzdist = min(zdist)
        print(minzdist)
        for i in range(len(zdist)):
            zdist[i] -= minzdist
        #conv3d[0,2] -= minzdist
        
        # Generate points for the transition from the cylinder to the shift
        for i in range(self.ngores * 2 + 1):
            self.verts.append(deepcopy(self.verts[i]))
            self.verts[-1].pos3d[0,2] += zdist[i]
            self.verts[-1].pos2d[0,1] += zdist[i]
        
        # Compute the 2D height of the convergence point
        dx = self.verts[1].pos2d[0,0] - self.gorewidth/2 
        dxy = math.sqrt((self.verts[1].pos3d[0,0] - conv3d[0,0])**2 
                       + (self.verts[1].pos3d[0,1] - conv3d[0,1])**2)
        conv2d = self.verts[self.ngores*2 + 2].pos2d[0,1] + math.sqrt(dxy**2 - dx**2)
        #print('conv2d',conv2d, self.startradius)
        #for i in range(len(self.verts)):
        #    print(self.verts[i].pos2d[0,0], self.verts[i].pos2d[0,1])
        
        # Add vertices surrounding convergence point 
        for i in range(0, self.ngores * 2 + 1, 2):
            i1 = i - 1
            i3 = i + 1 if i + 1 < len(zdist) else 1
            x13 = self.verts[i1].pos3d[0,0] - self.verts[i3].pos3d[0,0]
            y13 = self.verts[i1].pos3d[0,1] - self.verts[i3].pos3d[0,1]
            z13 = zdist[i1] - zdist[i3]
            d13 = math.sqrt(x13**2 + y13**2 + z13**2)
            if self.cwrot:
                xd = (self.gorewidth/2) * x13 / d13
                yd = (self.gorewidth/2) * y13 / d13
                zd = (self.gorewidth/2) * z13 / d13
            else:
                xd = -(self.gorewidth/2) * x13 / d13
                yd = -(self.gorewidth/2) * y13 / d13
                zd = -(self.gorewidth/2) * z13 / d13
            self.verts.append(el.Vertex([self.gorewidth*i/2, conv2d],
                                        [conv3d[0,0] + xd, conv3d[0,1] + yd, conv3d[0,2] + zd]))
            self.verts.append(el.Vertex([self.gorewidth*(i+1)/2 , conv2d],
                                        [conv3d[0,0], conv3d[0,1], conv3d[0,2]]))
        self.verts.pop(-1)
        '''
        '''
        # Generate vertices for the top half of the diagonal shift
        # Points for the transition from the twist to the cylinder should be on a 
        # straight line from points from the bottom edge of the twist
        self.cwrot = not self.cwrot
        dfrac = self.endradius/self.startradius
        #print(dfrac, self.endradius, self.startradius)
        zmax = -100.0
        ymax = -100.0
        for i in range(0, self.ngores * 2, 2):
            # Add a placeholder for the outer vertex
            self.verts.append(el.Vertex([i*self.gorewidth/2,conv2d],conv3d[0,:]))
            # Position of the inner vertex follows a straight line through the convergence point
            i1 = self.ngores*2 + 1 + i + 1
            i2 = (self.ngores*2 + 1)*2 + i + 1
            x21 = self.verts[i2].pos3d[0,0] - self.verts[i1].pos3d[0,0]
            y21 = self.verts[i2].pos3d[0,1] - self.verts[i1].pos3d[0,1]
            z21 = self.verts[i2].pos3d[0,2] - self.verts[i1].pos3d[0,2]
            xa21 = self.verts[i2].pos2d[0,0] - self.verts[i1].pos2d[0,0]
            ya21 = self.verts[i2].pos2d[0,1] - self.verts[i1].pos2d[0,1]
            self.verts.append(el.Vertex([self.verts[i2].pos2d[0,0] + dfrac * xa21,
                                         self.verts[i2].pos2d[0,1] + dfrac * ya21],
                                        [self.verts[i2].pos3d[0,0] + dfrac * x21,
                                         self.verts[i2].pos3d[0,1] + dfrac * y21,
                                         self.verts[i2].pos3d[0,2] + dfrac * z21]))
            zmax = max(zmax, self.verts[i2].pos3d[0,2] + dfrac * z21)
            ymax = max(ymax, self.verts[i2].pos2d[0,1] + dfrac * ya21)
        # Add final placeholder outer vertex
        self.verts.append(el.Vertex([self.ngores*self.gorewidth,conv2d],conv3d[0,:]))
        
        # Correct the 2D and 3D positions of the outer vertices
        dfrac = self.verts[(self.ngores*2 + 1)*3 + 1].pos2d[0,0]/self.gorewidth
        if dfrac < 0.5:
            scale = dfrac/(1 - dfrac*2)
        else:
            scale = (1 - dfrac)/(1 - (1 - dfrac)*2)
        print('dfrac',dfrac, scale, self.verts[(self.ngores*2 + 1)*3 + 1].pos2d)
        for i in range(0, self.ngores*2+1, 2):
            i1 = (self.ngores*2 + 1)*3 + i - 1 if i > 0 else (self.ngores*2 + 1)*4 + i - 1
            i2 = (self.ngores*2 + 1)*3 + i
            i3 = (self.ngores*2 + 1)*3 + i + 1 if i < self.ngores*2 else (self.ngores*2 + 1)*2 + i + 1
            x31 = self.verts[i3].pos3d[0,0] - self.verts[i1].pos3d[0,0]
            y31 = self.verts[i3].pos3d[0,1] - self.verts[i1].pos3d[0,1]
            z31 = self.verts[i3].pos3d[0,2] - self.verts[i1].pos3d[0,2]
            xa31 = self.verts[i3].pos2d[0,0] - self.verts[i1].pos2d[0,0]
            ya31 = self.verts[i3].pos2d[0,1] - self.verts[i1].pos2d[0,1]
            
            if dfrac < 0.5:
                self.verts[i2] = el.Vertex([self.verts[i3].pos2d[0,0] + scale * xa21,
                                            self.verts[i3].pos2d[0,1] + scale * ya21],
                                           [self.verts[i3].pos3d[0,0] + scale * x21,
                                            self.verts[i3].pos3d[0,1] + scale * y21,
                                            self.verts[i3].pos3d[0,2] + scale * z21])
                print(self.verts[i3].pos2d[0,0], self.verts[i1].pos2d[0,0], xa21, self.verts[i3].pos2d[0,0] + scale * xa21)
            else:
                pass
        
        
        for i in range(0, self.ngores * 2 + 1):
            self.verts.append(deepcopy(self.verts[-2 * self.ngores - 1]))
            self.verts[-1].pos3d[0,2] = zmax
            self.verts[-1].pos2d[0,1] = ymax
        '''
        '''
        ######
        # 2. Generate lower transition from cylinder to shift
        ######
        # Compute the vertical distance the shift plane deviates from each point, 
        # then use the lowest point to determine vertical positions of the points
        tiltaxis = np.array([math.cos(self.tiltrot), math.sin(self.tiltrot), 0])
        #print('tiltaxis',tiltaxis)
        vertdist = []
        for vertex in self.verts:
            # Project (x,y) 3D position onto tilt axis
            #print(vertex.pos3d, tiltaxis)
            #planexydist = vertex.pos3d[0,0] * tiltaxis[0] + vertex.pos3d[0,1] * tiltaxis[1]
            planexydist = np.dot(vertex.pos3d[0,:], tiltaxis)
            # Compute z distance of point to tilt plane
            #print('dist', planexydist, vertex.pos3d[0,:], tiltaxis, planexydist * math.tan(self.tiltangle))
            vertdist.append(planexydist * math.tan(self.tiltangle))
        #print('vertdist', vertdist)
        minvertdist = min(vertdist)
        
        for i in range(0, len(self.verts)):
            self.verts.append(deepcopy(self.verts[i]))
            self.verts[-1].translate2d([0.0, vertdist[i]-minvertdist])
            self.verts[-1].translate3d([0.0, 0.0, vertdist[i]-minvertdist])
        
        # Find the convergence point
        # Set up vertical distances of only inside points
        vertinside = []
        for i in range(1,len(vertdist),2):
            vertinside.append(vertdist[i])
        lowindex = self.ngores*2 + 1 + vertdist.index(min(vertinside)) 
        highindex = lowindex + self.ngores if lowindex + self.ngores < len(self.verts) else lowindex - self.ngores
        lowest, highest = self.verts[lowindex], self.verts[highindex]
        print(lowindex, highindex, lowest.pos3d, highest.pos3d)
        
        dxy = math.sqrt((highest.pos3d[0,0] - lowest.pos3d[0,0])**2 
                        + (highest.pos3d[0,1] - lowest.pos3d[0,1])**2)
        dz = highest.pos3d[0,2] - lowest.pos3d[0,2]
        dxyz = math.sqrt(dxy**2 + dz**2)
        # Distance from lowest point to convergence point
        dlowconv = (dxyz + dz) / 2.0
        dfract = dlowconv / dxyz
        #conv3d = [lowest.pos3d[0,0] + dfract * (highest.pos3d[0,0] - lowest.pos3d[0,0]),
        #          lowest.pos3d[0,1] + dfract * (highest.pos3d[0,1] - lowest.pos3d[0,1]),
        #          lowest.pos3d[0,2] + dfract * (highest.pos3d[0,2] - lowest.pos3d[0,2])]
        convheight2d = lowest.pos2d[0,1] + dlowconv
        #print('high/low', highindex, lowindex, highest.pos3d, lowest.pos3d)
        #print('conv3d', conv3d)
        conv3d = lowest.pos3d + dfract * (highest.pos3d - lowest.pos3d)
        conv3d = conv3d[0]
        print('conv3d', conv3d)
        self.verts.append(el.Vertex([0,convheight2d],conv3d))
        '''
        '''
        # Generate top circle end points - start off with points centered around axis at start height, then shift later to new cylinder center and end height
        temptopedge = []
        tempinner3d = np.array([-self.endradius, -self.endradius * math.tan(math.pi/self.ngores), self.startheight3d])
        tempouter3d = np.array([-self.endradius, self.gorewidth/2, self.startheight3d])
        
        temptopedge.append(el.Vertex([0.0, self.startheight2d], tempouter3d))
        temptopedge.append(el.Vertex([self.gorewidth/2 + -self.endradius * math.tan(math.pi/self.ngores), self.startheight2d], tempinner3d))
        
        # Replicate vertices for each gore, with the correct rotations and translations
        for i in range(1, self.ngores):
            temptopedge.append(deepcopy(temptopedge[-2]))
            temptopedge[-1].translate2d([self.gorewidth, 0.0])
            temptopedge[-1].rotate(math.pi * 2 / self.ngores, 'z')
            
            temptopedge.append(deepcopy(temptopedge[-2]))
            temptopedge[-1].translate2d([self.gorewidth, 0.0])
            temptopedge[-1].rotate(math.pi * 2 / self.ngores, 'z')
            
        # Add final edge vertex
        temptopedge.append(deepcopy(temptopedge[-2]))
        temptopedge[-1].translate2d([self.gorewidth, 0.0])
        temptopedge[-1].rotate(math.pi * 2 / self.ngores, 'z')
        
        # Using this top edge, compute the shift end points similar to before
        # Compute the vertical distance the shift plane deviates from each point, then use the lowest point to determine vertical positions of the points
        vertdisttop = []
        for vertex in temptopedge:
            # Project (x,y) 3D position onto tilt axis
            planexydist = vertex.pos3d[0,0] * tiltaxis[0] + vertex.pos3d[0,1] * tiltaxis[1]
            # Compute z distance of point to tilt plane
            vertdisttop.append(planexydist * math.tan(self.tiltangle))
        maxvertdisttop = max(vertdisttop)
        
        temptopshift = []
        for i in range(0, len(temptopedge)):
            temptopshift.append(deepcopy(temptopedge[i]))
            temptopshift[-1].translate2d([0.0, vertdisttop[i]+maxvertdisttop])
            temptopshift[-1].translate3d([0.0, 0.0, vertdisttop[i]-maxvertdisttop])
        
        # Find the convergence point
        # Set up vertical distances of only inside points
        vertinsidetop = []
        for i in range(1,len(vertdisttop),2):
            vertinsidetop.append(vertdisttop[i])
        lowindextop = vertinsidetop.index(min(vertinsidetop)) * 2 + 1
        highindextop = lowindextop + self.ngores if lowindextop + self.ngores < len(temptopshift) else lowindextop - self.ngores
        lowesttop, highesttop = temptopshift[lowindextop], temptopshift[highindextop]
        
        dxytop = math.sqrt((highesttop.pos3d[0,0] - lowesttop.pos3d[0,0])**2 
                           + (highesttop.pos3d[0,1] - lowesttop.pos3d[0,1])**2)
        dztop = highesttop.pos3d[0,2] - lowesttop.pos3d[0,2]
        dxyztop = math.sqrt(dxytop**2 + dztop**2)
        # Distance from lowest point to convergence point
        dlowconvtop = (dxyztop + dztop) / 2.0
        dfracttop = dlowconvtop / dxyztop
        conv3dtop = highesttop.pos3d - dfracttop * (highesttop.pos3d - lowesttop.pos3d)
        #conv3dtop = [highesttop.pos3d[0,0] - dfracttop * (highesttop.pos3d[0,0] - lowesttop.pos3d[0,0]),
        #             highesttop.pos3d[0,1] - dfracttop * (highesttop.pos3d[0,1] - lowesttop.pos3d[0,1]),
        #             highesttop.pos3d[0,2] - dfracttop * (highesttop.pos3d[0,2] - lowesttop.pos3d[0,2])]
        convheight2dtop = lowesttop.pos2d[0,1] + dlowconvtop
        #print(vertdisttop, vertinsidetop)
        #print('top high/low', highindextop, lowindextop, highesttop.pos3d, lowesttop.pos3d)
        #print('conv3dtop', conv3dtop)
        
        # Translate upper half points so the convergence point matches the lower half and add points to master list
        #conv3dshift = [-conv3dtop[0] + conv3d[0], -conv3dtop[1] + conv3d[1], -conv3dtop[2] + conv3d[2]]
        conv3dshift = conv3d - conv3dtop
        #print(conv3dshift)
        convheight2dshift = -convheight2dtop + convheight2d
        for vertex in temptopshift:
            vertex.translate2d([0, convheight2dshift])
            vertex.translate3d(conv3dshift)
            self.verts.append(vertex)
        for vertex in temptopedge:
            vertex.translate2d([0, convheight2dshift])
            vertex.translate3d(conv3dshift)
            self.verts.append(vertex)
        #for vertex in self.vertices:
        #    print(vertex.pos3d[0,0], vertex.pos3d[0,1], vertex.pos3d[0,2])
        
        # Store final top edge info
        #print('center',self.center, conv3d, self.startradius, self.endradius)
        self.topcenter = self.center + conv3d[:,:2]*((self.startradius + self.endradius)/self.startradius)
        '''
        
        '''
        #minz = min(zpos[:self.ngores-1])
        maxz, maxy = 0, 0
        #print('z', minz, maxz)
        # Function to compute error in 2D vs. 3D fold lengths for all mountains
        # and valleys in the twist. Use to optimize the vertex heights.
        # Define all distances based only on inner points. Outer point positions
        # will be computed later based on these.
        # Input vector h has indices [gi (starting at 1), hi (starting at 0), dh]
        def heights(h, dz=0):
            # Separate h into useful quantities
            # hb = 3D heights of bottom vertices (first fixed at 0)
            # ht = 3D heights of top vertices (will be offset for planarity later)
            # offset2d = extra length needed in 2d crease pattern
            hb = [0] + np.ndarray.tolist(h[:self.ngores-1]) + [0]
            ht = [0] + np.ndarray.tolist(h[self.ngores-1:self.ngores*2-2]) + [0]
            print(hb,ht)
            offset2d = h[-1]
            error = []
            
            # First, compute offset of 3D heights
            # Set such that the first 3D vertex is dz above the plane of the 
            # first three lower vertices
            vec1 = self.verts[3].pos3d - self.verts[1].pos3d
            vec2 = self.verts[5].pos3d - self.verts[1].pos3d
            normal = np.cross(vec1, vec2)
            normal /= np.linalg.norm(normal)
            intersect = (self.verts[1].pos3d - 
                         normal * np.tensordot(self.verts[1].pos3d - topedge[1].pos3d, normal)/
                         np.tensordot(np.array([[0,0,1]]), normal))
            #print(intersect)
            dz += intersect[0,2]
            #print(dz)
            
            # Compute error in length of mountain folds
            # Mountain folds connect outer points
            # If radii are less than max, z positions must be derived based on adjacent inner points
            for i in range(0, self.ngores):
                x3d = self.verts[i*2].pos3d[0,0] - topedge[i*2].pos3d[0,0]
                y3d = self.verts[i*2].pos3d[0,1] - topedge[i*2].pos3d[0,1]
                z3d = hb[i] - ht[i]
                #if i == 0:
                #    z3d = 0 - h[self.ngores + i - 1]
                #else:
                #    z3d = h[i-1] - h[self.ngores + i - 1] 
                #    
                x2d = self.verts[i*2].pos2d[0,0] - topedge[i*2].pos2d[0,0]
                y2d = z3d - offset2d
                print(x3d, y3d, z3d, x2d, y2d)
                error.append(math.sqrt(x3d**2 + y3d**2 + z3d**2) - math.sqrt(x2d**2 + y2d**2))
            # Compute error in length of valley folds
            for i in range(0, self.ngores-1):
                #print(self.verts[i*2+1].pos3d, topedge[i*2+1].pos3d)
                #print(self.verts[i*2+1].pos2d, topedge[i*2+1].pos2d)
                x3d = self.verts[i*2+1].pos3d[0,0] - topedge[i*2+1].pos3d[0,0]
                y3d = self.verts[i*2+1].pos3d[0,1] - topedge[i*2+1].pos3d[0,1]
                z3d = hb[i] - ht[i-1]
                #if i == 0:
                #    z3d = h[i] - h[self.ngores + i - 1]
                #elif i == self.ngores - 1:
                #    z3d = 0 - h[self.ngores + i - 1]
                #else: 
                #    z3d = h[i] - h[self.ngores + i - 1]
                x2d = self.verts[i*2+1].pos2d[0,0] - topedge[i*2+1].pos2d[0,0]
                y2d = z3d - offset2d
                #print(x3d, y3d, z3d, x2d, y2d)
                error.append(math.sqrt(x3d**2 + y3d**2 + z3d**2) - math.sqrt(x2d**2 + y2d**2))
            # Add error in final valley into last point
            i = self.ngores-1
            x3d = self.verts[i*2+1].pos3d[0,0] - topedge[i*2+1].pos3d[0,0]
            y3d = self.verts[i*2+1].pos3d[0,1] - topedge[i*2+1].pos3d[0,1]
            z3d = hb[i] - ht[i-1]
            x2d = self.verts[i*2+1].pos2d[0,0] - topedge[i*2+1].pos2d[0,0]
            y2d = z3d - offset2d
            error[-1] += math.sqrt(x3d**2 + y3d**2 + z3d**2) - math.sqrt(x2d**2 + y2d**2)
            #print(error)
            return error
        
        # Solve heights function to get the actual heights of the points
        init = np.array([0] * (self.ngores - 1) + [0] * (self.ngores - 1) + [3])
        #init = [1] * (self.ngores * 2)
        #print(init)
        print(heights(init))
        zpos = init
        zpos = fsolve(heights, init, args=(1))
        #print('zpos',zpos)
        #print('error',heights(zpos))
        # Add new points for transition from cylinder to shift
        for i in range(self.ngores*2+1):
            self.verts.append(deepcopy(self.verts[i]))
            if i == 0 or i >= self.ngores*2 - 1:
                self.verts[-1].pos3d[0,2] = - minz
                self.verts[-1].pos2d[0,1] = - minz
            else:
                self.verts[-1].pos3d[0,2] = - minz + zpos[math.floor((i-1)/2)]
                self.verts[-1].pos2d[0,1] = - minz + zpos[math.floor((i-1)/2)]
        # Add new points for top edge of shift
        for i in range(self.ngores*2+1):
            self.verts.append(deepcopy(topedge[i]))
            if i <= 1 or i > self.ngores*2 - 1:
                self.verts[-1].pos3d[0,2] = - minz + zpos[self.ngores - 1]
                self.verts[-1].pos2d[0,1] = - minz + zpos[self.ngores - 1] + zpos[-1]
            else:
                self.verts[-1].pos3d[0,2] = - minz + zpos[self.ngores + math.floor((i-2)/2)]
                self.verts[-1].pos2d[0,1] = - minz + zpos[self.ngores + math.floor((i-2)/2)] + zpos[-1]
            maxz = max(maxz, self.verts[-1].pos3d[0,2])
            maxy = max(maxy, self.verts[-1].pos2d[0,1])
        # Add new points for straight top edge
        for i in range(self.ngores*2+1):
            self.verts.append(deepcopy(topedge[i]))
            self.verts[-1].pos3d[0,2] = maxz
            self.verts[-1].pos2d[0,1] = maxy
        '''
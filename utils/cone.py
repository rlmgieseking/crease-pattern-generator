# -*- coding: utf-8 -*-

from utils.molecule import Molecule
from copy import deepcopy
import math
import utils.elements as el

"""
Cone class
Inherits from Molecule

New attributes:
    endradius: Radius of the cone at its final edge
    diffheight3d: Difference between height at the start edge and height at the
        end edge in 3D. Positive number means the cone goes up.
    endheight3d: Height at the end edge in 3D
    angle: Angle of the faces of the cone relative to its central axis
    anglediff: Difference in angle between the faces of this cone and the faces
        of the previous molecule, if one exists

Functions: 
    generateverts: Generate all vertices of the cone in 2D and 3D. Generates
        vertices along the start edge if it does not exist, and always 
        generates vertices along the end edge.
    generateedges: Generate all edges by connecting vertices and assign fold 
        directions.
    generatefaces: Generate all faces based on vertices for that face and the
        full list of edges in the cone
"""
class Cone(Molecule):
    def __init__(self, ngores, gorewidth, cwrot, 
                 startradius, startheight2d, startheight3d, startverts, startedges,
                 incomingangle, xyrot, axis, center, startedgetype, 
                 endradius, diffheight3d):
        Molecule.__init__(self, ngores, gorewidth, cwrot, 
                          startradius, startheight2d, startheight3d, startverts, startedges,
                          incomingangle, xyrot, axis, center, startedgetype)

        #print('Start cone',datetime.now())
        self.endradius = endradius
        self.diffheight3d = diffheight3d
        self.endheight3d = startheight3d + diffheight3d
        self.angle = math.atan2(diffheight3d, (endradius - startradius))
        #print('angle', self.angle)
        if self.incomingangle != None:
            self.anglediff = self.angle - self.incomingangle
        else:
            self.anglediff = None
        
        self.generateverts()
        self.generateedges()
        self.generatefaces()

        #print('End   cone',datetime.now())
    
    # Set up all vertices for the cone
    # Start edge and end edge each have 2 * ngores + 1 vertices
    # If the start vertices already exist, only the end vertices need to be generated
    def generateverts(self):
        # Assume that the first side is centered on the +x axis in 3D space
        # Relative 2D and 3D orientations are set up such that 2D shows the inside of the cone (assuming cone goes up)
        # Place cone such that base is centered at (0, 0, 0); translation and rotation come later
        
        # Generate points along the start edge, in order from left to right on page (cw around model)
        # Place first vertex pair in 3D space
        visshift = 0.01
        if self.startverts and len(self.startverts) == self.ngores * 2 + 1:
            self.verts = self.startverts
        else:
            # Generate start edge vertices
            if self.cwrot:
                startouter3d = [self.startradius, -self.gorewidth/2, 0]
                self.verts.append(el.Vertex([0.0, self.startheight2d], startouter3d))
                self.verts[-1].rotate(-math.pi * 2 / self.ngores, 'z') 
                startinner3d = [self.startradius, self.startradius * math.tan(math.pi/self.ngores), 0]
                startinner3dvis = [self.startradius - visshift, self.startradius * math.tan(math.pi/self.ngores), 0]
                self.verts.append(el.Vertex([self.gorewidth/2 - self.startradius * math.tan(math.pi/self.ngores), 
                                             self.startheight2d], 
                                            startinner3d,
                                            startinner3dvis))
            else: 
                startouter3d = [self.startradius, self.gorewidth/2, 0]
                self.verts.append(el.Vertex([0.0, self.startheight2d], startouter3d))
                startinner3d = [self.startradius, -self.startradius * math.tan(math.pi/self.ngores), 0]
                startinner3dvis = [self.startradius - visshift, -self.startradius * math.tan(math.pi/self.ngores), 0]
                self.verts.append(el.Vertex([self.gorewidth/2 + self.startradius * math.tan(math.pi/self.ngores), 
                                             self.startheight2d], 
                                            startinner3d,
                                            startinner3dvis))
    
            # Replicate verts for each gore, with the correct rotations and translations
            for i in range(1, self.ngores):
                self.verts.append(deepcopy(self.verts[-2]))
                self.verts[-1].translate2d([self.gorewidth, 0.0])
                self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
                
                self.verts.append(deepcopy(self.verts[-2]))
                self.verts[-1].translate2d([self.gorewidth, 0.0])
                self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
                
            # Add final edge vertex
            self.verts.append(deepcopy(self.verts[-2]))
            self.verts[-1].translate2d([self.gorewidth, 0.0])
            self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
        
        # Generate points along the end edge
        self.endheight2d = self.startheight2d + math.sqrt((self.endradius - self.startradius)**2 + (self.endheight3d - self.startheight3d)**2)
        # Place first vertex pair in 3D space
        if self.cwrot:
            endouter3d = [self.endradius, -self.gorewidth/2, self.diffheight3d]
            self.verts.append(el.Vertex([0.0, self.endheight2d], endouter3d))
            self.verts[-1].rotate(-math.pi * 2 / self.ngores, 'z') 
            endinner3d = [self.endradius, self.endradius * math.tan(math.pi/self.ngores), self.diffheight3d]
            endinner3dvis = [self.endradius - visshift, self.endradius * math.tan(math.pi/self.ngores), self.diffheight3d]
            self.verts.append(el.Vertex([self.gorewidth/2 - self.endradius * math.tan(math.pi/self.ngores), 
                                         self.endheight2d], 
                                        endinner3d,
                                        endinner3dvis))
        else: 
            endouter3d = [self.endradius, self.gorewidth/2, self.diffheight3d]
            self.verts.append(el.Vertex([0.0, self.endheight2d], endouter3d))
            endinner3d = [self.endradius, -self.endradius * math.tan(math.pi/self.ngores), self.diffheight3d]
            endinner3dvis = [self.endradius - visshift, -self.endradius * math.tan(math.pi/self.ngores), self.diffheight3d]
            self.verts.append(el.Vertex([self.gorewidth/2 + self.endradius * math.tan(math.pi/self.ngores), 
                                         self.endheight2d], 
                                        endinner3d,
                                        endinner3dvis))
        
        # Replicate verts for each gore, with the correct rotations and translations
        for i in range(1, self.ngores):
            self.verts.append(deepcopy(self.verts[-2]))
            self.verts[-1].translate2d([self.gorewidth, 0.0])
            self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
            
            self.verts.append(deepcopy(self.verts[-2]))
            self.verts[-1].translate2d([self.gorewidth, 0.0])
            self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
            
        # Add final edge vertex
        self.verts.append(deepcopy(self.verts[-2]))
        self.verts[-1].translate2d([self.gorewidth, 0.0])
        self.verts[-1].rotate(math.pi * 2 / self.ngores, 'z')
        
        # Correct the final 3D positions of vertices
        self.rotate_xy()
        self.rotate_axis()
        self.translate3d()
        # Select end vertices
        self.endverts = self.verts[2*self.ngores+1:]
        
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
            elif (self.anglediff > 0 and not self.cwrot) or (self.anglediff < 0 and self.cwrot):
                self.edges[i].direction   = 'V'
                self.edges[i+1].direction = 'M'
            else:
                self.edges[i].direction   = 'M'
                self.edges[i+1].direction = 'V'
        # Connect start to end edge
        for i in range(0, self.ngores*2, 2):
            self.edges.append(el.Edge(self.verts[i],                   self.verts[i+1  + self.ngores*2], 'V'))
            self.edges.append(el.Edge(self.verts[i+1],                 self.verts[i+2  + self.ngores*2], 'M'))
        # Final start to end edge
        self.edges.append(el.Edge(self.verts[self.ngores*2], self.verts[self.ngores*4+1], 'B'))
        # Correct first side edge to be a border
        self.edges[self.ngores*2].direction = 'B'
        # Connect end edge
        for i in range(self.ngores*2, self.ngores*4, 2):
            self.edges.append(el.Edge(self.verts[i+1], self.verts[i+2], 'B'))
            self.edges.append(el.Edge(self.verts[i+2], self.verts[i+3], 'B'))
        self.endedges = self.edges[-self.ngores*2:]
        
    
    def generatefaces(self):
        # Add faces
        for i in range(self.ngores*2):
            self.faces.append(el.Face(verts=[self.verts[i], 
                                             self.verts[i + self.ngores*2 + 1],
                                             self.verts[i + self.ngores*2 + 2],
                                             self.verts[i+1]],
                                      edgelist = self.edges))
        #print('faces',len(self.faces))

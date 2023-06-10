# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:34:02 2023

@author: gieseking
"""

import math
from utils.cone import Cone
from utils.diagshift import DiagonalShift

# Contains master list of all elements for the full model
class Model:
    def __init__(self, ngores=8, gorewidth=1.0, cwrot=False, currentradius=0.0, relativedim=True):
        self.ngores = ngores
        self.gorewidth = gorewidth
        self.relativedim = relativedim
        # Define maximum radius for a polygon of ngores sides of length gorewidth
        self.maxradius = gorewidth / 2 / math.tan(math.pi/ngores)
        
        # Direction of flanges from top view. From side view, cw = left, ccw = right
        self.cwrot = cwrot
        if self.relativedim:
            self.currentradius = currentradius * self.maxradius
            if abs(self.currentradius) > self.maxradius:
                print("Error: Relative dimensions are turned on")
                print("All heights and widths are defined as a fraction of the maximum possible radius")
                print("Relative radius cannot exceed 1")
                print("Relative radius of ", currentradius, "is impossible")
        else:
            self.currentradius = currentradius
            if abs(self.currentradius) > self.maxradius:
                print("Error: Radius cannot exceed the maximum possible radius of ", str(self.maxradius))
                print("Radius of ", str(currentradius), " is impossible for ", str(ngores),
                      " sides of width", str(gorewidth))
        self.currentheight2d = 0.0
        self.currentheight3d = 0.0
        self.currentcenter = [0.0, 0.0]
        self.currentaxis = [0.0, 0.0, 1.0]
        # Angle of all facets relative to current axis
        self.currentangle = None
        self.xyrot = 0.0
        
        self.verts = []
        self.edges = []
        self.faces = []
        self.molecules = []

    
    def add_cone(self, endradius, diffheight3d, startedgetype='auto', scaleradius=True):
        if self.relativedim and scaleradius:
            endradius *= self.maxradius
            diffheight3d *= self.maxradius
        if abs(endradius) > self.maxradius:
            print("Error: Radius of ", str(endradius), 
                  "exceeds the maximum possible radius of ", str(self.maxradius))
            endradius = self.maxradius * 0.99999
        elif abs(endradius) == self.maxradius:
            endradius = self.maxradius * 0.99999
        if len(self.molecules) > 0:
            startverts, startedges = self.molecules[-1].endverts, self.molecules[-1].endedges
        else:
            startverts, startedges = None, None
        self.molecules.append(Cone(self.ngores, self.gorewidth, self.cwrot, 
                                   self.currentradius, self.currentheight2d, self.currentheight3d, 
                                   startverts, startedges,
                                   self.currentangle, self.xyrot, self.currentaxis, self.currentcenter,
                                   startedgetype, endradius, diffheight3d))
        
        # Update current radius and height to match the end of the cone
        self.currentradius = endradius
        self.currentheight2d = self.molecules[-1].endheight2d
        if self.currentaxis == [0,0,1]:
            self.currentheight3d += diffheight3d
        else:
            phi = math.atan2(self.currentaxis[1], self.currentaxis[0])
            theta = math.acos(self.currentaxis[2]/math.sqrt(self.currentaxis[0]**2 + self.currentaxis[1]**2 + self.currentaxis[2]**2))
            self.currentheight3d  += diffheight3d * math.cos(theta)
            self.currentcenter[0] += diffheight3d * math.sin(theta) * math.cos(phi)
            self.currentcenter[1] -= diffheight3d * math.sin(theta) * math.sin(phi)
        self.currentangle = self.molecules[-1].angle
        #print('Update mod',datetime.now())
        
        # Add vertices and edges not in existing lists
        #print(len(self.molecules[-1].verts), len(self.verts), len(startverts))
        
        self.add_verts(self.molecules[-1])
        #print('Add  verts',datetime.now())
        self.add_edges(self.molecules[-1])
        #print('Add  edges',datetime.now())
        self.add_faces(self.molecules[-1])
        #print("# edges = ", len(self.edges))
        #print('Add faces',datetime.now())
        
    def add_curvedcone(self, endradius, diffheight3d, vertexradius, vertexheight, 
                       nsegments, startedgetype='auto'):
        if self.relativedim:
            endradius *= self.maxradius
            diffheight3d *= self.maxradius
            vertexradius *= self.maxradius
            vertexheight *= self.maxradius
        if abs(endradius) > self.maxradius:
            print("Error: Curved cone end radius of ", str(endradius), 
                  "exceeds the maximum possible radius of ", str(self.maxradius))
        if abs(vertexradius) > self.maxradius:
            print("Error: Curved cone vertex radius of ", str(vertexradius), 
                  "exceeds the maximum possible radius of ", str(self.maxradius))
        # A curved cone is a series of cones with start edges that are not folded (type=none)
        # Each half of the cone is a parabola
        # The two parabolas share a vertex with a given radius and height
        # Each half of the parabola is divided into ceil(nsegments/2) sub-cones
        
        # Define start and end points in a plane passing through the axis of the form
        startpt = [self.currentradius, 0.0]
        endpt   = [endradius, diffheight3d]
        vertex  = [vertexradius, vertexheight]
        angle   = math.atan2(endpt[0]-startpt[0], endpt[1]-startpt[1]) + math.pi/2
        #print(startpt, endpt, vertex, angle)
        # Translate points such that vertex is at the origin
        startpt[0] -= vertex[0]
        startpt[1] -= vertex[1]
        endpt[0] -= vertex[0]
        endpt[1] -= vertex[1]
        #print(startpt, endpt)
        # Rotate end points to parallel with x axis
        rotstart = [startpt[0]*math.cos(-angle) + startpt[1]*math.sin(-angle),
                    startpt[1]*math.cos(-angle) - startpt[0]*math.sin(-angle)]
        rotend   = [endpt[0]*math.cos(-angle) + endpt[1]*math.sin(-angle),
                    endpt[1]*math.cos(-angle) - endpt[0]*math.sin(-angle)]
        #print(rotstart, rotend)
        # Generate new points for sub-cones between start and end points
        newpts = []
        nsegs = math.ceil(nsegments/2)
        start_a = rotstart[1] / rotstart[0]**2
        start_gap = rotstart[0] / nsegs
        for i in range(nsegs):
            x = rotstart[0] - i * start_gap
            newpts.append([x, start_a * x**2])
        newpts.append([0.0, 0.0])
        end_a = rotend[1] / rotend[0]**2
        end_gap = rotend[0] / nsegs
        for i in range(nsegs):
            x = (i+1) * end_gap
            newpts.append([x, end_a * x**2])
        #print(newpts)
        # Rotate the new points back to the correct angle and un-translate based on vertex
        rotpts = []
        for pt in newpts:
            rotpts.append([pt[0]*math.cos(angle) + pt[1]*math.sin(angle) + vertex[0],
                           pt[1]*math.cos(angle) - pt[0]*math.sin(angle) + vertex[1]])
        #print(rotpts)
        # Use all of these points to generate cones
        for i in range(1, len(rotpts)):
            if i == 1:
                startedge = startedgetype
            else:
                startedge = 'none'
            self.add_cone(rotpts[i][0], rotpts[i][1]-rotpts[i-1][1], startedge, scaleradius=False)
        
    def add_diagshift(self, endradius, offsetfract, tiltrot, startedgetype='auto'):
        if len(self.molecules) > 0:
            startverts, startedges = self.molecules[-1].endverts, self.molecules[-1].endedges
        else:
            startverts, startedges = None, None
        if self.relativedim:
            endradius *= self.maxradius
        self.molecules.append(DiagonalShift(self.ngores, self.gorewidth, self.cwrot, 
                                            self.currentradius, self.currentheight2d, self.currentheight3d, 
                                            startverts, startedges,
                                            self.currentangle, self.xyrot, self.currentaxis, self.currentcenter,
                                            startedgetype, endradius, offsetfract, tiltrot))
        # Update current radius and height to match the end of the shift
        
        diffheight3d = self.molecules[-1].verts[-1].pos3d[0,2] - self.currentheight3d
        self.cwrot = not self.cwrot
        self.currentradius = endradius
        self.currentheight2d = self.molecules[-1].endheight2d
        if self.currentaxis == [0,0,1]:
            self.currentheight3d += diffheight3d
        else:
            phi = math.atan2(self.currentaxis[1], self.currentaxis[0])
            theta = math.acos(self.currentaxis[2]/math.sqrt(self.currentaxis[0]**2 + self.currentaxis[1]**2 + self.currentaxis[2]**2))
            self.currentheight3d  += diffheight3d * math.cos(theta)
            self.currentcenter[0] -= diffheight3d * math.sin(theta) * math.cos(phi)
            self.currentcenter[1] += diffheight3d * math.sin(theta) * math.sin(phi)
        self.currentcenter = self.molecules[-1].center
        #print('center',self.currentcenter)
        self.currentangle = self.molecules[-1].angle
        if self.cwrot:
            self.xyrot -= 2 * math.pi * (self.ngores/2 - 1) / self.ngores
        else:
            self.xyrot += 2 * math.pi * (self.ngores/2 - 1) / self.ngores
        #print('Update mod',datetime.now())
        
        # Add vertices and edges not in existing lists
        self.add_verts(self.molecules[-1])
        #print('Add  verts',datetime.now())
        self.add_edges(self.molecules[-1])
        #print('Add  edges',datetime.now())
        self.add_faces(self.molecules[-1])
        #print("# edges = ", len(self.edges))
        #print('Add faces',datetime.now())
        '''
        # Update current radius and height to match the end of the cone
        self.currentradius = endradius
        self.currentheight2d = self.molecules[-1].endheight2d
        if self.currentaxis == [0,0,1]:
            self.currentheight3d = self.molecules[-1].endheight3d
        else:
            diffheight3d = self.molecules[-1].endheight3d - self.currentheight3d
            phi = math.atan2(self.currentaxis[1], self.currentaxis[0])
            theta = math.acos(self.currentaxis[2]/math.sqrt(self.currentaxis[0]**2 + self.currentaxis[1]**2 + self.currentaxis[2]**2))
            self.currentheight3d  += diffheight3d * math.cos(theta)
            self.currentcenter[0] -= diffheight3d * math.sin(theta) * math.cos(phi)
            self.currentcenter[1] += diffheight3d * math.sin(theta) * math.sin(phi)
        self.currentangle = self.molecules[-1].angle
        
        # Add vertices and edges not in existing lists
        self.add_verts(self.molecules[-1])
        self.add_edges(self.molecules[-1])
        self.add_faces(self.molecules[-1])
        '''
    
    # Add new vertices from a molecule to the master list
    #
    def add_verts(self, molecule):
        self.verts.extend(molecule.verts[molecule.nstartverts:])

    '''
    # Check whether a vertex is already in the master list
    def match_vertex(self, vertex):
        #print(vertex.pos2d, vertex.pos3d)
        for i in range(0, len(self.verts)):
            if (np.allclose(self.verts[i].pos2d, vertex.pos2d) and 
                np.allclose(self.verts[i].pos3d, vertex.pos3d)):
            #if self.verts[i].pos2d == vertex.pos2d and self.verts[i].pos3d == vertex.pos3d:
                #print('Match found', vertex.pos2d, vertex.pos3d)
                return i
        return -1
    '''
    
    # Add new edges
    def add_edges(self, molecule):
        # Correct directions of start edges
        curredges = len(self.edges)
        for i in range(molecule.nstartedges):
            self.edges[curredges - molecule.nstartedges + i].direction = molecule.startedges[i].direction
        # Add new edges
        self.edges.extend(molecule.edges[molecule.nstartedges:])
        '''
        for edge in molecule.edges:
            # Find match for each endpoint
            end1num = self.match_vertex(edge.end1)
            end2num = self.match_vertex(edge.end2)
            if end1num == -1 or end2num == -1:
                print("No vertices found for edge")
            else:
                match = self.match_edge(edge)
                if match == -1:
                    self.edges.append(el.Edge(self.verts[end1num], self.verts[end2num], edge.direction))
                else: 
                    # Update edge direction to new value
                    self.edges[match].direction = edge.direction
        '''
    '''
    def match_edge(self, edge):
        #print('Trying to match edge ', edge.end1, edge.end2)
        for i in range(0, len(self.edges)):
            #print(self.edges[i].end1, self.edges[i].end2)
            if (np.allclose(self.edges[i].end1.pos2d, edge.end1.pos2d) and
                np.allclose(self.edges[i].end1.pos3d, edge.end1.pos3d) and
                np.allclose(self.edges[i].end2.pos2d, edge.end2.pos2d) and
                np.allclose(self.edges[i].end2.pos3d, edge.end2.pos3d)):
            #if self.edges[i].end1.pos2d == edge.end1.pos2d and self.edges[i].end2.pos2d == edge.end2.pos2d and self.edges[i].end1.pos3d == edge.end1.pos3d and self.edges[i].end2.pos3d == edge.end2.pos3d:
                #print('Match found', edge.end1.pos2d, edge.end2.pos2d)
                return i
            if (np.allclose(self.edges[i].end2.pos2d, edge.end1.pos2d) and
                np.allclose(self.edges[i].end2.pos3d, edge.end1.pos3d) and
                np.allclose(self.edges[i].end1.pos2d, edge.end2.pos2d) and
                np.allclose(self.edges[i].end1.pos3d, edge.end2.pos3d)):
            #if self.edges[i].end2.pos2d == edge.end1.pos2d and self.edges[i].end1.pos2d == edge.end2.pos2d and self.edges[i].end2.pos3d == edge.end1.pos3d and self.edges[i].end1.pos3d == edge.end2.pos3d:
                #print('Match found', edge.end1.pos2d, edge.end2.pos2d)
                return i
        return -1
    '''
    def add_faces(self, molecule):
        for face in molecule.faces:
            self.faces.append(face)

        
        
        
        
        
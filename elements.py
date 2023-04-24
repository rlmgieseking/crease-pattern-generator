# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 08:49:19 2023
Elements are the basic components of folded models:
    Vertex
    Edge
    Face (to be implemented)
Each has its 2D and 3D positional information

@author: gieseking
"""

import math
import numpy as np

"""
Vertex attributes:
    2D position (x, y) 
    3D position (x, y, z)
    
Functions:
    check_valid_position_format: Make sure the 2D and 3D positions have
        the correct number of values and that all values are numbers
"""
class Vertex:
    def __init__(self, pos2d, pos3d):
        self.pos2d = np.array([pos2d]).astype(float)
        self.pos3d = np.array([pos3d]).astype(float)
        #if len(self.pos2d) != 2 or len(self.pos3d) != 3:
        #    print(print("Error: Position ", self.pos2d, self.pos3d, " is invalid"))
        self.check_valid_position_format()
        return
    
    # Check that the 2D position contains 2 values and the 3D position contains 3 values
    def check_valid_position_format(self):
        if self.pos2d.shape != (1, 2):
            print("Error: 2D position ", self.pos2d, " is invalid")
        if self.pos3d.shape != (1, 3):
            print("Error: 3D position ", self.pos3d, " is invalid")
        #print(self.pos2d, self.pos3d)
        return
    
    def translate2d(self, distance):
        distance = np.array(distance)
        self.pos2d += distance
            
    def translate3d(self, distance):
        distance = np.array(distance)
        self.pos3d += distance
    
    # Rotates clockwise around the specified axis. Angle in radians.
    def rotate(self, theta, axis):
        if axis.lower() == 'x':
            rmat = np.matrix([[ 1, 0              , 0              ],
                              [ 0, math.cos(theta), math.sin(theta)],
                              [ 0,-math.sin(theta), math.cos(theta)]])
        elif axis.lower() == 'y': 
            rmat = np.matrix([[ math.cos(theta), 0, math.sin(theta)],
                              [ 0              , 1, 0              ],
                              [-math.sin(theta), 0, math.cos(theta)]])
        elif axis.lower() == 'z': 
            rmat = np.matrix([[ math.cos(theta), math.sin(theta), 0 ],
                              [-math.sin(theta), math.cos(theta), 0 ],
                              [ 0              , 0              , 1 ]])
        self.pos3d = np.array((rmat * self.pos3d.reshape((3,1))).reshape(3))
    
"""
Edge attributes:
    Endpoints (2 Vertex objects)
    Direction ('B', 'M', 'V', '0' for border, mountain, valley, or placeholder with no line)

Functions:
    check_valid_endpoint_format: Make sure the 2D and 3D positions have
        the correct number of values and that all values are numbers
    check_valid_distance: Make sure the 2D and 3D distances between endpoints are the same
    check_valid_direction: Make sure the direction is valid
"""
class Edge:
    def __init__(self, end1, end2, direction):
        self.end1 = end1
        self.end2 = end2
        self.direction = direction
        
        self.check_valid_endpoint_format()
        self.check_valid_distance()
        
    def check_valid_endpoint_format(self):
        end1_valid = isinstance(self.end1, Vertex)
        end2_valid = isinstance(self.end2, Vertex)
        if not end1_valid or not end2_valid:
            print('Edge has an invalid end point', self.end1, self.end2)
        
    def check_valid_distance(self):
        allowed_error = 0.02
        
        dist2d = np.linalg.norm(self.end1.pos2d - self.end2.pos2d)
        dist3d = np.linalg.norm(self.end1.pos3d - self.end2.pos3d)
                
        if ((dist2d > dist3d * (1.0 + allowed_error) or dist2d < dist3d / (1.0 + allowed_error)) and
             abs(dist2d - dist3d) > allowed_error/100):
            print("Error: 2D distance and 3D distances have a significant difference")
            print("2D distance", dist2d, self.end1.pos2d, self.end2.pos2d)
            print("3D distance", dist3d, self.end1.pos3d, self.end2.pos3d)
            pass
        else:
            #print("Distances valid")
            pass
    
"""
Face attributes:
    List of corners (at least 3 Vertex objects)
    List of edges   (at least 3 Edge   objects; must be same length as vertices)
    Color (to be used later)
"""
class Face:
    def __init__(self, verts = [], edges = [], edgelist=[], color=(1.0, 1.0, 1.0)):
        self.verts = verts
        self.edges = edges
        self.color = color
        
        #self.checkValid(edgelist)
        return
    
    # If face is defined based only on a list of verts, 
    # extract edges from full list of edges
    def getEdges(self, edgelist):
        for edge in edgelist:
            if edge.end1 in self.verts and edge.end2 in self.verts:
                self.edges.append(edge)
        return
    
    # Generate list of verts from provided edges
    def getVerts(self):
        for edge in self.edges:
            if edge.end1 not in self.verts:
                self.verts.append(edge.end1)
            if edge.end2 not in self.verts:
                self.verts.append(edge.end2)
        return
    
    # Check the validity of the given face
    def checkValid(self, edgelist):
        if len(self.verts) == 0 and len(self.edges) == 0:
            print("Error: Face has no edges or vertices")
        elif len(self.verts) == 0:
            self.getVerts()
        elif len(self.edges) == 0:
            self.getEdges(edgelist)
        
        if len(self.verts) != len(self.edges):
            print("Error: Face has mismatched numbers of edges and verts")
            print(len(self.verts), len(self.edges))
            for vertex in self.verts:
                print('Vertex', vertex.pos2d)
            for edge in self.edges:
                print('Edge', edge.end1.pos2d, edge.end2.pos2d)
        if len(self.verts) < 3:
            print("Error: Face has fewer than 3 vertices")
            print(len(self.verts), len(self.edges))
        
        # Compute surface normal
        if len(self.verts) >= 3:
            vec1 = self.verts[1].pos3d - self.verts[0].pos3d
            vec2 = self.verts[2].pos3d - self.verts[0].pos3d
            self.normal = np.cross(vec1, vec2)
            self.normal /= np.linalg.norm(self.normal)
            #print(self.normal)
        # Check coplanarity of vertices
        self.coplanar = True
        for i in range(3,len(self.verts)):
            vec = self.verts[i].pos3d - self.verts[0].pos3d
            vec = vec.reshape((3,1))
            #print(abs(np.dot(self.normal, vec)))
            if abs(np.dot(self.normal, vec)) > 0.01:
                self.coplanar = False
                print('Error: Face is not planar')
                print('Normal = ', self.normal)
                print('Vec = ', vec)
                for j in range(len(self.verts)):
                    print(self.verts[j].pos3d)
        
    
    '''
    def sortVertices(self):
        if len(self.vertices) <= 3:
            return
        for i in range(1, len(self.vertices)):
            validside = False
            if (el.Edge(self.vertices[i-1], self.vertices[i]) in self.edges or
                el.Edge(self.vertices[i], self.vertices[i-1]) in self.edges):
                validside = True
    '''
    
    
    
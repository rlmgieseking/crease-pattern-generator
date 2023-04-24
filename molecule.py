# -*- coding: utf-8 -*-
"""
General Molecule class
Molecules are collections of elements for a segment of the final form
"""
import math
import numpy as np
import elements as el
from copy import deepcopy
from scipy.optimize import fsolve

class Molecule:
    def __init__(self, ngores, gorewidth, cwrot, 
                 startradius, startheight2d, startheight3d, startverts, startedges,
                 incomingangle, xyrot, axis, center, startedgetype):
        self.ngores = ngores
        self.gorewidth = gorewidth
        self.cwrot = cwrot
        self.startradius = startradius
        self.startheight2d = startheight2d
        self.startheight3d = startheight3d
        self.startverts = startverts
        if startverts:
            self.nstartverts = len(startverts)
        else: 
            self.nstartverts = 0
        self.startedges = startedges
        if startverts:
            self.nstartedges = len(startedges)
        else: 
            self.nstartedges = 0
        self.incomingangle = incomingangle
        self.xyrot = xyrot
        self.axis = axis
        self.center = center
        self.startedgetype = startedgetype
        
        self.verts = []
        self.edges = []
        self.faces = []
        return

    # Rotate all vertices in the xy plane
    def rotate_xy(self):
        for vertex in self.verts[self.nstartverts:]:
            vertex.rotate(self.xyrot, 'z')
    
    # Rotate around (0,0,0) such that z is aligned with given axis
    def rotate_axis(self):
        # phi = angle in xy plane from +x, theta = angle from +z
        phi = math.atan2(self.axis[1], self.axis[0])
        theta = math.acos(self.axis[2]/math.sqrt(self.axis[0]**2 + self.axis[1]**2 + self.axis[2]**2))
        
        for vertex in self.verts[self.nstartverts:]:
            vertex.rotate(-phi, 'z')
            vertex.rotate(theta, 'y')
            vertex.rotate(phi, 'z')
    
    # Translate to the correct offset for the final height
    def translate3d(self):
        # If startverts exist, they are already positioned correctly
        # Position only vertices after startverts
        for vertex in self.verts[self.nstartverts:]:
            vertex.translate3d([self.center[0], self.center[1], self.startheight3d])
    
    # Function to compute heights of vertices given convergence point xyz coords.
    # Height of first inner vertex is defined as 0
    # Input vector h has indices [conv3d, hi (starting at 1)]
    # conv = xy coordinates of convergence point
    # toptwist = True if paper goes up to convergence point, False if down
    def heights(self, h, verts, conv=[0,0], toptwist = True):
        fac = 1e-4
        error = []
            
        # Check deviation from coplanarity for convergence point
        conv3d = np.array([[conv[0], conv[1], h[0]]])
        #print('conv',conv3d)
        # Choose i, j, and k to be the first 3 inner vertices
        
        i, j, k = 1, 3, 5
        vec1 = verts[j].pos3d - verts[i].pos3d
        #print('jk', j, k, (j+1)//2, (k+1)//2)
        vec1[0,2] = h[(j-1)//2]
        vec2 = verts[k].pos3d - verts[i].pos3d
        vec2[0,2] = h[(k-1)//2]
        normal = np.cross(vec1, vec2)
        normal /= np.linalg.norm(normal)
        #print('pos', verts[i].pos3d, verts[j].pos3d, h[(j+1)//2], verts[k].pos3d, h[(k+1)//2])
        #print('normal',normal, vec1, vec2)
        vec = conv3d - verts[i].pos3d
        #print('vec', conv3d, verts[i].pos3d, vec)
        #print('dot', float(np.tensordot(normal, vec)))
        error.append(float(np.tensordot(normal, vec)) + fac*abs(h[i]))
            
        # Compute target paper length based on first inner vertex
        dist3d = math.sqrt((verts[1].pos3d[0,0] - conv[0])**2 +
                           (verts[1].pos3d[0,1] - conv[1])**2 +
                           (0 - h[0])**2)
        #print('dist', dist3d)
        dx2d = verts[1].pos2d[0,0] - self.gorewidth/2
        target = math.sqrt(dist3d**2 - dx2d**2)
            
        # Compute error in paper lengths for subsequent inner vertices
        for i in range(1,self.ngores):
            dist3d = math.sqrt((verts[i*2+1].pos3d[0,0] - conv[0])**2 +
                               (verts[i*2+1].pos3d[0,1] - conv[1])**2 +
                               (h[i] - h[0])**2)
            #print('dist', dist3d, target)
            dy2d = math.sqrt(dist3d**2 - dx2d**2)
            if toptwist:
                error.append(target - dy2d - h[i] + fac*abs(h[i]))
            else:
                error.append(target - dy2d + h[i] + fac*abs(h[i]))
        
        #print(error)
        return error
    
    # Find error in 3D location of a point, given its 2D distance from the convergence point
    # 3D location should be coplanar with the 3 points given in pts, with distance
    # dy from pts[0] and dx from pts[1]
    def locpt(self, pos3d, dx, dy, verts, vertind):
        error = []
        #print(pos3d, pos3d.shape)
        pos3d = pos3d.reshape((1,3))
        #print('pos3d',pos3d)
        
        # Error in height
        i, j, k = vertind[0], vertind[1], vertind[2]
        vec1 = verts[j].pos3d - verts[i].pos3d
        vec2 = verts[k].pos3d - verts[i].pos3d
        normal = np.cross(vec1, vec2)
        normal /= np.linalg.norm(normal)
        #print('normal',normal, vec1, vec2)
        intersect = (pos3d - 
                     np.array([[0,0,1]]) * np.tensordot(pos3d - verts[i].pos3d, normal)/
                     np.tensordot(np.array([[0,0,1]]), normal))
        #print('intersect',intersect)
        error.append(intersect[0,2] - pos3d[0,2])
        
        # Error in x distance
        xdist = np.linalg.norm(verts[j].pos3d - pos3d)
        error.append(xdist - dx)
        
        # Error in y distance
        ydist = np.linalg.norm(verts[i].pos3d - pos3d)
        error.append(ydist - dy)
        
        return error

    def generatehalftwist(self, verts, conv3d, toptwist=True):
        # Solve heights function to get the actual heights of the points
        if toptwist:
            init = np.array([1] + [0] * (self.ngores - 1))
        else:
            init = np.array([-1] + [0] * (self.ngores - 1))
        zpos = fsolve(self.heights, init, args=(verts, conv3d[0,0:2], toptwist))
        minz = min(0, min(zpos[1:])) if toptwist else max(0,max(zpos[1:]))
        print('zpos',zpos)
        #print('error',self.heights(zpos, verts, conv3d[0,0:2], True))
        conv3d[0,2] = zpos[0] - minz
        print('conv3d',conv3d)
        if toptwist:
            conv2d = - minz + math.sqrt((verts[1].pos3d[0,0] - conv3d[0,0])**2 +
                                        (verts[1].pos3d[0,1] - conv3d[0,1])**2 +
                                        (0 - zpos[0])**2 - 
                                        (verts[1].pos2d[0,0] - self.gorewidth/2)**2)
        else:
            conv2d = - minz - math.sqrt((verts[1].pos3d[0,0] - conv3d[0,0])**2 +
                                        (verts[1].pos3d[0,1] - conv3d[0,1])**2 +
                                        (0 - zpos[0])**2 - 
                                        (verts[1].pos2d[0,0] - self.gorewidth/2)**2)
        #print('conv',conv3d, conv2d)
        dfrac = (min(self.gorewidth - verts[1].pos2d[0,0], verts[1].pos2d[0,0])/
                 (2*abs(verts[1].pos2d[0,0] - self.gorewidth/2)))
        zpos[0] = 0.0
        #print('dfrac',dfrac)
        newmin = 0
        
        # Add new points for transition from cylinder to shift
        # Outer vertices get placeholder heights, corrected in next step
        for i in range(self.ngores*2+1):
            verts.append(deepcopy(verts[i]))
            if i % 2 == 1:
                verts[-1].pos3d[0,2] = - minz + zpos[(i-1)//2]
                verts[-1].pos2d[0,1] = - minz + zpos[(i-1)//2]
            else:
                if self.cwrot:
                    if i == 0:
                        h = - minz + zpos[i//2] + dfrac * (zpos[i//2] - zpos[-1])
                    elif i == self.ngores*2:
                        h = - minz + zpos[0] + dfrac * (zpos[0] - zpos[(i-2)//2])
                    else:
                        h = - minz + zpos[i//2] + dfrac * (zpos[(i)//2] - zpos[(i-2)//2])
                else:
                    if i == 0:
                        h = - minz + zpos[-1] - dfrac * (zpos[i//2] - zpos[-1])
                    elif i == self.ngores*2:
                        h = - minz + zpos[(i-2)//2] - dfrac * (zpos[0] - zpos[(i-2)//2])
                    else:
                        h = - minz + zpos[(i-2)//2] - dfrac * (zpos[(i)//2] - zpos[(i-2)//2])
                verts[-1].pos3d[0,2] = h
                verts[-1].pos2d[0,1] = h
            newmin = min(newmin, verts[-1].pos3d[0,2]) if toptwist else max(newmin, verts[-1].pos3d[0,2])
        # If the new lowest vertex is below 0, correct all the vertices
        if newmin != 0:
            for i in range(self.ngores*2+1):
                verts[i + self.ngores*2+1].pos3d[0,2] -= newmin
                verts[i + self.ngores*2+1].pos2d[0,1] -= newmin
            conv2d -= newmin
            conv3d[0,2] -= newmin
        print('conv3d',conv3d)
        # Add vertices surrounding convergence point 
        # Inner vertices get exact positions
        # Outer vertices get placeholder 3D positions to be fixed in the next step
        for i in range(0, self.ngores * 2 + 1, 2):
            i1 = i + self.ngores*2     if i > 0 else self.ngores*4
            i3 = i + self.ngores*2 + 2 if i < self.ngores*2 else self.ngores*2 + 2
            x13 = verts[i1].pos3d[0,0] - verts[i3].pos3d[0,0]
            y13 = verts[i1].pos3d[0,1] - verts[i3].pos3d[0,1]
            z13 = verts[i1].pos3d[0,2] - verts[i3].pos3d[0,2]
            d13 = math.sqrt(x13**2 + y13**2 + z13**2)
            if self.cwrot:
                xd = (self.gorewidth/2) * x13 / d13
                yd = (self.gorewidth/2) * y13 / d13
                zd = (self.gorewidth/2) * z13 / d13
            else:
                xd = -(self.gorewidth/2) * x13 / d13
                yd = -(self.gorewidth/2) * y13 / d13
                zd = -(self.gorewidth/2) * z13 / d13
            verts.append(el.Vertex([self.gorewidth*i/2, conv2d],
                                        [conv3d[0,0] + xd, conv3d[0,1] + yd, conv3d[0,2] + zd]))
            verts.append(el.Vertex([self.gorewidth*(i+1)/2 , conv2d],
                                        [conv3d[0,0], conv3d[0,1], conv3d[0,2]]))
        verts.pop(-1)
        
        # Fix positions of outer vertices surrounding the convergence point
        vertind = [0,0,0]
        for i in range(self.ngores+1):
            vertind[0] = i*2 + self.ngores*2 + 1
            vertind[1] = i*2 + self.ngores*4 + 3 if i < self.ngores else i*2 + self.ngores*4 + 1
            vertind[2] = vertind[0] + 2 if i < self.ngores else vertind[0] + 2 - self.ngores*2
            v = i*2 + self.ngores*4 + 2
            #print('vertind', i, v, vertind)
            dx = abs(verts[vertind[1]].pos2d[0,0] - verts[v].pos2d[0,0])
            dy = abs(verts[vertind[0]].pos2d[0,1] - verts[v].pos2d[0,1])
            #print(verts[v].pos3d)
            verts[v].pos3d = np.array([fsolve(self.locpt,verts[v].pos3d, 
                                         args=(dx, dy, verts, vertind))])
            #print(self.locpt(verts[v].pos3d, dx, dy, verts, vertind))
        
        return verts, conv3d, conv2d
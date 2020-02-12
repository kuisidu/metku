# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 13:58:54 2018

Class for modelling cold-formed sections, according to 

EN 1993-1-3, Annex C

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt

class OpenProf:
    """
    OpenProf Class for general open cold-formed profiles
    Main purpose is to calculate cross-sectional properties as
    described in EN 1993-1-3, Annex C.
    """

    def __init__(self,nodes,t):
        """ Constructor
            input:
                nodes .. list of nodal coordinates
                    nodes[i] = [y_i, z_i]
                t .. array of wall thicknesses or a single value
        """        
        n = len(nodes)
        self.nodes = np.array(nodes)
        if isinstance(t,float):
            self.t = t*np.ones(n-1)
        else:
            self.t = np.array(t)
            
        self.n = n
        
        self.segments = []
        
        for i in range(1,n):                      
            self.add_segment(self.nodes[i],self.nodes[i-1],self.t[i-1])
            #self.add_segment(self.nodes[i,:],self.nodes[i-1,:],self.t[i-1])
    
    def y_coord(self):
        """ Return y coordinates """
        return self.nodes[:,0]
    
    def z_coord(self):
        """ Return z coordinates """
        return self.nodes[:,1]
    
    def add_segment(self,n1,n2,t):
        """ Adds new line segment to cross section """
        new_segment = Segment(n1,n2,t)
        self.segments.append(new_segment)
                  
    def draw(self):
        """ draw nodes """
        for node in self.nodes:
            plt.plot(node[0],node[1],'ro')            
        
        for i in range(self.n):
            plt.text(self.nodes[i][0],self.nodes[i][1],str(i))
            
        
        """ draw line segments """
        for s in self.segments:            
            X = s.nodes
            plt.plot(X[:,0],X[:,1],'b')
            Xmid = X[0,:]+0.5*(X[1,:]-X[0,:])
            plt.text(Xmid[0],Xmid[1],str(i))
                
    
        ygc, zgc, A = self.centroid()        
        ysc, zsc = self.shear_center()
        plt.plot(ygc,zgc,'+r')
        plt.plot(ysc,zsc,'or')
        
        plt.axis('equal')
        plt.show()        

    def area(self):        
        """ Cross-sectional area """
        A = 0.0
            
        for s in self.segments:
            A += s.area()        
        return A
                
    def first_moments(self):
        """ First moments of area  """                
        Sy0 = 0.0        
        Sz0 = 0.0
        for s in self.segments:            
            Sy0 += s.first_moment_y()
            Sz0 += s.first_moment_z()
        return Sy0, Sz0
            
    def centroid(self):
        """ Centroid of area """
        Sy0, Sz0 = self.first_moments()
        A = self.area()
        zgc = Sy0/A
        ygc = Sz0/A
            
        return ygc, zgc, A
    
    def second_moments(self):
        """ Second moments of area """                
        Iy0 = 0.0        
        Iz0 = 0.0
        for s in self.segments:            
            Iy0 += s.second_moment_y()
            Iz0 += s.second_moment_z()
        return Iy0, Iz0
    
    def second_moments_gc(self):
        """ Second moments of area through the center of gravity """
        ygc, zgc, A = self.centroid()
        Iy0, Iz0 = self.second_moments()
        Iy = Iy0-A*zgc**2
        Iz = Iz0-A*ygc**2
        return Iy, Iz
    
    def product_moment(self):
        """ Product moment with respect to original y- and z-axes """
        Iyz0 = 0.0
        for s in self.segments:            
            Iyz0 += s.product_moment()
            
        return Iyz0
    
    def product_moment_gc(self):
        """ Product moment with respect to center of gravity """
        A = self.area()
        Sy0, Sz0 = self.first_moments()
        Iyz0 = self.product_moment()
        return Iyz0-Sy0*Sz0/A
    
    def sectorial_coordinates(self):
        """ Determine sectorial coordinates """
        n = self.n
        y = self.y_coord()
        z = self.z_coord()
        w = np.zeros(n)
        for i in range(n-1):
            w[i+1] = w[i] + y[i]*z[i+1]-y[i+1]*z[i]
            
        return w
    
    def principal_axes(self):
        """ Principal axes """
        Iy, Iz = self.second_moments_gc()
        Iyz = self.product_moment_gc()
        d = Iz-Iy
        if abs(d) > 1e-8:
            a = .5*np.atan(2*Iyz/d)
        else:
            a = 0
            
        I1 = 0.5*(Iy+Iz+np.sqrt(d^2+4*Iyz^2))
        I2 = .5*(Iy+Iz-np.sqrt(d^2+4*Iyz^2))
        
        return [I1,I2], a
    
    def mean_sectorial_moment(self):
        n = self.n
        dA = [s.area() for s in self.segments]
        Iw = 0
        w = self.sectorial_coordinates()
        
        for i in range(n-1):
            Iw += 0.5*dA[i]*(w[i+1]+w[i])
            
        return Iw
    
        
    def mean_sectorial_coordinate(self):
        Iw = self.mean_sectorial_moment()
        A = self.area()
        return Iw/A
    
    def sectorial_constants(self):
        dA = [s.area() for s in self.segments]
        Iyw0 = 0.0
        Izw0 = 0.0
        Iww0 = 0.0
        
        Iw = self.mean_sectorial_moment()
        A = self.area()
        Sy0, Sz0 = self.first_moments()
        
        y = self.y_coord()
        z = self.z_coord()
        w = self.sectorial_coordinates()
        
        for i in range(self.n-1):
            Iyw0 += dA[i]*(2*y[i]*w[i]+2*y[i+1]*w[i+1]+y[i]*w[i+1]+y[i+1]*w[i])/6
            Izw0 += dA[i]*(2*z[i]*w[i]+2*z[i+1]*w[i+1]+z[i]*w[i+1]+z[i+1]*w[i])/6
            Iww0 += dA[i]*(w[i+1]**2+w[i]**2+w[i+1]*w[i])/3                    
    
        Iyw = Iyw0 - Sz0*Iw/A
        Izw = Izw0 - Sy0*Iw/A
        Iww = Iww0 - Iw**2/A
        
        return Iyw, Izw, Iww
    
    def shear_center(self):
        
        Iyw, Izw, Iww = self.sectorial_constants()
        Iy, Iz = self.second_moments_gc()
        Iyz = self.product_moment_gc()
        
        d = Iy*Iz-Iyz**2
        
        if d > 1e-8:
            y_sc = (Izw*Iz-Iyw*Iyz)/d
            z_sc = (-Iyw*Iy+Izw*Iyz)/d
        else:
            y_sc = 0.0
            z_sc = 0.0
        
        return y_sc, z_sc
    
    def warping_constant(self):
        
        y_sc, z_sc = self.shear_center()        
        Iyw, Izw, Iww = self.sectorial_constants()
        
        return Iww + z_sc*Iyw - y_sc*Izw
    
    def torsion_constant(self):
        dA = np.array([s.area() for s in self.segments])
        
        return sum(dA*np.array(self.t)**2/3)
    
    def torsion_modulus(self):
        
        return self.torsion_constant()/min(self.t)
"""        
        # adds point x to the section and splits the
        # line segment that contains x
        def self = Split(self,x)
            n = find(self.y == x(1))
            if isempty(n)
                # find points with matching z-coordinate
                n = find(self.z==x(2))
                # of these points, find the ones that have lower
                # value of y
                nlower = n(self.y(n) < x(1))
                n1 = nlower(1)-1
                n2 = nlower(1)
            
            n1
            n2
            self.y = [self.y(1:n1)x(1)self.y(n2:)]
            self.z = [self.z(1:n1)x(2)self.z(n2:)]
            self.t = [self.t(1:n1)self.t(n1)self.t(n1+1:)]
"""
        
class Segment:
    """ Line segment that is a part of a cross section """
    
    def __init__(self,n1,n2,t):
        """ Constructor
            input:
                n1 .. node 1 (list of coordinates)
                n2 .. node 2 (list of coordinates)
                t .. thickness
        """
        
        """ nodes is a numpy array with nodal coordinates
            nodes[:,0] .. y coordinates
            nodes[:,1] .. z coordinates
        """
        self._n1 = n1
        self._n2 = n2
        self.nodes = np.array([n1,n2])
        self.t = t
        
    def length(self):
        """ Length """
        return np.sqrt(sum((self.nodes[1,:]-self.nodes[0,:])**2))

    def area(self):
        """ Cross-sectional area """
        return self.t*np.sqrt(sum((self.nodes[1,:]-self.nodes[0,:])**2))
    
    def first_moment_y(self):
        """ First moment of area with respect to y axis"""
        dA = self.area()
        return sum(self.nodes[:,1])*dA*0.5
    
    def second_moment_y(self):
        """ Second moment of area with respect to y axis"""
        dA = self.area()
        return (sum(self.nodes[:,1]**2)+self.nodes[0,1]*self.nodes[1,1])*dA/3
    
    def first_moment_z(self):
        """ First moment of area with respect to z axis"""
        dA = self.area()
        return sum(self.nodes[:,0])*dA*0.5
    
    def second_moment_z(self):
        """ Second moment of area with respect to z axis"""
        dA = self.area()
        return (sum(self.nodes[:,0]**2)+self.nodes[0,0]*self.nodes[1,0])*dA/3
    
    def product_moment(self):
        """ Product moment of areaq with respect to original y- and z-axes """
        dA = self.area()
        Iyz0 = (sum(2*self.nodes.prod(1))+self.nodes[0,0]*self.nodes[1,1]+\
           self.nodes[0,1]*self.nodes[1,0])*dA/6
        return Iyz0

if __name__ == "__main__":
    
    n = []
    n.append([20,3])
    n.append([20,0])
    n.append([0,0])
    n.append([0,60])
    n.append([20,60])
    n.append([20,57])
    
    C = OpenProf(n,1.5)
    C.draw()
    
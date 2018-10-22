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
                nodes .. array of nodal coordinates
                    nodes[i][0] .. y-coordinates
                    nodes[i][1] .. z-coordinates
                t .. array of wall thicknesses or a single value
        """        
        n = len(nodes)
        self.nodes = np.array(nodes)
        if len(t) == 1:
            self.t = t[0]*np.ones(n-1)
        else:
            self.t = np.array(t)
            
        self.n = n
        
        self.segments = []
        
        for i in range(1,n):                      
            self.add_segment(self.nodes[i,:],self.nodes[i-1,:],self.t[i-1])
    
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
            plt.text(self.nodes[i,0],self.nodes[i,1],str(i))
            
        
        """ draw line segments """
        for s in self.segments:            
            X = s.nodes
            plt.plot(X[:,0],X[:,1],'b')
            Xmid = X[0,:]+0.5*(X[1,:]-X[0,:])
            plt.text(Xmid[0],Xmid[1],str(i))
                
    
        ygc, zgc = self.centroid()        
        #sc = self.ShearCenter
        plt.plot(ygc,zgc,'+r')
        #plot(sc(1),sc(2),'or')
        
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
        w = np.zeros(n,1)
        for i in range(n):
            w[i+1] = w[i] + y[i]*z[i+1]-y[i+1]*z[i]
        
        return w
        
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
        
"""                
        
            
        
        
        def Iw = Imean(self)
           dA = self.Areas
           Iw = 0
           w = self.SectorialCoordinates
           for i = 1:self.n-1
               Iw = Iw+0.5*dA(i)*(w(i+1)+w(i))
            
        
        
        def wMean = SectorialMean(self)
            Iw = self.Imean
            A = self.Area
            wMean = Iw/A
        
        
        def Iw0 = SectorialConstants0(self)
            dA = self.Areas
            Iw0 = [000]
            y = self.y
            z = self.z
            w = self.SectorialCoordinates
            for i = 1:self.n-1
                Iw0 = Iw0 + dA(i)*[...
                    (2*y(i)*w(i)+2*y(i+1)*w(i+1)+y(i)*w(i+1)+y(i+1)*w(i))/6...
                    (2*z(i)*w(i)+2*z(i+1)*w(i+1)+z(i)*w(i+1)+z(i+1)*w(i))/6...
                    (w(i+1)^2+w(i)^2+w(i+1)*w(i))/3
                    ]
            
        
        
        # Iww = [Iyw,Izw,Iww] (w = omega)
        def Iww = SectorialConstants(self)
            A = self.Area
            S0 = self.MomentArea
            Iw = self.Imean
            Iw0 = self.SectorialConstants0
            Iww = Iw0 - [S0([2,1])Iw]*Iw/A
        
        
        def [sc,varargout] = ShearCenter(self)
            I2 = self.SecondMomentArea
            Iy = I2(1)
            Iz = I2(2)
            Iyz = I2(3)
            denom = Iy*Iz-Iyz^2
            if abs(denom) > 1e-5
                Iww = self.SectorialConstants
                Iyw = Iww(1)
                Izw = Iww(2)
                sc = [Izw*Iz-Iyw*Iyz-Iyw*Iy+Izw*Iyz]/denom
                if nargout > 1
                    varargout{1} = Iww
                
            
        
        
        def Iw = WarpingConstant(self)
            [sc,Iww] = self.ShearCenter
            Iw = Iww(3)+sc(2)*Iww(1)-sc(1)*Iww(2)
        
        
        def [It,Wt] = TorsionConstant(self)
            t = self.t
            dA = self.Areas
            It = dot(dA,t.^2/3)
            Wt = It/min(t)
        
        
        def [Ip,a] = PrincipalAxis(self)
            I2 = self.SecondMomentArea
            Iy = I2(1)
            Iz = I2(2)
            Iyz = I2(3)
            d = Iz-Iy
            if abs(d) > 1e-3
                a = .5*atan(2*Iyz/d)
                if abs(a) < 1e-8
                    a = 0
                                    
            else
                a = 0
            
            Ip(1) = .5*(Iy+Iz+sqrt(d^2+4*Iyz^2))
            Ip(2) = .5*(Iy+Iz-sqrt(d^2+4*Iyz^2))            
 """       
 
class Segment:
    """ Line segment that is a part of a cross section """
    
    def __init__(self,n1,n2,t):
        """ Constructor
            input:
                n1 .. node 1 (coordinates)
                n2 .. node 2 (coordinates)
                t .. thickness
        """
        
        """ nodes is a numpy array with nodal coordinates
            nodes[:,0] .. y coordinates
            nodes[:,1] .. z coordinates
        """
        self.nodes = np.array([n1,n2])
        self.t = t

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
    
    C = OpenProf(n,[1.5])
    C.draw()
    
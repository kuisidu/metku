# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 13:58:54 2018

Class for modelling cold-formed sections, according to 

EN 1993-1-3, Annex C

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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
        
        self.nodes = []
        for node in nodes:
            self.nodes.append(Node(y=node[0],z=node[1]))
        
        #self.nodes = np.array(nodes)
        if isinstance(t,float):
            self.t = t*np.ones(n-1)
        else:
            self.t = np.array(t)
                        
        self.segments = []
        
        for i in range(1,n):                      
            self.add_segment(self.nodes[i-1],self.nodes[i],self.t[i-1])
            #self.add_segment(self.nodes[i,:],self.nodes[i-1,:],self.t[i-1])

    @property
    def n(self):
        """ Number of nodes """
        return len(self.nodes)
    
    def y_coord(self):
        """ Return y coordinates """
        return np.array([node.y for node in self.nodes])
        #return self.nodes[:,0]
    
    def z_coord(self):
        """ Return z coordinates """
        return np.array([node.z for node in self.nodes])
        #return self.nodes[:,1]
    
    def add_segment(self,n1,n2,t):
        """ Adds new line segment to cross section """
        new_segment = Segment(n1,n2,t)
        self.segments.append(new_segment)
                  
    def draw(self,node_labels=False,seg_labels=False,axes_on=True,
             sc_on=True,gc_on=True,coord_axes=False):
        
        fig, ax = plt.subplots(1)         
        """ draw nodes """
        for i, node in enumerate(self.nodes):
            if node_labels:
                node.draw(str(i),ax)
            else:
                node.draw(axes=ax)
            #plt.plot(node[0],node[1],'ro')            
            #plt.text(node[0],node[1],str(i))
        
        #for i in range(self.n):
        #    plt.text(self.nodes[i][0],self.nodes[i][1],str(i))
            
        
        """ draw line segments """
        for i, s in enumerate(self.segments):
            if seg_labels:
                s.draw(str(i),axes=ax)
            else:
                s.draw(axes=ax)
            #X = s.nodes
            #plt.plot(X[:,0],X[:,1],'b')
            #Xmid = X[0,:]+0.5*(X[1,:]-X[0,:])
            #plt.text(Xmid[0],Xmid[1],str(i))
                
    
        if gc_on:
            ygc, zgc, A = self.centroid()        
            ax.plot(ygc,zgc,'+r')
        
        if sc_on:
            ysc, zsc = self.shear_center()        
            ax.plot(ysc,zsc,'or')
        
        ax.axis('equal')
        
        if axes_on is False:
            plt.axis('off')
        #ax.show()        

        if coord_axes:
            x = self.y_coord()
            y = self.z_coord()
            dx = 0.5*(np.max(x)-np.min(x))
            dy = 0.5*(np.max(y)-np.min(y))
            hw = 3
            ax.arrow(0,0,min(dx,dy),0,head_width=hw,facecolor='k')
            ax.text(min(dx,dy),0.05*dy,'y',verticalalignment='bottom')
            ax.arrow(0,0,0,min(dx,dy),head_width=hw,facecolor='k')
            ax.text(-0.05*dy,min(dx,dy),'z',horizontalalignment='right')

        return fig

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
    
    def Iy0(self):
        """ Second moment of area with respect to Y-axis """
        Iy0 = 0.0        
        for s in self.segments:            
            Iy0 += s.second_moment_y()           
        return Iy0
    
    def Iz0(self):
        """ Second moment of area with respect to Z-axis """
        Iz0 = 0.0        
        for s in self.segments:            
            Iz0 += s.second_moment_z()           
        return Iz0
    
    def Iy(self):
        """ Second moment of area with respect to centroid """
        ygc, zgc, A = self.centroid()
        Iy0 = self.Iy0()
        return Iy0-A*zgc**2
    
    def Iz(self):
        """ Second moment of area with respect to centroid """
        ygc, zgc, A = self.centroid()
        Iz0 = self.Iz0()
        return Iz0-A*ygc**2
    
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
        #print(y,z)
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
    
    def split_segment(self,n,s):
        """ Splits a segment into two by adding a point 
            input:
                n .. index of the segment
                s .. local coordinate along the segment, where
                     the new point is to be added        
        """
        
        seg = self.segments[n]
        
        # Get nodal coordinates
        nodal_coords = seg.coord        
        new_point = nodal_coords[0] + s*(nodal_coords[1]-nodal_coords[0])    
        new_node = Node(new_point[0],new_point[1])
        
        #print(seg.nodes[0].y,seg.nodes[0].z)
        
        new_seg1 = Segment(seg.nodes[0],new_node,seg.t)
        new_seg2 = Segment(new_node,seg.nodes[1],seg.t)
        
        
        self.nodes.insert(n+1,new_node)
        self.segments.remove(seg)
        self.segments.insert(n,new_seg1)
        self.segments.insert(n+1,new_seg2)
    
    def export_cufsm(self):
        """ Write data to be exported to cufsm in Matlab """
        
        """ Write nodes """
        print("Nodes:")
        for i, node in enumerate(self.nodes):        
            print('{0:3g} {1:6.4f} {2:6.4f} 1 1 1 1 1.000'.format(i+1,node.y,node.z))

        """ Write elements """
        print("Elements:")
        for i, seg in enumerate(self.segments):        
            print('{0:3g} {1:6.4f} {2:6.4f} {3:6.5f} 100'.format(i + 1, i + 1, i + 2, seg.khp))

class Node:
    """ Node of an open section profile """
    
    def __init__(self,y,z):
        """ Constructor 
            y is the horizontal axis and
            z is the vertical axis
        """
        
        self.y = y
        self.z = z
    
    def __repr__(self):
        
        return "y = {0:4.3f}, z = {1:4.3f}".format(self.y,self.z)
    
    @property
    def coord(self):
        """ Return numpy array of coordinates """
        return np.array([self.y,self.z])
    
    def draw(self,label=None,axes=None):        
        
        if axes is None:
            plt.plot(self.y,self.z,'ok',markersize=4)
            if label is not None:
                plt.text(self.y,self.z,label)
        else:
            axes.plot(self.y,self.z,'ok',markersize=4)
            if label is not None:
                axes.text(self.y,self.z,label)
        
        
        
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
        
    def __repr__(self):
        
        n1 = 'n1 = [{0:4.2f}, {1:4.2f} '.format(self._n1.y,self._n1.z)
        n2 = 'n2 = [{0:4.2f}, {1:4.2f}] '.format(self._n2.y,self._n2.z)
        
        return n1 + n2 + " t = {0:4.2f}".format(self.t)
    
    @property
    def coord(self):
        """ Nodal coordinates """
        return np.array([self._n1.coord, self._n2.coord])
                   
    @property
    def y(self):
        """ Y-coordinates of the segment """
        return np.array([self._n1.y,self._n2.y])
    
    @property
    def z(self):
        """ Z-coordinates of the segment """
        return np.array([self._n1.z,self._n2.z])
    
    def mid_point(self):
        """ coordinates of the middle point of the segment """
        
        return self._n1.coord + 0.5*(self._n2.coord-self._n1.coord)
    
    
    def length(self):
        """ Length """
        return np.linalg.norm(self._n2.coord-self._n1.coord)
        #return np.sqrt(sum((self.nodes[1,:]-self.nodes[0,:])**2))

    def area(self):
        """ Cross-sectional area """
        return self.t*self.length()
        #return self.t*np.sqrt(sum((self.nodes[1,:]-self.nodes[0,:])**2))
    
    def first_moment_y(self):
        """ First moment of area with respect to y axis"""
        dA = self.area()
        return sum(self.z)*dA*0.5
        #return sum(self.nodes[:,1])*dA*0.5
    
    def second_moment_y(self):
        """ Second moment of area with respect to y axis"""
        dA = self.area()
        z = self.z
        return (z[1]**2+z[0]**2+z[1]*z[0])*dA/3
        #return (sum(self.nodes[:,1]**2)+self.nodes[0,1]*self.nodes[1,1])*dA/3
    
    def first_moment_z(self):
        """ First moment of area with respect to z axis"""
        dA = self.area()
        return sum(self.y)*dA*0.5
        #return sum(self.nodes[:,0])*dA*0.5
    
    def second_moment_z(self):
        """ Second moment of area with respect to z axis"""
        dA = self.area()
        y = self.y
        return (y[1]**2+y[0]**2+y[1]*y[0])*dA/3
        #return (sum(self.nodes[:,0]**2)+self.nodes[0,0]*self.nodes[1,0])*dA/3
    
    def product_moment(self):
        """ Product moment of areaq with respect to original y- and z-axes """
        dA = self.area()
        y = self.y
        z = self.z
        Iyz0 = (2*y[0]*z[0] + 2*y[1]*z[1] + y[0]*z[1] + y[1]*z[0])*dA/6
        #Iyz0 = (sum(2*self.nodes.prod(1))+self.nodes[0,0]*self.nodes[1,1]+\
        #   self.nodes[0,1]*self.nodes[1,0])*dA/6
        return Iyz0
    
    def draw(self,label=None,axes=None):
        """ Draw segment """
        
        if axes is None:
            plt.plot(self.y,self.z,'b')
            if label is not None:            
                mid = self.mid_point()
                plt.text(mid[0],mid[1],label)
        else:
            lw0 = 3
            
            axes.plot(self.y,self.z,linewidth=self.t*lw0,color='b')
            
            #xy = (self.y[0],self.z[0])        
            #seg1 = patches.Rectangle(xy,width=self.length(),height=0.5*self.t,fill=False)            
            #seg2 = patches.Rectangle(xy,width=self.length(),height=0.5*self.t,fill=False)            
            #axes.add_patch(seg)
        
            if label is not None:     
                if self.t > 0.0:
                    mid = self.mid_point()
                    axes.text(mid[0],mid[1],label)        


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
    
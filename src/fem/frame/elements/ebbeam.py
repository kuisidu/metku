# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:14:54 2018

Euler-Bernoulli beam element in 2D

@author: kmela
"""

import numpy as np
from fem.frame.frame_fem import Element

class EBBeam(Element):
    """ Euler-Bernoulli beam element 
        
        Attributes:
            bending_moment -- bending moment in the end nodes [list]
            shear_force -- shear force in the end nodes [list]
    """
    
    def __init__(self,n1,n2,section,material):
        Element.__init__(self,n1,n2,section,material)
        
        self.bending_moment = [0.0,0.0]
        self.shear_force = [0.0,0.0]        
    
    def transformation_matrix(self):
        """ From local to global coordinates """
        c = self.direction_cosines()
        L = np.zeros((6,6))
        L[0,0:2] = c[0:2]
        L[1,0] = -L[0,1]
        L[1,1] = L[0,0]
        L[2,2] = 1.0
        L[3:6,3:6] = L[0:3,0:3]
        return L
    
    def global_dofs(self):
        """ Get numbering of element globall degrees of freedom """        
        return  self.nodes[0].dofs + self.nodes[1].dofs
    
    def local_stiffness_matrix(self):
        """ Stiffness matrix in local coordinates """
        E = self.material.young
        A = self.section.area
        I1 = self.section.Iy
        Le = self.length()
        
        rodc = E*A/Le
        EI1 = E*I1
        k0 = np.zeros((6,6))
        #
        k0[[0,3],[0,3]] = rodc
        k0[[0,3],[3,0]] = -rodc
        #
        bend1 = 12*EI1/Le**3
        k0[[1,4],[1,4]] = bend1
        k0[[1,4],[4,1]] = -k0[1,1]
        #
        bend2 = 6*EI1/Le**2
        k0[[1,2,1,5],[2,1,5,1]] = bend2
        #
        k0[[2,4,4,5],[4,2,5,4]] = -bend2
        #
        k0[[2,5],[5,2]] = 2*EI1/Le
        #
        k0[[2,5],[2,5]] = 4*EI1/Le
          
        return k0
    
    def stiffness_matrix(self):
        """ Compute the stiffness matrix """
        
        k0 = self.local_stiffness_matrix()                
        
        #k0 = CheckReleases[fem.elem[N],k0];
                        
        # local-global matrix
        L = self.transformation_matrix()
        
        # globaali elementin matriisi
        ke = L.transpose().dot(k0.dot(L))
        
        return ke

    def local_geometric_stiffness_matrix(self):
        """ Geometric stiffness matrix in local coordinates
            From: Cook et. al 1989, Section 14.2
        """
        P = self.axial_force[1]
        Le = self.length()
        
        g0 = P/30/Le
        keg0 = g0*np.array([[36,3*Le,-36,3*Le],
                            [3*Le,4*Le**2,-3*Le,-Le**2],
                            [-36,-3*Le,36,-3*Le],
                            [3*Le,-Le**2,-3*Le,4*Le**2]])
        kg0 = np.zeros([6,6])
        q = [1,2,4,5]
        kg0[np.ix_(q,q)] = keg0
        
        return kg0
    
    def geometric_stiffness_matrix(self):
        """ Geometric stiffness matrix in global coordinates
            From: Cook et. al (1989), Section 14.2
        """       
        
        # local-global transformation matrix
        L = self.transformation_matrix()
        
        kG0 = self.local_geometric_stiffness_matrix()
        
        kG = L.transpose().dot(kG0.dot(L))
        
        return kG

    def equivalent_nodal_loads(self,q):
        """ Equivalent nodal loads for load in vector q
        """
                
        floc = np.zeros(6)
        
        T = self.transformation_matrix()
        L = self.length()
        
        # Load vector transformed into element local coordinate system
        qloc = T[:3,:3].dot(q)
        
        
        # Construct nodal load in local coordinates
        # qloc[0,1,2] = axial, shear, moment of node 1
        # qloc[3,4,5] = axial, shear, moment ofnode 2
        # Axial force
        floc[[0,3]] = 0.5*qloc[0]*L
        # Shear force
        floc[[1,4]] = 0.5*qloc[1]*L
        # Moment
        floc[2] = qloc[1]*L**2/12.0
        floc[5] = -floc[2]
        
        # Store local loads
        self.floc = floc
        
        #print(self.floc)
        
        fglob = T.transpose().dot(floc)
        return fglob

    def shape_fun(self,x):
        """ Evaluate Hermitean shape functions at local coordinate x """
        L = self.length()
        N = np.zeros(4)
        N[0] = 1-3*x**2/L**2 + 2*x**3/L**3
        N[1] = x-2*x**2/L+x**3/L**2
        N[2] = 3*x**2/L**2-2*x**3/L**3
        N[3] = -x**2/L**2 + x**3/L**2
         
        return N


# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:12:04 2018

Rod element, for pin-jointed members with only axial forces

@author: kmela
"""

import numpy as np

try:
    from src.framefem import Element
except:
    from framefem import Element
        


class Rod(Element):
    """ Rod element, carries only axial loads """

    def __init__(self, n1, n2, section, material):
        Element.__init__(self, n1, n2, section, material)
        
        # Dimension of the element
        self.dim = len(n1.coord)

    def transformation_matrix(self):
        """ From local to global coordinates """
        c = self.direction_cosines()
        L1 = np.hstack((c, np.zeros(len(c))))
        L2 = np.hstack((np.zeros(len(c)), c))
        L = np.vstack((L1, L2))
        return L

    def global_dofs(self):
        """ Get numbering of element global degrees of freedom """
        d = self.dim            
        return self.nodes[0].dofs[:d] + self.nodes[1].dofs[:d]


    def init_dofs(self):
        """ For truss members there are no rotational degrees of freedom
            so these can be set to be neglected (=1)
        """
        for n in self.nodes:            
            if self.dim == 2:
                n.dofs[2] = 1
            else:
                n.dofs[3] =  1
                n.dofs[4] =  1
                n.dofs[5] =  1
            
            

    def local_stiffness_matrix(self):
        """ Stiffness matrix in local coordinates """
        E = self.material.young
        A = self.section.area
        Le = self.length()

        """ Stiffness matrix in local coordinates """
        return E * A / Le * np.array([[1, -1], [-1, 1]])

    def stiffness_matrix(self):
        """ Compute the stiffness matrix """
        k0 = self.local_stiffness_matrix()

        """ Transformation matrix """
        L = self.transformation_matrix()

        """ Transform stiffness matrix to global coordinates """
        ke = L.transpose().dot(k0.dot(L))

        return ke

    def local_geometric_stiffness_matrix(self):
        """ Geometric stiffness matrix in local coordinates
            From: Cook et. al 1989, Section 14.2
            
            TODO: Modify for 3D!
        """
        P = self.axial_force[1]
        Le = self.length()

        return P / Le * np.array([[0, 0, 0, 0], [0, 1, 0, -1], [0, 0, 0, 0], [0, -1, 0, 1]])

    def geometric_stiffness_matrix(self):
        """ Geometric stiffness matrix in global coordinates
            From: Cook et. al (1989), Section 14.2
        """
        """ Generate first the transformation matrix
            that takes into account two DOFs per node
            From: Cook et. al (1989), Section 7.5
            
            TODO: Modify for 3D!
        """
        c = self.direction_cosines()
        L = np.array([[c, 0, 0], [-c[1], c[0], 0, 0], [0, 0, c], [0, 0, -c[1], c[0]]])

        kG0 = self.local_geometric_stiffness_matrix()

        kG = L.transpose().dot(kG0.dot(L))

        return kG

    def equivalent_nodal_loads(self, q):
        """ Equivalent nodal loads for load in vector q
            Ignore bending moment
            
            TODO: Modify for 3D!
        """

        fglob = np.zeros(4)

        dc = self.direction_cosines()
        L = self.length()

        fglob[[0, 2]] = 0.5 * dc[1] * q[0] * L
        fglob[[1, 3]] = 0.5 * dc[0] * q[1] * L

        return fglob


    def internal_forces(self):
        """ Calculates internal forces (axial force)
        """
        q = self.local_displacements()
        E = self.material.young
        A = self.section.area
        L = self.length()
        self.axial_force[0] = E * A / L * (q[1] - q[0])


    def nodal_displacements(self):
        """ Get nodal displacements of an element
            Requires previously performed structural analysis such
            that nodal displacements are available.
        """

        d = self.dim
        return np.hstack((self.nodes[0].u[:d], self.nodes[1].u[:d]))
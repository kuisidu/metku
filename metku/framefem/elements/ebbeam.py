# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:14:54 2018

Euler-Bernoulli beam element in 2D

@author: kmela
"""

import numpy as np
from functools import lru_cache

CACHE_BOUND = 2**10

try:
    from metku.framefem import Element
except:
    from framefem import Element


class EBBeam(Element):
    """ Euler-Bernoulli beam element

        Parameters:
        -----------
        :param n1: start node
        :param n2: end node
        :param section: element's section properties
        :param material: element's material properties

        :type n1: FEMNode
        :type n2: FEMNode
        :type section: Section
        :type material: Material

        Variables:
        ----------
        :ivar bending_moment: bending moment values at element's ends
        :ivar shear_force: shear force values at element's ends

        :vartype bending_moment: list
        :vartype shear_force: list
    """

    def __init__(self, n1, n2, section, material):
        Element.__init__(self, n1, n2, section, material)

        self.bending_moment = [0.0, 0.0]
        self.shear_force = [0.0, 0.0]


    def transformation_matrix(self):
        """ Calculates transformation matrix from local to global coordinates

            Returns:
            --------
            :return: element's transformation matrix
            :rtype: np.array

         """
        c = self.direction_cosines()
        L = np.zeros((6, 6))
        L[0, 0:2] = c[0:2]
        L[1, 0] = -L[0, 1]
        L[1, 1] = L[0, 0]
        L[2, 2] = 1.0
        L[3:6, 3:6] = L[0:3, 0:3]
        return L

    def global_dofs(self):
        """ Get numbering of element's global degrees of freedom

            Returns:
            --------
            :return: list of end nodes dofs added together
            :rtype: list
        """
        return self.nodes[0].dofs + self.nodes[1].dofs

    @lru_cache(CACHE_BOUND)
    def local_stiffness_matrix(self, E, A, I1, Le):
        """ Generates stiffness matrix in local coordinates


            Returns:
            --------
            :return: element's stiffness matrix in local coordinates
            :rtype: np.array

        """
        rodc = E * A / Le
        EI1 = E * I1
        k0 = np.zeros((6, 6))
        #
        k0[[0, 3], [0, 3]] = rodc
        k0[[0, 3], [3, 0]] = -rodc
        #
        bend1 = 12 * EI1 / Le ** 3
        k0[[1, 4], [1, 4]] = bend1
        k0[[1, 4], [4, 1]] = -k0[1, 1]
        #
        bend2 = 6 * EI1 / Le ** 2
        k0[[1, 2, 1, 5], [2, 1, 5, 1]] = bend2
        #
        k0[[2, 4, 4, 5], [4, 2, 5, 4]] = -bend2
        #
        k0[[2, 5], [5, 2]] = 2 * EI1 / Le
        #
        k0[[2, 5], [2, 5]] = 4 * EI1 / Le

        return k0

    def stiffness_matrix(self):
        """ Computes element's global stiffness matrix

            Returns:
            ---------
            :return: element's stiffness matrix in global coordinates
            :rtype: np.array

        """
        E = self.material.young
        A = self.section.A
        I1 = self.section.I[0]
        Le = self.length()
        k0 = self.local_stiffness_matrix(E, A, I1, Le)

        # k0 = CheckReleases[fem.elem[N],k0];

        # local-global matrix
        L = self.transformation_matrix()

        # globaali elementin matriisi
        ke = L.transpose().dot(k0.dot(L))

        return ke

    def local_geometric_stiffness_matrix(self):
        """ Geometric stiffness matrix in local coordinates
            source: Cook et. al 1989, Section 14.2

            Returns:
            ---------
            :return: element's geometric stiffness matrix in local coordinates
            :rtype: np.array

        """
        P = self.axial_force[1]
        Le = self.length()

        g0 = P / 30 / Le
        keg0 = g0 * np.array([[36, 3 * Le, -36, 3 * Le],
                              [3 * Le, 4 * Le ** 2, -3 * Le, -Le ** 2],
                              [-36, -3 * Le, 36, -3 * Le],
                              [3 * Le, -Le ** 2, -3 * Le, 4 * Le ** 2]])
        kg0 = np.zeros([6, 6])
        q = [1, 2, 4, 5]
        kg0[np.ix_(q, q)] = keg0

        return kg0

    def geometric_stiffness_matrix(self):
        """ Geometric stiffness matrix in global coordinates \n
            From: Cook et. al (1989), Section 14.2

            Returns:
            ---------
            :return: element's geometric stiffness matrix in global coordinates
            :rtype: np.array

        """

        # local-global transformation matrix
        L = self.transformation_matrix()

        kG0 = self.local_geometric_stiffness_matrix()

        kG = L.transpose().dot(kG0.dot(L))

        return kG

    def equivalent_nodal_loads(self, q):
        """ Equivalent nodal loads for load in vector q

            Returns:
            ---------
            :return: Equivalent nodal loads for load in vector q
            :rtype: np.array
        """

        floc = np.zeros(6)

        T = self.transformation_matrix()
        L = self.length()

        # Load vector transformed into element local coordinate system
        qloc = T[:3, :3].dot(q)

        # Construct nodal load in local coordinates
        # qloc[0,1,2] = axial, shear, moment of node 1
        # qloc[3,4,5] = axial, shear, moment ofnode 2
        # Axial force
        floc[[0, 3]] = 0.5 * qloc[0] * L
        # Shear force
        floc[[1, 4]] = 0.5 * qloc[1] * L
        # Moment
        floc[2] = qloc[1] * L ** 2 / 12.0
        floc[5] = -floc[2]

        # Store local loads
        self.floc = floc

        # print(self.floc)

        fglob = T.transpose().dot(floc)
        return fglob

    def shape_fun(self, x):
        """ Evaluate Hermitean shape functions at local coordinate x

            Returns:
            ---------
            :return: Hermitean shape functions at local coordinate x
            :rtype: np.array


        """
        L = self.length()
        N = np.zeros(4)
        N[0] = 1 - 3 * x ** 2 / L ** 2 + 2 * x ** 3 / L ** 3
        N[1] = x - 2 * x ** 2 / L + x ** 3 / L ** 2
        N[2] = 3 * x ** 2 / L ** 2 - 2 * x ** 3 / L ** 3
        N[3] = -x ** 2 / L ** 2 + x ** 3 / L ** 2

        return N

    def internal_forces(self):
        """ Calculate internal forces
            NOTE: these internal forces do not take
            loads along the element into account!

            Works only for a single load case!
        """

        """ Get nodal displacements in local coordinates
            and multiply them with local element stiffness matrix
            to get internal forces in member's local coordinate system.
        """
        q = self.local_displacements()
        # ke = self.k0
        E = self.material.young
        A = self.section.A
        I1 = self.section.I[0]
        Le = self.length()
        ke = self.local_stiffness_matrix(E, A, I1, Le)

        if self.floc.size > 0:
            R = ke.dot(q) - self.floc
        else:
            R = ke.dot(q)

        """ Any load on the element not acting on a node must be
            taken into account here. This requires the following steps:
            1. Identify loads acting on element
            2. Compute their equivalent nodal loads
            3. Modify internal forces

            Probably an easy solution is to include an
            Attribute floc for equivalent nodal loads. The nodal loads are
            saved there, when they are first calculated for analysis
        """

        self.axial_force[:2] = [-R[0], R[3]]
        self.bending_moment[:2] = [-R[2], R[5]]
        self.shear_force[:2] = [R[1], -R[4]]
        
class EBBeam3D(Element):
    """ Euler-Bernoulli beam element in 3D

        Parameters:
        -----------
        :param n1: start node
        :param n2: end node
        :param nref: reference node for orienting the element
        :param section: element's section properties
        :param material: element's material properties

        :type n1: FEMNode
        :type n2: FEMNode
        :type section: Section
        :type material: Material

        Variables:
        ----------
        :ivar bending_moment: bending moment values at element's ends
        :ivar shear_force: shear force values at element's ends

        :vartype bending_moment: list
        :vartype shear_force: list
        
        TODO:
            equivalent_nodal_loads
            internal_forces
            test all functions
    """
    
    def __init__(self, n1, n2, nref, section, material):
        Element.__init__(self, n1, n2, section, material)

        # Bending moment with respect to z axis
        self.bending_moment_z = [0.0, 0.0]
        # Shear force with respect to z-axis
        self.shear_force_z = [0.0, 0.0]
        # Torque moment
        self.torque = 0.0
        
        """ Reference node for orienting the member
            This node must not be on the line passing through nodes n1 and n2.
        """
        self.nref = nref
        
    def direction_cosines(self):
        """ Calculates direction cosines for a 3D-beam element """
        
        l0 = np.zeros((3,3))
        
        X = self.coord()
        L = self.length()
        dX = X[1, :] - X[0, :]        
        # direction cosines along the local x axis
        l0[0,:] = dX / L
        
        # calculate direcion cosines along the local z axis
        # coordinates of the reference node
        X3 = self.nref.coord
        dX13 = X3-X[0,:]
        l13 = np.linalg.norm(dX13)
        
        # unit vector between point 1 and reference point
        V13 = dX13/l13
        VZ0 = np.cross(l0[0,:],V13)
        l0[2,:] = VZ0/np.linalg.norm(VZ0)
        
        l0[1,:] = np.cross(l0[2,:],l0[1,:])
        
        return l0
        
    def transformation_matrix(self):
        """ Calculates transformation matrix from local to global coordinates

            Returns:
            --------
            :return: element's transformation matrix
            :rtype: np.array

         """
        l0 = self.direction_cosines()
        L = np.zeros((12, 12))
        L[0:3, 0:3] = l0
        L[3:6, 3:6] = l0
        L[6:9, 6:9] = l0
        L[9:12, 9:12] = l0
        return L

    def global_dofs(self):
        """ Get numbering of element's global degrees of freedom

            Returns:
            --------
            :return: list of end nodes dofs added together
            :rtype: list
        """
        return self.nodes[0].dofs + self.nodes[1].dofs
    
    
    def local_stiffness_matrix(self):
        """ Generates stiffness matrix in local coordinates


            Returns:
            --------
            :return: element's stiffness matrix in local coordinates
            :rtype: np.array

        """
        E = self.material.young
        G = self.material.shear_modulus
        A = self.section.area
        I1 = self.section.Iy
        I2 = self.section.Iz
        J = self.section.J
        
        Le = self.length()

        rodc = E * A / Le
        EI1 = E * I1
        EI2 = E * I2
        torc = G*J/Le
        k0 = np.zeros((12, 12))
        # Axial stiffness
        k0[[0, 5], [0, 5]] = rodc
        k0[[0, 5], [5, 0]] = -rodc
        # Torsional stiffness
        k0[[3, 9], [3, 9]] = torc
        k0[[3, 9], [9, 3]] = -torc
        # Bending with respect to y axis
        bend1y = 12 * EI1 / Le ** 3
        k0[[1, 7], [1, 7]] = bend1y
        k0[[1, 7], [7, 1]] = -bend1y
        #
        bend2y = 6 * EI1 / Le ** 2
        k0[[1, 1, 5, 11], [5, 11, 1, 1]] = bend2y
        #
        k0[[5, 7, 7, 11], [7, 5, 11, 7]] = -bend2y
        #
        k0[[5, 11], [11, 5]] = 2 * EI1 / Le
        #
        k0[[5, 11], [5, 11]] = 4 * EI1 / Le
        
        bend1z = 12 * EI2 / Le**3
        k0[[2, 8], [2, 8]] = bend1z
        k0[[8, 2], [2, 8]] = -bend1z

        bend2z = 6 * EI2 / Le ** 2
        k0[[8, 4, 8, 10], [4, 8, 10, 8]] = bend2z
        #
        k0[[2, 10, 2, 4], [10, 2, 4, 2]] = -bend2z

        k0[[4, 10], [10, 4]] = 2 * EI2 / Le
        #
        k0[[4, 10], [4, 10]] = 4 * EI2 / Le

        return k0
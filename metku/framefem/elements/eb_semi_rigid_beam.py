# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:16:56 2018

Euler-Bernoulli beam in 2D with rotational stiffness at nodes

@author: kmela
"""

import numpy as np
from functools import lru_cache

CACHE_BOUND = 2**10
from .ebbeam import EBBeam

# Constants
kRigid = 1e20  # articifial large stiffness for rigid joints
kHinge = 1e-4  # artificial small stiffness for hinges


class EBSemiRigidBeam(EBBeam):
    """ Euler-Bernoulli beam element with possibility to include
        semi-rigid joints. Based on Jalkanen (2004):
        JOUSTAVASTI TUETTU TASOPALKKIELEMENTTI, Journal of Structural Mechanics


        Parameters:
        -----------

        :param n1: element's first node
        :param n2: element's second node
        :param section: element's section
        :param material: element's material
        :param rot_stiff: element's ends rotational stiffness [n1 , n2] [kNm/rad]

        :type n1: FEMNode
        :type n2:  FEMNode
        :type section: Section
        :type material: Material
        :type rot_stiff: list

        Variables:
        ----------
        
        :ivar rot_stiff:
            
        Check EBBeam for more info

    """

    def __init__(self, n1, n2, section, material, rot_stiff=[np.inf, np.inf]):
        """

        """
        EBBeam.__init__(self, n1, n2, section, material)

        self.bending_moment = [0.0, 0.0]
        self.shear_force = [0.0, 0.0]
        for i in range(2):
            if rot_stiff[i] <= kHinge:
                rot_stiff[i] = kHinge
            elif rot_stiff[i] >= kRigid:
                rot_stiff[i] = kRigid

        # print(rot_stiff)
        self.rot_stiff = rot_stiff

    def transformation_matrix(self):
        """ From local to global coordinates """
        c = self.direction_cosines()
        L = np.zeros((6, 6))
        L[0, 0:2] = c[0:2]
        L[1, 0] = -L[0, 1]
        L[1, 1] = L[0, 0]
        L[2, 2] = 1.0
        L[3:6, 3:6] = L[0:3, 0:3]
        return L

    def gamma(self):
        """ Eq. (12) of Jalkanen (2004) """

        gamma = []

        E = self.material.young
        I1 = self.section.I[0]
        L = self.length()
        g0 = E * I1 / L

        for k in self.rot_stiff:
            """ k = 0 means a hinged connection """
            if k == 0.0:
                gamma.append(np.inf)
            elif k == kRigid:
                gamma.append(0.0)
            else:
                gamma.append(g0 / k)
        return gamma

    def delta(self):
        """ Eq. (15) of Jalkanen (2004) """
        gamma = self.gamma()
        delta = (1 + 4 * gamma[0]) * (1 + 4 * gamma[1]) - 4 * gamma[0] * gamma[1]
        return delta

    def matrix_S(self):
        """ Eq. (16) of Jalkanen (2004) """
        S = np.zeros((4, 4))
        L = self.length()
        gamma = self.gamma()

        S[1, 0:4] = [-6 * gamma[0] / L * (1 + 2 * gamma[1]), -4 * gamma[0] * (1 + 3 * gamma[1]),
                     6 * gamma[0] / L * (1 + 2 * gamma[1]), -2 * gamma[0]]
        S[3, 0:4] = [-6 * gamma[1] / L * (1 + 2 * gamma[0]), -2 * gamma[1], 6 * gamma[1] / L * (1 + 2 * gamma[0]),
                     -4 * gamma[1] * (1 + 3 * gamma[0])]

        return S

    def global_dofs(self):
        """ Get numbering of element global degrees of freedom """
        return self.nodes[0].dofs + self.nodes[1].dofs

    @lru_cache(CACHE_BOUND)
    def local_stiffness_matrix(self, E, A, I1, Le):
        """ Stiffness matrix in local coordinates """

        rodc = E * A / Le
        EI1 = E * I1
        bend1 = 12 * EI1 / Le ** 3
        bend2 = 6 * EI1 / Le ** 2
        bend3 = 4 * EI1 / Le
        bend4 = 2 * EI1 / Le

        # k00 = np.zeros((4, 4))
        k0 = np.zeros((6, 6))

        k00 = np.array([[bend1, bend2, -bend1, bend2],
                        [bend2, bend3, -bend2, bend4],
                        [-bend1, -bend2, bend1, -bend2],
                        [bend2, bend4, -bend2, bend3]])


        delta = self.delta()
        S = self.matrix_S()

        k1 = 1 / delta * k00.dot(S) + \
             1 / delta * S.transpose().dot(k00) + \
             1 / delta ** 2 * S.transpose().dot(k00).dot(S)

        # s2 = S[1,:]
        # s4 = S[3,:]
        k2 = np.zeros((4, 4))
        for i in range(0, 2):
            if self.rot_stiff[i] < kRigid:
                s = S[2 * i + 1, :]
                # print(s)
                k2 += self.rot_stiff[i] / delta ** 2 * np.outer(s, s)

        # k2 = self.rot_stiff[0]/delta**2*np.outer(s2,s2) + self.rot_stiff[1]/delta**2*np.outer(s4,s4)


        # print(k00)
        # print(k1)
        # print(k2)

        k0[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = k00 + k1 + k2
        k0[[0, 3], [0, 3]] = rodc
        k0[[0, 3], [3, 0]] = -rodc

        return k0

    def stiffness_matrix(self):
        """ Compute the stiffness matrix """

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
            From: Cook et. al 1989, Section 14.2
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
        """ Geometric stiffness matrix in global coordinates
            From: Cook et. al (1989), Section 14.2
        """

        # local-global transformation matrix
        L = self.transformation_matrix()

        kG0 = self.local_geometric_stiffness_matrix()

        kG = L.transpose().dot(kG0.dot(L))

        return kG

    def equivalent_nodal_loads(self, q):
        """ Equivalent nodal loads for load in vector q
        """

        floc0 = np.zeros(6)
        floc = np.zeros(6)

        T = self.transformation_matrix()
        L = self.length()

        # Load vector transformed into element local coordinate system
        # qloc[0] = qx, qloc[1] = qy
        qloc = T[:3, :3].dot(q)

        delta = self.delta()
        S = self.matrix_S()

        # Construct nodal load in local coordinates
        floc0[[0, 3]] = 0.5 * qloc[0] * L
        floc0[[1, 4]] = 0.5 * qloc[1] * L
        floc0[2] = qloc[1] * L ** 2 / 12.0
        floc0[5] = -floc0[2]

        floc[[0, 3]] = 0.5 * qloc[0] * L  # axial nodal loads
        floc[1] = floc0[1] + S[1, 0] / delta * floc0[2] + S[3, 0] / delta * floc0[5]
        floc[2] = (1 + S[1, 1] / delta) * floc0[2] + S[3, 1] / delta * floc0[5]
        floc[4] = S[1, 2] / delta * floc0[2] + floc0[4] + S[3, 2] / delta * floc0[5]
        floc[5] = S[1, 3] / delta * floc0[2] + (1 + S[3, 3] / delta) * floc0[5]

        self.floc = floc

        # print(self.floc)

        fglob = T.transpose().dot(floc)
        return fglob


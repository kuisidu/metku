# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
class PointLoad(Load):
    """ Class for point loads (forces and moments)

        Parameters:
        ------------
        :param sid:  load id
        :param node: node subjected to the load
        :param v: load vector [Fx, Fy, Mz] [kN, kN, kNm]
        :param f: scaling factor

        :type sid: int
        :type node: FEMNode
        :type v: list
        :type f: float

        Variables:
        --------------
        :ivar sid:  load id
        :ivar node: node subjected to the load
        :ivar v: load vector [Fx, Fy, Mz] [kN, kN, kNm]
        :ivar f: scaling factor

        :vartype sid: int
        :vartype node: FEMNode
        :vartype v: list
        :vartype f: float

    """

    def __init__(self, sid, node, v, f):
        Load.__init__(self, sid)

        self.node = node
        """ node subjected to the load"""
        self.v = v
        """ load vector"""
        self.f = f
        """ scaling factor"""

    def load_and_dofs(self):
        """ Returns the loads to be inserted to the global load vector
            and the corresponding global degrees of freedom

            Returns:
            ---------

            :return: load to be inserted to the global load vector and
                    the corresponding global degrees of freedom
            :rtype: (int, int)

        """

        """ Load vector """
        F = self.f * np.array(self.v)

        """ Degrees of freedom """
        dofs = np.array(self.node.dofs)

        """ Find free degrees of freedom (those with numbering
                                          greater than -1)
        """
        nzero_dofs = dofs >= 0

        """ Get free degrees of freedom and the corresponding part
            of the load vector
        """
        d = dofs[nzero_dofs]
        Fnz = F[nzero_dofs]

        return Fnz, d

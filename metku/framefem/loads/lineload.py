# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
import numpy as np

class LineLoad(Load):
    """ Class for 2D line loads (forces and moments)
    """

    def __init__(self, sid, eid, xloc, qval, direction, coords=1):
        """ Input:
            sid -- load id
            eid -- element subjected to the load
            xloc -- starting and ending locations of the load (local coordinates)
                    xloc[0] .. staring coordinate 0 <= xloc[0] <= 1
                    xloc[1] .. ending coordinate 0 <= xloc[0] < xloc[1] <= 1
            qval -- value of the load at the coordinates xloc
            direction -- 0 .. x-axis
                         1 .. y-axis
            coords .. 0 .. local
                      1 .. global
        """
        Load.__init__(self, sid)

        self.elem = eid
        self.xloc = xloc
        self.qval = qval
        self.dir = direction
        self.coords = coords

    def load_and_dofs(self):
        """ Returns the loads to be inserted to the global load vector
            and the corresponding global degrees of freedom
        """

        """ Get load vector """
        if self.dir == "x":
            q = np.array([self.qval[0], 0, 0])
        else:
            q = np.array([0, self.qval[0], 0])

        """ Get equivalent nodal loads """
        F = self.elem.equivalent_nodal_loads(q)

        """ Number of degrees of freedom """
        dofs = np.array(self.elem.global_dofs())

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

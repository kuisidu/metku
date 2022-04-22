# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela


class PointLoad:

    def __init__(self, coord, F, load_id=0, scale=1):

        self.load_id = load_id
        self.node = None
        self._coord = coord
        self.F = F
        self.scale = scale

    @property
    def Fx(self):
        return self.F[0]

    @Fx.setter
    def Fx(self, val):
        self.F[0] = val

    @property
    def Fy(self):
        return self.F[1]

    @Fy.setter
    def Fy(self, val):
        self.F[1] = val

    @property
    def Fz(self):
        return self.F[2]

    @Fz.setter
    def Fz(self, val):
        self.F[2] = val

    @property
    def coord(self):
        return self.node.coord

    @coord.setter
    def coord(self, val):
        # 2D
        if len(val) == 2:
            x, y = val
            self.node.x = x
            self.node.y = y
        # 3D
        elif len(val) == 3:
            x, y, z = val
            self.node.x = x
            self.node.y = y
            self.node.z = z


if __name__ == '__main__':
    pl = PointLoad([0, 0], [0, 100, 0])
    print(pl)
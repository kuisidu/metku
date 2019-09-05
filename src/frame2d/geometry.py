import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import *

class GeometryMember(Segment):

    def __new__(cls, p1, p2, **kwargs):
        dim = len(p1)
        if dim == 2:
            return GeometryMember2D(p1, p2, **kwargs)
        elif dim == 3:
            return GeometryMember3D(p1, p2, **kwargs)

    def __init__(self, start_coord, end_coord):
        super().__init__()
        self.c0 = start_coord
        self.c1 = end_coord

    @property
    def start_coord(self):
        return np.asarray(self.c0)

    @start_coord.setter
    def start_coord(self, val):
        self.c0 = val


    @property
    def end_coord(self):
        return np.asarray(self.c1)

    @end_coord.setter
    def end_coord(self, val):
        self.c1 = val

    @property
    def x0(self):
        return self.start_coord[0]

    @property
    def y0(self):
        return self.start_coord[1]

    @property
    def x1(self):
        return self.end_coord[0]

    @property
    def y1(self):
        return self.end_coord[1]

    @property
    def coordinates(self):
        return np.asarray([self.start_coord, self.end_coord])



    @property
    def angle(self):
        """
        Returns member's angle in radians
        """
        if (self.x1 - self.x0) == 0:
            angle = math.radians(90)
        else:
            angle = np.arctan((self.y1 - self.y0) / (self.x1 - self.x0))
        return angle


    def to_global(self, loc_coord):
        """
        Returns local coordinate as global coordinates

        Parameters:
        ------------
        :param loc_coord: local coordinate
        :return: global coordinate

        :type coord: float
        :rtype: np.array
        """

        dist = loc_coord * self.length.evalf()
        coord = self.start_coord + dist * np.asarray(self.direction.unit)
        return coord

    def intersection(self, mem):
        """
        Calculates coordinate where two members intersect

        Parameters:
        -----------
        :param mem: member
        :return: coordinates for intersection

        source: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
        """
        point = super().intersection(mem)
        if len(point):
            return np.asarray(point)[0]
        return point

class GeometryMember2D(GeometryMember, Segment2D):

    def __new__(cls, p1, p2, **kwargs):
        return Segment2D.__new__(cls, p1, p2, **kwargs)

    def __init__(self, start_coord, end_coord):
        super().__init__(start_coord, end_coord)

    def plot(self, show=True):
        """
        Draws member
        """
        plt.plot([self.x0, self.x1], [self.y0, self.y1])
        if show:
            plt.show()

    @property
    def angle(self):
        """
        Returns member's angle in radians
        """
        if (self.x1 - self.x0) == 0:
            angle = math.radians(90)
        else:
            angle = np.arctan((self.y1 - self.y0) / (self.x1 - self.x0))
        return angle

    @property
    def v(self):
        """
        Unit vector
        :return: unit vector
        """
        v = np.array([np.cos(self.angle), np.sin(self.angle)])
        return v

    @property
    def u(self):
        """
        Vector perpendicular to member
        :return:
        """
        u = np.array([-np.sin(self.angle), np.cos(self.angle)])
        return u


class GeometryMember3D(GeometryMember, Segment3D):

    def __new__(cls, p1, p2, **kwargs):
        return Segment3D.__new__(cls, p1, p2, **kwargs)

    def __init__(self, start_coord, end_coord):
        super().__init__(start_coord, end_coord)


    @property
    def z0(self):
        return self.start_coord[2]

    @property
    def z1(self):
        return self.end_coord[2]


    @property
    def v(self):
        """
        Unit vector
        :return: unit vector
        """
        v = (self.end_coord - self.start_coord) / self.length

        return v

    @property
    def u(self):
        """
        Vector perpendicular to member
        :return:
        """
        u = np.array([-np.sin(self.angle), np.cos(self.angle)])
        return u


    def plot(self, show=True):
        """
        Plots member
        :param show:
        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot([self.x0, self.x1],
                [self.y0, self.y1],
                [self.z0, self.z1])
        if show:
            plt.show()


if __name__ == '__main__':

    mems = []
    tc = GeometryMember([0, 100], [1000, 500])
    bc = GeometryMember([0, 0], [1000, 0])
    mems.extend([tc, bc])
    b = 0
    t = 0
    for i in range(10):
        if i % 2:
            b = i / 9
        else:
            t = i / 9
        w = GeometryMember(bc.to_global(b), tc.to_global(t))
        mems.append(w)
    w = GeometryMember(bc.to_global(1), tc.to_global(1))
    mems.append(w)

    for mem in mems:
        mem.plot(False)
    plt.show()


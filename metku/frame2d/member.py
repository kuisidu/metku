# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
from metku.frame2d.frame2d import *
from metku.frame2d.geometry import *
from metku.frame2d.materials import MATERIALS
from metku.framefem.elements import *
from metku.framefem.framefem import *
from metku.sections.steel import *


class FrameMember(GeometryMember):
    """
    CLASS DEFINITION
    """
    num_members = 0

    def __new__(cls, p1, p2, **kwargs):
        """
        Checks whether the dimension of the coordinates is 2 or 3
        :param p1: start coordinate
        :param p2: end coordinate
        :param kwargs: optional keyword arguments
        :return: 2D or 3D FrameMember instance
        """
        dim = len(p1)
        if dim == 2:
            return FrameMember2D(p1, p2, **kwargs)
        elif dim == 3:
            return FrameMember3D(p1, p2, **kwargs)

    def __init__(self, start_coord, end_coord, **kwargs):
        super().__init__(start_coord, end_coord)

        # Default values
        values = dict(
            profile='IPE 100',
            num_elements=4,
            Sj1=np.inf,
            Sj2=np.inf,
            material='S355'
        )
        # Update values with given kwargs
        values.update(kwargs)

        # Assign unique id
        self.mem_id = self.num_members
        FrameMember.num_members += 1

        # Initialize nodes and elements
        self.nodes = {}
        self.elements = {}
        self.nels = values['num_elements']

        # Initialize cross-sectional properties
        self.cross_section = None
        self.steel_member = None  # TODO: CHANGE THIS NAME TO A GENERAL ONE
        self.material = MATERIALS[values['material']]
        self.profile = values['profile']

        # Initialize rotational stiffness
        self.Sj1 = values['Sj1']
        self.Sj2 = values['Sj2']

        # List of intersecting objects and nodes
        self.intersecting_objs = []
        self.intersecting_nodes = {}

        # Create nodes and elements
        self.generate_nodes()
        self.generate_elements()

    @property
    def num_elements(self):
        return self.nels

    @num_elements.setter
    def num_elements(self, val):
        """
        Change number of elements
        TODO!!
        """
        self.update()

    @property
    def nodal_coordinates(self):
        """
        Returns nodal coordinates as a list
        :return: nodal coordinate list
        """
        coords = [list(node.coord) for node in self.nodes.values()]
        intersect = [list(node.coord) for node in self.intersecting_nodes.values()]
        if len(intersect):
            for coord in intersect:
                coords.append(coord)
        return coords

    @property
    def fy(self):
        return self.material.fy

    @property
    def profile(self):
        return str(self.cross_section)

    @profile.setter
    def profile(self, val):
        """
        Changes member's profile to given value.
        Calculates new cross-sectional properties and sets these new values
        to member's elements.

        :param val: string, profile name e.g. 'IPE 100'
        """
        if isinstance(val, list) and len(val) == 5:
            h, b, tf, tw, r = val
            if isinstance(self.cross_section, CustomISection):
                self.cross_section.H = h
                self.cross_section.B = b
                self.cross_section.tf = tf
                self.cross_section.tw = tw
                self.cross_section.r = r
                self.cross_section.cross_section_properties()
            else:
                self.cross_section = CustomISection(h, b, tf, tw, r)
        else:
            val = val.upper()
            splitted_val = val.split(" ")
            profile_type = splitted_val[0]
            if profile_type == 'IPE' or profile_type == 'HE':
                try:
                    H = int(splitted_val[1])
                    catalogue = True
                except ValueError:
                    H = float(splitted_val[1])
                    catalogue = False

            elif profile_type == 'WI':
                vals = splitted_val[1].replace('X', '-').split('-')
                H = float(vals[0])
                b1 = float(vals[3])
                b2 = float(vals[5])
                tw = float(vals[1])
                tf1 = float(vals[2])
                tf2 = float(vals[4])

            else:
                vals = splitted_val[1].split('X')
                if len(vals) == 2:
                    H = float(vals[0])
                    T = float(vals[1])
                elif len(vals) == 3:
                    H = float(vals[0])
                    B = float(vals[1])
                    T = float(vals[2])
                else:
                    raise TypeError(f'{splitted_val[1]} is not valid profile')

            if profile_type == 'IPE':
                self.cross_section = IPE(H, self.fy)

            elif profile_type == 'WI':
                self.cross_section = WISection(H, tw, [b1, b2], [tf1, tf2],
                                               self.fy)

            elif profile_type == 'HE':
                if splitted_val[2] == 'A':
                    self.cross_section = HEA(H, self.fy)
                elif splitted_val[2] == 'AA':
                    self.cross_section = HEAA(H, self.fy)
                elif splitted_val[2] == 'B':
                    self.cross_section = HEB(H, self.fy)
                elif splitted_val[2] == 'C':
                    self.cross_section = HEC(H, self.fy)
                elif splitted_val[2] == 'M':
                    self.cross_section = HEM(H, self.fy)

            elif profile_type == 'CHS':
                self.cross_section = CHS(H, T, self.fy)

            elif profile_type == 'RHS':
                self.cross_section = RHS(H, B, T, self.fy)

            elif profile_type == 'SHS':
                self.cross_section = SHS(H, T, self.fy)
            else:
                raise ValueError(
                    '{} is not valid profile type!'.format(profile_type))

        # Change steel_member objects properties
        if self.steel_member:
            self.steel_member.profile = self.cross_section

    @GeometryMember.start_coord.setter
    def start_coord(self, val):
        self.c0 = val
        self.update()


    @GeometryMember.end_coord.setter
    def end_coord(self, val):
        self.c1 = val
        self.update()

    def update(self):
        """
        Updates nodal coordinates
        """
        self.update_coordinates()
        self.update_elements()

    @property
    def locs(self):
        """
        List of nodes' local coordinates
        """
        if len(self.nodes):
            locs = [self.p1.distance(node.coord) / self.length for node in self.nodes.values()]
            intersect_locs = [self.p1.distance(node.coord) / self.length for node in self.intersecting_nodes.values()]
            if len(intersect_locs):
                locs.extend(intersect_locs)
        else:
            locs = list(np.arange(0, 1, 1 / self.num_elements))
            locs.append(1)

        return sorted(set(locs))

    def generate_elements(self):
        """
        Generates elements
        """
        for i in range(len(self.locs) -1):
            n0 = self.nodes[i]
            n1 = self.nodes[i + 1]
            elem = EBSemiRigidBeam(n0, n1, self.cross_section, self.material)
            self.elements[i] = elem

    def update_elements(self):
        """
        Updates elements
        """
        elements = self.elements.values()
        nodes = self.nodes
        nodes.update(self.intersecting_nodes)
        nodes = list(nodes.values())
        nodes.sort(key=lambda n: self.p1.distance(n.coord))
        for i in range(len(self.locs) -1):
            try:
                elem = elements[i]
                elem.n1 = nodes[i]
                elem.n2 = nodes[i + 1]
            except:
                n0 = nodes[i]
                n1 = nodes[i + 1]
                elem = EBSemiRigidBeam(n0, n1, self.cross_section,
                                       self.material)
                self.elements[i] = elem


    def generate_nodes(self):
        """
        Generates nodes
        """
        for i, loc in enumerate(self.locs):
            x, y = self.to_global(loc)
            nid = self.mem_id * 100 + i
            node = FEMNode(nid, x, y)
            self.nodes[i] = node

    def plot_elements(self, show=True):
        """
        Plots elements and nodes
        """
        X = []
        Y = []
        for node in self.nodes.values():
            X.append(node.x)
            Y.append(node.y)
            plt.scatter(node.x, node.y)
            plt.text(node.x, node.y, node.nid)
        plt.plot(X, Y)
        if show:
            plt.show()


class FrameMember2D(FrameMember, GeometryMember2D):

    def __new__(cls, p1, p2, **kwargs):
        return GeometryMember2D.__new__(cls, p1, p2, **kwargs)

    def __init__(self, start_coord, end_coord, **kwargs):
        super().__init__(start_coord, end_coord, **kwargs)

    def update_coordinates(self):
        """
        Updates nodes' coordinates
        """
        all_nodes = self.nodes
        all_nodes.update(self.intersecting_nodes)
        all_nodes = list(all_nodes.values())
        all_nodes.sort(key=lambda n: self.p1.distance(n.coord))
        for i, loc in enumerate(self.locs):
            x, y = self.to_global(loc)
            node = all_nodes[i]
            node.x = x
            node.y = y



class FrameMember3D(FrameMember, GeometryMember3D):

    def __new__(cls, p1, p2, **kwargs):
        return GeometryMember3D.__new__(cls, p1, p2, **kwargs)

    def __init__(self, start_coord, end_coord, **kwargs):
        super().__init__(start_coord, end_coord, **kwargs)

    def generate_nodes(self):
        """
        Generates nodes
        """
        locs = list(np.arange(0, 1, 1 / self.num_elements))
        locs.append(1)

        for i, loc in enumerate(locs):
            x, y, z = self.to_global(loc)
            nid = self.mem_id * 100 + i
            node = FEMNode(nid, x, y, z)
            self.nodes[i] = node

    def update_coordinates(self):
        """
        Updates nodes' coordinates
        """
        for i, loc in enumerate(self.locs):
            x, y, z = self.to_global(loc)
            node = self.nodes[i]
            node.x = x
            node.y = y
            node.z = z

    def plot_elements(self, show=True):
        """
        Plots elements and nodes
        """
        X = []
        Y = []
        Z = []

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        for elem in self.elements.values():
            n1, n2 = elem.nodes
            X.append(n1.x)
            Y.append(n1.y)
            Z.append(n1.z)
            ax.scatter(n1.x, n1.y, n1.z)
            ax.text(n1.x, n1.y, n1.z, n1.nid)
        X.append(n2.x)
        Y.append(n2.y)
        Z.append(n2.z)
        ax.scatter(n2.x, n2.y, n2.z)
        ax.text(n2.x, n2.y, n2.z, n2.nid)

        ax.plot(X, Y, Z)
        if show:
            plt.show()


if __name__ == '__main__':
    mem = FrameMember([0, 0], [0, 1000], num_elements=5, number_elements=3)
    print(mem.plot())

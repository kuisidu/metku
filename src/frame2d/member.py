from src.frame2d.frame2d import *
from src.frame2d.geometry import *
from src.framefem.elements import *
from src.framefem.framefem import *
from src.frame2d.materials import MATERIALS
from src.sections.steel import *
from mpl_toolkits.mplot3d import Axes3D


class FrameMember(GeometryMember):
    """
    CLASS DEFINITION
    """
    num_members = 0

    def __new__(cls, p1, p2, **kwargs):
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
        self.nels= values['num_elements']

        # Initialize cross-sectional properties
        self.cross_section = None
        self.steel_member = None # TODO: CHANGE THIS NAME TO A GENERAL ONE
        self.material = MATERIALS[values['material']]
        self.profile = values['profile']

        # Initialize rotational stiffness
        self.Sj1 = values['Sj1']
        self.Sj2 = values['Sj2']


    @property
    def num_elements(self):
        return self.nels

    @num_elements.setter
    def num_elements(self, val):
        """
        Change number of elements
        TODO!!
        """
        self.nodes = {}
        self.elements = {}
        self.nels = val

        self.generate_nodes()
        self.generate_elements()

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
        self.update_coordinates()

    @GeometryMember.end_coord.setter
    def end_coord(self, val):
        self.c1 = val
        self.update_coordinates()


    def update_coordinates(self):
        """
        Updates nodes' coordinates
        """
        locs = list(np.arange(0, 1, 1 / self.num_elements))
        locs.append(1)

        for i, loc in enumerate(locs):
            x, y = self.to_global(loc)
            node = self.nodes[i]
            node.x = x
            node.y = y

    def generate_elements(self):
        """
        Generates elements
        """
        for i in range(self.num_elements):
            n0 = self.nodes[i]
            n1 = self.nodes[i + 1]
            elem = EBSemiRigidBeam(n0, n1, self.cross_section, self.material)
            self.elements[i] = elem

    def generate_nodes(self):
        """
        Generates nodes
        """
        locs = list(np.arange(0, 1, 1 / self.num_elements))
        locs.append(1)

        for i, loc in enumerate(locs):
            x, y = self.to_global(loc)
            nid = self.mem_id * 100 + i
            node = FEMNode(nid, x, y)
            self.nodes[i] = node

    def plot_elements(self):
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
        plt.show()


class FrameMember2D(FrameMember, GeometryMember2D):

    def __new__(cls, p1, p2, **kwargs):
        return GeometryMember2D.__new__(cls, p1, p2, **kwargs)

    def __init__(self, start_coord, end_coord, **kwargs):

        super().__init__(start_coord, end_coord, **kwargs)

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

    def plot_elements(self):
        """
        Plots elements and nodes
        """
        X = []
        Y = []
        Z = []

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        for node in self.nodes.values():
            X.append(node.x)
            Y.append(node.y)
            Z.append(node.z)
            ax.scatter(node.x, node.y, node.z)
            ax.text(node.x, node.y, node.z, node.nid)
        ax.plot(X, Y, Z)
        plt.show()

if __name__ == '__main__':
    mem = FrameMember([0, 1000], [1000, 0])
    print(mem.coordinates)
    mem.generate_nodes()
    mem.plot()

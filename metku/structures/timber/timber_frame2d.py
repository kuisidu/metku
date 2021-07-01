"""
Timber Frame 2D

@author Viktor Haimi
"""

from structures.timber.timber_member import TimberMember
from sections.timber.timber_section import TimberSection, GypsumPlasterboardF
import framefem.framefem as fem
import frame2d.frame2d as f2d
import numpy as np
import matplotlib.pyplot as plt

class TimberFrame2D():
    def __init__(self, R=0, supports=None):
        self.f = fem.FrameFEM()
        self.members = []
        self.nodal_coordinates =[]
        self.point_loads = {}
        self.line_loads = {}
        self.supports = {}
        self.support_nodes = []
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.r = []
        self.is_calculated = False
        self.is_generated = False
        self.truss = []
        self.joints = {}
        self.nodes = []
        self.R = R

        if supports:
            self.generate_supports(supports)

    @property
    def L(self):
        """ Width/span of the frame """
        x_coordinates = [mem.coordinates[0][0] for mem in
                         self.members]
        x_coordinates.extend(
            [mem.coordinates[1][0] for mem in self.members])
        x_min = min(x_coordinates)
        x_max = max(x_coordinates)
        L = x_max - x_min
        return L

    @property
    def H(self):
        """ Height of the frame """
        y_coordinates = [mem.coordinates[0][1] for mem in
                         self.members]
        y_coordinates.extend(
            [mem.coordinates[1][1] for mem in self.members])
        y_min = min(y_coordinates)
        y_max = max(y_coordinates)
        H = y_max - y_min
        return H

    def add(self, this):

        if isinstance(this, TimberMember):
            this.mem_id = len(self.members)
            self.members.append(this)
            this.R = self.R
            this.calc_nodal_coordinates()

            for mem in self.members:
                coord = this.line_intersection(this.coordinates)
                if isinstance(coord, list):
                    this.add_node_coord(coord)
                    mem.add_node_coord(coord)

            # TODO iffittely että ei lisätä samoja nodeja
            self.nodal_coordinates.extend(this.nodal_coordinates)
            self.nodal_coordinates.sort()

        # POINTLOADS
        elif isinstance(this, f2d.PointLoad):
            """ If a point load with same 'name' is already included
                in the frame, add a number to the end of the name
                of the new load such that each point load has a unique
                name.
    
                NOTE: using 'len' to provide new id number might fail,
                if point loads are deleted such that the number of point loads
                changes along the way.
            """
            if this.name in self.point_loads.keys():
                this.name += str(len(self.point_loads))
            self.point_loads[this.name] = this

            """ If the location of the point load is in between
                member end nodes, add a node to the corresponding member
            """
            for member in self.members:
                if member.point_intersection(this.coordinate):
                    member.add_node_coord(this.coordinate)

        # LINELOADS
        elif isinstance(this, f2d.LineLoad):

            if this.name in self.line_loads.keys():
                this.name += str(len(self.line_loads))
            self.line_loads[this.name] = this

        # SUPPORTS
        elif isinstance(this, f2d.Support):
            # this.supp_id = len(self.supports)
            supp_label = len(self.supports)
            self.supports[supp_label] = this
            # self.supports[this.supp_id] = this
            self.support_nodes.append(this.coordinate)

            """ If the support is located between end nodes of a member,
                add a node to that member
            """
            for member in self.members:
                #coord = member.point_intersection(this.coordinate)
                if member.point_intersection(this.coordinate):
                    member.add_node_coord(this.coordinate)

        # WRONG TYPE
        else:
            print(type(this), " is not supported.")
            raise TypeError

    def add_materials_and_sections(self):
        """ Adds members' materials and sections to the fem model
        """
        # Iterate through all frame's members and calculate cross-sections
        # properties and add those and member's material properties to
        # calculation model
        for member in self.members:
            member.add_material(self.f)
            member.add_section(self.f)

    def calc_nodal_forces(self):
        """ Calculates nodal forces and saves values to
            self.nodal_forces - dict
        """
        # Iterate through every frame's member's node and calculate forces
        # acting on that node
        for member in self.members:
            member.calc_nodal_forces()
            for node in member.nodal_forces:
                self.nodal_forces[node] = member.nodal_forces[node]

    def calc_nodal_displacements(self):
        """ Calculates nodal displacements and saves values to
            self.nodal_displacements -dict
        """
        # Iterate through every frame's member's node and calculate
        # displacements for that node
        for member in self.members:
            member.calc_nodal_displacements(self.f)
            for node in member.nodal_displacements.keys():
                self.nodal_displacements[node] = member.nodal_displacements[
                    node]

    def calculate(self, load_id=2, support_method='ZERO'):
        """ Calculates forces and displacements

            Parameters
            ----------
            :param load_id: Id of the loads to be added in the calculation model

            :type load_id: int / str
        """

        lcase_ids = [lc.load for lc in self.f.loadcases.values()]
        print(lcase_ids)
        if self.is_calculated == False:
            self.is_calculated = True
            """ Support ID is always 1! """
            self.f.nodal_dofs()

        # If load_id == 'ALL' calculates all load cases
        # calls recursively itself for each case
        if str(load_id).upper() == 'ALL':
            for lid in lcase_ids:
                self.calculate(load_id=lid,
                               support_method=support_method)
        else:
            self.f.linear_statics(support_method=support_method,
                                  lcase=load_id)
            self.calc_nodal_forces()
            self.calc_nodal_displacements()
            self.assign_forces()
            self.check_members_strength()
            # self.alpha_cr, _ = self.f.linear_buckling(k=4)

    def assign_forces(self):

        for member in self.members:
            member.assign_forces()

    def check_members_strength(self):
        """ Checks if members can bear their loads
        """
        self.r.clear()
        for member in self.members:
            member.check_cross_section()
            self.r.append(member.r)

    def generate_supports(self, supp_type):
        supp_type = supp_type.upper()
        if supp_type == 'FIXED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(f2d.FixedSupport(coord))

        elif supp_type == 'YHINGED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(f2d.YHingedSupport(coord))

        elif supp_type == 'XHINGED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(f2d.XHingedSupport(coord))

        elif supp_type == 'XYHINGED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(f2d.XYHingedSupport(coord))
        else:
            raise Exception('Väärän tyyppinen tuki')

    def generate(self):
        """ Generates the frame and truss FEM model
        """
        # If the frame hasn't been created, create members and calculate nodal
        # coordinates, otherwise initialize frame.
        if not self.is_generated:
            self.is_generated = True
        for member in self.members:
            member.calc_nodal_coordinates() # Generate FEM nodes for a member
            member.round_coordinates()
            member.generate(self.f)
            for coord in member.nodal_coordinates:
                if coord not in self.nodal_coordinates:
                    self.nodal_coordinates.append(coord)
            self.nodes.extend(list(member.nodes.values()))

        # Generate TrussJoints
        if len(self.truss):
            for truss in self.truss:
                for joint in truss.joints.values():
                    joint.generate(self.f)

        # Generate eccentricity elements
        for member in self.members:
            member.generate_eccentricity_elements(self.f)

        # Remove duplicate nodes
        self.nodes = list(set(self.nodes))

        # Add supports
        for support in self.supports.values():
            support.add_support(self.f)

        # Add point loads (if any)
        for pointLoad in self.point_loads.values():
            pointLoad.add_load(self.f)
            # Creates new loadcase if one with same load_id does not exist
            lcase_ids = [lc.load for lc in self.f.loadcases.values()]
            if pointLoad.load_id not in lcase_ids:
                self.f.add_loadcase(supp_id=1,
                                    load_id=pointLoad.load_id)

        # Add line loads (if any)
        for lineLoad in self.line_loads.values():
            member = lineLoad.member
            lineLoad.element_ids = member.lineload_elements(
                lineLoad.coordinates)
            lineLoad.add_load(self.f)
            # Creates new loadcase if one with same load_id does not exist
            lcase_ids = [lc.load for lc in self.f.loadcases.values()]
            if lineLoad.load_id not in lcase_ids:
                self.f.add_loadcase(supp_id=1,
                                    load_id=lineLoad.load_id)

    def plot(self, print_text=True, show=True,
             loads=True, color=False, axes=None):
        """ Plots the frame

            Parameters
            ----------
            :param print_text: Set true to print member's profiles and names (default: True)
            :param show: Set true to show the plot (default: True)
            :param loads: Set true to show loads (default: True)
            :param color: Set true to show members' utilization ratio (default: False)

            :type print_text : bool
            :type show: bool
            :type loads: bool
            :type color: bool

            Colors' meaning:

                blue -- member has load
                green -- member can bear its loads
                red -- member breaks under its loads
                black -- member is added, but not designed
        """
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes

        if self.is_calculated and color:
            color = True

        # Plot members
        for member in self.members:
            member.plot(print_text, color, ax)

        # Plot joints
        for joint in self.joints.values():
            joint.plot(color=color)

        # Plot supports
        for support in self.supports.values():
            node_coord = support.coordinate
            if support.dofs == [-1, -1, -1] or support.dofs == [1, 1, 1]:
                marker = 's'
            elif support.dofs == [1]:
                if node_coord[0] == 0:
                    marker = '3'
                else:
                    marker = '4'
            else:
                marker = '2'
            ax.scatter(node_coord[0], node_coord[1], s=50, c='k',
                       marker=marker)

        # if self.truss:
        #    self.truss.plot(show=False, print_text=print_text, color=color)
        if loads:
            self.plot_loads()
        ax.axis('equal')
        if show:
            plt.axis('equal')
            plt.show()

    def plot_loads(self):
        """ Plots loads

            Point loads are blue arrows
            Starting point is the given coordinate

            Line loads are red arrows
            Plots loads above given member

            TODO: Load value plotting
        """

        for load in self.point_loads.values():
            x, y = load.coordinate
            scl = max(self.L, self.H) * 1e-1
            dx, dy, _ = load.v
            plt.arrow(x, y, np.sign(dx) * scl, np.sign(dy) * scl,
                      head_width=scl * 1e-1, ec='b')
            # lt.scatter(x, y, c='b', marker='*')

        for lineload in self.line_loads.values():
            c0, c1 = lineload.coordinates
            q1, q2 = lineload.values
            dq = (q2 - q1) / 10
            x1, y1 = c0
            x2, y2 = c1
            num_arrows = 10
            dx = (x2 - x1) / num_arrows
            dy = (y2 - y1) / num_arrows
            if dx:
                X = np.arange(x1, x2 + dx, dx)
            else:
                X = np.ones(11) * x1
            if dy:
                Y = np.arange(y1, y2 + dy, dy)
            else:
                Y = np.ones(11) * y1
            scl = max(self.L, self.H) * 8e-2
            for i, (x, y) in enumerate(zip(X, Y)):
                q = q1 + dq * i
                q_scl = q / max(abs(q2), abs(q1))
                if lineload.direction == 'y':
                    # Moves arrows above member
                    y += (np.sign(q1) - 2e-1) * scl * q_scl
                    dx = 0
                    dy = -np.sign(q1) * scl * q_scl
                    plt.arrow(x, y, dx, dy,
                              head_width=scl * 1e-1, ec='r',
                              head_starts_at_zero=False)
                else:
                    x -= (np.sign(q1) + 0.25) * scl * q_scl
                    dx = np.sign(q1) * scl * q_scl
                    dy = 0
                    plt.arrow(x, y, dx, dy,
                              head_width=scl * 1e-1, ec='y')

            # plt.plot([c0[0], c1[0]], [c0[1], c1[1]], c='b')

            # x0, y0 = self.f.elements[lineload.element_ids[0]].nodes[0].x
            # x1, y1 = self.f.elements[lineload.element_ids[-1]].nodes[1].x
            # plt.plot([x0, x1], [y0, y1], c='b')

    def plot_normal_force(self, show=True):
        """ Plots normal force and utilization ratio

            Parameters
            ------------
            :param show: Shows the plotted diagram, allows to plot multiple
                        diagrams to be plotted on same picture
            :type show: bool


        """
        for member in self.members:
            member.plot_normal_force()
        # if self.truss:
        #    for member in self.truss.members.values():
        #        member.plot_normal_force()
        if show:
            plt.show()

    def plot_deflection(self, scale=1, prec=4, show=True, load_id=2):
        """ Draws deflected shape of the frame

            Parameters
            ----------
            :param scale: Scaling factor
            :param prec: Precision for displacement value
            :param show: Shows the plotted diagram, allows to plot multiple
                        diagrams to be plotted on same picture

            :type scale : float
            :type prec : int
            :type show: bool

        """
        # if self.truss:
        #    self.truss.plot_deflection(scale, show=False)

        self.plot(print_text=False, show=False)
        self.calc_nodal_displacements()
        max_locs = []
        for member in self.members:
            X = []
            Y = []

            member.calc_nodal_displacements(self.f)
            """ For columns,  highlight maximum horizontal displacement (max_x)
                For beams, highlight maximum vertical displacements (max_y)
            """
            max_x = 0
            max_y = 0

            """ Calculate the deflected location of each node of the member """
            for i, node in enumerate(member.nodes.values()):
                x0, y0 = node.coord
                x1 = node.u[load_id][0]
                y1 = node.u[load_id][1]
                x = x0 + x1 * (scale)
                y = y0 + y1 * (scale)

                """ Update value and location of the maximum displacmeent """
                if abs(x1) >= abs(max_x) and member.mtype == "column":
                    max_x = abs(x1)
                    loc_max_x = x
                    loc_max_y = y
                if abs(y1) >= abs(max_y) and member.mtype != "column":
                    max_y = abs(y1)
                    loc_max_x = x
                    loc_max_y = y

                """ Store deflected locations to X and Y """
                X.append(x)
                Y.append(y)

            """ Plot deflected locations """
            plt.plot(X, Y, color='gray')
            if (loc_max_x, loc_max_y) not in max_locs:

                max_locs.append((loc_max_x, loc_max_y))

                if member.mtype != "column":
                    plt.plot(loc_max_x, loc_max_y, 'ro')
                    plt.text(loc_max_x, loc_max_y,
                             "{0:5.{1}g} mm".format(max_y, prec))

                else:
                    plt.plot(loc_max_x, loc_max_y, 'ro')
                    plt.text(loc_max_x, loc_max_y,
                             "{0:5.{1}g} mm".format(max_x, prec))

        if show:
            plt.show()

    def bmd(self, scale=1):
        """ Draws bending moment diagram

            Parameters
            ----------
            :param scale: Scaling factor
            :type scale : int
        """
        self.plot(print_text=False, show=False, color=False)
        for member in self.members:
            member.bmd_test(scale)
        # for truss in self.truss:
        #     truss.bmd(scale)
        plt.show()

    @property
    def weight(self):
        """ Calculates frame's weight
        """
        weight = 0
        for member in self.members:
            weight += member.weight
        return weight

if __name__ == '__main__':
    from materials.timber_data import T
    fr = TimberFrame2D()
    # coord1 = [[0, 0], [0, 5000]]
    # coord2 = [[0, 5000], [10000, 5000]]
    # coord3 = [[10000, 0], [10000, 5000]]
    # coord4 = [[0, 5000], [0, 10000]]
    # coord5 = [[0, 10000], [10000, 10000]]
    # coord6 = [[10000, 5000], [10000, 10000]]
    # fr.add(TimberMember(T.GL32c, 200, 200, coord1, 'instantaneous', 1))
    # fr.add(TimberMember(T.GL32c, 400, 200, coord2, 'medium_term', 1, num_elements=8))
    # fr.add(TimberMember(T.GL32c, 200, 200, coord3, 'medium_term', 1))
    # fr.add(TimberMember(T.GL32c, 200, 200, coord4, 'instantaneous', 1))
    # fr.add(TimberMember(T.GL32c, 400, 200, coord5, 'medium_term', 1, num_elements=8))
    # fr.add(TimberMember(T.GL32c, 200, 200, coord6, 'medium_term', 1))
    # #fr.add(f2d.PointLoad([2000.0, 5000], [0, -100, 0]))
    # fr.add(f2d.LineLoad(fr.members[1], [-10, -10], 'y'))
    # fr.add(f2d.LineLoad(fr.members[4], [-8, -8], 'y'))
    # fr.add(f2d.LineLoad(fr.members[0], [2, 3], 'x'))
    # fr.add(f2d.LineLoad(fr.members[3], [3, 4], 'x'))
    # fr.add(f2d.FixedSupport([0,0]))
    # fr.add(f2d.FixedSupport([10000,0]))
    # fr.generate()
    # fr.calculate()
    # #print(fr.members[0].nodal_displacements)
    # fr.plot()
    # fr.bmd(10)
    #
    # #print(fr.members[0].Ned, fr.members[1].Ved, fr.members[2].Med)
    # print(' Norm     ', 'Leikk     ', 'Mom   ', 'Vään  ', 'TaivEto  ', 'TaivPur  ', 'MV')
    # print(fr.members[0].check_section())
    # print(fr.members[1].check_section())
    # print(fr.members[2].check_section())
    # print(fr.members[3].check_section())
    # print(fr.members[4].check_section())
    # print(fr.members[5].check_section())
    # #print(fr.members[3].ned, ' ned')
    # print(fr.members[4].myed,' myed')
    # #print(fr.members[0].mzed)
    # #print(fr.members[0].vyed)
    # print(fr.members[4].vzed, ' vzed')
    # #fr.plot()


    coord1 = [[0, 0], [0, 5000]]
    coord3 = [[0, 5000], [6000, 5000]]
    coord2 = [[6000, 5000], [6000, 0]]
    fr.add(TimberMember(T.GL32c, 200, 200, coord1, 'instantaneous', 1))
    fr.members[0].section = TimberSection(fr.members[0].material, fr.members[0].B, fr.members[0].H, fire_protection_generic=GypsumPlasterboardF(1))
    fr.add(TimberMember(T.GL32c, 400, 200, coord3, 'medium_term', 1, num_elements=8, Sj1=np.inf, Sj2=np.inf))
    fr.add(TimberMember(T.GL32c, 200, 200, coord2, 'instantaneous', 1))
    fr.members[2].R = 60
    #fr.add(TimberMember(T.Kerto_Q_21_24))
    #fr.add(f2d.PointLoad([50.0, 2500.0], [30000, 0, 0]))
    fr.add(f2d.LineLoad(fr.members[1], [-10, -10], 'y', load_id=1))
    fr.add(f2d.LineLoad(fr.members[0], [10, 10], 'x'))
    fr.add(f2d.XYHingedSupport([0, 0]))
    fr.add(f2d.XYHingedSupport([6000, 0]))
    fr.generate()
    fr.calculate()
    #print(fr.members[0].check_section())
    fr.members[1].check_cross_section()
    print(fr.members[1].r)
    fr.bmd(1)
    fr2 = TimberFrame2D()


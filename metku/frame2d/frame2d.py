"""
@author: Jaakko Huusko
"""

import os
import sys
import time
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.path as mpath
from loadIDs import LoadIDs

try:
    import metku.framefem.framefem as fem
    from metku.framefem.elements import EBBeam, EBSemiRigidBeam, Rod
    from metku.sections.steel import *
    from metku.sections.steel.catalogue import *
    from metku.frame2d.materials import MATERIALS
    from metku.structures.steel.steel_member import SteelMember
    from metku.sections.timber.timber_section import TimberSection, TaperedSection, SmallElementSection, PulpettiPalkki, HarjaPalkki, KaarevaHarjaPalkki, KaarevaPalkki, MahaPalkki
    from metku.structures.timber.timber_member import TimberMember
    from metku.materials.timber_data import T, Timber
    from metku.sections.steel.catalogue import ipe_profiles, h_profiles, rhs_profiles, shs_profiles, chs_profiles
    from metku.sections.steel.CHS import CHS
    from metku.sections.steel.RHS import RHS, SHS
    from metku.sections.steel.ISection import IPE, HEA, HEB
except:
    import framefem.framefem as fem
    from framefem.elements import EBBeam, EBSemiRigidBeam, Rod
    from materials.steel_data import Steel
    from sections.steel import *
    from frame2d.materials import MATERIALS
    from structures.steel.steel_member import SteelMember
    from sections.timber.timber_section import TimberSection, TaperedSection, SmallElementSection, PulpettiPalkki, HarjaPalkki, KaarevaHarjaPalkki, KaarevaPalkki, MahaPalkki
    from structures.timber.timber_member import TimberMember
    from materials.timber_data import T, Timber
    from sections.steel.catalogue import ipe_profiles, h_profiles, rhs_profiles, shs_profiles, chs_profiles
    from sections.steel.CHS import CHS
    from sections.steel.RHS import RHS, SHS
    from sections.steel.ISection import IPE, HEA, HEB

from time import time
# Rounding precision (number of decimals)
PREC = 3

if sys.version_info[0] == 3 and sys.version_info[1] < 9:
    from typing import List, Tuple, Dict
    list_array = List[np.ndarray]
    #list_out = List[float]
    llist_out = List[List[float]]
    lllist_out = List[List[List[float]]]
    
    dict_str_float = Dict[str, float]
    
    list_out = List
    tuple_out = Tuple
else:
    list_array = list[np.ndarray]
    #list_out = list[float]
    llist_out = list[list[float]]
    lllist_out = list[list[list[float]]]
    
    dict_str_float = dict[str: float]
    
    list_out = list
    tuple_out = tuple


# -----------------------------------------------------------------------------
class Frame2D:
    """ Class for 2D frame
    
        Methods:
        ------------
        add
        generate
        calculate
        plot

        Parameters:
        ------------
        :param simple: creates a simple frame with given list
                    [storeys, bays, storey height, bay length]
        :param num_elements: number of elements per member
        :param supports: Generates supports to simple frame
                    'FIXED', 'XHINGED', 'YHINGED', 'XYHINGED'
        :param fem: FrameFEM -instance, used to share fem model
                    with other Frame2D -instances
                    default: creates a new  FrameFEM -instance

        :type simple: list
        :type num_elements: int
        :type supports: string
        :type fem: FrameFEM


        Variables:
        ----------
        :ivar nodes: list of nodes in the FEM model
        :ivar nodal_coordinates: list of nodal coordinates
        :ivar num_elements: number of elements per member
        :ivar f: Finite element model of the frame (FrameFEM)
        :ivar alpha_cr: critical load factor from linear stability analysis (list)
        :ivar members: dict of members (FrameMember)
        :ivar support_nodes: list of supported nodes
        :ivar point_loads: dict of point loads. key: name of the load (string) 
        :ivar line_loads: dict of line loads. key: name of the load (string)
        :ivar nodal_forces:
        :ivar nodal_displacements:
        :ivar joints:
        :ivar r: member utilization ratios (list)
        :ivar is_generated:
        :ivar is_calculated:
        :ivar self_weight:
        :ivar truss:
        :ivar simple:
    """

    def __init__(self, simple=None, num_elements=None, supports=None,
                 fem_model=None, create_beams=True):

        if fem_model:
            self.f = fem_model
        else:
            self.f = fem.FrameFEM()
        self.alpha_cr = []
        self.members = {}
        self.support_nodes = []
        self.num_elements = num_elements
        self.load_ids = []  # List of load case ID numbers
        self.point_loads = {}
        self.line_loads = {}
        self.supports = {}
        self.nodes = []
        self.nodal_coordinates = []
        #self.nodal_forces = {}          # For storing nodal forces in various load cases
        #self.nodal_displacements = {}   # For storing nodal displacements in various load cases
        self.joints = {}
        self.r = []
        self.is_generated = False
        self.is_calculated = False
        self.self_weight = False
        self.truss = []
        self.simple = simple
        self.create_beams = create_beams

        if simple:
            self.storeys = simple[0]
            self.bays = simple[1]
            self.storey_height = simple[2]
            self.bay_length = simple[3]
            self.generate_members()

        if supports:
            self.generate_supports(supports)

    @property
    def R(self):
        return self.members[0].R
    @R.setter
    def R(self, r):
        for mem in self.members.values():
            mem.R = r

    @property
    def L(self):
        """ Width/span of the frame """
        x_coordinates = [mem.coordinates[0][0] for mem in
                         self.members.values()]
        x_coordinates.extend(
            [mem.coordinates[1][0] for mem in self.members.values()])
        x_min = min(x_coordinates)
        x_max = max(x_coordinates)
        L = x_max - x_min
        return L

    @property
    def H(self):
        """ Height of the frame """
        y_coordinates = [mem.coordinates[0][1] for mem in
                         self.members.values()]
        y_coordinates.extend(
            [mem.coordinates[1][1] for mem in self.members.values()])
        y_min = min(y_coordinates)
        y_max = max(y_coordinates)
        H = y_max - y_min
        return H

    @property
    def beams(self):
        """ Returns a list of beam members """
        return [mem for mem in self.members.values() if mem.mtype == 'beam']

    @property
    def columns(self):
        """ Returns a list of column members """
        return [mem for mem in self.members.values() if mem.mtype == 'column']

    @property
    def load_cases(self):
        """ Returns a list of load case labels """
        return self.load_ids
        #return list(set([lc.load for lc in self.f.loadcases.values()]))

    def add(self, this):
        """ Adds given item to the frame.

            Parameters:
            -----------

            :param this:  item to be added to frame

            :type this: FrameMember, Support, PointLoad, LineLoad, Truss2D
        """
        # MEMBERS
        if isinstance(this, FrameMember):
            # Give member an unique id
            # id is used in creating elements with wanted cross-sectional properties
            this.mem_id = int(len(self.members))

            this.frame2d = self

            self.members[this.mem_id] = this
            if self.num_elements:
                this.calc_nodal_coordinates(self.num_elements)
            else:
                this.calc_nodal_coordinates()

            # Check if members intersect
            # If coordinate is outside of member's coordinates,
            # it's rejected in member's add_node_coord function
            for member in self.members.values():
                coord = this.line_intersection(member.coordinates)
                if isinstance(coord, list):
                    this.add_node_coord(coord)
                    member.add_node_coord(coord)

            if str(this.nodal_coordinates) not in [str(c) for c in
                                                   self.nodal_coordinates]:
                self.nodal_coordinates.extend(this.nodal_coordinates)
                self.nodal_coordinates.sort()




        # POINTLOADS
        elif isinstance(this, PointLoad):

            """ If a point load with same 'name' is already included
                in the frame, add a number to the end of the name
                of the new load such that each point load has a unique
                name.
                
                NOTE: using 'len' to provide new id number might fail,
                if point loads are deleted such that the number of point loads
                changes along the way.
            """
            if not (this.load_id in self.load_ids):
                self.load_ids.append(this.load_id)
            
            if this.name in self.point_loads.keys():
                this.name += str(len(self.point_loads))
            self.point_loads[this.name] = this

            """ If the location of the point load is in between
                member end nodes, add a node to the corresponding member
            """
            for member in self.members.values():
                if member.point_intersection(this.coordinate):
                    member.add_node_coord(this.coordinate)


        # LINELOADS
        elif isinstance(this, LineLoad):

            if not (this.load_id in self.load_ids):
                self.load_ids.append(this.load_id)
            
            if this.name in self.line_loads.keys():
                this.name += str(len(self.line_loads))
            self.line_loads[this.name] = this


        # SUPPORTS
        elif isinstance(this, Support):
            # this.supp_id = len(self.supports)
            supp_label = len(self.supports)
            self.supports[supp_label] = this
            # self.supports[this.supp_id] = this
            self.support_nodes.append(this.coordinate)

            """ If the support is located between end nodes of a member,
                add a node to that member
            """
            for member in self.members.values():
                #coord = member.point_intersection(this.coordinate)
                if member.point_intersection(this.coordinate):
                    member.add_node_coord(this.coordinate)

        # TRUSS
        elif type(this).__name__ == 'Truss2D':
            if self.f != this.f:
                this.f = self.f
            self.truss.append(this)

            # Bottom chord coordinates
            for bchord in this.bottom_chords:
                c0, c1 = bchord.coordinates
                for col in self.columns:
                    coord = col.line_intersection(bchord.coordinates)
                    # If bottom chord and column intersect
                    if isinstance(coord, list):
                        col.add_node_coord(coord)
                        c0, c1 = bchord.coordinates
                        coord = np.asarray(coord)
                        # Distance between coordinates
                        if np.linalg.norm(np.asarray(coord) - c0) <= 1e-2:
                            c0 = coord + col.cross_section.h / 2 * bchord.unit
                            bchord.coordinates = [c0, c1]
                            col.ecc_coordinates.append([coord, c0])
                        elif np.linalg.norm(np.asarray(coord) - c1) <= 1e-2:
                            c1 = coord - col.cross_section.h / 2 * bchord.unit
                            bchord.coordinates = [c0, c1]
                            col.ecc_coordinates.append([coord, c1])
                        bchord.calc_nodal_coordinates()
                        for joint in bchord.joints:
                            if joint.loc == 0:
                                joint.loc = (col.cross_section.h / 2 + joint.g1) / bchord.length
                            elif joint.loc == 1:
                                joint.loc = 1 - (col.cross_section.h / 2 + joint.g1) / bchord.length
                            else:
                                joint.calc_nodal_coordinates()

                bchord.nodal_coordinates.sort(
                    key=lambda c: np.linalg.norm(
                        np.asarray(c0) - np.asarray(c)))

            # Top chord coordinates
            for tchord in this.top_chords:
                c0, c1 = tchord.coordinates
                for member in self.members.values():
                    coord = member.line_intersection(tchord.coordinates)
                    if isinstance(coord, list):
                        member.add_node_coord(coord)
                        # Create eccentricity to connection
                        if member.mtype == "column":

                            joint0 = [j for j in tchord.joints if
                                      j.loc <= 0.02]
                            joint1 = [j for j in tchord.joints if
                                      j.loc >= 0.98]

                            # Chord's position vector
                            v = np.array([math.cos(tchord.angle),
                                          math.sin(tchord.angle)])
                            # Vector perpendicular to v
                            u = np.array([-math.sin(tchord.angle),
                                          math.cos(tchord.angle)])

                            # Start coordinate
                            if coord == c0 and joint0:
                                joint0 = joint0[0]
                                web = list(joint0.webs.values())[0]
                                theta = abs(joint0.chord.angle - web.angle)
                                # ecc_x = member.h / 2 + abs(joint0.g1*math.cos(joint0.chord.angle) +
                                #        web.h/2 * math.sin(theta))
                                ecc_x = member.h / 2 + joint0.g1

                                coord0 = np.asarray(joint0.coordinate)
                                # Change joint's coordinate
                                joint0.coordinate = list(
                                    c0 + ecc_x * tchord.unit)
                            if coord == c0:
                                # Calculate chord's new start coordinate
                                # c0 = [c0[0] + tchord.h / 2000 * u[0], c0[1] + tchord.h / 2000 * u[1]]
                                ecc_y = -tchord.h / 2 * tchord.perpendicular[1]
                                # c0 = [c0[0] + member.h/2000, c0[1]]
                                # Add eccentricity elements coordinates to list
                                coord[1] -= ecc_y
                                member.ecc_coordinates.append([coord, c0])
                                tchord.columns[c0[0]] = member
                                temp_coords = sorted(member.coordinates,
                                                     key=lambda x: x[1])
                                temp_coords[1] = coord
                                member.coordinates = temp_coords
                                member.calc_nodal_coordinates()

                            # End coordinate
                            if coord == c1 and joint1:
                                joint1 = joint1[0]
                                web = list(joint1.webs.values())[0]
                                theta = abs(joint1.chord.angle - web.angle)
                                # ecc_x = member.h / 2 + abs(joint1.g1*1000*math.cos(joint1.chord.angle) +
                                #        web.h/2 * math.sin(theta))  
                                ecc_x = member.h / 2 + joint1.g1
                                coord1 = np.asarray(joint1.coordinate)
                                joint1.coordinate = list(
                                    c1 - ecc_x * tchord.unit)
                            if coord == c1:
                                # Calculate chord's new end point
                                # c1 = [c1[0] + tchord.h/2000 * u[0], c1[1] + tchord.h/2000 * u[1]]
                                # c1 = [c1[0] - member.h/2000, c1[1]]
                                # member.ecc_coordinates.append([coord, c1])
                                # tchord.columns[c1[0]] = member
                                ecc_y = -tchord.h / 2 * tchord.perpendicular[1]
                                # c0 = [c0[0] + member.h/2000, c0[1]]
                                # Add eccentricity elements coordinates to list
                                coord[1] -= ecc_y
                                member.ecc_coordinates.append([coord, c1])
                                tchord.columns[c1[0]] = member
                                temp_coords = sorted(member.coordinates,
                                                     key=lambda x: x[1])
                                temp_coords[1] = coord
                                member.coordinates = temp_coords
                                member.calc_nodal_coordinates()

                # Change chord's ends' coordinates
                tchord.coordinates = [c0, c1]
                tchord.calc_nodal_coordinates()

            # Add truss's members to frame's members dict
            for key in this.members.keys():
                new_id = len(self.members)
                this.members[key].mem_id = new_id
                self.members[new_id] = this.members[key]

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
        for member in self.members.values():
            member.add_material(self.f)
            member.add_section(self.f)

    # This method is not utilized later in code ?!?!?!
    def calc_nodal_coordinates(self):
        """ Calculates nodal coordinates and saves values to a list.
            These values are later used for creating the nodes.
        """
        # Iterate through every frame's member and calculate nodal coordinates.
        # If the x,y -coordinate is not in the list add that to the list
        for member in self.members.values():
            member.calc_nodal_coordinates(self.num_elements)
            for coordinate in member.nodal_coordinates:
                if coordinate not in self.nodal_coordinates:
                    self.nodal_coordinates.append(coordinate)

    def calc_nodal_forces(self,lcase=LoadIDs.ULS):
        """ Calculates nodal forces and saves values to
            self.nodal_forces - dict       
        """
        # Iterate through every frame's member's node and calculate forces
        # acting on that node
        for member in self.members.values():
            member.calc_nodal_forces(lcase)
            #for node in member.nodal_forces:
            #    self.nodal_forces[node] = member.nodal_forces[node]

    def calc_nodal_displacements(self,lcase=LoadIDs.ULS):
        """ Calculates nodal displacements and saves values to
            self.nodal_displacements -dict
        """
        # Iterate through every frame's member's node and calculate 
        # displacements for that node
        for member in self.members.values():
            member.calc_nodal_displacements(self.f,lcase)
            #for node in member.nodal_displacements.keys():
            #    self.nodal_displacements[node] = member.nodal_displacements[
            #        node]


    def update_model(self):
        """
        Updates frame's geometry
        :return:
        """
        for mem in self.members.values():
            mem.calc_nodal_coordinates()
        if len(self.truss):
            for t in self.truss:
                for j in t.joints.values():
                    j.calc_nodal_coordinates()

    def calculate(self, load_id=LoadIDs.ULS, support_method='ZERO'):
        """ Calculates forces and displacements
            
            Parameters
            ----------
            :param load_id: Id of the loads to be added in the calculation model
            
            :type load_id: int / str
        """
        # print(f'calculation happening with load id {load_id}')
        if self.is_calculated == False:
            self.is_calculated = True
            """ Support ID is always 1! """
            self.f.nodal_dofs()

        # If load_id == 'ALL' calculates all load cases
        # calls recursively itself for each case
        if str(load_id).upper() == 'ALL':
            #print("Calculate all cases")
            lcase_ids = self.load_cases #[lc.load for lc in self.f.loadcases.values()]            
            for lid in lcase_ids:                
            #for lid in self.f.loadcases.values():
                self.calculate(load_id=lid,
                               support_method=support_method)
        else:
            #print('Calculate case:' + str(load_id))

            self.f.linear_statics(support_method=support_method,
                                  lcase=load_id)
            self.calc_nodal_forces(load_id)
            self.calc_nodal_displacements(load_id)
            self.assign_forces(load_id)
            self.design_members(load_id)
            #self.check_members_strength()
            # self.alpha_cr, _ = self.f.linear_buckling(k=4)

    def optimize_members(self, prof_type="CURRENT"):
        """ Find minimum smallest profiles for frame members
            for given loads.
        """
        
        PROFILES_CHANGED = True
        kmax = 10
        k = 0
        
        while PROFILES_CHANGED and k < kmax:            
            self.calculate('all','REM')
            for member in self.members.values():
                PROFILES_CHANGED = member.optimum_design(prof_type)

            k += 1
        
        
    def design_members(self,load_id = LoadIDs.ULS):
        """ Designs frame members (checks resistance)
        """
        
        for member in self.members.values():
            member.design(load_id)

        #self.check_members_strength()
        
    def check_members_strength(self):
        """ Checks if members can bear their loads
        """
        self.r.clear()
        for member in self.members.values():
            member.check_cross_section()
            self.r.append(member.r)

    def delete_member(self, id):
        """ Removes member from frame
        TODO!
        """

        member = self.members[id]
        member.delete()

    def generate_nodes(self):
        """ Generates nodes to coordinates in self.nodal_coordinates
        """
        # Start node indexing from zero
        index = 0
        # Iterate through all x,y-coordinates in list and create FrameFEM.Node-
        # objects to those coordinates
        for coordinate in self.nodal_coordinates:
            x = coordinate[0]
            y = coordinate[1]
            self.f.add_node(x, y)
            self.nodes.append(index)
            # if the y-coordinate is zero, save that coordinate as support
            # node coordinate
            if y == 0:
                self.support_nodes.append(self.nodal_coordinates.index([x, y]))
            index += 1
        # Iterate through all frame's members and add created node-objects
        # to their member's own list
        for member in self.members.values():
            member.add_nodes(self.nodal_coordinates)

    def generate_columns(self):
        """ Generates columns and adds them to self.members -dict

        """
        # x-coordinate
        x = 0
        # y-coordinates 1: start, 2: end
        y1 = 0
        y2 = self.storey_height
        # number of columns
        num_cols = (self.bays + 1) * self.storeys
        # create columns as FrameMember object
        for i in range(num_cols):
            coords = [[x, y1], [x, y2]]
            col = SteelColumn(coords)
            self.add(col)
            # If the y2-coordinate equals frame's height, move the x-coordinate
            # by one bay length and initialize y-coordinates
            if y2 == self.storeys * self.storey_height:
                x += self.bay_length
                y1 = 0
                y2 = self.storey_height
            else:
                y1 = y2
                y2 += self.storey_height

    def generate_beams(self):
        """ Generates beams and adds them to self.members -dict """
        # x-coordinate, 1:start, 2: end
        x1 = 0
        x2 = self.bay_length
        # y-coordinate
        y = self.storey_height
        # number of beams
        num_beams = self.bays * self.storeys
        # create beams as FrameMember object
        for i in range(num_beams):
            coords = [[x1, y], [x2, y]]
            beam = SteelBeam(coords)
            self.add(beam)
            # If y-coordinate equals frame's height, initialize y-coordinate
            # and change x-coordinates
            if y == self.storeys * self.storey_height:
                y = self.storey_height
                x1 = x2
                x2 += self.bay_length
            else:
                y += self.storey_height

    def generate_members(self):
        """ Generates columns and beams
            member_id keeps track for running numbering for id's
        """
        self.generate_columns()
        if self.create_beams:
            self.generate_beams()

    def generate_supports(self, supp_type):
        supp_type = supp_type.upper()
        if supp_type == 'FIXED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(FixedSupport(coord))

        elif supp_type == 'YHINGED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(YHingedSupport(coord))

        elif supp_type == 'XHINGED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(XHingedSupport(coord))

        elif supp_type == 'XYHINGED':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(XYHingedSupport(coord))

    def generate(self, self_weight=False):
        """ Generates the frame and truss FEM model
        """
        # If the frame hasn't been created, create members and calculate nodal
        # coordinates, otherwise initialize frame.
        if not self.is_generated:
            self.is_generated = True
        for member in self.members.values():
            member.calc_nodal_coordinates() # Generate FEM nodes for a member
            member.round_coordinates()
            member.generate(self.f)
            for coord in member.nodal_coordinates:
                if coord not in self.nodal_coordinates:
                    self.nodal_coordinates.append(coord)
            self.nodes.extend(list(member.nodes.values()))

        for member in self.members.values():
            if hasattr(member.member, 'add_neighbours_support'):
                member.member.add_neighbours_support()

        # Generate TrussJoints
        if len(self.truss):
            for truss in self.truss:
                for joint in truss.joints.values():
                    joint.generate(self.f)

        # Generate eccentricity elements
        for member in self.members.values():
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

        if self_weight:
            self.add_self_weight()

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


    def to_csv(self, filename, num_frames=0, s=1):

        with open(filename + '.csv', 'w') as f:
            for member in self.members.values():
                for y in range(num_frames):
                    y = y * s
                    member.round_coordinates()
                    start, end = member.coordinates
                    x0, z0 = start
                    x1, z1 = end
                    profile = member.profile
                    profile = profile.replace(" ", "")
                    profile = profile.replace("X", "*")
                    profile = profile.strip()

                    f.write(f'{x0}, {y}, {z0}, {x1}, {y}, {z1},{profile} \n')
        print(f'{filename}.csv created to: \n{os.getcwd()}')

    def to_robot(self, filename, num_frames=1, s=1,
                 brace_profile="SHS 50x50x2"):
        """  Creates an Autodesk Robot Structural Analysis .str -file
        
            Parameters
            ----------
            :param filename: Name of the created file
            :param num_frames: Number of frames 
            :param s: Spacing between created frames
            :param brace_profile: Profile for vertical braces
                
            :type filename: string
            :type num_frames: int
            :type s: float
            :type brace_profile: string
        """
        nodes = []
        elements = []
        profiles = []
        material = []
        releases = {}
        pointloads = []
        vert_braces = []
        for i in range(num_frames):
            for member in self.members.values():
                n1, n2 = member.coordinates
                if num_frames != 1:
                    n1 = [n1[0], i * s, n1[1]]
                    n2 = [n2[0], i * s, n2[1]]
                if n1 not in nodes:
                    nodes.append(n1)
                if n2 not in nodes:
                    nodes.append(n2)
                elem = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                elements.append(elem)
                splitted_val = member.profile.split(" ")
                profile_type = splitted_val[0]
                # Add columns end node to vertical braces list
                if num_frames != 1 and member.mtype == 'top_chord':
                    if nodes.index(n1) + 1 not in vert_braces:
                        vert_braces.append(nodes.index(n1) + 1)

                    if nodes.index(n2) + 1 not in vert_braces:
                        vert_braces.append(nodes.index(n2) + 1)

                if profile_type == 'RHS':
                    profile = 'RRHS ' + splitted_val[1]
                elif profile_type == 'HE':
                    profile = splitted_val[0] + splitted_val[2] + splitted_val[
                        1]
                elif profile_type == "SHS":
                    dims = splitted_val[1].split("X")
                    profile = 'RRHS ' + str(dims[0]) + 'X' + str(
                        dims[0]) + 'X' + str(dims[1])
                else:
                    profile = member.profile
                # Webs
                if member.mtype == "web" or member.mtype == "beam" or \
                        member.mtype == "top_chord" or member.mtype == "bottom_chord":
                    Sj1 = min(1e10, member.Sj1)
                    Sj2 = min(1e10, member.Sj2)
                    releases[elements.index(
                        elem) + 1] = f'ORIgin RY Hy={Sj1}  END RY Hy={Sj2}'

                profiles.append(profile)
                material.append(member.material)

                # Eccentricity elements
                if len(member.ecc_coordinates):
                    for coords in member.ecc_coordinates:
                        n1, n2 = coords
                        if num_frames != 1:
                            n1 = [n1[0], i * s, n1[1]]
                            n2 = [n2[0], i * s, n2[1]]
                        if n1 not in nodes:
                            nodes.append(n1)
                        if n2 not in nodes:
                            nodes.append(n2)
                        elem = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                        elements.append(elem)
                        profiles.append("HEA 1000")
                        material.append("S355")
                        # if n1 < n2:
                        #    releases[elements.index(elem) +1] = 'END RY'
                        # else:
                        #    releases[elements.index(elem) +1] = 'ORIgin RY'

            if self.truss:
                trusses = self.truss
            else:
                trusses = [self]

            for truss in trusses:
                for joint in truss.joints.values():

                    # Y joint
                    if len(joint.nodal_coordinates) == 2:
                        for n, c in joint.nodal_coordinates.items():
                            joint.nodal_coordinates[n] = list(c)
                        n1, n2 = sorted(joint.nodal_coordinates.values())
                        if num_frames != 1:
                            n1 = [n1[0], i * s, n1[1]]
                            n2 = [n2[0], i * s, n2[1]]
                        if n1 not in nodes:
                            nodes.append(n1)
                        if n2 not in nodes:
                            nodes.append(n2)
                        # Eccentricity element
                        elem = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                        elements.append(elem)
                        profiles.append("HEA 1000")
                        material.append("S355")
                    # K joint
                    elif joint.joint_type == "K":
                        for n, c in joint.nodal_coordinates.items():
                            joint.nodal_coordinates[n] = list(c)
                        coords = sorted(joint.nodal_coordinates.values())
                        n1, n2, n3, n4, n5 = coords
                        if num_frames != 1:
                            n1 = [n1[0], i * s, n1[1]]
                            n2 = [n2[0], i * s, n2[1]]
                            n3 = [n3[0], i * s, n3[1]]
                            n4 = [n4[0], i * s, n4[1]]
                            n5 = [n5[0], i * s, n5[1]]
                        if n1 not in nodes:
                            nodes.append(n1)
                        if n2 not in nodes:
                            nodes.append(n2)
                        if n4 not in nodes:
                            nodes.append(n4)
                        if n5 not in nodes:
                            nodes.append(n5)
                        # Eccentricity elements
                        elem1 = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                        elem2 = [nodes.index(n4) + 1, nodes.index(n5) + 1]
                        elements.append(elem1)
                        elements.append(elem2)
                        profiles.append("HEA 1000")
                        material.append("S355")
                        profiles.append("HEA 1000")
                        material.append("S355")
                        # releases[elements.index(elem1)+1] = "END RY"
                        # releases[elements.index(elem2)+1] = "END RY"
                    else:
                        pass

            # except:
            #    pass

            for pointload in self.point_loads.values():
                try:
                    n1 = pointload.coordinate
                    if num_frames != 1:
                        n1 = [n1[0], i * s, n1[1]]
                    idx = nodes.index(n1) + 1
                except ValueError:
                    n1 = pointload.coordinate
                    if num_frames != 1:
                        n1 = [n1[0], i * s, n1[1]]
                    n1 = [val for val in n1]
                    nodes.append(n1)
                    idx = nodes.index(n1) + 1
                FX, FZ, MY = pointload.v
                if FX:
                    pointloads.append(f'    {idx} FX={FX}')
                elif FX and FZ:
                    pointloads.append(f'    {idx} FX={FX}  FZ={FZ}')
                elif FX and MY:
                    pointloads.append(f'    {idx} FX={FX}  MY={MY}')
                elif FZ:
                    pointloads.append(f'    {idx} FZ={FZ}')
                elif FZ and MY:
                    pointloads.append(f'    {idx} FZ={FZ}  MY={MY}')
                elif MY:
                    pointloads.append(f'    {idx} MY={MY}')
                elif FX and FZ and MY:
                    pointloads.append(f'    {idx} FX={FX}  FZ={FZ}  MY={MY}')

        # Vertical braces
        for i in range(len(vert_braces) - 3):
            n1 = vert_braces[i]
            n2 = vert_braces[i + 3]
            # Eccentricity element
            elem = [n1, n2]
            elements.append(elem)
            splitted_val = brace_profile.split(" ")
            splitted_val = [val.upper() for val in splitted_val]
            profile_type = splitted_val[0]
            if profile_type == 'RHS':
                profile = 'RRHS ' + splitted_val[1]
            elif profile_type == 'HE':
                profile = splitted_val[0] + splitted_val[2] + splitted_val[1]
            elif profile_type == "SHS":
                dims = splitted_val[1].split("X")
                profile = 'RRHS ' + str(dims[0]) + 'X' + str(
                    dims[0]) + 'X' + str(dims[2])
            else:
                profile = brace_profile
            profiles.append(profile)
            material.append("S355")
            releases[elements.index(elem) + 1] = f'ORIgin RY  END RY'

        with  open(filename + '.str', 'w') as f:
            f.write("ROBOT97 \n")
            if num_frames != 1:
                f.write("FRAme SPAce \n")
            f.write("NUMbering DIScontinuous \n")
            f.write(f'NODes {len(nodes)}  ELEments {len(elements)} \n')
            f.write("UNIts \n")
            f.write("LENgth=mm	Force=N \n")
            f.write("NODes \n")
            for i, node in enumerate(nodes):
                if num_frames != 1:
                    f.write(
                        f'{i + 1}   {node[0]}    {node[1]}    {node[2]} \n')
                else:
                    f.write(f'{i + 1}   {node[0]}    {node[1]} \n')

            f.write('ELEments \n')
            for i, element in enumerate(elements):
                f.write(f' {i + 1}  {element[0]}    {element[1]} \n')

            f.write('PROperties \n')
            for i in range(len(elements)):
                f.write(f' "{material[i]}" \n')
                f.write(f' {i + 1}  ')
                f.write(f' {profiles[i]}    \n')

            f.write("SUPports \n")
            if self.is_calculated:
                for sup in self.supports.values():
                    for i in range(num_frames):
                        n1 = sup.coordinate
                        if num_frames != 1:
                            n1 = [n1[0], i * s, n1[1]]
                        idx = nodes.index(n1) + 1
                        if sup.dofs == [0]:
                            f.write(f'{idx}  UX  \n')
                        elif sup.dofs == [0, 1]:
                            f.write(f'{idx}  UX UZ \n')
                        elif sup.dofs == [0, 1, 2]:
                            f.write(f'{idx}  \n')
                        elif sup.dofs == [1]:
                            f.write(f'{idx}  UZ  \n')
                        elif sup.dofs == [1, 2]:
                            f.write(f'{idx}  UZ RY  \n')
                        else:
                            f.write(f'{idx}  RY \n')
            else:
                for sup in self.supports.values():
                    for i in range(num_frames):
                        n1 = sup.coordinate
                        if num_frames != 1:
                            n1 = [n1[0], i * s, n1[1]]
                        idx = nodes.index(n1) + 1
                        if sup.dofs == [1, 0, 0]:
                            f.write(f'{idx}  UX  \n')
                        elif sup.dofs == [1, 1, 0]:
                            f.write(f'{idx}  UX UZ \n')
                        elif sup.dofs == [1, 1, 1]:
                            f.write(f'{idx}  \n')
                        elif sup.dofs == [0, 1, 0]:
                            f.write(f'{idx}  UZ  \n')
                        elif sup.dofs == [0, 1, 1]:
                            f.write(f'{idx}  UZ RY  \n')
                        else:
                            f.write(f'{idx}  RY \n')

            f.write("RELeases \n")
            for elem in releases.keys():
                f.write(f'ELEments{elem} {releases[elem]} \n')

            f.write("LOAds \n")
            f.write("CASe # 1 LC1 \n")
            if len(self.line_loads):
                f.write("ELEments \n")
                for i in range(num_frames):
                    for lineload in self.line_loads.values():
                        n1, n2 = lineload.member.coordinates
                        if num_frames != 1:
                            n1 = [n1[0], i * s, n1[1]]
                            n2 = [n2[0], i * s, n2[1]]
                        n1_idx = nodes.index(n1) + 1
                        n2_idx = nodes.index(n2) + 1
                        idx = elements.index([n1_idx, n2_idx]) + 1
                        if lineload.direction == 'y':
                            dir = 'PZ'
                        else:
                            dir = 'PX'
                        q0, q1 = lineload.values
                        if q0 != q1:
                            f.write(f' {idx} X=0.0 {dir}={q0} TILl  ')
                            f.write(f'X=1.000  {dir}={q1}      RElative \n')

                        else:
                            f.write(f' {idx} {dir}={q0} \n')
            if len(pointloads):
                f.write("NODes \n")
                for val in pointloads:
                    f.write(val + "\n")

            f.write("END")

        print(f'{filename}.str created to: \n{os.getcwd()}')

    def plot(self, print_text=True, show=True,
             loads=True, color=False, axes=None, save=False):
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
        for member in self.members.values():            
            member.plot(print_text, color, ax)

        # Plot joints
        for joint in self.joints.values():
            joint.plot(color=color)

        # Plot supports
        for support in self.supports.values():
            node_coord = support.coordinate
            if support.dofs == [-1, -1, -1] or support.dofs == [1, 1, 1] or support.dofs == [0, 1, 2]:
                marker = 's'
            elif support.dofs == [1]:
                marker = '^'
            elif support.dofs == [0]:
                marker = '>'
            else:
                marker = 'D'
            ax.scatter(node_coord[0], node_coord[1], s=50, c='k',
                        marker=marker)

        # if self.truss:
        #    self.truss.plot(show=False, print_text=print_text, color=color)
        if loads:
            self.plot_loads()
        ax.axis('equal')
        if save:
            plt.savefig('default.svg', format='svg')
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

    def plot_deflection(self, scale=1, prec=4, show=True, save=False, load_id=LoadIDs.ULS):
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
        self.calc_nodal_displacements(lcase=load_id)
        max_locs = []
        for member in self.members.values():
            X = []
            Y = []

            member.calc_nodal_displacements(self.f, lcase=load_id)
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
                             "{0:5.{1}g} mm".format(max_y,prec))

                else:
                    plt.plot(loc_max_x, loc_max_y, 'ro')
                    plt.text(loc_max_x, loc_max_y,
                             "{0:5.{1}g} mm".format(max_x,prec))

        if save:
            plt.savefig('deflections.svg', format='svg')
        if show:
            plt.show()

    def plot_buckling(self, scale=1, k=4, show=True, load_id=LoadIDs.ULS, axes=None):
        """ Draws buckling shapes of the frame
            
            Parameters
            ----------
            :param scale: Scaling factor
            :param k: number of buckling modes
            :param show: Shows the plotted diagram, allows to plot multiple
                        diagrams to be plotted on same picture
            
            :type scale : float, optional
            :type k: int
            :type show: bool
        """

        """ Perform linear buckling analysis 
            output:
                w .. list of eigenvalues
                v .. list of eigenmodes (vectors)
        """
        w, v = self.f.linear_buckling(k=k)

        """ Find the buckling mode corresponding to the smallest eigenvalue """
        nd = np.argmin(w)
        vcrit = v[nd]
        acrit = np.min(w).real
        print("w=",acrit)
        print("v=",vcrit)
        self.alpha_cr = acrit

        sorted_nd = np.argsort(w)

        #for j in range(v.shape[1]):
        for j in sorted_nd[:1]:
            if w[j] < 1e5:
                if axes is None:
                    fig, ax = plt.subplots(1)
                else:
                    ax = axes
                if show:
                    self.plot(print_text=False, show=False, loads=False,
                              color=False,axes=ax)
                for member in self.members.values():
                    X = []
                    Y = []
                    for i, node in enumerate(member.nodes.values()):
                        x0, y0 = node.coord
                        if not isinstance(node.v[0], int):
                            x1 = node.v[j][0]
                            y1 = node.v[j][1]
                        else:
                            x1 = node.v[0]
                            y1 = node.v[1]

                        print(x1,y1)
                        x = x0 + x1 * (scale)
                        y = y0 + y1 * (scale)
                        X.append(x)
                        Y.append(y)

                    ax.plot(X, Y, color='m')
                ax.set_title(f'Buckling shape {j + 1},  ' r"$\alpha_{cr}$" f' = {w[j].real:.2f}')
                # if self.truss:
                #    self.truss.plot_buckling(show=False)

                plt.show()

    def bmd(self, scale=1, save=False, show=True, load_id=LoadIDs.ULS):
        """ Draws bending moment diagram

            Parameters
            ----------
            :param scale: Scaling factor
            :type scale : int
        """
        self.plot(print_text=False, show=False, color=False)
        for member in self.members.values():
            member.bmd(scale, load_id=load_id)
        # for truss in self.truss:
        #     truss.bmd(scale)
        if save:
            plt.savefig('bending moment diagram.svg', format='svg')
        if show:
            plt.show()

    def plot_normal_force(self, show=True, save=False):
        """ Plots normal force and utilization ratio

            Parameters
            ------------
            :param show: Shows the plotted diagram, allows to plot multiple
                        diagrams to be plotted on same picture
            :type show: bool


        """
        for member in self.members.values():
            member.plot_normal_force()
        # if self.truss:
        #    for member in self.truss.members.values():
        #        member.plot_normal_force()
        if save:
            plt.savefig('normal forces.svg', format='svg')
        if show:
            plt.show()

    @property
    def weight(self):
        """ Calculates frame's weight
        """
        weight = 0
        for member in self.members.values():
            weight += member.weight
        return weight

    def add_self_weight(self):
        """ Adds self weight to frame's members
        """
        if self.self_weight is False:
            self.self_weight = True
            for member in self.members.values():
                member.add_self_weight(self)
            # self.add_loads('self_weight')
            # self.generate_frame()
        else:
            self.remove_self_weight()
            self.add_self_weight()

    def remove_self_weight(self):
        """ Removes self weight from frame's members
        """
        self.self_weight = False
        del (self.line_loads["self_weight"])
        for member in self.members.values():
            member.remove_self_weight()
        self.generate_frame()

    def assign_forces(self, load_id=LoadIDs.ULS):

        for member in self.members.values():
            member.assign_forces(load_id=load_id)

    def rigid_joints(self):
        """ Set all beam-column joints to rigid
        """
        for member in self.members.values():
            if member.mtype == "beam":
                member.alpha1 = np.inf
                member.alpha2 = np.inf
                member.Sj1 = np.inf
                member.Sj2 = np.inf

    def hinge_joints(self):
        """ Set all beam-column joints to hinges
        """
        MIN_VAL = 0
        for member in self.members.values():
            if member.mtype != "column":
                member.alpha1 = MIN_VAL
                member.alpha2 = MIN_VAL
                member.Sj1 = MIN_VAL
                member.Sj2 = MIN_VAL

    def print_info(self, prec, R=0, beam_characteristic_divider=300,
                   beam_quasi_permanent_divider=200, column_divider=150, apxspt_kaava6_55=False):

        if LoadIDs.ACC in self.load_ids:
            self.R = R
            self.calculate(LoadIDs.ACC)
            self.R = 0

        for mem in self.members.values():
            if isinstance(mem, TimberFrameMember):
                mem.print_info(self.load_ids, self, prec, firetime=R, beam_characteristic_divider=beam_characteristic_divider,
                               beam_quasi_permanent_divider=beam_quasi_permanent_divider, column_divider=column_divider,
                               apxspt_kaava6_55=apxspt_kaava6_55)

# -----------------------------------------------------------------------------
class FrameMember:
    """ Class for frame member
    
        Written by Jaakko Huusko
        
        Attributes:
        ========   
        f -- FrameFEM class
        elements -- dictionary of elements, key is element id e.g. 0
        nodes -- list of nodes of the member
        nodal_coordinates -- list of lists consisting nodal coordinates
        loc -- list local coordinates, 0 <= loc <= 1
        coordinates -- list of members start and end coordinates
        profile -- name of member's profile
        material -- material code of used steel
        length -- length of member in meters
        mtype -- type of member, "beam" or "column"
        mem_id -- member id
        nodal_forces -- dict of nocal forces, key: node id
        nodal_displacements -- dict of nocal displacements, key: node id
        Sj1 -- rigidness of the left end connection, only for beams
        Sj2 -- rigidness of the right end connection, only for beams
        ned -- max normal force acting on member in kilo Newtons
        ved -- max shear force acting on member in kilo Newtons
        med -- max bending moment acing on member in kilo Newton meters
        has_load -- boolean value if member has load
        is_designed -- boolean value if member is designed
        is_strong_enough -- boolean value if member is strong enough
        r -- member's maximum utilization ratio
        self_weight -- boolean value if self weight is added
            
        Methods:
        ========   
        calc_nodal_coordinates -- calculates coordinates of nodes
        add_nodes -- adds nodes to SemiRigigFrame nodes list
        generate_elements -- generates elements for member
        calc_nodal_forces -- calculates force values on nodes
        calc_nodal_displacement -- calculates nodes' displacements
        bmd -- draws bending moment diagram
        add_line_load -- adds line load for all elements in member
        add_self_weight -- adds self-weight as line load for all elements
        remove_self_weight -- removes self-weight load from f.loads list
        cross_section_properties -- get's cross-section properties
        add_material -- adds member's material to f.materials
        add_section -- adds membetr's cross-section to f.sections
        check_cross_section -- checks if current cross-section is strong enough
        design_members -- designs member for current forces
            
        Uses classes:
        ========   
        Frame2D
        FrameFEM
        I_sections
        CrossSection
        SteelMember
            
    """

    def __init__(self, coordinates, mem_id="", profile="IPE 100",
                 material="S355", num_elements=6,
                 Sj1=np.inf, Sj2=np.inf, mtype="", LT_buckling=False,
                 reverse=False):

        self.active = True          # If 'active' is False, this memer is not included in analysis or plotting.
        self.frame2d = None
        self.element_ids = []
        self.elements = {}
        self.ecc_elements = {}
        self.ecc_coordinates = []
        self.ecc_element_nodes = {}
        self.member = None
        # start node, FEMNode object
        self.n1 = None
        # end node, FEMNode object
        self.n2 = None
        self.nodes = {}
        self.added_coordinates = []
        self.nodal_coordinates = []
        self.loc = []
        self.__coordinates = [[round(c, PREC) for c in coords] for coords in coordinates]
        # self.coordinates = coordinates
        self.cross_section = None
        self.__profile = profile
        self.material = material
        self.profile = profile
        # self.profile_idx = PROFILES.index(self.profile)
        # self.length = self.calc_length()
        self.mem_id = mem_id
        self.mtype = mtype
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.num_elements = num_elements
        #self.loads = {}
        self.loads = []
        self.alpha1 = 25
        self.alpha2 = 25
        self.__Sj1 = None
        self.__Sj2 = None
        self.Sj1 = Sj1
        self.Sj2 = Sj2
        self.ned = 0
        self.ved = 0
        self.med = 0
        self.reverse = reverse
        self.has_load = False        
        self.is_designed = False
        self.is_strong_enough = False
        #self.r = np.ones(7)
        self.r = {}
        self.self_weight = False
        self.is_generated = False
        # Lineload values
        self.q0 = 0
        self.q1 = 0
        # Create SteelMember object
        self.member = SteelMember(self.cross_section, self.length,
                                  Lcr=[1.0, 1.0], mtype=self.mtype)

        self.members = [self.member]
        self.hinges = [False,False]

    def copy(self):
        return type(self)(self.coordinates)

    def __repr__(self):
        return f"{type(self).__name__}({self.coordinates})"

    @property
    def Y(self):
        """ Returns y coordinates of the member ends """
        (x1, y1), (x2, y2) = self.coordinates
        return np.asarray([y1, y2])

    @Y.setter
    def Y(self, val):
        """ Sets y coordinates of the member ends """
        (x1, y1), (x2, y2) = self.coordinates
        self.coordinates = [[x1, val], [x2, val]]

    @property
    def start_coord(self):
        return self.coordinates[0]

    @property
    def end_coord(self):
        return self.coordinates[1]

    def to_global(self, loc):
        """
        Returns local coordinate in global coordinates
        :param loc: local coordinate
        :return: global coordinate
        """
        Lx = self.length * loc
        return self.coordinates[0] + Lx * self.unit

    @property
    def perpendicular(self):
        """
        Reuturns vector perpendicular to member
        :return:
        """
        unit = self.unit
        if self.reverse:
            return np.array([-unit[1], unit[0]])

        return np.array([unit[1], -unit[0]])

    @property
    def unit(self):
        """
        Returns direction vector of the member.
        :return:
        """
        start, end = np.asarray(self.coordinates)
        unit = np.array(end - start) / self.length
        return unit

    @property
    def A(self):
        return self.cross_section.A

    @A.setter
    def A(self, val):
        self.cross_section.A = val

    @property
    def Iy(self):
        return self.cross_section.Iy

    @Iy.setter
    def Iy(self, val):
        self.cross_section.Iy = val

    @property
    def E(self):
        return self.material.young

    @E.setter
    def E(self, val):

        self.material.young = val

    @property
    def fy(self):
        return self.material.fy

    @fy.setter
    def fy(self,val):
        self.material.fy = val

    @property
    def nu(self):
        return self.material.nu

    @property
    def rho(self):
        return self.material.density

    @rho.setter
    def rho(self, val):
        self.material.density = val

    @property
    def angle(self):
        """
        Returns member's angle in radians
        """
        start_node, end_node = self.coordinates
        x0, y0 = start_node
        x1, y1 = end_node
        if abs(x1 - x0) < 1e-6:
            # 90 degrees
            angle = 0.5*np.pi # np.radians(90)
        else:
            angle = np.arctan((y1 - y0) / (x1 - x0))
            if angle < 0:
                return np.pi + angle #np.radians(180) + angle
        return angle

    @property
    def coordinates(self):
        if self.n1 and self.n2:
            self.__coordinates = [list(self.n1.coord), list(self.n2.coord)]
            # self.calc_nodal_coordinates(self.num_elements)
        return self.__coordinates

    @coordinates.setter
    def coordinates(self, val):
        if not self.n1:
            self.__coordinates = val
            self.calc_nodal_coordinates()
        else:
            self.__coordinates = val
            x1, y1 = val[0]
            x2, y2 = val[1]
            self.n1.coord = np.array([x1, y1])
            self.n2.coord = np.array([x2, y2])
            for mem in self.n1.parents:
                mem.calc_nodal_coordinates()
            for mem in self.n2.parents:
                mem.calc_nodal_coordinates()

    @property
    def dx(self):
        vals = {}
        for load_id, displacements in self.nodal_displacements.items():
            dx_max = 0
            for node_id, disp_vals in displacements.items():
                dx, dy, rz = disp_vals
                if abs(dx) > abs(dx_max):
                    vals[load_id] = dx
                    dx_max = dx
            if load_id not in vals:
                vals[load_id] = 0

        return vals

    @property
    def dy(self):
        vals = {}
        for load_id, displacements in self.nodal_displacements.items():
            dy_max = 0
            for node_id, disp_vals in displacements.items():
                dx, dy, rz = disp_vals
                if abs(dy) > abs(dy_max):
                    vals[load_id] = dy
                    dy_max = dy
            if load_id not in vals:
                vals[load_id] = 0

        return vals
    
    @property
    def load_ids(self):
        """ Returns a list of load id's """
        n = list(self.elements.keys())[0]
        return [val for val in self.elements[n].fint.keys()]


    @property
    def MRd(self):
        """ Returns bending resistance
            Units: Nmm
        """
        return self.cross_section.MRd

    @property
    def MbRd(self):
        """ Returns Lateral-torsional buckling resistance
            Units: Nmm
        """
        return self.member.MbRd

    @property
    def NRd(self):
        """ Returns axial force resistance
            Units: N
        """
        return self.cross_section.NRd

    @property
    def NbRd(self):
        """ Returns buckling resistance
            Units: N
        """
        return self.member.NbRd

    @property
    def VRd(self):
        """ Returns shear force resistance
            Units: N
        """
        return self.cross_section.VRd

    @property
    def weight(self):
        """
        Function that calculates member's weight
        :return weight: float, weight of the member
        """

        weight = self.A * self.length * self.rho
        return weight

    @property
    def h(self):
        """
        Property, returns member's cross-section's height
        """
        try:
            return self.cross_section.h
        except AttributeError:
            return self.cross_section.H

    @h.setter
    def h(self, val):
        """
        Sets cross-section's height to given value.
        Changes profile and sets new cross-sectional properties to member's elements
        """

        splitted = self.profile.split()
        if len(splitted) == 2:
            profile_type = self.profile.split()[0].upper()
            if profile_type == "IPE":
                self.profile = profile_type + " " + str(val)
            elif profile_type == "SHS":
                val = max(val, 25)
                self.profile = profile_type + " " + str(val) + 'X' + str(
                    val) + 'X' + str(max(round(val / 20), 2))
        elif len(splitted) == 3:
            self.profile = splitted[0] + " " + str(val) + " " + splitted[2]

    @property
    def material(self):
        """
        Returns member's material e.g. S355
        """
        return self.__material

    @material.setter
    def material(self, val):
        """
        Changes member's matrial and gets new material properties
        from tables_and_tuples.mat -dict
        """
        self.__material = MATERIALS[val]

        # self.fy = MATERIALS[self.material]['f_y']
        # self.E = MATERIALS[self.material]['E']
        # self.fu = MATERIALS[self.material]['f_u']
        # self.nu = MATERIALS[self.material]['v']
        # self.rho = MATERIALS[self.material]['rho']

    @property
    def profile(self):
        """
        Returns member's profile e.g. IPE 200, HE 120 A
        """
        # if self.cross_section is not None:
        #     return str(self.cross_section)
        return self.__profile

    @profile.setter
    def profile(self, val):
        """
        Changes member's profile to given value.
        Calculates new cross-sectional properties and sets these new values
        to member's elements.
        
        :param val: string, profile name e.g. 'IPE 100'
        """

        #print(val, isinstance(val, SteelSection))

        if isinstance(val, SteelSection):
            self.cross_section = val
            self.__profile = str(val)
        elif isinstance(val, list) and len(val) == 5:
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
            self.__profile = val.upper()
            splitted_val = self.profile.split(" ")
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

        # Change member's elements' properties
        if len(self.elements) > 0:
            for element in self.elements.values():
                element.section = self.cross_section
                #element.section.A = self.cross_section.A
                #element.section.Iy = self.cross_section.I[0]
        # Change steel_member objects properties
        if self.member:
            self.member.profile = self.cross_section

    @property
    def Sj1(self):
        return self.__Sj1

    @Sj1.setter
    def Sj1(self, val):

        if val < 1e-4:
            val = 1e-4
        elif val > 1e20:
            val = 1e20
        if len(self.element_ids) > 0:
            idx = self.element_ids[0]
            ele = self.elements[idx]
            ele.rot_stiff[0] = val

        self.__Sj1 = val

    @property
    def Sj2(self):
        return self.__Sj2

    @Sj2.setter
    def Sj2(self, val):

        if val < 1e-4:
            val = 1e-4
        elif val > 1e20:
            val = 1e20
        if len(self.element_ids) > 0:
            idx = self.element_ids[-1]
            ele = self.elements[idx]
            ele.rot_stiff[1] = val

        self.__Sj2 = val

    def shape(self, x):
        """
        Returns chord's y value corresponding to x
        """
        start_node, end_node = self.coordinates
        x0, y0 = start_node
        Lx = end_node[0] - start_node[0]
        Ly = end_node[1] - start_node[1]
        try:
            k = Ly / Lx
        except:
            k = 0

        return k * (x - x0) + y0

    @property
    def length(self):
        start_node, end_node = np.asarray(self.coordinates)
        return np.linalg.norm(end_node - start_node)

    def alpha_to_Sj(self):
        """ Change alpha value to Sj
        """
        # SFS-EN 1993-1-8 page 60, hinge 0.5 < alpha < 25 rigid
        self.Sj1 = self.alpha1 * self.E * self.I_y * 1e-12 / self.length
        self.Sj2 = self.alpha2 * self.E * self.I_y * 1e-12 / self.length

    def Sj_to_alpha(self):
        """ Change Sj value to alpha
        """
        # SFS-EN 1993-1-8 page 60, hinge 0.5 < alpha < 25 rigid
        self.alpha1 = self.Sj1 * self.length / (self.E * self.I_y * 1e-12)
        self.alpha2 = self.Sj2 * self.length / (self.E * self.I_y * 1e-12)

    def generate(self, fem_model):
        """ Generates the member
        """
        # If the frame hasn't been created, create members and calculate nodal
        # coordinates, otherwise initialize frame.

        if not self.is_generated:
            self.add_nodes(fem_model)
            self.add_material(fem_model)
            self.add_section(fem_model)
            self.generate_elements(fem_model)
            self.is_generated = True

    def calc_nodal_coordinates(self, num_elements=0):
        """
            Calculates node locations along member.
            Coordinates used are global coordinates.
            Adds locations to a list where they can be accessed later.
            
            :param num_elements: int, number of elements to be created
        """
        # if num_elements != 0:
        #     self.num_elements = num_elements
        # self.nodal_coordinates = []
        # self.nodal_coordinates.extend(self.added_coordinates)
        # start_node, end_node = self.coordinates
        # x0, y0 = start_node
        # Lx = end_node[0] - start_node[0]
        # Ly = end_node[1] - start_node[1]
        # dx = Lx / self.num_elements
        # dy = Ly / self.num_elements
        # for i in range(self.num_elements):
        #     x = i * dx + x0
        #     y = i * dy + y0
        #     loc = i / self.num_elements
        #     node = [x, y]
        #     if node not in self.nodal_coordinates:
        #         self.nodal_coordinates.append(node)
        #     if loc not in self.loc:
        #         self.loc.append(loc)
        # if [end_node[0], end_node[1]] not in self.nodal_coordinates:
        #     self.nodal_coordinates.append([end_node[0], end_node[1]])
        # if 1 not in self.loc:
        #     self.loc.append(1)
        #
        # if self.is_generated:
        #     for j, node in enumerate(self.nodes.values()):
        #         # print(self.nodal_coordinates)
        #         # print(j, len(self.nodal_coordinates))
        #         node.coord = np.array(self.nodal_coordinates[j])

        # clear list of nodal coordinates
        self.nodal_coordinates.clear()
        """ add optional coordinates that are derived from point loads and supports
            etc. that are located between end points of the member.
        """
        self.nodal_coordinates.extend(self.added_coordinates)
        
        """ By default, 5 elements are created, but this can be changed by
            'num_elements'
        """
        if num_elements > 0:
            self.num_elements = num_elements

                
        dloc = 1 / self.num_elements # step length in local coordinate system
        #for loc in np.arange(0, 1 + dloc, dloc):
        for loc in np.linspace(0, 1, self.num_elements + 1):
            
            loc = round(loc, PREC) # round local coordinate to PREC decimals
            
            """ WHY IS THIS if LOOP HERE? """
            if loc not in self.loc:
                self.loc.append(loc)

            coord = list(self.to_global(loc))
            self.nodal_coordinates.append(coord)

        """ Sort coordinates according to their distance from the
            first end node
        """
        c0, c1 = self.coordinates
        self.nodal_coordinates.sort(
            key=lambda c: np.linalg.norm(np.asarray(c0) - np.asarray(c)))

        if self.is_generated:
            for j, node in enumerate(self.nodes.values()):
                node.coord = np.array(self.nodal_coordinates[j])

    def add_node_coord(self, coord):
        """
        Adds coordinate to member's coordinates and
        calculates new coordinate's local location
        If coordinate is not between member's coordinates, reject it
        :param coord: array of two float values, node's coordinates
        """
        # if coord not in self.nodal_coordinates and self.point_intersection(coord):
        coord = list(coord)
        if coord not in self.added_coordinates and self.point_intersection(
                coord):
            self.added_coordinates.append(coord)
            self.nodal_coordinates.append(coord)
            c0, c1 = self.coordinates
            self.nodal_coordinates.sort(
                key=lambda c: np.linalg.norm(np.asarray(c0) - np.asarray(c)))
            start_node, end_node = self.coordinates
            x, y = coord
            x1, y1 = start_node
            dx = x - x1
            dy = y - y1
            dz = math.sqrt((dx ** 2 + dy ** 2))
            z = dz / self.length

            if x < x1:
                z *= -1
            z = round(z, PREC)
            if z not in self.loc:
                self.loc.append(z)
                self.loc.sort()

    def add_nodes(self, fem_model):
        """ Creates nodes to previously calculated locations
        """
        idx = None
        for coordinate in self.nodal_coordinates:
            x, y = coordinate
            if coordinate in fem_model.nodal_coords:
                idx = fem_model.nodal_coords.index(coordinate)
                self.nodes[idx] = fem_model.nodes[idx]
                self.nodes[idx].parents.append(self)
            else:
                idx = fem_model.nnodes()
                fem_model.add_node(x, y)
                self.nodes[idx] = fem_model.nodes[idx]
                self.nodes[idx].parents.append(self)
            if not self.n1:
                self.n1 = fem_model.nodes[idx]
                # Add last node
        if not self.n2 and idx:
            self.n2 = fem_model.nodes[idx]

    def generate_elements(self, fem_model):
        """ Generates elements between nodes
            For beam member's first and last element, creates semi-rigid end
            elements. Elements are added to a list
            :param fem_model: FrameFEM -object
        """
        index = fem_model.nels()
        node_ids = list(self.nodes.keys())
        # EBBeam -elements
        if self.mtype == "column":

            for i in range(len(self.nodes) - 1):
                n1 = node_ids[i]
                n2 = node_ids[i + 1]

                if isinstance(self.cross_section, TaperedSection):
                    self.elements[index] = \
                        EBBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                               self.cross_section.small_elements[i],
                               self.material)
                else:
                    self.elements[index] = \
                        EBBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                               self.cross_section,
                               self.material)

                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
                index += 1
        elif self.mtype == 'bar':
            """ Simple bar element that carries only axial forces """
            n1 = node_ids[0]
            n2 = node_ids[1]
            self.elements[index] = \
                Rod(fem_model.nodes[n1], fem_model.nodes[n2],
                        self.cross_section,
                        self.material)
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
        elif self.mtype == "beam":
            if len(self.nodes) == 2:
                n1 = node_ids[0]
                n2 = node_ids[1]
                
                """ EBSemiRigidBeam is a special element that has
                    rotational stiffness at the ends.
                """
                if isinstance(self.cross_section, TaperedSection):
                    self.elements[index] = \
                        EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                                        self.cross_section.small_elements[0],
                                        self.material, rot_stiff=[np.inf, self.Sj2])
                else:
                    self.elements[index] = \
                        EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                                        self.cross_section,
                                        self.material,
                                        rot_stiff=[self.Sj1, self.Sj2])
                
                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
            else:
                for i in range(len(self.nodes) - 1):
                    n1 = node_ids[i]
                    n2 = node_ids[i + 1]

                    # EBSemiRigid -element, first element
                    if i == 0:
                        if isinstance(self.cross_section, TaperedSection):
                            self.elements[index] = \
                                EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                                       self.cross_section.small_elements[i],
                                       self.material, rot_stiff=[self.Sj1, np.inf])
                        else:
                            self.elements[index] = \
                                EBSemiRigidBeam(fem_model.nodes[n1],
                                            fem_model.nodes[n2],
                                            self.cross_section,
                                            self.material,
                                            rot_stiff=[self.Sj1, np.inf])
                    # EBSemiRigid -element, last element
                    elif i == (len(self.nodes) - 2):
                        if isinstance(self.cross_section, TaperedSection):
                            self.elements[index] = \
                                EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                                       self.cross_section.small_elements[i],
                                       self.material, rot_stiff=[np.inf, self.Sj2])
                        else:
                            self.elements[index] = \
                                EBSemiRigidBeam(fem_model.nodes[n1],
                                            fem_model.nodes[n2],
                                            self.cross_section,
                                            self.material,
                                            rot_stiff=[np.inf, self.Sj2])
                    # EBBeam -elements
                    else:
                        if isinstance(self.cross_section, TaperedSection):
                            self.elements[index] = \
                                EBBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                                       self.cross_section.small_elements[i],
                                       self.material)
                        else:
                            self.elements[index] = \
                                EBBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                                   self.cross_section,
                                   self.material)

                    fem_model.add_element(self.elements[index])
                    self.element_ids.append(index)
                    index += 1
        else:
            """ For other kinds of elements, use regular
                Euler-Bernoulli beam elements.
            """
            N = len(self.nodes)
            for i in range(N - 1):
                n1 = node_ids[i]
                n2 = node_ids[i + 1]

                if isinstance(self.cross_section, TaperedSection):
                    self.elements[index] = \
                        EBBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                               self.cross_section.small_elements[i],
                               self.material)
                else:
                    self.elements[index] = \
                        EBBeam(fem_model.nodes[n1], fem_model.nodes[n2],
                               self.cross_section,
                               self.material)

                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
                
                if i == 0 and self.hinges[0]:
                    """ There is a hinge in the first end """
                    fem_model.add_release(index,[2])
                
                if i == N-2 and self.hinges[1]:
                    """ There is a hinge in the first end """
                    fem_model.add_release(index,[5])
                
                index += 1

    def generate_eccentricity_elements(self, fem_model):
        """ Generates eccentricity elements to connections
        """
        # Material id
        mat_id = len(fem_model.materials)
        # Add material
        # Material(E, nu, rho)
        fem_model.add_material(210e3, 0.3, 7850e-9)
        # Add section
        # BeamSection(A, Iy)
        # HE 1000 A
        # A = 34680
        # Iy = 5538000000
        sect_id = len(fem_model.sections)
        sect = HEA(1000)
        fem_model.add_section(sect)
        for coordinates in self.ecc_coordinates:
            for i, coord in enumerate(coordinates):
                coordinates[i] = [round(c, PREC) for c in coord]
            idx = fem_model.nels()
            # n1 is the node connected to the column
            # n2 is hinged connection to chord
            n1 = fem_model.nodes[fem_model.nodal_coords.index(coordinates[0])]
            n2 = fem_model.nodes[fem_model.nodal_coords.index(coordinates[1])]
            self.ecc_elements[idx] = EBBeam(n1, n2,
                                            fem_model.sections[sect_id],
                                            fem_model.materials[mat_id])
            # rot_stiff=[np.inf, 0])
            fem_model.add_element(self.ecc_elements[idx])

    def calc_nodal_forces(self,load_id=LoadIDs.ULS):
        """ Gets calculated nodal forces on member's elements and
            adds them to dict where node is the key and froces are value
            self.nodal_forces[node_id] = [Fx, Fy, Mz]
        """

        def absmax(vals):
            min_val = min(vals)
            max_val = max(vals)
            abs_max = max(abs(min_val), abs(max_val))
            if abs_max == abs(min_val):
                return min_val
            else:
                return max_val

        """ Store also maximum bending moment, shear force and axial force """
        max_med = 0
        max_ved = 0
        max_ned = 0

        """ i is node index """
        i = 0        
        """ node_ids is a list of node indices for the member """
        node_ids = list(self.nodes.keys())
        node_forces = {}
        for element in self.elements.values():
            node_id = node_ids[i]
            #axial_force = absmax(element.axial_force)
            axial_force = absmax(element.fint[load_id]['fx'])
            if abs(axial_force) > abs(max_ned):
                max_ned = axial_force

            #shear_force = absmax(element.shear_force)
            shear_force = absmax(element.fint[load_id]['fz'])
            if abs(shear_force) > abs(max_ved):
                max_ved = shear_force

            #bending_moment = absmax(element.bending_moment)
            bending_moment = absmax(element.fint[load_id]['my'])
            if abs(bending_moment) > abs(max_med):
                max_med = bending_moment

            """
            self.nodal_forces[node_id] = [axial_force,
                                          shear_force,
                                          bending_moment]
            """
            node_forces[node_id] = [axial_force,
                                    shear_force,
                                    bending_moment]
            i += 1

            if i == len(self.elements):

                node_id = node_ids[i]
                #axial_force = absmax(element.axial_force)
                axial_force = absmax(element.fint[load_id]['fx'])
                if abs(axial_force) > abs(max_ned):
                    max_ned = axial_force

                #shear_force = absmax(element.shear_force)
                shear_force = absmax(element.fint[load_id]['fz'])
                if abs(shear_force) > abs(max_ved):
                    max_ved = shear_force

                #bending_moment = absmax(element.bending_moment)
                bending_moment = absmax(element.fint[load_id]['my'])
                if abs(bending_moment) > abs(max_med):
                    max_med = bending_moment

                """
                self.nodal_forces[node_id] = [axial_force,
                                              shear_force,
                                              bending_moment]
                """
                node_forces[node_id] = [axial_force,
                                    shear_force,
                                    bending_moment]

        self.nodal_forces[load_id] = node_forces
        self.med = max_med
        self.ved = max_ved
        self.ned = max_ned

    def calc_nodal_displacements(self, fem_model, lcase=LoadIDs.ULS):
        """ Calculates nodal displacements and saves them to a dict
            :param fem_model: FrameFEM -object
        """
        #for node in self.nodes:
        #    self.nodal_displacements[node] = fem_model.nodes[node].u[lcase]
            #self.nodal_displacements[lcase] = fem_model.nodes[node].u[lcase]
        
        U = {}
        
        for n, node in self.nodes.items():
            #print(n,node)
            U[n] = node.u[lcase]
        
        self.nodal_displacements[lcase] = U

    def add_self_weight(self, frame):
        """ Adds self-weight to the member's loads
        """
        if not self.self_weight:
            self.self_weight = True
            load_id = "self_weight"
            # self.weight is kg's, multiplier changes it to N's
            multiplier = 10
            value = -1 * multiplier * self.weight / self.length            
            frame.add(LineLoad(self, [value, value], 'y'))

    def remove_self_weight(self):
        """ Removes self-weight from loads
        """
        if self.self_weight:
            self.self_weight = False
            del (self.loads["self_weight"])

        if len(self.loads) == 0:
            self.has_load = False

    def add_material(self, fem_model):
        """ Adds member's material information to calculation model
            Young's modulus MPa
            nu is dimensionless
            Density kg/m^3
        """
        fem_model.add_material(self.E, self.nu, self.rho)

    def add_section(self, fem_model):
        """ Units: mm A and I_y are in mm**2 and mm**4, respectively"""
        s = fem.BeamSection(self.cross_section.A, self.cross_section.I[0])
        fem_model.add_section(s)

    @property
    def locs(self):
        """ Locations of member nodes in local coordinates 
            loc = 0 is for node 'start' and loc = 1 for node 'end'
        """
        locs = []
        start_coord, end_coord = np.asarray(self.coordinates)
        for node in self.nodes.values():
            ncoord = np.asarray(node.coord)            
            loc = np.linalg.norm(start_coord - ncoord) / self.length
            locs.append(round(loc, PREC))
        return locs


    def assign_forces(self,load_id=LoadIDs.ULS):
        """
        Assigns forces to SteelMember -object.
        Creates new sections if none are created,
        otherwise updates old values
        """
        self.member.length = self.length
        self.member.profile = self.cross_section
        self.member.clear_sections()
        for i, element in enumerate(self.elements.values()):
            #M1, M2 = element.bending_moment
            #N1, N2 = element.axial_force
            #V1, V2 = element.shear_force
            M1, M2 = element.fint[load_id]['my']
            N1, N2 = element.fint[load_id]['fx']
            V1, V2 = element.fint[load_id]['fz']
            loc = self.locs[i]
            self.member.add_section(ned=N1, vzed=V1,
                                    myed=M1, loc=loc)
        loc = self.locs[i + 1]
        self.member.add_section(ned=N2, vzed=V2,
                                myed=M2, loc=loc)

        #
        # if not self.steel_member.nsect():
        #     for i, node in enumerate(self.nodal_forces.keys()):
        #         forces = self.nodal_forces[node]
        #         ned = forces[0]
        #         ved = forces[1]
        #         med = forces[2]
        #         loc = self.locs[i]
        #         self.steel_member.add_section(ned=ned, vzed=ved,
        #                                       myed=med, loc=loc)
        #
        # else:
        #     for i, node in enumerate(self.nodal_forces.keys()):
        #         Fx, Fy, Mz = self.nodal_forces[node]
        #         self.steel_member.ned[i] = Fx
        #         self.steel_member.vzed[i] = Fy
        #         self.steel_member.myed[i] = Mz


    def check_cross_section(self):
        """ Checks if cross-section is strong enough.
            Gives boolean value for self.is_strong_enough
            Saves list of stress ratios to self.r
            
            TODO!
            This has to be modified for multiple load cases.
        """

        self.r = np.zeros_like(self.r)
        # Cross-sectional stress ratios in member's nodes' locations
        self.r[:4] = self.member.check_sections()
        # Buckling about y - and z - axis
        buckling_r = self.member.check_buckling()
        self.r[4] = buckling_r[0]
        self.r[5] = buckling_r[1]
        # Lateral-Torsional buckling
        self.r[6] = self.member.check_LT_buckling()
        # Set cross-section's maximum forces for getting results later
        self.cross_section.Ned = self.ned
        self.cross_section.Ved = self.ved
        self.cross_section.Med = self.med

        if max(self.r) > 1.0:
            self.is_strong_enough = False
        else:
            self.is_strong_enough = True

    def design(self,lcase=LoadIDs.ULS):
        """ Designs the member (check resistance) 
            This has to be modified for multiple load cases.
        """
        
        """ Here, the internal forces of the member at sections must
            somehow be inserted. Note that SteelMember does not yet
            provide several load cases.
        """

        rmax = self.member.design()
        
        if rmax > 1.0:
            self.is_strong_enough = False
        else:
            self.is_strong_enough = True
            
        self.r[lcase] = rmax
        
        #self.r = max(max(r1), r2)

    def optimum_design(self, prof_type='CURRENT'):
        """ Goes through all profiles in list and stops iterating when 
            cross-section can bear it's loads.
            If member is previously designed, iterating starts from 3 profiles
            before currrent profile.
        """
        
        prof_type = prof_type.upper()
        
        if prof_type == "CURRENT":
            """ Use the cross section family of the current profile """
            if isinstance(self.cross_section,IPE):
                prof_type = "IPE"
            elif isinstance(self.cross_section,HEA):
                prof_type = "H"
            elif isinstance(self.cross_section,HEB):
                prof_type = "H"
            elif isinstance(self.cross_section,RHS):
                prof_type = "RHS"
            elif isinstance(self.cross_section,SHS):
                prof_type = "SHS"
            elif isinstance(self.cross_section,CHS):
                prof_type = "CHS"
            else:
                prof_type = "IPE"
        
        if prof_type == "IPE":            
            profiles = ipe_profiles.keys()
        elif prof_type == "H":
            profiles = h_profiles.keys()
        elif prof_type == "RHS":
            profiles = rhs_profiles.keys()
        elif prof_type == "CHS":
            profiles = chs_profiles.keys()
        elif prof_type == "SHS":
            profiles = shs_profiles.keys()
        else:
            raise ValueError(f'"{prof_type}" is not a valid profile type!')

        initial_profile = self.profile

        for profile in profiles:
            self.profile = profile
            
            for load_id in self.load_ids:            
                self.design(load_id)
                            
            if all(r<=1.0 for r in self.r.values()):# <= 1.0:
                break
            """
            self.profile = profile
            self.check_cross_section()
            for smem in self.members:
                r1 = smem.check_sections()

                r2 = max(smem.check_beamcolumn())
                self.r = max(max(r1), r2)
            if self.r <= 1.0:
                break
            """
        
        """ If the profile changed during iteration, return 'True'. """
        if self.profile != initial_profile:
            return True
        else:
            return False
        
    def line_intersection(self, coordinates):
        """
        Calculates coordinate where two members intersect
        :param coordinates: array of two arrays of two float values, [[x1,y1],[x2,y2]]
        :return [px, py]: array of two float values, coordinates for intersection
        source: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
        """
        if isinstance(coordinates, FrameMember):
            coordinates = coordinates.coordinates
        start_node, end_node = self.coordinates
        x1, y1 = start_node
        x2, y2 = end_node
        x3, y3 = coordinates[0]
        x4, y4 = coordinates[1]

        if ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)) == 0:
            return None
        else:
            px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (
                        x3 * y4 - y3 * x4)) / (
                         (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
            py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (
                        x3 * y4 - y3 * x4)) / (
                         (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))

            return [round(px, 5), round(py, 5)]

    def point_intersection(self, coordinate):
        """
        Calculates if member and point intersect
        :param coordinate: array of two float values
        :return bool: true if point and member intersect
        """

        start_node, end_node = self.coordinates
        x1, y1 = start_node
        x2, y2 = end_node
        x, y = coordinate
        # Coordinate is between member's coordinates
        if x1 <= x <= x2 and y1 <= y <= y2 or 'chord' in self.mtype:
            if x2 - x1 == 0:
                y1, y2 = sorted([y1, y2])
                return x == x1 and y1 <= y <= y2
            else:
                k = (y2 - y1) / (x2 - x1)
                b = y1 - k * x1
                return math.isclose(y, k * x + b,
                                    rel_tol=1e-3)  # y == k * x + b

        return False

    def plot(self, print_text=True, c='k', axes=None):
        
        if self.active:
            if axes is None:
                fig, ax = plt.subplots(1)
            else:
                ax = axes
    
            X = self.coordinates
            if c:
                if self.is_strong_enough:
                    color = 'green'
                else:
                    color = 'red'
            else:
                color = 'k'
            # Plot members
            ax.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color)
            # Plot text
    
            if self.mtype == 'beam':
                horzalign = 'center'
                vertalign = 'bottom'
    
            elif self.mtype == 'column':
                horzalign = 'right'
                vertalign = 'center'
    
            else:
                horzalign = 'center'
                vertalign = 'center'
    
            x, y = self.to_global(0.3) - self.perpendicular * 50
            rot = np.degrees(self.angle)
    
            if print_text:
                ax.text(x, y, str(self.mem_id) + ": " + str(self.cross_section),
                         rotation=rot, horizontalalignment=horzalign,
                         verticalalignment=vertalign)
            else:
                ax.text(x, y, str(self.mem_id),
                         rotation=rot, horizontalalignment=horzalign,
                         verticalalignment=vertalign)

    def plot_results(self):
        """ Plots utilization ratios for different resistances """
        
        UN, UV, UM, UMN = self.cross_section.section_resistance(
            return_list=True)
        B = max(self.member.check_buckling())
        LT = self.member.check_LT_buckling()

        x = np.arange(6)
        Y = [UN, UV, UM, UMN, B, LT]
        c = []
        for val in Y:
            if val < 0.9:
                c.append('lightgreen')
            elif val > 1.0:
                c.append('r')
            else:
                c.append('ly')
        plt.bar(x, Y, color=c, width=0.75)
        X = ['N', 'V', 'M', 'MN', 'B', 'LT']
        plt.xticks(x, X)
        plt.title(f'Member {self.mem_id} stress ratios')
        line_x = np.arange(-0.5, len(Y))
        line_y = [1] * len(line_x)
        plt.plot(line_x, line_y, '--', c='k')
        plt.show()

    def bmd(self, scale=1, load_id = LoadIDs.ULS):
        """ Plots bending moment diagram """

        # Scales Nmm to kNm
        unit_scaler = 1e-6

        # Member's position vector
        v = self.unit
        # Vector perpendicular to v
        u = -self.perpendicular
        X = []
        Y = []
        start_coord, end_coord = self.coordinates
        x01, y01 = start_coord
        x02, y02 = end_coord
        X.append(x01)
        Y.append(y01)
        self.calc_nodal_forces(load_id=load_id)

        moment_values = [x[2] for x in self.nodal_forces[load_id].values()]

        for elem in self.elements.values():
            x0, y0 = elem.nodes[0].coord
            #bending_moment = elem.bending_moment[0]
            bending_moment = elem.fint[load_id]['my'][0]
            val = bending_moment * unit_scaler * scale
            x, y = np.array([x0, y0]) - val * u
            X.append(x)
            Y.append(y)
            horzalign = 'center'
            vertalign = 'center'
            if bending_moment == max(moment_values) or bending_moment == min(
                    moment_values):
                plt.text(x, y, f'{bending_moment * unit_scaler:.2f} kNm',
                         horizontalalignment=horzalign)
        x0, y0 = elem.nodes[1].coord
        #bending_moment = elem.bending_moment[1]
        bending_moment = elem.fint[load_id]['my'][1]
        val = bending_moment * unit_scaler * scale
        x, y = np.array([x0, y0]) - val * u
        X.append(x)
        Y.append(y)
        X.append(x02)
        Y.append(y02)
        plt.text(x, y, f'{bending_moment * unit_scaler:.2f} kNm',
                 horizontalalignment=horzalign)
        plt.plot(X, Y, color='gray')

    def sfd(self, scale):
        """
        Plots shear force diagram
        :param scale: float, scaling factor
        """

        X = []
        Y = []
        self.calc_nodal_forces()
        X.append(self.nodal_coordinates[0][0])
        Y.append(self.nodal_coordinates[0][1])
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            if self.mtype == "beam":
                x0 = self.nodal_coordinates[i][0]
                y0 = self.nodal_coordinates[i][1]
                y1 = self.nodal_forces[node][1] / (1000 / scale)
                x = x0
                y = y0 - y1
                X.append(x)
                Y.append(y)
            elif self.mtype == "column":
                x0 = self.nodal_coordinates[i][0]
                y0 = self.nodal_coordinates[i][1]
                x1 = self.nodal_forces[node][1] / (1000 / scale)
                x = x0 + x1
                y = y0
                X.append(x)
                Y.append(y)
        X.append(self.nodal_coordinates[i][0])
        Y.append(self.nodal_coordinates[i][1])
        plt.plot(X, Y, color='gray')

    def plot_normal_force(self):

        # Scales N to kN
        unit_scaler = 1e-3

        X = self.coordinates
        if self.ned < 0:
            r = min(1, max(self.r[4:5]))
            alpha = max(0.1, r)
            color = (1, 0.5 - r / 2, 0, alpha)
        else:
            r = min(1, self.r[0])
            alpha = max(0.1, r)
            color = (0, 0.5 - r / 2, 1, alpha)
        # Plot members
        plt.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color=color,
                 linewidth=2)
        # Plot text
        # Calculate text location
        delta_y = X[1][1] - X[0][1]
        delta_x = X[1][0] - X[0][0]

        if delta_x == 0:
            rot = 90
        elif delta_y == 0:
            rot = 0
        else:
            rot = math.degrees(math.atan(delta_y / delta_x))

        x = (X[0][0] + X[1][0]) / 2
        y = (X[1][1] + X[0][1]) / 2
        horzalign = 'center'
        vertalign = 'center'

        if self.ned > 0:
            NED = max([f[0] for f in self.nodal_forces.values()])
        else:
            NED = min([f[0] for f in self.nodal_forces.values()])

        plt.text(x, y,
                 f'{NED * unit_scaler:.2f}  kN,\nr: {r * 100:.2f} %',
                 rotation=rot, horizontalalignment=horzalign,
                 verticalalignment=vertalign)

    def lineload_elements(self, coords):

        # Returns the elements that are under line load
        #:param coords: 2D-array of float, line load's start and end coordinates
        #:return element_ids: array of int, elements id's

        start, end = coords
        start = [round(c, PREC) for c in start]
        end = [round(c, PREC) for c in end]
        start_idx = self.nodal_coordinates.index(start)
        end_idx = self.nodal_coordinates.index(end)
        element_ids = self.element_ids[start_idx:end_idx]
        return element_ids

    def round_coordinates(self, prec=PREC):
        """ Rounds the nodal coordinates with the precision of 'prec'
            (number of decimals)
        """

        for i, coord in enumerate(self.nodal_coordinates):
            self.nodal_coordinates[i] = [round(c, prec) for c in coord]

    def add_hinge(self, loc=0):
        """ Adds a hinge to the member
        input:
            loc .. location along the member (should be 0 or 1)
        """
        if loc == 0:
            self.hinges[0] = True
        else:
            self.hinges[1] = True


class SteelFrameMember(FrameMember):
    def __init__(self, coordinates, mem_id='', profile='IPE 100', material="S355", num_elements=4, **kwargs):
        super().__init__(coordinates, mem_id, profile, material, num_elements, **kwargs)
        self.member = SteelMember(self.cross_section, self.length,
                                        Lcr=[1.0, 1.0], mtype=self.mtype)

        self.members = [self.member]


class TimberFrameMember(FrameMember):
    def __init__(self, coordinates: list_out[list_out[int]], profile: (str, TimberSection), mem_id: str='', num_elements: int=5,
                 varnished: bool=False, ldc: str='inst', sc: int=1, Sj1: float=np.inf, Sj2: float=np.inf, mtype: str="", lateral_support_y: list_out[float]=None,
                 lateral_support_z: list_out[float]=None, edge_load: str='compression', beta: (float, list_out[list_out[float]])=None, k_vol: float=None):
        """

        @param coordinates:
        @param profile: TimberSection
        @param mem_id:
        @param num_elements: Elementtien määrä
        @param varnished: lakkaus
        @param ldc: Load duration class. options: 'perm' = permanent, 'lt' = long-term, 'mt' = medium-term,
                    'st' = short-term, 'inst' = instantaneous
        @param sc: Service class.
        @param Sj1: Sauvan alkupään jäykkyys
        @param Sj2: Sauvan loppupään jäykkyys
        @param mtype: member type: 'column', 'beam'
        @param lateral_support_y: locations of the supports in global y direction (This has nothing to do with FEM)
        @param lateral_support_z: locations of the supports in z direction (This has nothing to do with FEM)
        @param edge_load: 'compression' or 'tension'
        @param beta: Eulerin nurjahduskerroin
        @param k_vol: pakotettu k_vol
        """
        self.__profile = None
        if isinstance(profile, TimberSection):
            super().__init__(coordinates, mem_id, profile, profile.material, num_elements, Sj1, Sj2, mtype)
        else:
            super().__init__(coordinates, mem_id, profile, T.GL30c, num_elements, Sj1, Sj2, mtype)

        if isinstance(self.cross_section, TaperedSection):
            if isinstance(self.cross_section, (KaarevaHarjaPalkki, KaarevaPalkki)):
                self.cross_section.parent = self
                self.cross_section.update_circles()
            elif isinstance(self.cross_section, MahaPalkki):
                self.cross_section.parent = self
            self.cross_section.num_elements = num_elements
            self.cross_section.create_elements()
        if lateral_support_z is None:
            lateral_support_z = []
        if lateral_support_y is None:
            lateral_support_y = []
        self.member = TimberMember(self.cross_section, self.material, self.length,
                                   ldc, sc, mtype=mtype, varnished=varnished, Sj1=self.Sj1, Sj2=self.Sj2,
                                   nodes=self.nodes, parent=self, lateral_support_y=lateral_support_y,
                                   lateral_support_z=lateral_support_z, edge_load=edge_load, beta=beta, k_vol=k_vol)
        self.members = [self.member]
        self.cross_section.timber_member = self.member

    @property
    def R(self):
        return self.cross_section.R

    @R.setter
    def R(self, val):
        self.cross_section.R = val

    @property
    def profile(self):
        return str(self.__profile)

    @profile.setter
    def profile(self, val):
        if isinstance(val, TimberSection):
            self.__profile = val
            self.cross_section = val
        elif isinstance(val, str):
            splitted_val = val.split(" ")
            profile_type = splitted_val[0]
            vals = splitted_val[1].split('X')
            if self.__profile is None:
                if profile_type == 'TS':
                    b = float(vals[0])
                    h = float(vals[1])
                    self.__profile = TimberSection(b, h)
                elif profile_type == 'HP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    hap = float(vals[2])
                    self.__profile = HarjaPalkki(b, h0, hap)
                elif profile_type == 'PP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    hap = float(vals[2])
                    self.__profile = PulpettiPalkki(b, h0, hap)
                elif profile_type == 'MP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    alph = float(vals[2])
                    self.__profile = MahaPalkki(b, h0, alph)
                elif profile_type == 'KP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    alph = float(vals[2])
                    beta = float(vals[3])
                    rin = float(vals[4])
                    self.__profile = KaarevaPalkki(b, h0, alph, beta, rin)
                elif profile_type == 'KHP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    alph = float(vals[2])
                    beta = float(vals[3])
                    rin = float(vals[4])
                    self.__profile = KaarevaHarjaPalkki(b, h0, alph, beta, rin)
                self.cross_section = self.__profile
            else:
                if profile_type == 'TS':
                    b = float(vals[0])
                    h = float(vals[1])
                    self.cross_section.B = b
                    self.cross_section.H = h
                if profile_type == 'HP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    hap = float(vals[2])
                    self.cross_section.B = b
                    self.cross_section.H0 = h0
                    self.cross_section.Hap = hap
                elif profile_type == 'PP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    hap = float(vals[2])
                    self.cross_section.B = b
                    self.cross_section.H0 = h0
                    self.cross_section.Hap = hap
                elif profile_type == 'MP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    alph = float(vals[2])
                    self.cross_section.B = b
                    self.cross_section.H0 = h0
                    self.cross_section.Alpha = alph
                elif profile_type == 'KP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    alph = float(vals[2])
                    beta = float(vals[3])
                    rin = float(vals[4])
                    self.cross_section.B = b
                    self.cross_section.H0 = h0
                    self.cross_section.Alpha = alph
                    self.cross_section.Beta = beta
                    self.cross_section.R_in = rin
                elif profile_type == 'KHP':
                    b = float(vals[0])
                    h0 = float(vals[1])
                    alph = float(vals[2])
                    beta = float(vals[3])
                    rin = float(vals[4])
                    self.cross_section.B = b
                    self.cross_section.H0 = h0
                    self.cross_section.Alpha = alph
                    self.cross_section.Beta = beta
                    self.cross_section.R_in = rin

    @property
    def material(self):
        """
        Returns member's material e.g. GL30c
        """
        return self.cross_section.material

    @material.setter
    def material(self, val):
        """
        Changes member's matrial and gets new material properties
        """
        self.__material = val

    def add_material(self, fem_model):
        """ Adds member's material information to calculation model
            Young's modulus MPa
            TODO tarkista G parametri (sen suuntaus ja arvo)
            Density kg/m^3
        """
        fem_model.add_material(self.member.E0mean, None, self.member.material.rhomean, G=self.member.material.Gmean)

    def add_self_weight(self, frame, load_id=LoadIDs.ULS):
        """ Adds self-weight to the member's loads
        """
        if not self.self_weight:
            self.self_weight = True

            if isinstance(self.cross_section, PulpettiPalkki):
                multiplier = 10
                value1 = -1 * multiplier * self.cross_section.get_A(0) * self.rho * 1e-9
                value2 = -1 * multiplier * self.cross_section.get_A(1) * self.rho * 1e-9
                direction = 'y'
                frame.add(LineLoad(self, [value1, value2], direction, load_id=load_id))
            elif isinstance(self.cross_section, HarjaPalkki):
                print('Self weight cannot be added to Harjapalkki')
            elif isinstance(self.cross_section, KaarevaHarjaPalkki):
                print('Self weight cannot be added to KaarevaHarjapalkki')
            elif isinstance(self.cross_section, KaarevaPalkki):
                print('Self weight cannot be added to KaarevaPalkki')
            elif isinstance(self.cross_section, MahaPalkki):
                print('Self weight cannot be added to MahaPalkki')
            elif isinstance(self.cross_section, TimberSection):
                # self.weight is kg's, multiplier changes it to N's
                multiplier = 10
                value = -1 * multiplier * self.A * self.rho * 1e-9
                direction = 'y'
                frame.add(LineLoad(self, [value, value], direction, load_id=load_id))

    @property
    def weight(self):
        if isinstance(self.cross_section, (KaarevaPalkki, KaarevaHarjaPalkki, MahaPalkki, HarjaPalkki)):
            return self.member.V_tot() * self.material.rhomean * 1e-9
        elif isinstance(self.cross_section, PulpettiPalkki):
            return self.cross_section.get_A(0.5) * self.length * self.rho * 1e-9
        elif isinstance(self.cross_section, TimberSection):
            return self.A * self.length * self.rho * 1e-9

    def check_cross_section(self):
        """ Checks if cross-section is strong enough.
            Gives boolean value for self.is_strong_enough
            Saves list of stress ratios to self.r
        """
        self.r = np.zeros_like(self.r) # [0, 1, 2, .. 8]
        # Cross-sectional stress ratios in member's nodes' locations
        self.r[:6] = self.member.check_sections()
        # Buckling about y - and z - axis
        buckling_r = self.member.check_buckling()
        self.r[6] = buckling_r[0]
        self.r[7] = buckling_r[1]
        # Lateral-Torsional buckling
        self.r[8] = self.member.check_LT_buckling()

        if max(self.r) > 1.0:
            self.is_strong_enough = False
        else:
            self.is_strong_enough = True

    def print_utility_factor(self, apxspt_kaava6_55=False):
        END = '\033[0m'
        print(f'{self.cross_section.__class__.__name__} käyttöasteet')
        print('______________________________________')
        ratios = self.member.check_sections()
        COLOR = self.get_printing_color(ratios[0])
        print(f'{COLOR}EC5(6.1)\t  UN\t\t {round(ratios[0] * 100, PREC)} %{END}')
        COLOR = self.get_printing_color(ratios[1])
        print(f'{COLOR}EC5(6.13)\t  UV\t\t {round(ratios[1] * 100, PREC)} %{END}')
        COLOR = self.get_printing_color(ratios[2])
        print(f'{COLOR}EC5(6.11)\t  UM\t\t {round(ratios[2] * 100, PREC)} %{END}')
        # COLOR = self.get_printing_color(ratios[3])
        # print(f'{COLOR}EC5(6.14)\t  UT\t\t {round(ratios[3] * 100, PREC)} %{END}')
        COLOR = self.get_printing_color(ratios[4])
        print(f'{COLOR}EC5(6.17)\t  UBT\t\t {round(ratios[4] * 100, PREC)} %{END}')
        COLOR = self.get_printing_color(ratios[5])
        print(f'{COLOR}EC5(6.19)\t  UBC\t\t {round(ratios[5] * 100, PREC)} %{END}')

        buck = self.member.check_buckling()
        COLOR = self.get_printing_color(buck[0])
        print(f'{COLOR}EC5(6.23)  buckling y\t {round(buck[0] * 100, PREC)} %{END}')
        COLOR = self.get_printing_color(buck[1])
        print(f'{COLOR}EC5(6.24)  buckling z\t {round(buck[1] * 100, PREC)} %{END}')
        ms = self.member.get_max_sigma_per_segment()
        # print(f'abs max sigma per y segment: {ms[0]}')
        # print(f'abs max sigma per z segment: {ms[1]}')
        # print('sauvan hoikkuus', u'\u03BBy', '=\t', round(self.member.lamda()[0][ms[0].index(max(ms[0]))], PREC))
        # print('sauvan hoikkuus', u'\u03BBz', '=\t', round(self.member.lamda()[1][ms[1].index(max(ms[1]))], PREC))

        if self.mtype in ('beam', 'rigid-beam'):
            lt = self.member.check_LT_buckling()
            COLOR = self.get_printing_color(lt)
            print(f'{COLOR}EC5(6.35/33) LT buckling {round(lt * 100, PREC)} %{END}')
            UNcp = self.member.check_perpendicular_compression(location=True)
            COLOR = self.get_printing_color(UNcp[0])
            print(f'{COLOR}EC5(6.2)\t  UNcp\t\t {round(UNcp[0] * 100, PREC)} %{END}')
            if UNcp[1] is not None:
                print('leimapaineen solmu: x:', UNcp[1].x, 'y:', UNcp[1].y)

            if isinstance(self.cross_section, (KaarevaPalkki, KaarevaHarjaPalkki, HarjaPalkki, MahaPalkki)):
                apbs = self.member.check_apex_bending_stress()
                COLOR = self.get_printing_color(apbs)
                print(f'{COLOR}EC5(6.41)\t APXBS\t\t {round(apbs * 100, PREC)} %{END}')
                show = True
                if self.material.type == 'lvl':
                    if self.material.direction == 'flat':
                        show = False
                if show:
                    appt = self.member.check_apex_perpendicular_tension()
                    COLOR = self.get_printing_color(appt)
                    print(f'{COLOR}EC5(6.50)\t APXPT\t\t {round(appt * 100, PREC)} %{END}')
                    apspt = self.member.check_apex_shear_perpendicular_tension_combined(kaava6_55=apxspt_kaava6_55)
                    COLOR = self.get_printing_color(apspt)
                    print(f'{COLOR}EC5(6.53)\t APXSPT\t\t {round(apspt * 100, PREC)} %{END}')
                else:
                    print()
            brC, brFd = self.member.bracing()
            print('Z-suunnan stabiloivan tuen vaadittu jousijäykkyys')
            print(f'EC5(9.34)\t\tC\t\t {round(brC, 1)} N/mm')
            print('Z-suunnan stabiloivan tuen voima (1. muoto)')
            print(f'EC5(9.35)\t\tFd\t\t {round(brFd/1000, 2)} kN')
        print('--------------------------------------')

    def get_printing_color(self, val):
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        if val > 1:
            return RED
        elif 0.8 < val <= 1:
            return YELLOW
        else:
            return GREEN

    def print_info(self, ids, parent, prec: int = 1, firetime: int = 0, beam_characteristic_divider=300,
                   beam_quasi_permanent_divider=200, column_divider=150, apxspt_kaava6_55=False):
        print(f'\n{self.mem_id}: {self.cross_section}')
        print('------------------------------------------------')
        if LoadIDs.ULS in ids:
            print('ULS')
            parent.assign_forces(LoadIDs.ULS)
            self.print_utility_factor(apxspt_kaava6_55)
        if LoadIDs.ACC in ids:
            print('ACC')
            print(f'Paloaika {firetime}')
            parent.R = firetime
            parent.assign_forces(LoadIDs.ACC)
            self.print_utility_factor(apxspt_kaava6_55)
            parent.R = 0

        i_list = self.member.radius_of_gyration(get_locs=True)
        iy_paikka = i_list[1][0][list(i_list[0][0]).index(max(i_list[0][0]))]
        print(f'Laskentaparametrit on laskettu suurimman jännityksen\n'
              f'kohdasta suhteellisella etäisyydellä {iy_paikka} sauvan alkupäästä.\n')

        print(f'A: {round(self.cross_section.get_A(iy_paikka) / 1e4, prec)} dm^2')
        print(f'Iy: {round(self.cross_section.get_I(iy_paikka)[0] / 1e8, prec)} dm^4')
        print(f'Iz: {round(self.cross_section.get_I(iy_paikka)[1] / 1e8, prec)} dm^4')
        print(f'Wy: {round(self.cross_section.get_Wel(iy_paikka)[0] / 1e6, prec)} dm^3')
        print(f'Wz: {round(self.cross_section.get_Wel(iy_paikka)[1] / 1e6, prec)} dm^3')
        print(f'Sauvan pituus: {round(self.length, prec)} mm')
        print(f'Sauvan tilavuus: {round(self.member.V_tot() / 1e9, prec)} m^3')

        print('\nMatriaali info')
        print(f'käyttöluokka: {self.member.sc}')
        print(f'aikaluokka: {self.member.ldc}')
        print(f'Materiaali: {self.material}')
        if self.material.type == 'lvl':
            if self.material.direction == 'edge':
                print(f'sauva asetettu syrjälleen')
            else:
                print(f'sauva asetettu lappeeleen')

        print(f'gammaM: {self.member.gammaM}')
        print(f'kdef: {self.member.get_kdef()}')
        print(f'kmod: {self.member.kmod}')
        print(f'k_L: {self.member.k_L}')

        print(f'fmd: {round(self.member.fmd, prec)} N/mm^2')
        print(f'ft0d: {round(self.member.ft0d, prec)} N/mm^2')
        try:
            print(f'ft90d: {round(self.member.ft90d, prec)} N/mm^2')
        except:
            pass
        print(f'fc0d: {round(self.member.fc0d, prec)} N/mm^2')
        print(f'fc90d: {round(self.member.fc90d, prec)} N/mm^2')
        print(f'fvd: {round(self.member.fvd, prec)} N/mm^2')
        try:
            print(f'ft90edged: {round(self.member.ft90edged, prec)} N/mm^2')
            print(f'fr0d: {round(self.member.fr0d, prec)} N/mm^2')
        except:
            pass

        try:
            print(f'Eulerin beta-arvo y-suunta {self.member.beta[0]}')
            print(f'Eulerin beta-arvo z-suunta {self.member.beta[1]}')
        except:
            print(f'Eulerin beta-arvo {self.member.beta}')
        print(f'\niy: {round(max(i_list[0][0]), prec)}')
        print(f'iz: {round(max(i_list[0][1]), prec)}')
        ln = self.member.Ln()
        print(f'Lny: {round(max(ln[0]), prec)}')
        print(f'Lnz: {round(max(ln[1]), prec)}')

        lamda = self.member.lamda()
        print(f'lambda y: {round(max(lamda[0]), prec)}')
        print(f'lambda z: {round(max(lamda[1]), prec)}')
        l_rel = self.member.lamda_rel()
        print(f'lambda rel y: {round(max(l_rel[0]), prec)}')
        print(f'lambda rel z: {round(max(l_rel[1]), prec)}')
        sigma = self.member.sigma()
        print(f'sigma_crit y: {round(max(sigma[0]), prec)} N/mm^2')
        print(f'sigma_crit z: {round(max(sigma[1]), prec)} N/mm^2')
        print(f'beta_c: {self.member.beta_c()}')
        print(f'ky: {round(max(self.member.k()[0]), prec)}')
        print(f'kz: {round(max(self.member.k()[1]), prec)}')
        print(f'kcy: {round(max(self.member.k_c()[0]), prec)}')
        print(f'kcz: {round(max(self.member.k_c()[1]), prec)}')

        try:
            print(f'kr: {self.member.k_r()}')
            print(f'V: {round(self.member.V() / 1e6, prec)} dm^3')
            print(f'V tot: {round(self.member.V_tot() / 1e6, prec)} dm^3')
            print(f'k vol: {round(self.member.k_vol(), prec)}')
            print(f'k dis: {round(self.member.k_dis(), prec)}')
            print(f'k_p: {round(self.member.k_p(), prec)}')
            print(f'k_l: {round(self.member.k_l(), prec)}')
        except:
            pass

        if self.mtype in ('beam', 'rigid-beam'):
            print('\nLT Buckling')
            lt = self.member.check_LT_buckling(True)
            print(f'LT buckling laskettu paikassa {lt["loc"]} kaavalla 6.{lt["kaava"]}')
            print(f'l ef: {round(self.member.l_ef(lt["i"], 1, lt["loc"]), prec)} mm')
            print(f'sigma mcrit: {round(self.member.sigma_mcrit(lt["i"], 1, lt["loc"]), prec)} MPa')
            print(f'lambda relm: {round(self.member.lamda_relm(lt["i"], 1, lt["loc"]), prec)}')
            print(f'k crit: {round(self.member.k_crit(lt["i"], 1, lt["loc"]), prec)}')

        print('\nSisäiset voimasuureet kuormitusyhdistelmittäin')

        for id in ids:
            print(f'\nid: {id}')
            self.assign_forces(id)
            print(f'Ned: {round(max(self.member.ned, key=abs) / 1e3, prec)} kN ({self.member.loc[self.member.ned.index(max(self.member.ned, key=abs))]})')
            print(f'Myed: {round(max(self.member.myed, key=abs) / 1e6, prec)} kNm ({self.member.loc[self.member.myed.index(max(self.member.myed, key=abs))]})')
            print(f'Mzed: {round(max(self.member.mzed, key=abs) / 1e6, prec)} kNm ({self.member.loc[self.member.mzed.index(max(self.member.mzed, key=abs))]})')
            print(f'Vyed: {round(max(self.member.vyed, key=abs) / 1e3, prec)} kN ({self.member.loc[self.member.vyed.index(max(self.member.vyed, key=abs))]})')
            print(f'Vzed: {round(max(self.member.vzed, key=abs) / 1e3, prec)} kN ({self.member.loc[self.member.vzed.index(max(self.member.vzed, key=abs))]})')
            stod_list = [self.member.sigma_t0d(n) for n in range(self.member.nsect())]
            print(f'Sigma t0d: {round(max(stod_list, key=abs), prec)} N/mm^2 ({self.member.loc[stod_list.index(max(stod_list, key=abs))]})')
            print(f'Sigma c0d: {round(max(stod_list, key=abs), prec)} N/mm^2 ({self.member.loc[stod_list.index(max(stod_list, key=abs))]})')
            try:
                print(f'sigma t90d: {round(self.member.sigma_t90d(), prec)} N/mm^2')
                print(f'km alpha: {round(self.member.k_m_alpha(), prec)}')
            except:
                pass
            myd_list = [self.member.sigma_md(n)[0] for n in range(self.member.nsect())]
            mzd_list = [self.member.sigma_md(n)[1] for n in range(self.member.nsect())]
            print(f'sigma myd: {round(max(myd_list, key=abs), prec)} N/mm^2 ({self.member.loc[myd_list.index(max(myd_list, key=abs))]})')
            print(f'sigma mzd: {round(max(mzd_list, key=abs), prec)} N/mm^2 ({self.member.loc[mzd_list.index(max(mzd_list, key=abs))]})')
            ty = [self.member.tau_d(n)[0] for n in range(self.member.nsect())]
            tz = [self.member.tau_d(n)[1] for n in range(self.member.nsect())]
            print(f'tau yd: {round(max(ty, key=abs), prec)} N/mm^2 ({self.member.loc[ty.index(max(ty, key=abs))]})')
            print(f'tau zd: {round(max(tz, key=abs), prec)} N/mm^2 ({self.member.loc[tz.index(max(tz, key=abs))]})')

            if id in (LoadIDs.SLS_Characteristic, LoadIDs.SLS_Quasi_permanent):
                print('Tähän kuormitusyhdistelyyn liittyvä taipuma/siirtymä')
                if self.mtype in ('beam', 'rigid-beam'):
                    print(f'siirtymä x: {round(abs(self.dx[id]), prec)} mm')
                    print(f'taipuma y: {round(abs(self.dy[id]), prec)} mm')
                    if id == LoadIDs.SLS_Characteristic:
                        print(f'suurin sallittu taipuma L/{beam_characteristic_divider}: {round(self.length / beam_characteristic_divider, prec)} mm')
                    elif id == LoadIDs.SLS_Quasi_permanent:
                        print(f'suurin sallittu taipuma L/{beam_quasi_permanent_divider}: {round(self.length / beam_quasi_permanent_divider, prec)} mm')
                else:
                    print(f'siirtymä x: {round(abs(self.dx[id]), prec)} mm')
                    print(f'siirtymä y: {round(abs(self.dy[id]), prec)} mm')
                    print(f'suurin sallittu siirtymä L/{column_divider}: {round(self.length / column_divider, prec)} mm')
        print(f'')


class SteelBeam(SteelFrameMember):
    def __init__(self, coordinates, alpha1=25, alpha2=25, mem_id="",
                 profile="IPE 100", material="S355", num_elements=4, **kwargs):
        super().__init__(coordinates, mem_id, profile, material, num_elements,
                         **kwargs)

        self.mtype = 'beam'
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.check_mtype()

    def check_mtype(self):
        start_node, end_node = self.coordinates
        x1, y1 = start_node
        x2, y2 = end_node

        try:
            angle = abs(y2 - y1) / abs(x2 - x1)

            if angle >= 0.8:
                self.mtype = 'column'
        except ZeroDivisionError:
            self.mtype = 'column'


class SteelColumn(SteelFrameMember):
    def __init__(self, coordinates, mem_id="", profile="IPE 100",
                 material="S355", num_elements=4, LT_buckling=False):
        super().__init__(coordinates, mem_id, profile, material, num_elements,
                         LT_buckling=LT_buckling)

        self.mtype = 'column'
        self.check_mtype()

    def check_mtype(self):
        start_node, end_node = self.coordinates
        x1, y1 = start_node
        x2, y2 = end_node

        try:
            angle = abs(y2 - y1) / abs(x2 - x1)
            if angle < 0.8:
                self.mtype = 'beam'
        except ZeroDivisionError:
            None


# --------------------------- LOAD CLASSES ----------------------------

class LoadCombination:

    def __init__(self, name="", desc=""):
        self.name = name
        self.desc = desc
        self.idx = None
        self.loadcases = []


    def add(self, loadcase):


        if isinstance(loadcase, LoadCase):
            self.loadcases.append(loadcase)
        else:
            raise ValueError(f"Can't add {type(loadcase)},"
                             f"added load must be a LoadCase -object!")



class LoadCase:
    """
    Class for constructing load cases.
    LoadCase can have multiple different loads acting on members

    :param name: name of the load case
    :param desc: description
    """

    def __init__(self, name="", desc="", ltype="A"):
        """


        """
        self.name = name
        self.desc = desc
        self.idx = None
        self.ltype = ltype
        self.loads = []


    def add(self, load):
        """
        Adds load to load case
        :param load: load to be added

        """
        if isinstance(load, Load):
            self.loads.append(load)

        else:
            raise ValueError(f"Can't add {type(load)},"
                             f"added load must be a Load -object!")


class Load:
    """ General class for loads """

    def __init__(self, load_id, ltype='live', name='Load', f=1.0):
        """ Constructor """
        self.load_id = load_id
        self.name = name
        self.ltype = ltype


class PointLoad(Load):
    """ Class for point loads (forces and moments) """

    def __init__(self, coordinate, v, load_id=LoadIDs.ULS, ltype='live',
                 name='PointLoad',
                 f=1.0):
        super().__init__(load_id, ltype, name, f)
        """ Input:
            sid -- load id
            nid -- node subjected to the load
            v -- load vector: v[0] = Fx, v[1] = Fy, v[2] = Mz
            f -- scaling factor
        """

        self.__coordinate = coordinate
        self.v = np.asarray(v)
        self.f = f
        self.node = None

    @property
    def coordinate(self):
        if self.node:
            return self.node.coord
        else:
            return self.__coordinate

    def add_load(self, fem_model):

        if self.coordinate in fem_model.nodal_coords:
            idx = fem_model.nodal_coords.index(self.coordinate)
            self.node = fem_model.nodes[idx]
            fem_model.add_load(
                fem.PointLoad(self.load_id, self.node, self.v, self.f))


class LineLoad(Load):
    def __init__(self, coordinates, values, direction, load_id=LoadIDs.ULS, f=1.0,
                 ltype='live', name='LineLoad',coord_sys="global"):
        """ Class for line loads. Line load can be uniform or trapezoidal
        
            input:
                coordinates .. FrameMember type object that is the subject of the load
                values .. list type object. values[0] and values[1] are the magnitudes
                          of the load at the ends of the member
                direction .. 'x' or 'y'
                load_id .. integer type load identifier
                f .. scaling factor
                ltype .. type of load
                name .. name of load
                coord_sys .. "global" or 'local'. If 'local' the load is given in the
                             local coordinate system of the member.
        """
        super().__init__(load_id, ltype, name, f)
        if isinstance(coordinates, FrameMember):
            self.member = coordinates
            self.member.has_load = True
            self.member.q0 = values[0]
            self.member.q1 = values[1]            
            self.member.loads.append(self)
        else:
            self.mem_id = None
            self.member = None

        self.values = np.asarray(values)
        self.direction = direction
        self.f = f
        self.element_ids = None
        self.coord_sys = coord_sys

    @property
    def coordinates(self):
        return self.member.coordinates

    def calc_k(self):
        """ Evaluate the change in the value of a line load
            per element on a member
        """
        v0, v1 = self.values
        k = (v1 - v0) / len(self.member.elements)
        return k

    def add_load(self, fem_model):
        """ Add the load to the FEM model """
        #k = self.calc_k()
        v0, v1 = self.values
        qvals = np.linspace(v0,v1,len(self.member.elements)+1)        
        for i, elem_id in enumerate(self.member.elements.keys()):
            #v1 = (i * k) + v0
            load = fem.LineLoad(self.load_id,
                                fem_model.elements[elem_id],
                                [0.0, 1.0],
                                [qvals[i],qvals[i+1]],
                                #[v0, v1],
                                self.direction,
                                self.coord_sys)
            fem_model.add_load(load)
            v0 = v1


# --------------------------- SUPPORT CLASSES ----------------------------
class Support:
    def __init__(self, coordinate, dofs, supp_id=1):
        """ Constructor
            coordinate .. nodal coordinates [list]
            dofs .. list of supported degrees of freedom
            supp_id .. each support has an ID that relates it to a given
                    support system. NOTE! In 'calculate' method, supp_id = 1 is always assumed!
        """
        self.node = None
        self.node_id = None
        self.coordinate = [round(c, PREC) for c in coordinate]
        self.supp_id = supp_id
        self.dofs = dofs

    @property
    def coordinate(self):
        if self.node and list(self.node.coord) != self.__coordinate:
            self.__coordinate = list(self.node.coord)
        return self.__coordinate

    @coordinate.setter
    def coordinate(self, val):
        self.__coordinate = val
        if self.node:
            self.node.coord = np.array(val)

    def add_support(self, fem_model):
        if self.coordinate in fem_model.nodal_coords:
            idx = fem_model.nodal_coords.index(self.coordinate)
            self.node_id = idx
            self.node = fem_model.nodes[idx]
            fem_model.add_support(self.supp_id, self.node_id, self.dofs,
                                  val=0.0)
        else:
            pass


class FixedSupport(Support):
    '''
    In ROBOT this setting means that UX, UZ and RY check marks are present in the support dialog
    '''
    def __init__(self, coordinate, supp_id=1):
        # super().__init__(coordinate, [1, 1, 1], supp_id)
        super().__init__(coordinate, [0, 1, 2], supp_id)


class XHingedSupport(Support):
    '''
    In ROBOT this setting means that UX direction check marks is present in the support dialog, but UZ and RY are empty
    '''
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0], supp_id)


class YHingedSupport(Support):
    '''
    In ROBOT this setting means that UZ direction check mark is present in the support dialog, but UX and RY are empty
    '''
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [1], supp_id)


class XYHingedSupport(Support):
    '''
    In ROBOT this setting means that UX and UZ direction check marks are present in the support dialog, but RY is empty
    '''
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0, 1], supp_id)


class Hinge(Support):
    '''
    In ROBOT this setting means all of the check marks removed from the support dialog
    '''
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [], supp_id)



def test():
    frame = Frame2D()
    mem = SteelBeam([[0, 0], [2000, 0]], num_elements=10)
    frame.add(FixedSupport([0, 0]))
    frame.add(FixedSupport([2000, 0]))
    frame.add(LineLoad(mem, [-10, -10], 'y'))
    frame.add(mem)
    frame.generate()
    frame.calculate()
    frame.bmd(100)
    
    return frame

def test_portal():
    frame = Frame2D()
    H = 6000
    L = 8000
    X = [[0,0],[0,H],[L,H],[L,0]]
    
    col1 = SteelColumn([X[0],X[1]], profile="HE 240 A")
    col2 = SteelColumn([X[3],X[2]], profile="HE 240 A")
    beam = SteelBeam([X[1],X[2]], profile="IPE 300")
    frame.add(FixedSupport(X[0]))
    frame.add(FixedSupport(X[3]))
    frame.add(LineLoad(beam, [-20, -20], 'y'))
    frame.add(col1)
    frame.add(col2)
    frame.add(beam)
    frame.generate()
    frame.calculate()
    #frame.bmd(10)
    
    return frame


if __name__ == '__main__':
    fr = Frame2D()
    x0 = [0,0]
    x1 = [8000,0]
    fr.add(FrameMember([x0,[0,6000]]))
    fr.add(FrameMember([[0,6000],[8000,6000]],mtype='rigid_beam'))
    #fr.members[1].add_hinge(0)
    #fr.members[1].add_hinge(1)
    fr.add(FrameMember([[8000,6000],x1]))
                        
    fr.add(LineLoad(fr.members[1],[-20,-20],'y',coord_sys='local',load_id=LoadIDs.ULS))
    fr.add(LineLoad(fr.members[0],[8,8],'x',coord_sys='global',load_id=LoadIDs.SLS_Characteristic))
    fr.add(FixedSupport(x0))
    fr.add(FixedSupport(x1))

    # coord1 = [[0,0], [0, 5000]]
    # coord2 = [[0,5000], [10000, 7000]]
    # coord3 = [[10000, 0], [10000, 7000]]
    #coord1 = [[0, 0], [10000, 0]]
    #fr.add(FrameMember(coord1))
#    fr.add(FrameMember(coord2))
#    fr.add(FrameMember(coord3))
    #fr.add(PointLoad([5000, 0], [0, -2000, 0]))
    #fr.add(XYHingedSupport(coord1[0]))
    #fr.add(XYHingedSupport(coord1[1]))
    fr.generate()
    fr.calculate(load_id='all',support_method="REM")
    fr.optimize_members("CURRENT")
    fr.bmd(6)

    #frame.plot_buckling(scale=1000,k=4)

    """
    for mem in frame.members.values():
        for key, value in mem.nodal_forces.items():
            print(key, value)
    """
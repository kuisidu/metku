"""
Created on Fri Jan 19 10:39:56 2018

@author: huuskoj
"""
import sys

# path needs to be added to use tables_and_tuples and I_sections
sys.path.append('.\End-plate')

import numpy as np
import math
import matplotlib.pyplot as plt

import fem.frame.frame_fem as fem
from fem.frame.elements.ebbeam import EBBeam
from fem.frame.elements.eb_semi_rigid_beam import EBSemiRigidBeam

from sections.steel.ISection import IPE, HEA, HEAA, HEB, HEC, HEM
from sections.steel.RHS import RHS
from sections.steel.CHS import CHS

from steel_member import SteelMember

import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename

import tables_and_tuples

# profiles from lightest to heaviest
PROFILES = [
    'IPE 80', 'IPE 100', 'IPE 120', 'HE 100 AA', 'IPE 140',
    'HE 120 AA', 'IPE 160', 'HE 100 A', 'HE 140 AA', 'IPE 180',
    'HE 120 A', 'HE 100 B', 'IPE 200', 'HE 160 AA', 'HE 140 A',
    'IPE 220', 'HE 120 B', 'HE 180 AA', 'HE 160 A', 'IPE 240',
    'HE 100 C', 'HE 140 B', 'HE 200 AA', 'HE 180 A', 'IPE 270',
    'HE 120 C', 'HE 220 AA', 'HE 100 M', 'IPE 300', 'HE 200 A',
    'HE 160 B', 'HE 240 AA', 'HE 140 C', 'IPE 330', 'HE 220 A',
    'HE 180 B', 'HE 120 M', 'HE 260 AA', 'IPE 360', 'HE 160 C',
    'HE 240 A', 'HE 280 AA', 'HE 200 B', 'HE 140 M', 'IPE 400',
    'HE 260 A', 'HE 300 AA', 'HE 180 C', 'HE 220 B', 'HE 320 AA',
    'HE 160 M', 'HE 280 A', 'IPE 450', 'HE 340 AA', 'HE 200 C',
    'HE 240 B', 'HE 360 AA', 'HE 300 A', 'HE 180 M', 'IPE 500',
    'HE 400 AA', 'HE 260 B', 'HE 220 C', 'HE 320 A', 'HE 450 AA',
    'HE 200 M', 'HE 280 B', 'HE 340 A', 'IPE 550', 'HE 500 AA',
    'HE 360 A', 'HE 300 B', 'HE 220 M', 'HE 240 C', 'HE 550 AA',
    'IPE 600', 'HE 400 A', 'HE 320 B', 'HE 600 AA', 'HE 260 C',
    'HE 340 B', 'HE 650 AA', 'HE 450 A', 'HE 360 B', 'HE 280 C',
    'HE 700 AA', 'HE 500 A', 'HE 400 B', 'HE 240 M', 'HE 550 A',
    'HE 450 B', 'HE 800 AA', 'HE 260 M', 'HE 300 C', 'HE 600 A',
    'HE 320 C', 'HE 500 B', 'HE 280 M', 'HE 650 A', 'HE 900 AA',
    'HE 550 B', 'HE 700 A', 'HE 600 B', 'HE 1000 AA', 'HE 800 A',
    'HE 650 B', 'HE 300 M', 'HE 700 B', 'HE 320 M', 'HE 340 M',
    'HE 360 M', 'HE 900 A', 'HE 400 M', 'HE 800 B', 'HE 450 M',
    'HE 500 M', 'HE 1000 A', 'HE 550 M', 'HE 600 M', 'HE 900 B',
    'HE 650 M', 'HE 700 M', 'HE 1000 B', 'HE 800 M', 'HE 900 M',
    'HE 1000 M', 'IPE 750']

IPE_PROFILES = ['IPE 80', 'IPE 100', 'IPE 120', 'IPE 140', 'IPE 160',
                'IPE 180', 'IPE 200', 'IPE 220', 'IPE 240', 'IPE 270',
                'IPE 300', 'IPE 330', 'IPE 360', 'IPE 400', 'IPE 450',
                'IPE 500', 'IPE 550', 'IPE 600', 'IPE 750']


# -----------------------------------------------------------------------------
class Frame2D:
    """ Class for semi-rigid frame
    
        Written by Jaakko Huusko
        
        Attributes:
        ========
        f -- FrameFEM-class
        storeys -- number of storeys
        bays -- number of bays
        storey_height -- storey height in meters
        bay_length -- bay length in meters
        members -- dict of members of the frame, key: member name 
        num_member -- number of members in the frame
        support_nodes -- node id's of supports
        num_elements -- number of elements in one member
        line_loads -- dict containing all line loads, key: load id
        nodes -- nodes in list, list index is node id
        nodal_coordinates -- coordinates of every node in a list,
                             list index is node id
        nodal_forces -- dict of nodal forces, key: node id
        nodal_displacements -- dict of nocal displacements, key: node id
        weight -- frame's weight in kg's
        r -- list of every members' utilization ratios
        is_generated -- boolean value if frame is generated
        is_designed -- boolean value if frame is designed
        is_calculated -- boolean value if frame is calculated
        is_strong_enough -- boolean value if frame is strong enough
        penalty_val -- penalty value used in differential evolution
        optimizer -- optimization algorithm used for optimization
        optimize_joints -- boolean value if joints are optimized
        
        
        Methods:
        ========
        add --
        add_materials_and_sections -- 
        calc_nodal_coordinates --
        calc_nodal_forces --
        calc_nodal_displacements --
        generate_columns -- 
        generate_beams -- 
        generate_members --
        generate_nodes -- 
        generate -- 
        calculate --
        design_members --
        frame_weight --
        add_self_weight --
        remove_self_weight --
        design_frame --
        check_members_strength -- 
        rigid_joints --
        hinge_joints --
       
        Uses classes:
        ========   
        FrameFEM
        FrameMember
        
        Examples:
        ========
        TODO
    
    """

    def __init__(self, simple=None, num_elements=5, supports=None):

        self.f = fem.FrameFEM()
        self.members = {}
        self.num_members = 0
        self.support_nodes = []
        self.num_elements = num_elements
        self.point_loads = {}
        self.line_loads = {}
        self.supports = {}
        self.nodes = []
        self.nodal_coordinates = []
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.r = []
        self.is_generated = False
        self.is_designed = False
        self.is_calculated = False
        self.is_strong_enough = False
        self.is_optimized = False
        self.self_weight = False
        self.optimize_joints = True
        self.load_robot = False
        
        #I added this
        self.simple=simple
   
        

        self.penalty_val = 1000


        if simple:
            self.storeys = simple[0]
            self.bays = simple[1]
            self.storey_height = simple[2]
            self.bay_length = simple[3]
            self.generate_members()

        if supports:
            self.generate_supports(supports)

    def add(self, this):
 
        # MEMBERS
        if isinstance(this, FrameMember):
            # Give member an unique id
            # id is used in creating elements with wanted cross-sectional properties
            this.mem_id = int(len(self.members))

            self.members[this.mem_id] = this
            this.calc_nodal_coordinates(self.num_elements)

            # Check if members intersect
            for member in self.members.values():
                coord = this.line_intersection(member.coordinates)
                if isinstance(coord, list) and \
                        this.coordinates[0][0] <= coord[0] <= this.coordinates[1][0] and \
                        this.coordinates[0][1] <= coord[1] <= this.coordinates[1][1]:
                    this.add_node_coord(coord)
                    if member.coordinates[0][0] <= coord[0] <= member.coordinates[1][0] and \
                        member.coordinates[0][1] <= coord[1] <= member.coordinates[1][1]:
                            member.add_node_coord(coord)

            self.nodal_coordinates.extend(this.nodal_coordinates)
            self.nodal_coordinates.sort()

        # POINTLOADS
        elif isinstance(this, PointLoad):

            if this.name in self.point_loads.keys():
                this.name += str(len(self.point_loads))
            self.point_loads[this.name] = this

            for member in self.members.values():
                coord = member.point_intersection(this.coordinate)
                if coord:
                    member.add_node_coord(this.coordinate)


        # LINELOADS
        elif isinstance(this, LineLoad):

            if this.name in self.line_loads.keys():
                this.name += str(len(self.line_loads))
            self.line_loads[this.name] = this

        # SUPPORTS
        elif isinstance(this, Support):
            this.supp_id = len(self.supports)
            self.supports[this.supp_id] = this
            self.support_nodes.append(this.coordinate)

            for member in self.members.values():
                coord = member.point_intersection(this.coordinate)
                if coord:
                    member.add_node_coord(this.coordinate)

        else:
            print(type(this), " is not supported.")
            raise TypeError
            
    def add_materials_and_sections(self):
        """ Adds materials and sections to a list in fem.FrameFEM()
        """
        # Iterate through all frame's members and calculate cross-sections
        # properties and add those and member's material properties to 
        # calculation model
        for member in self.members.values():
            member.add_material(self.f)
            member.add_section(self.f)
            
            
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

    def calc_nodal_forces(self):
        """ Calculates nodal forces and saves values to
            self.nodal_forces - dict       
        """
        # Iterate through every frame's member's node and calculate forces
        # acting on that node
        for member in self.members.values():
            member.calc_nodal_forces()
            for node in member.nodal_forces:
                self.nodal_forces[node] = member.nodal_forces[node]

    def calc_nodal_displacements(self):
        """ Calculates nodal displacements and saves values to
            self.nodal_displacements -dict
        """
        # Iterate through every frame's member's node and calculate 
        # displacements for that node
        for member in self.members.values():
            member.calc_nodal_displacements(self.f)
            for node in member.nodal_displacements.keys():
                self.nodal_displacements[node] = member.nodal_displacements[node]
                
    def calculate(self, load_id=2):
        """ Calculates forces and displacements
            
            Parameters
            ----------
            load_id : int
                Load id
        """
        if self.is_calculated == False:
            self.is_calculated = True
            self.f.add_loadcase(supp_id=1, load_id=load_id)
            self.f.nodal_dofs()

        self.f.linear_statics()
        self.calc_nodal_forces()
        self.calc_nodal_displacements()
        self.check_members_strength()

    def design_members(self, symmetry=True):
        """ Desgins frame's members
            
            Parameters
            ----------
            symmetry : bool
                Set true to design members symmetrically
        """
        self.is_designed = True
        for member in self.members.values():
            member.design_member()

        self.check_members_strength()
                
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
        num_cols = (self.bays+1) * self.storeys
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
        """ Generates beams and adds them to self.members -dict

            Parameters
            ----------
            member_id: int
                member id of the last created column
        """
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
        self.generate_beams()

    def generate_supports(self, supp_type):
        if supp_type == 'fixed':
            for coord in self.nodal_coordinates:
                if coord[1] == 0:
                    self.add(FixedSupport(coord))

        elif supp_type == 'hingedY':
            pass



   

    def generate(self):
        """ Generates the frame
        """
        # If the frame hasn't been created, create members and calculate nodal
        # coordinates, otherwise initialize frame.

        for member in self.members.values():
            member.generate(self.f)

            for coord in member.nodal_coordinates:
                if coord not in self.nodal_coordinates:
                    self.nodal_coordinates.append(coord)

        for support in self.supports.values():
            support.add_support(self.f)

        for pointLoad in self.point_loads.values():
            pointLoad.add_load(self.f)

        for lineLoad in self.line_loads.values():
            member = lineLoad.member
            lineLoad.element_ids = member.lineload_elements(lineLoad.coordinates)
            lineLoad.add_load(self.f)

    def load_robot_csv(self):
        self.load_robot = True
        tk = Tk()
        tk.withdraw()
        # NODES AND SUPPORTS
        f_nodes = askopenfilename(title="Select Autodesk Robot csv-file for nodes",
                                  filetype=(("csv files", "*.csv"), ("All files", "*.*")))
        df_nodes = pd.read_csv(f_nodes, delimiter=';', encoding='utf-16')
        nodes = {}
        for idx, row in df_nodes.iterrows():
            try:
                x = float(str(row['X (m)']).replace(',', '.'))
                z = float(str(row['Z (m)']).replace(',', '.'))
                node = [x, z]
                idx = int(row['Node'])
                nodes[idx] = node

                if row[3] != 'N/A':
                    dofs = []
                    for i, value in enumerate(row[3:6]):
                        if value == 'fixed':
                            dofs.append(i)
                    if dofs:
                        self.supports[len(self.supports)] = Support(node,
                                                                    dofs,
                                                                    supp_id=idx)
            except ValueError:
                pass

        # MEMBERS
        f_members = askopenfilename(title="Select Autodesk Robot csv-file for members",
                                    filetype=(("csv files", "*.csv"), ("All files", "*.*")))

        df_members = pd.read_csv(f_members, delimiter=';', encoding='utf-16')
        members = {}
        for i, row in df_members.iterrows():
            try:
                n1 = int(row['Node 1'])
                n2 = int(row['Node 2'])
                coordinates = [nodes[n1], nodes[n2]]
                profile = row['Section']
                material = row['Material']
                members[int(row['Bar'])] = [coordinates, profile, material, i]
            except ValueError:
                pass

        for mem_id in members.keys():
            coordinates, profile, material, idx = members[mem_id]
            member = SteelBeam(coordinates, profile=profile,
                               material=material, mem_id=idx)
            self.members[mem_id] = member
            member.calc_nodal_coordinates(self.num_elements)

        self.calc_nodal_coordinates()
        f_loads = askopenfilename(title="Select Autodesk Robot csv-file for loads",
                                  filetype=(("csv files", "*.csv"), ("All files", "*.*")))
        df_loads = pd.read_csv(f_loads, delimiter=';', encoding='utf-16')

        for idx, row in df_loads.iterrows():
            try:
                loads = row['Load values'].split('=')
                direction = loads[0]
                load_val, unit = loads[1].split('(')
                load_val = float(load_val.replace(',', '.'))
                unit = unit[:-1]
                load_type = row['Load type']
                if load_type == 'nodal force':
                    nodes_list = row['List'].split(' ')
                    for node in nodes_list:
                        node = int(node)
                        coord = nodes[node]
                        if direction == 'FX':
                            value = [load_val, 0, 0]
                        elif direction == 'FZ':
                            value = [0, load_val, 0]
                        elif direction == 'CY':
                            value = [0, 0, load_val]

                        load = PointLoad(coord, value)
                        self.add(load)
                elif load_type == 'uniform load' and load_val != 0:
                    mem_ids = row['List'].split(' ')
                    for mem_id in mem_ids:
                        mem_id = int(mem_id)
                        coordinates = self.members[mem_id].coordinates
                        if direction == 'PZ':
                            value = [load_val, load_val]
                            load_dir = 'y'

                        load = LineLoad(coordinates, value, load_dir)
                        self.add(load)

            except AttributeError:
                pass

        tk.destroy()
        self.generate()

    def plot(self, print_text=True, show=True, loads=True):
        """ Plots the frame
            
            Parameters
            ----------
            print_text : boolean
                Set true to print member's profiles and names
            Color's meaning:
                blue -- member has load
                green -- member can bear its loads
                red -- member breaks under its loads
                black -- member is added, but not designed
        """
        for member in self.members.values():
            member.plot(print_text)

        # Plot supports
        for support in self.supports.values():
            node_coord = support.coordinate
            if support.dofs == [0, 1, 2]:
                marker = 's'
            elif support.dofs == [0]:
                if node_coord[0] == 0:
                    marker = '3'
                else:
                    marker = '4'
            else:
                marker = '2'
            plt.scatter(node_coord[0], node_coord[1], s=50, c='k', marker=marker)
        if loads:
            self.plot_loads()
        plt.axis('equal')
        if show:
            plt.show()

    def plot_loads(self):

        for load in self.point_loads.values():
            x, y = load.coordinate
            plt.scatter(x, y, c='b', marker='*')

        for lineload in self.line_loads.values():
            x0, y0 = lineload.coordinates[0]
            x1, y1 = lineload.coordinates[1]
            plt.plot([x0, x1], [y0, y1], c='b')

    def plot_deflection(self, scale=1, prec=4):
        """ Draws deflected shape of the frame
            
            Parameters
            ----------
            scale : float, optional
                Scaling factor
            prec : int, optional
                Precision for displacement value
        """
        self.plot(print_text=False, show=False)
        self.calc_nodal_displacements()
        for member in self.members.values():
            X = []
            Y = []

            member.calc_nodal_displacements(self.f)
            max_x = 0
            max_y = 0
            for i in range(len(member.nodes)):
                node = member.nodes[i]
                x0 = member.nodal_coordinates[i][0]
                y0 = member.nodal_coordinates[i][1]
                x1 = member.nodal_displacements[node][0]
                y1 = member.nodal_displacements[node][1]
                x = x0 + x1 * (scale) / 1000
                y = y0 + y1 * (scale) / 1000
                if abs(x1) >= abs(max_x) and member.mtype == "column":
                    max_x = abs(x1)
                    loc_max_x = x
                    loc_max_y = y
                if abs(y1) >= abs(max_y) and member.mtype == "beam":
                    max_y = abs(y1)
                    loc_max_x = x
                    loc_max_y = y
                X.append(x)
                Y.append(y)
            plt.plot(X, Y, color='gray')

            if member.mtype == "beam":
                plt.plot(loc_max_x, loc_max_y, 'ro')
                plt.text(loc_max_x, loc_max_y,
                         (str(max_y)[0:prec + 1] + " mm"))

            if member.mtype == "column":
                plt.plot(loc_max_x, loc_max_y, 'ro')
                plt.text(loc_max_x, loc_max_y,
                         (str(max_x)[0:prec + 1] + " mm"))

        plt.show()

    def bmd(self, scale=1):
        """ Draws bending moment diagram
            
            Parameters
            ----------
            scale : int, optional
                Scaling factor
        """
        self.plot(print_text=False, show=False)
        for member in self.members.values():
            member.bmd(scale)
        plt.show()

    def smd(self, scale=1):
        """ Draws bending moment diagram

            Parameters
            ----------
            scale : int, optional
                Scaling factor
        """
        self.plot(print_text=False, show=False)
        for member in self.members.values():
            member.smd(scale)
        plt.show()



    @property
    def weight(self):
        """ Calculates frame's weight and saves it to self.weight
        """
        weight = 0
        for member in self.members.values():
            weight += member.weight
        return weight

    def add_self_weight(self):
        """ Adds self weight to frame's members
        """
        if self.self_weight == False:
            self.self_weight = True
            for member in self.members.values():
                member.add_self_weight()
            self.add_loads('self_weight')
            self.generate_frame()
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

    def check_members_strength(self):
        """ Checks if members can bear their loads
        """
        self.is_strong_enough = True
        self.r.clear()
        for member in self.members.values():
            member.check_cross_section()
            self.r.append(member.r)
            if not member.is_strong_enough:
                self.is_strong_enough = False

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
        MIN_VAL = 1e-20
        for member in self.members.values():
            if member.mtype == "beam":
                member.alpha1 = MIN_VAL
                member.alpha2 = MIN_VAL
                member.Sj1 = MIN_VAL
                member.Sj2 = MIN_VAL


# -----------------------------------------------------------------------------
class SimpleFrame:

    def __init__(self, vals, num_elements=5, supports='fixed'):

        self.storeys = vals[0]
        self.bays = vals[1]
        self.storey_height = vals[2]
        self.bay_length = vals[3]
        self.supp_type = supports




# -----------------------------------------------------------------------------
class FrameMember(SteelMember):
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
        SemiRigidFrame
        FrameFEM
        I_sections
        CrossSection
        SteelMember
            
    """

    def __init__(self, coordinates, mem_id="",
                 profile="IPE 100", material="S355",
                 Sj1=np.inf, Sj2=np.inf):

        self.element_ids = []
        self.elements = {}
        self.cross_section = None
        self.nodes = []
        self.nodal_coordinates = []
        self.loc = []
        self.coordinates = coordinates
        self.material = material
        self.cross_section = None
        self.__profile = profile
        self.profile = profile
        self.profile_idx = PROFILES.index(self.profile)
        self.length = self.calc_length()
        self.mem_id = mem_id
        self.mtype = ""
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.loads = {}
        self.alpha1 = 25
        self.alpha2 = 25
        self.__Sj1 = None
        self.__Sj2 = None
        self.Sj1 = Sj1
        self.Sj2 = Sj2
        self.ned = 0
        self.ved = 0
        self.med = 0
        self.has_load = False
        self.is_designed = False
        self.is_strong_enough = False
        self.r = 1
        self.self_weight = False
        self.is_generated = False
        # Lineload values
        self.q0 = 0
        self.q1 = 0


    @property
    def weight(self):
        """
        Function that calculates member's weight
        :return weight: float, weight of the member
        """
        
        weight = self.cross_section.A * self.length * self.rho * 1000
        return weight

    @property
    def h(self):
        """
        Property, returns member's cross-section's height
        """
        return self.__h

    @h.setter
    def h(self, val):
        """
        Sets cross-section's height to givel value.
        Changes profile and sets new cross-sectional properties to member's elements
        """
        self.__h = val
        splitted = self.profile.split()
        if len(splitted) == 2:
            profile_type = self.profile.split()[0]
            self.profile = profile_type + " " + str(val)
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
        self.__material = val
        self.fy = tables_and_tuples.mat[self.material]['f_y']
        self.E = tables_and_tuples.mat[self.material]['E']
        self.fu = tables_and_tuples.mat[self.material]['f_u']
        self.nu = tables_and_tuples.mat[self.material]['v']
        self.rho = tables_and_tuples.mat[self.material]['rho']

    @property
    def profile(self):
        """
        Returns member's profile e.g. IPE 200, HE 120 A
        """
        return self.__profile

    @profile.setter
    def profile(self, val):
        """
        Changes member's profile to given value.
        Calculates new cross-sectional properties and sets these new values
        to member's elements.
        
        :param val: string, profile name e.g. 'IPE 100'
        """
        self.__profile = val.upper()
        splitted_val = self.profile.split(" ")
        profile_type = splitted_val[0]
        
        print ("profile", self.profile)
        print ("splitted val",splitted_val)
        
        
        try:
            kwargs = tables_and_tuples.profile[self.profile]
            H = int(splitted_val[1])
            catalogue = True
        except KeyError:
            H = float(splitted_val[1])
            catalogue = False
        
        if profile_type == 'IPE':
            self.cross_section = IPE(H, self.fy, catalogue=False)
            
        elif profile_type == 'RHS':
            vals = splitted_val[1].split('X')
            H = vals[0]
            B = vals[2]
            T = vals[4]
            self.cross_section = RHS(H, B, T, self.fy)

            pass
        elif profile_type == 'CHS':
            self.cross_section = CHS(H, T, self.fy)
            
        elif profile_type == 'HE':
            print ("splitted val",splitted_val)
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
            else:
                raise ValueError('{} is not valid profile!'.format(val))
        elif profile_type == 'SHS':
            pass
        
        else:
            raise ValueError('{} is not valid profile type!'.format(profile_type))

        if len(self.elements) > 0:
            for element in self.elements.values():
                element.section.A = self.cross_section.A * 1e-6
                element.section.Iy = self.cross_section.I[0] * 1e-12

    @property
    def Sj1(self):
        return self.__Sj1

    @Sj1.setter
    def Sj1(self, val):

        if len(self.element_ids) > 0:
            idx = self.element_ids[0]
            ele = self.elements[idx]
            ele.rot_stiff = [val, np.inf]

        self.__Sj1 = val

    @property
    def Sj2(self):
        return self.__Sj2

    @Sj2.setter
    def Sj2(self, val):

        if len(self.element_ids) > 0:
            idx = self.element_ids[-1]
            ele = self.elements[idx]
            ele.rot_stiff = [np.inf, val]

        self.__Sj2 = val

    def calc_length(self):
        x0, y0 = self.coordinates[0]
        x1, y1 = self.coordinates[1]
        x = (x1 - x0) ** 2
        y = (y1 - y0) ** 2
        L = np.sqrt(x + y)
        return L

    def alpha_to_Sj(self):
        """ Change alpha value to Sj
        """
        # SFS-EN 1993-1-8 page 60, hinge 0.5 < alpha < 25 rigid
        self.Sj1 = self.alpha1 * self.E * self.I_y * 1e-12 / self.length
        self.Sj2 = self.alpha2 * self.E * self.I_y * 1e-12 / self.length

    def Sj_to_alpha(self):
        """ Change Sj value to alpha
        """
        # SFS-EN 1993-1-8 alpha
        self.alpha1 = self.Sj1 * self.length / (self.E * self.I_y * 1e-12)
        self.alpha2 = self.Sj2 * self.length / (self.E * self.I_y * 1e-12)

    def generate(self, fem_model):
        """ Generates the member
        """
        # If the frame hasn't been created, create members and calculate nodal
        # coordinates, otherwise initialize frame.
        if not self.is_generated:
            self.is_generated = True
            # self.calc_nodal_coordinates(num_elems)
            self.add_nodes(fem_model)
            self.add_material(fem_model)
            self.add_section(fem_model)
            self.generate_elements(fem_model)

    def calc_nodal_coordinates(self, num_elements):
        """
            Calculates node locations along member
            Coordinates used are global coordinates
            Adds locations to a list where they can be accessed later
            :param num_elements: int, number of elements to be created
        """
        start_node = self.coordinates[0]
        end_node = self.coordinates[1]

        x0 = start_node[0]
        y0 = start_node[1]
        Lx = end_node[0] - start_node[0]
        Ly = end_node[1] - start_node[1]
        dx = Lx / num_elements
        dy = Ly / num_elements
        for i in range(num_elements):
            x = (i) * dx + x0
            y = (i) * dy + y0
            loc = i / num_elements
            node = [x, y]
            if node not in self.nodal_coordinates:
                self.nodal_coordinates.append(node)
            if loc not in self.loc:
                self.loc.append(loc)
        if [end_node[0], end_node[1]] not in self.nodal_coordinates:
            self.nodal_coordinates.append([end_node[0], end_node[1]])
        if 1 not in self.loc:
            self.loc.append(1)

    def add_node_coord(self, coord):
        """
        Adds coordinate to memeber's coordinates and
        calculates new coordinates local location
        :param coord: array of two float values, node's coordinates
        """
        if coord not in self.nodal_coordinates:
            self.nodal_coordinates.append(coord)
            self.nodal_coordinates = sorted(self.nodal_coordinates)
            x, y = coord
            x1, y1 = self.coordinates[0]
            dx = x - x1
            dy = y - y1
            dz = math.sqrt((dx ** 2 + dy ** 2))
            z = dz / self.length
            self.loc.append(z)
            self.loc.sort()

    def add_nodes(self, fem_model):
        """ Creates nodes to previously calculated locations
        """
        for coordinate in self.nodal_coordinates:
            x, y = coordinate
            if coordinate in fem_model.nodal_coords:
                self.nodes.append(fem_model.nodal_coords.index(coordinate))
            else:
                self.nodes.append(fem_model.nnodes())
                fem_model.add_node(x, y)

    def generate_elements(self, fem_model):
        """ Generates elements between nodes
            For beam member's first and last element, creates semi-rigid end
            elements. Elements are added to a list
            :param fem_model: FrameFEM -object
        """
        index = fem_model.nels()

        # EBBeam -elements
        if self.mtype == "column":
            
            for i in range(len(self.nodes) - 1):
                n1 = self.nodes[i]
                n2 = self.nodes[i + 1]
                self.elements[index] = \
                    EBBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                           fem_model.sections[self.mem_id], \
                           fem_model.materials[self.mem_id])

                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
                index += 1

        if self.mtype == "beam":
            if len(self.nodes) == 2:
                n1 = self.nodes[0]
                n2 = self.nodes[1]
                self.elements[index] = \
                        EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                                        fem_model.sections[self.mem_id], \
                                        fem_model.materials[self.mem_id], \
                                        rot_stiff=[self.Sj1, self.Sj2])
                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
            else:
                for i in range(len(self.nodes) - 1):
                    n1 = self.nodes[i]
                    n2 = self.nodes[i + 1]
                    # EBSemiRigid -element, first element
                    if i == 0:
                        self.elements[index] = \
                            EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                                            fem_model.sections[self.mem_id], \
                                            fem_model.materials[self.mem_id], \
                                            rot_stiff=[self.Sj1, np.inf])
                    # EBSemiRigid -element, last element
                    elif i == (len(self.nodes) - 2):
                        self.elements[index] = \
                            EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                                            fem_model.sections[self.mem_id], \
                                            fem_model.materials[self.mem_id], \
                                            rot_stiff=[np.inf, self.Sj2])
                    # EBBeam -elements
                    else:
                        self.elements[index] = \
                            EBBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                                   fem_model.sections[self.mem_id], \
                                   fem_model.materials[self.mem_id])
    
                    fem_model.add_element(self.elements[index])
                    self.element_ids.append(index)
                    index += 1

    def calc_nodal_forces(self):
        """ Calculates nodal forces on member's nodes          
        """
        self.med = 0
        self.ved = 0
        self.ned = 0

        i = 0
        for element in self.elements.values():
            axial_force = element.axial_force[0]
            if abs(axial_force) > abs(self.ned):
                self.ned = axial_force

            shear_force = element.shear_force[0]
            if abs(shear_force) > abs(self.ved):
                self.ved = shear_force

            bending_moment = element.bending_moment[0]
            if abs(bending_moment) > abs(self.med):
                self.med = bending_moment

            self.nodal_forces[self.nodes[i]] = [axial_force,
                                                shear_force,
                                                bending_moment]
            i += 1

            if i == len(self.elements):
                axial_force = element.axial_force[1]
                if abs(axial_force) > abs(self.ned):
                    self.ned = axial_force

                shear_force = element.shear_force[1]
                if abs(shear_force) > abs(self.ved):
                    self.ved = shear_force

                bending_moment = element.bending_moment[1]
                if abs(bending_moment) > abs(self.med):
                    self.med = bending_moment

                self.nodal_forces[self.nodes[i]] = [axial_force,
                                                    shear_force,
                                                    bending_moment]

    def calc_nodal_displacements(self, fem_model):
        """ Calculates nodal displacements and saves them to a dict
            :param fem_model: FrameFEM -object
        """
        for node in self.nodes:
            self.nodal_displacements[node] = fem_model.nodes[node].u

    def testbmd(self, scale):
        """ Plots member's bending moment diagram
            Plots currently only horizontal diagrams
            :param scale: float, scales the diagram
        """
        try:
            x1 = []
            y1 = []
            i = 0
            for node in self.nodes:
                forces = self.nodal_forces[node]
                bending_moment = forces[2]
                x1.append(i)
                y1.append(bending_moment)
                i += self.length / len(self.elements)
            x2 = [0] * len(x1)
            max_val = max(y1)
            min_val = min(y1)
            if abs(max_val) > abs(min_val) and max_val != y1[0] \
                    and max_val != y1[-1]:
                val = min_val
            else:
                val = max_val
            max_loc = int(y1.index(val))
            plt.gca().invert_yaxis()
            plt.fill_between(x1, y1, color='lightgray')
            plt.plot(x1, y1)
            plt.plot(x1, x2, color='black')
            plt.text(x1[0], y1[0], str(y1[0]) + "kNm")
            plt.text(x1[-1], y1[-1], str(y1[-1]) + "kNm")
            plt.text(x1[max_loc], y1[max_loc], str(val) + "kNm")
        except ValueError:
            print("Error! Calculate results first.")

    def add_self_weight(self):
        """ Adds self-weight to the member's loads
        """
        if not self.self_weight:
            self.self_weight = True
            load_id = "self_weight"
            # self.weight is kg's, multiplier changes it to kN's
            multiplier = 1e-2
            value = -1 * multiplier * self.weight / self.length
            direction = 'y'
            self.add_line_load(load_id, value, direction)

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
            Young's modulus kN/m^2
            nu is dimensionless
            Density kg/m^3
        """
        fem_model.add_material(self.E, self.nu, self.rho)

    def add_section(self, fem_model):
        """ Units: m**2 A and I_y are in mm**2 and mm**4, respectively"""
        s = fem.BeamSection(self.cross_section.A * 1e-6, self.cross_section.I[0] * 1e-12)
        fem_model.add_section(s)

    def check_cross_section(self):
        """ Checks if cross-section is strong enough.
            Gives boolean value for self.is_strong_enough
            Saves list of stress ratios to self.r
        """
        self.r = []
        self.cross_section.Med = self.ned
        self.cross_section.Ved = self.ved
        self.cross_section.Ned = self.ned
        steel_member = SteelMember(self.cross_section, self.length,
                                   Lcr=[1.0, 1.0], mtype=self.mtype)
        i = 0
        for node in self.nodal_forces.keys():
            forces = self.nodal_forces[node]
            ned = forces[0] * 1e3  # kN to N
            ved = forces[1] * 1e3  # kN to N
            med = forces[2] * 1e6  # kNm to Nmm
            loc = self.loc[i]
            steel_member.add_section(ned=ned, vzed=ved, myed=med, loc=loc)
            i += 1
        """         
        self.r = max(max(steel_member.check_sections()),\
                     (max(steel_member.check_buckling())),\
                     steel_member.check_LT_buckling())
        """
        # Cross-sectional stress ratios in members' nodes' locations
        self.r = steel_member.check_sections()
        # Buckling about y - and z - axis 
        buckling_r = steel_member.check_buckling()
        self.r.append(buckling_r[0])
        self.r.append(buckling_r[1])
        # Lateral-Torsional buckling
        self.r.append(steel_member.check_LT_buckling())

        if max(self.r) > 1.0:
            self.is_strong_enough = False
        else:
            self.is_strong_enough = True

    def design_member(self):
        """ Goes through all profiles in list and stops iterating when 
            cross-section can bear it's loads.
            If member is previously designed, iterating starts from 5 profiles
            before currrent profile.
        """
        i = max(0, self.profile_idx - 5)
        self.is_strong_enough = False
        while not self.is_strong_enough:
            i += 1
            self.profile = PROFILES[i]
            self.check_cross_section()
            self.profile_idx = i

    def line_intersection(self, coordinates):
        """
        Calculates coordinate where two members intersect
        :param coordinates: array of two arrays of two float values, [[x1,y1],[x2,y2]]
        :return [px, py]: array of two float values, coordinates for intersection
        """
        # source: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
        x1, y1 = self.coordinates[0]
        x2, y2 = self.coordinates[1]
        x3, y3 = coordinates[0]
        x4, y4 = coordinates[1]

        if ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)) == 0:
            return None
        else:
            px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / (
                (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
            py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / (
                (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))

            #if  self.coordinates[0][0] < px < self.coordinates[1][0] and \
             #   self.coordinates[0][1] < py < self.coordinates[1][1]:

            return [px, py]
            #else:
            #    return None

    def point_intersection(self, coordinate):
        
        #Calculates point where member and point load intersect
        #:param coordinate: array of two float values
        
        x1, y1 = self.coordinates[0]
        x2, y2 = self.coordinates[1]
        x, y = coordinate
        try:
            k = (y2 - y1) / (x2 - x1)
            b = y1 - k * x1
            return y == k * x + b

        except ZeroDivisionError:
            return x == x1 and y1 <= y <= y2

    def plot(self, print_text):

        X = self.coordinates
        if self.is_strong_enough:
            color = 'green'
        else:
            if self.r == 1:
                color = 'black'
            else:
                color = 'red'
        # Plot members
        plt.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color)
        # Calculate text location
        if self.mtype == 'beam':
            x = X[0][0] + self.length / 2
            y = X[1][1]
            rot = 'horizontal'
            horzalign = 'center'
            vertalign = 'bottom'
        elif self.mtype == 'column':
            x = X[1][0]
            y = X[0][1] + self.length / 2
            rot = 'vertical'
            horzalign = 'right'
            vertalign = 'center'

        # Plot text
        if print_text:
            plt.text(x, y, str(self.mem_id) + ": " + self.profile,
                     rotation=rot, horizontalalignment=horzalign, verticalalignment=vertalign)

    def bmd(self, scale):
        """
        Plots bending moment diagram
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
                y1 = self.nodal_forces[node][2] / (1000 / scale)
                x = x0
                y = y0 - y1
                X.append(x)
                Y.append(y)
            elif self.mtype == "column":
                x0 = self.nodal_coordinates[i][0]
                y0 = self.nodal_coordinates[i][1]
                x1 = self.nodal_forces[node][2] / (1000 / scale)
                x = x0 + x1
                y = y0
                X.append(x)
                Y.append(y)
        X.append(self.nodal_coordinates[i][0])
        Y.append(self.nodal_coordinates[i][1])
        plt.plot(X, Y, color='gray')

    def bmd_test(self):

        X = []
        Y = []
        x0 = self.coordinates[0][0]
        y0 = self.coordinates[0][1]
        dx = self.length / 20
        V0 = self.nodal_forces[self.nodes[0]][1]
        M0 = self.nodal_forces[self.nodes[0]][2]

        X.append(self.coordinates[0][0])
        Y.append(self.coordinates[0][1])
        for i in range(21):
            if self.mtype == 'beam':
                x = (dx * i  + x0)
                y = -1 * (V0*x +0.5 * ((self.q1 - self.q0)/self.length + self.q0)*x**2 + M0)
                y = y/1e2 + y0
                plt.text(x, y, f'{y:.2f}')
            if self.mtype == 'column':
                y = dx * i + y0
                x = 1 * (V0*y +0.5 * ((self.q1 - self.q0)/self.length + self.q0)*y**2 + M0)
                x = x/1e2 + x0

            X.append(x)
            Y.append(y)

        X.append(self.coordinates[1][0])
        Y.append(self.coordinates[1][1])
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

    def lineload_elements(self, coords):
    
    #Returns the elements that are under line load
    #:param coords: 2D-array of float, line load's start and end coordinates
    #:return element_ids: array of int, elements id's
   
        start, end = coords
        start_idx = self.nodal_coordinates.index(start)
        end_idx = self.nodal_coordinates.index(end)
        element_ids = self.element_ids[start_idx:end_idx]
        return element_ids


class SteelBeam(FrameMember):
    def __init__(self, coordinates, alpha1=25, alpha2=25, mem_id="", profile="IPE 100", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'beam'
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.check_mtype()

    def check_mtype(self):
        x1, y1 = self.coordinates[0]
        x2, y2 = self.coordinates[1]

        try:
            angle = abs(y2 - y1) / abs(x2 - x1)

            if angle >= 0.8:
                self.mtype = 'column'
        except ZeroDivisionError:
            self.mtype = 'column'


class SteelColumn(FrameMember):
    def __init__(self, coordinates, mem_id="", profile="IPE 100", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'column'
        self.check_mtype()

    def check_mtype(self):
        x1, y1 = self.coordinates[0]
        x2, y2 = self.coordinates[1]

        try:
            angle = abs(y2 - y1) / abs(x2 - x1)
            if angle < 0.8:
                self.mtype = 'beam'
        except ZeroDivisionError:
            None


# --------------------------- LOAD CLASSES ----------------------------
class Load:
    """ General class for loads """

    def __init__(self, load_id, ltype='live', name='Load', f=1.0):
        """ Constructor """
        self.load_id = load_id
        self.name = name
        self.type = ''


class PointLoad(Load):
    """ Class for point loads (forces and moments) """

    def __init__(self, coordinate, v, load_id=2, ltype='live', name='PointLoad',
                 f=1.0):
        super().__init__(load_id, ltype, name, f)
        """ Input:
            sid -- load id
            nid -- node subjected to the load
            v -- load vector: v[0] = Fx, v[1] = Fy, v[2] = Mz
            f -- scaling factor
        """

        self.coordinate = coordinate
        self.v = v
        self.f = f
        self.node = None

    def add_load(self, fem_model):

        if self.coordinate in fem_model.nodal_coords:
            idx = fem_model.nodal_coords.index(self.coordinate)
            self.node = fem_model.nodes[idx]
            fem_model.add_load(fem.PointLoad(self.load_id, self.node, self.v, self.f))
        else:
            pass


class LineLoad(Load):
    def __init__(self, coordinates, values, direction, load_id=2,f=1.0,
                 ltype='live', name='LineLoad'):
        super().__init__(load_id, ltype, name, f)

        if isinstance(coordinates, FrameMember):
            self.member = coordinates
            self.member.has_load = True
            self.member.q0 = values[0]
            self.member.q1 = values[1]
            coordinates = self.member.coordinates
        else:
            self.mem_id = None
            self.member = None
        self.coordinates = coordinates
        self.values = values
        self.direction = direction
        self.f = f
        self.k = self.calc_k()
        self.element_ids = None

    def calc_k(self):
        x0, y0 = self.coordinates[0]
        x1, y1 = self.coordinates[1]
        x = (x1 - x0) ** 2
        y = (y1 - y0) ** 2
        L = np.sqrt(x + y)
        v1, v2 = self.values
        k = (v2 - v1) / L

        return k

    def add_load(self, fem_model):
        v0, v1 = self.values
        start = v0
        end = v1
        for elem_id in self.element_ids:
            load = fem.LineLoad(self.load_id,
                                fem_model.elements[elem_id],
                                [0.0, 1.0],
                                [v0, v1],
                                self.direction)
            fem_model.add_load(load)

            v0 = v1
            v1 = (1 + self.k) * v0


# --------------------------- SUPPORT CLASSES ----------------------------
class Support:
    def __init__(self, coordinate, dofs, supp_id=1):

        self.coordinate = coordinate
        self.supp_id = supp_id
        self.dofs = dofs
        self.node = None
        self.node_id = None

    def add_support(self, fem_model):

        if self.coordinate in fem_model.nodal_coords:
            idx = fem_model.nodal_coords.index(self.coordinate)
            self.node_id = idx
            self.node = fem_model.nodes[idx]
            fem_model.add_support(self.supp_id, self.node_id, self.dofs, val=0.0)
        else:
            pass


class FixedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0, 1, 2], supp_id)


class XHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0], supp_id)


class YHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [1], supp_id)


class XYHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0, 1], supp_id)


class Hinge(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [], supp_id)

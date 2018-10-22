"""
Created on Fri Jan 19 10:39:56 2018

@author: huuskoj
"""
import sys
# path needs to be added to use tables_and_tuples and I_sections
sys.path.append('S:\91202_METKU\Kristo\Python\src\End-plate')
sys.path.append('S:\91202_METKU\Kristo\Python\src\sections\steel')

from timeit import default_timer as timer

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fmin_slsqp, differential_evolution

from steel_section import SteelSection
from ISection import ISection, HEA, IPE


import frame_fem as fem
import I_sections
import anna_HEA_sections
from pso import pso
from steel_member import SteelMember

import tables_and_tuples
#from profile_tables import A_Iy_Iz_matrix

# profiles from lightest to heaviest
PROFILES = [
                'IPE 80', 'IPE 100', 'IPE 120', 'HE 100 AA', 'IPE 140',
                'HE 120 AA', 'IPE 160','HE 100 A', 'HE 140 AA', 'IPE 180',
                'HE 120 A', 'HE 100 B', 'IPE 200', 'HE 160 AA','HE 140 A',
                'IPE 220', 'HE 120 B', 'HE 180 AA', 'HE 160 A', 'IPE 240',
                'HE 100 C','HE 140 B', 'HE 200 AA', 'HE 180 A', 'IPE 270',
                'HE 120 C', 'HE 220 AA', 'HE 100 M','IPE 300', 'HE 200 A',
                'HE 160 B', 'HE 240 AA', 'HE 140 C', 'IPE 330', 'HE 220 A',
                'HE 180 B', 'HE 120 M', 'HE 260 AA', 'IPE 360', 'HE 160 C',
                'HE 240 A', 'HE 280 AA','HE 200 B', 'HE 140 M', 'IPE 400',
                'HE 260 A', 'HE 300 AA', 'HE 180 C', 'HE 220 B','HE 320 AA',
                'HE 160 M', 'HE 280 A', 'IPE 450', 'HE 340 AA', 'HE 200 C',
                'HE 240 B','HE 360 AA', 'HE 300 A', 'HE 180 M', 'IPE 500',
                'HE 400 AA', 'HE 260 B', 'HE 220 C','HE 320 A', 'HE 450 AA',
                'HE 200 M', 'HE 280 B', 'HE 340 A', 'IPE 550', 'HE 500 AA',
                'HE 360 A', 'HE 300 B', 'HE 220 M', 'HE 240 C', 'HE 550 AA',
                'IPE 600', 'HE 400 A','HE 320 B', 'HE 600 AA', 'HE 260 C',
                'HE 340 B', 'HE 650 AA', 'HE 450 A', 'HE 360 B','HE 280 C',
                'HE 700 AA', 'HE 500 A', 'HE 400 B', 'HE 240 M', 'HE 550 A',
                'HE 450 B','HE 800 AA', 'HE 260 M', 'HE 300 C', 'HE 600 A',
                'HE 320 C', 'HE 500 B', 'HE 280 M','HE 650 A', 'HE 900 AA',
                'HE 550 B', 'HE 700 A', 'HE 600 B', 'HE 1000 AA', 'HE 800 A',
                'HE 650 B', 'HE 300 M', 'HE 700 B', 'HE 320 M', 'HE 340 M',
                'HE 360 M', 'HE 900 A','HE 400 M', 'HE 800 B', 'HE 450 M',
                'HE 500 M', 'HE 1000 A', 'HE 550 M', 'HE 600 M','HE 900 B',
                'HE 650 M', 'HE 700 M', 'HE 1000 B', 'HE 800 M', 'HE 900 M',
                'HE 1000 M','IPE 750']

IPE_PROFILES = ['IPE 80', 'IPE 100', 'IPE 120', 'IPE 140', 'IPE 160',
                'IPE 180', 'IPE 200', 'IPE 220', 'IPE 240', 'IPE 270',
                'IPE 300', 'IPE 330', 'IPE 360', 'IPE 400', 'IPE 450',
                'IPE 500', 'IPE 550', 'IPE 600', 'IPE 750']
#                       A mm^2            Iy mm^4
IPE_A_I_matrix = [[764.3401836602552, 801376.1398006069],
                  [1032.3219599741, 1710119.3948450384],
                  [1321.0219599741001, 3177531.491612866],
                  [1642.6019599741, 5412237.4783326965],
                  [2009.130995059227, 8692922.280149484],
                  [2394.7309950592266, 13169581.740920702],
                  [2848.41065788307, 19431661.97530533],
                  [3337.0506578830696, 27718364.6843977],
                  [3911.6216529422964, 38916214.50755847],
                  [4594.5016529422965, 57897773.31813715],
                  [5381.201652942297, 83561026.91741039], 
                  [6260.623980236907, 117668927.27986754],
                  [7272.9239802369075, 162656174.3911142],
                  [8446.3576397669, 231283456.02658707],
                  [9882.077639766901, 337429141.0374633],
                  [11552.1576397669, 481985026.72959846],
                  [13441.60263153228, 671164648.8112286],
                  [15598.44263153228, 920833980.6682992],
                  [17060.07972311255, 1507046150.557715]]

HEA_PROFILES = ['HE 100 A','HE 120 A','HE 140 A','HE 160 A','HE 180 A',
                'HE 200 A','HE 220 A','HE 240 A','HE 260 A','HE 280 A',
                'HE 300 A','HE 320 A','HE 340 A','HE 360 A','HE 400 A',
                'HE 450 A','HE 500 A','HE 550 A','HE 600 A','HE 650 A',
                'HE 700 A','HE 800 A','HE 900 A','HE 1000 A']

#                       A mm^2            Iy mm^4
HEA_A_I_matrix =[[2120, 3492000],
                 [2530, 6062000],
                 [3140, 10330000],
                 [3880, 16730000],
                 [4530, 25100000],
                 [5380, 36920000],
                 [6430, 54100000],
                 [7680, 77630000],
                 [8680, 104500000],
                 [9730, 136700000],
                 [11250, 182600000],
                 [12440, 229300000],
                 [13350, 276900000],
                 [14280, 330900000],
                 [15900, 450700000],
                 [17800, 637200000],
                 [19750, 869700000],
                 [21180, 1119000000],
                 [22650, 1412000000],
                 [24160, 1752000000],
                 [26050, 2153000000],
                 [28580, 3034000000],
                 [32050, 4221000000],
                 [34680,5538000000]]



OPTIMIZERS = ["PSO", "SLSQP", "DIF_EVOL"]
#-----------------------------------------------------------------------------
class SemiRigidFrame:
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
        nodal_forces -- dict of nocal forces, key: node id
        nodal_displacements -- dict of nocal displacements, key: node id
        weight -- frame's weight in kg's
        r -- list of every members' utilization ratios
        is_generated -- boolean value if frame is generated
        is_designed -- boolean value if frame is designed
        is_calculated -- boolean value if frame is calculated
        is_strong_enough -- boolean value if frame is strong enough
        is_optimized -- boolean value if frame has been optimized
        penalty_val -- penalty value used in differential evolution
        optimizer -- optimization algorithm used for optimization
        optimize_joints -- boolean value if joints are optimized
        
        
        Methods:
        ========
        generate_columns -- 
        generate_beams -- 
        generate_members --
        generate_nodes -- 
        add_materials_and_sections -- 
        add_supports --
        generate_elements --
        calc_nodal_coordinates --
        calc_nodal_forces --
        calc_nodal_displacements --
        generate_frame --
        initialize_frame --
        draw_frame --
        draw_deflected_frame -- 
        draw_bmd -- 
        add_point_load -- 
        add_loads --
        calculate --
        design_members --
        frame_weight --
        add_self_weight --
        remove_self_weight --
        design_frame --
        check_members_strength -- 
        find_closest_profile --
        continuous_to_discrete --
        costf_weight --
        func --
        constraints --
        log_optimization --
        initial_x0 --
        optimize_frame --
        print_member_info --
        rigid_joints --
        hinge_joints --
        
        Uses classes:
        ========   
        FrameFEM
        FrameMember
        
        Examples:
        ========
        SRF_example1.py
        SRF_example2.py
    
    """
    
    def __init__(self, storeys=2, bays=2, storey_height=5,\
                 bay_length=5, num_elements=5):
        
        self.f = fem.FrameFEM()
        self.storeys = storeys
        self.bays = bays
        self.storey_height = storey_height      
        self.bay_length = bay_length            
        self.height = self.storey_height * self.storeys
        self.members = {}
        self.num_members = 0
        self.support_nodes = []
        self.num_elements = num_elements
        self.line_loads = {}
        self.point_loads = {}
        self.nodes = []
        self.nodal_coordinates = []
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.r = []
        self.weight = 0                         
        self.is_generated = False
        self.is_designed = False
        self.is_calculated = False
        self.is_strong_enough = False
        self.is_optimized = False
        self.self_weight = False
        self.penalty_val = 1000
        self.optimizer = "PSO"
        self.optimize_joints = True
        
        self.generate_frame()
  
      
    def generate_columns(self):
        """ Generates columns and adds them to self.members -dict
            
            Returns
            -------
            i : int
                the current member id number
        """
        # x-coordinate
        x = 0
        # y-coordinates 1: start, 2: end
        y1 = 0
        y2 = self.storey_height
        # number of columns
        num_cols = (self.bays+1)*self.storeys
        # create columns as FrameMember object
        for i in range(num_cols):
            self.members["Column"+str(i)] = FrameMember([[x, y1],
                                                         [x, y2]],  
                                                         mtype="column",
                                                         length=(y2-y1),
                                                         mem_id=i,
                                                         f=self.f)
            # If the y2-coordinate equals frame's height, move the x-coordinate
            # by one bay length and initialize y-coordinates
            if y2 == self.storeys * self.storey_height:            
                x += self.bay_length
                y1 = 0
                y2 = self.storey_height
            else:
                y1 = y2
                y2 += self.storey_height
        
        return i
    
    
    def generate_beams(self, member_id):
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
        num_beams = self.bays*self.storeys
        # create beams as FrameMember object
        for i in range(num_beams):
            member_id += 1
            self.members["Beam"+str(i)] = FrameMember([[x1,y],
                                                     [x2,y]],
                                                     mtype="beam",
                                                     length=(x2-x1),
                                                     mem_id=member_id,
                                                     f=self.f)
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
        member_id = self.generate_columns()
        self.generate_beams(member_id)
        self.num_members = len(self.members.keys())
      
        
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
                self.support_nodes.append(self.nodal_coordinates.index([x,y]))
            index += 1
        # Iterate through all frame's members and add created node-objects
        # to their member's own list
        for key in self.members.keys():
            self.members[key].add_nodes(self.nodal_coordinates)           
    
    
    def add_materials_and_sections(self):
        """ Adds materials and sections to a list in fem.FrameFEM()
        """
        # Iterate through all frame's members and calculate cross-sections
        # properties and add those and member's material properties to 
        # calculation model
        for key in self.members.keys():
            member = self.members[key]
            #if member.profile != 'Custom':
                #if member.mtype=="beam":
                    #member.steel_section_properties_IPE()
                #else:
                    #member.steel_section_properties_HEA()
                    
            member.add_material()
            member.add_section()
       
        
    def add_supports(self, supp_id=1):
        """ Adds supports          
            
            Parameters
            ----------
            supp_id : int
                number of support id
        """
        # Iterate through every support node and create supports
        for node in self.support_nodes:
            self.f.add_support(supp_id, node, [0,1,2], 0)
             
            
    def generate_elements(self):
        """ Generates elements for each frame's member
        """
        # Element indexing starts from zero
        index = 0
        # Iterate through every frame's member and generate elements for each
        # member
        for key in self.members.keys():
            self.members[key].generate_elements(index)
            index += self.num_elements
                    
            
    def calc_nodal_coordinates(self): 
        """ Calculates nodal coordinates and saves values to a list.
            These values are later used for creating the nodes.
        """
        # Iterate through every frame's member and calculate nodal coordinates.
        # If the x,y -coordinate is not in the list add that to the list
        for key in self.members.keys():
            self.members[key].calc_nodal_coordinates(self.num_elements)
            for coordinate in self.members[key].nodal_coordinates:
                if coordinate not in self.nodal_coordinates:
                    self.nodal_coordinates.append(coordinate)


    def calc_nodal_forces(self):
        """ Calculates nodal forces and saves values to
            self.nodal_forces - dict       
        """
        # Iterate through every frame's member's node and calculate forces
        # acting on that node
        for key in self.members.keys():
            self.members[key].calc_nodal_forces()
            for node in self.members[key].nodal_forces:
                self.nodal_forces[node] = self.members[key].nodal_forces[node]

   
    def calc_nodal_displacements(self):
        """ Calculates nodal displacements and saves values to
            self.nodal_displacements -dict
        """
        # Iterate through every frame's member's node and calculate 
        # displecements for that node
        for key in self.members.keys():
            member = self.members[key]
            member.calc_nodal_displacements()
            for node in member.nodal_displacements.keys():
                self.nodal_displacements[node]=member.nodal_displacements[node]

 
    def generate_frame(self):
        """ Generates the frame
        """
        # If the frame hasn't been created, create members and calculate nodal
        # coordinates, otherwise initialize frame.
        if self.is_generated == False:
            self.is_generated = True
            self.generate_members()
            self.calc_nodal_coordinates()
        else:
            self.initialize_frame()
        self.generate_nodes()  
        self.add_materials_and_sections()
        self.generate_elements()
        self.add_supports()
  
    
    def initialize_frame(self):
        """ Initializes the frame's properties
        """
        # Create a new FrameFEM -object and initialize frame's properties to
        # empty lists and dicts
        self.f = fem.FrameFEM()
        self.support_nodes = []
        self.nodes = []
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.r = []
        self.is_calculated = False
        # Initialize frame's members' properties
        for key in self.members.keys():
            member = self.members[key]
            member.elements = {}
            member.nodes = []
            member.nodal_forces = {}
            member.nodal_displacements = {}
            member.nodes = []
            member.f = self.f
            member.med = 0
            member.ved = 0
            member.ned = 0

    
    def draw_frame(self, print_text=True, draw=True):
        """ Draws the frame
            
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
        for key in self.members.keys():
            member = self.members[key]
            X = member.coordinates
            color = 'black'
            if member.has_load:
                color = 'blue'
            if self.is_calculated and print_text:
                if member.is_strong_enough:
                    color = 'green'
                else:
                    color = 'red'
            # Plot members
            plt.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color)
            # Calculate text location
            if member.mtype == "beam":
                y = X[0][1]
                x = member.length/5 + X[0][0]
            if member.mtype == "column":
                x = X[0][0]
                y = member.length/2 + X[0][1]
                node_id = member.nodes[-1]
                node_coord = member.nodal_coordinates[-1]
                # Plot nodes
                plt.plot(node_coord[0], node_coord[1], 'go')
                plt.text(node_coord[0], node_coord[1], str(node_id))
            # Plot text
            if print_text:
                plt.text(x, y, str(key) + ": " + member.profile)
        # Plot supports
        for node in self.support_nodes:
            node_coord = self.nodal_coordinates[node]
            plt.plot(node_coord[0], node_coord[1], 'rs')
            plt.text(node_coord[0], node_coord[1], str(node))
        if self.is_calculated == True:
            print("Frame weight: {} kg".format(self.weight))
        if draw:
            plt.show()
 
                  
    def draw_deflected_frame(self, scale=1, prec=4):
        """ Draws deflected shape of the frame
            
            Parameters
            ----------
            scale : float, optional
                Scaling factor
            prec : int, optional
                Precision for displacement value
        """
        self.draw_frame(print_text=False, draw=False)
        self.calc_nodal_displacements()
        for key in self.members.keys():
            X = []
            Y = []
            member = self.members[key]
            member.calc_nodal_displacements()
            max_x = 0
            max_y = 0
            for i in range(len(member.nodes)):
                node = member.nodes[i]
                x0 = member.nodal_coordinates[i][0]
                y0 = member.nodal_coordinates[i][1]               
                x1 = member.nodal_displacements[node][0]
                y1 = member.nodal_displacements[node][1]                
                x = x0 + x1*(scale)
                y = y0 + y1*(scale)                
                if abs(x1) > abs(max_x) and member.mtype == "column":
                    max_x = abs(x1)
                    loc_max_x = x
                    loc_max_y = y                    
                if abs(y1) > abs(max_y) and member.mtype == "beam":
                    max_y = abs(y1)
                    loc_max_x = x
                    loc_max_y = y                
                X.append(x)
                Y.append(y)
            plt.plot(X,Y,color='r')
            
            if member.mtype == "beam":
                plt.plot(loc_max_x, loc_max_y, 'ro')
                plt.text(loc_max_x, loc_max_y,
                         (str(max_y*1000)[0:prec+1] + " mm"))
            
            if member.mtype == "column":
                plt.plot(loc_max_x, loc_max_y, 'ro')
                plt.text(loc_max_x, loc_max_y,
                         (str(max_x*1000)[0:prec+1] + " mm"))

            
    def draw_bmd(self, scale=1):
        """ Draws deflected shape of the frame
            
            Parameters
            ----------
            scale : int, optional
                Scaling factor
        """
        self.draw_frame(print_text=False, draw=False)
        for key in self.members.keys():
            X = []
            Y = []
            member = self.members[key]
            member.calc_nodal_forces()
            X.append(member.nodal_coordinates[0][0])
            Y.append(member.nodal_coordinates[0][1])
            for i in range(len(member.nodes)):
                node = member.nodes[i]
                if member.mtype == "beam":
                    x0 = member.nodal_coordinates[i][0]
                    y0 = member.nodal_coordinates[i][1]                   
                    y1 = member.nodal_forces[node][2]/(1000/scale)                    
                    x = x0
                    y = y0 - y1                    
                    X.append(x)
                    Y.append(y)                    
                if member.mtype == "column":
                    x0 = member.nodal_coordinates[i][0]
                    y0 = member.nodal_coordinates[i][1]                    
                    x1 = member.nodal_forces[node][2]/(1000/scale)                    
                    x = x0 + x1
                    y = y0                    
                    X.append(x)
                    Y.append(y)
            X.append(member.nodal_coordinates[i][0])
            Y.append(member.nodal_coordinates[i][1])
            plt.plot(X,Y,color='gray')
 
    
    def add_point_load(self, load_id, node, load_vec, factor=1):
        """ Saves point load values to a dict
            
            Parameters
            ----------
            load_id : int
                Load id
            node : int
                Node subjected to load
            load_vec : 1D-array of float
                Load vector [Fx, Fy, Mz]
            factor : float
                Load scaling factor
        """        
        if load_id not in self.point_loads.keys():  
            self.point_loads[load_id] = []            
        load = [node, load_vec, factor]
        self.point_loads[load_id].append(load)

   
    def add_loads(self, load_id=2):
        """ Adds all loads to calculation model
            
            Parameters
            ----------
            load_id : int
                Load id
        """  
        # Add line loads to calculation model
        # If load id can't be found in self.loads, create a new key-value pair
        # and add loads with the load id to calculation model
        if load_id not in self.line_loads.keys():
            self.line_loads[load_id] = []
            for key in self.members.keys():
                member = self.members[key]
                if load_id != 'self_weight':
                    member.add_load(load_id)
                if load_id in member.loads.keys():
                    self.line_loads[load_id].append([key] +
                                                    member.loads[load_id])

        # Add loads from self.loads to calculation model
        else:
            for load in self.line_loads[load_id]:
                member = self.members[load[0]]
                member.add_load(load_id)
        # Add point loads to calculation model
        if load_id in self.point_loads.keys():
            for load in self.point_loads[load_id]:
                point_load = fem.PointLoad(load_id,
                                           self.f.nodes[load[0]],
                                           load[1],
                                           load[2])
                self.f.add_load(point_load)
        # load_id != 'self_weight' denies the possibility to add self weight
        # accidentally more than once to the calculation model
        if self.self_weight and load_id != 'self_weight':
            for array in self.line_loads['self_weight']:
                key = array[0]
                member = self.members[key]
                value = member.loads['self_weight'][0]
                direction = member.loads['self_weight'][1]
                for eid in member.elements.keys():
                    load = fem.LineLoad(load_id,                # Load id
                                        self.f.elements[eid],   # element
                                        [0.0,1.0],         # local coordinates
                                        [value[0], value[1]],   
                                        direction)
                    self.f.add_load(load)
                
               
    def calculate(self, load_id=2):
        """ Calculates forces and displacements
            
            Parameters
            ----------
            load_id : int
                Load id
        """
        if self.is_calculated == False:
            self.is_calculated = True
            self.add_loads(load_id)
            self.f.add_loadcase(supp_id=1, load_id=load_id)
            self.f.nodal_dofs()

        self.f.linear_statics()
        self.calc_nodal_forces()
        self.calc_nodal_displacements()
        self.check_members_strength()
        self.frame_weight()

        
    def design_members(self, symmetry=True):
        """ Desgins frame's members
            
            Parameters
            ----------
            symmetry : bool
                Set true to design members symmetrically
        """
        self.is_designed = True
        if symmetry:
            # Change beam's profile
            # value helps calculating the member inidices
            value = 0
            for _ in range(math.floor(self.bays/2)):   
                for i in range(self.storeys):
                    key1 = "Beam" + str(i + value)
                    key2 = "Beam" + str(i - value 
                                        + int(self.storeys*(self.bays-1)))
                    # 
                    member1 = self.members[key1]
                    member2 = self.members[key2]
                    # design members
                    member1.design_member()
                    member2.design_member()
                    # change weaker beam's profile to stonger beam's profile
                    if member1.profile_idx > member2.profile_idx:
                        member2.profile = member1.profile
                    else:
                        member1.profile = member2.profile                    
                value += self.storeys                     
            # If there're beams on the symmetry axis
            if self.bays % 2 != 0:
                for i in range(self.storeys):
                    key = "Beam" + str(i + value)
                    member = self.members[key]
                    member.design_member()
            # Change column's profile
            # value helps calculating the member inidices
            value = 0
            for j in range(math.ceil(self.bays/2)):     
                for i in range(self.storeys):
                    key1 = "Column" + str(i + value)
                    key2 = "Column" + str(i + self.storeys * (self.bays - j))
                    
                    member1 = self.members[key1]
                    member2 = self.members[key2]
                    # desing member
                    member1.design_member()
                    member2.design_member()
                    # change weaker column profile to stonger column profile
                    if member1.profile_idx > member2.profile_idx:
                        member2.profile = member1.profile
                    else:
                        member1.profile = member2.profile                       
                value += self.storeys                
            # If there're columns on the symmetry axis
            if self.bays % 2 == 0:
                for i in range(self.storeys):
                    key = "Column" + str(i + self.bays)
                    member = self.members[key]
                    member.design_member()
        # if symmetry = False
        else:
            for key in self.members.keys():
                member = self.members[key]
                member.design_member()
        
        self.check_members_strength()
        self.frame_weight()


    def frame_weight(self):
        """ Calculates frame's weight and saves it to self.weight
        """
        self.weight = 0
        for key in self.members.keys():
            member = self.members[key]
            #if member.profile != "Custom":
                #if member.mtype=="beam":
                    #member.steel_section_properties_IPE()
                #else:
                    #member.steel_section_properties_HEA()
                             
            self.weight += member.weight
            
            
    def add_self_weight(self):
        """ Adds self weight to frame's members
        """
        if self.self_weight == False:
            self.self_weight = True
            for key in self.members.keys():
                member = self.members[key]
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
        for key in self.members.keys():
            member = self.members[key]
            member.remove_self_weight()
        self.generate_frame()

    def design_frame(self, symmetry=True): 
        """ Designs frame's members
            Parameters
            ----------
            symmetry : bool
                Set true if frame is symmetrical
        """
        if not self.is_calculated:
            self.generate_frame()
            self.calculate()
        while not self.is_strong_enough:
            self.design_members(symmetry)
            self.generate_frame()   # generates the designed frame
            self.calculate()     # calculates new forces and displacements      
        
            
    def check_members_strength(self):
        """ Checks if members can bear their loads
        """
        self.is_strong_enough = True
        self.r.clear()
        for key in self.members.keys():
            member = self.members[key]
            
            if member.mtype=="beam":
                member.check_steel_section_IPE()
                self.r.append(member.r)
                if not member.is_strong_enough:
                    self.is_strong_enough = False
            else:
                member.check_steel_section_HEA()
                self.r.append(member.r)
                if not member.is_strong_enough:
                    self.is_strong_enough = False
                
    def find_closest_profile_IPE(self, AI):
        """ Finds closest matching profile using euclidian distance
            
            Parameters
            ----------
            AI : 1-D-array of float
                Contains optimized member's area (A) and second moment of 
                inertia about y-axis (Iy)
        """
        # HEA, HEB and IPE profiles
        #M = np.asarray(A_Iy_Iz_matrix)
        # IPE profiles only
        
    
        
        M = np.asarray(IPE_A_I_matrix)
        M = np.sqrt(np.einsum('ij,ij->i', M-AI, M-AI))
        index = np.argmin(M)
        #profile = PROFILES[index]
        profile = IPE_PROFILES[index]
    
        return profile
    
    def find_closest_profile_HEA(self, AI):
        """ Finds closest matching profile using euclidian distance
            
            Parameters
            ----------
            AI : 1-D-array of float
                Contains optimized member's area (A) and second moment of 
                inertia about y-axis (Iy)
        """
        # HEA, HEB and IPE profiles
        #M = np.asarray(A_Iy_Iz_matrix)
        # HEA profiles only
 
        M = np.asarray(HEA_A_I_matrix)
        M = np.sqrt(np.einsum('ij,ij->i', M-AI, M-AI))
        index = np.argmin(M)
        #profile = PROFILES[index]
        profile = HEA_PROFILES[index]
            
        return profile
        
    def continuous_to_discrete(self):
        """ Changes custom profile measurements to discrete profile
        """
        """
        # HEA, HEB and IPE profiles
        for key in self.members.keys():
            member = self.members[key]
            AIyIz = [member.A, member.I_y, member.I_z]
            member.profile = self.find_closest_profile(AIyIz)
        """
        # IPE profiles only
        for key in self.members.keys():
            member = self.members[key]
            AI = [member.cross_section.A, member.cross_section.I[0]]
            if member.mtype == 'beam':
                member.profile = self.find_closest_profile_IPE(AI)
                
            else:
                member.profile = self.find_closest_profile_HEA(AI)
                
        self.generate_frame()
        self.calculate()
        
        
        
        
        
        
        
        
    def costf_weight(self, X=[]):
        """ Cost function used for optimization
            
            Parameters
            ----------
            X : 1-D-array of float, optional
                List containing rotational stiffness values [Sj1, Sj2]
                and h and b values in order
                [alpha11, alpha12,..., alphan2, h1, b1, h2,...,bn]
            Return:
                weight : float
                    Weight of the members
        """ 
        # Total number of beams
        num_beams = self.bays * self.storeys
        # Number of symmetrical beams
        sym_beams = int(self.storeys*(self.bays -self.bays%2)/2)
        # Number of beams to be designed
        beams = num_beams - sym_beams
        # Number of connections
        joints = num_beams - 2*sym_beams
        # If list X is empty, create a list x with current frame's dimensions
        if len(X) == 0:
            X = self.initial_x0()
        # If joints aren't optimized, no need to change the list
        if not self.optimize_joints:
            x = X      
        else:           
            # remove joint rigidities from list
            x = X[joints:]
        
        
        rho = 7850e-9
        # IPE values
        K = [0.0256, 3.2831, 0.0155, 2.4921,  0.3486, 30.178, 0.0155, 2.4921]
       
        # HEA values
        Q = [0.0294, 5.7651, 0.014, 4.2949,  1.0323, 2.5368, 0.0781, 2,9206]
       
        
        weight = 0
        index = 0
      
        # weight of the beams
        for i in range(beams):
            #crec que esta malament perque aixi s agafa com h les b tambe
            
            
            h = x[i*2]
            bvect = x[(i*2)+1]
            
            t_f = K[0] * h + K[1]
            t_w = K[2] * h + K[3]
            b = K[4] * h + K[5]
            r = K[6] * h + K[7]
            
           
          
            # Area of I-profile
            A = 2.0*t_f*b + (h - 2.0*t_f)*t_w + (4.0 - math.pi)*r**2.0
            weight += A * rho * self.bay_length * 1000  # kg's
        x = x[beams+2:]
         
        # weight of the columns
        for i in range(int(len(x)/2)):
            h = x[i*2]
            bvect = x[(i*2)+1]
            t_f = Q[0] * h + Q[1]
            t_w = Q[2] * h + Q[3]
            
            if h<=270:
                b = Q[4] * h + Q[5]
                r = Q[6] * h + Q[7]
            else:
                b=300
                if h<=690:
                    r=27
                else:
                    r=30
           
            
            # Area of HEA-profile
            A = 2.0*t_f*b + (h - 2.0*t_f)*t_w + (4.0 - math.pi)*r**2.0
            weight += A * rho * self.storey_height * 1000 # kg's
        
        
        # Add penalty for every non-satisfactory constraint
        if self.optimizer == "DIF_EVOL":
            value = self.penalty_val
            a = np.asarray(self.constraints(X)) < 0
            a = sum(a)
            penalty = a * value
            return weight + penalty
        
        return weight
    
    
    def func2(self, x):
        print(len(x))
        for i, member in enumerate(self.members.values()):
            member.h = x[i*2]
            print(member.h)
        
        #self.generate_frame()
        self.calculate()
       
    def func(self, x):
        """ Function, that assigns x's-values to frame's members
            
            Parameters
            ----------
            x : 1-D-array of float
                Array containing rotational stifness values [alpha1, alpha2]
                and h and b values
        """ 
        num_beams = self.bays * self.storeys
        sym_beams = int(self.storeys*(self.bays -self.bays%2)/2)
        index = 0
        if self.optimize_joints:
            # Change beam's joint rigidities
            # Symmetrical beams
            i = 0
            for i in range(sym_beams):
                key = "Beam" + str(i)
                member = self.members[key]
                member.alpha1 = x[index]
                member.alpha2 = x[index+1]
                index += 2
            # Beams on symmetry axis  
            for j in range(num_beams - 2*sym_beams):
                key = "Beam" + str(i+j)
                member = self.members[key]
                member.alpha1 = x[index]
                member.alpha2 = x[index]
                index += 1            
        # Change beam's profile
        # Symmetrical beams
        value = 0
        for _ in range(math.floor(self.bays/2)):
            for i in range(self.storeys):
                key1 = "Beam" + str(i + value)
                key2 = "Beam" + str(i - value 
                                    + int(self.storeys*(self.bays-1)))              
                member1 = self.members[key1]
                member2 = self.members[key2]                
                #member1.profile = "Custom"
                #member2.profile = "Custom"  
                member1.h = x[index]
                member2.h = x[index]
                #member1.custom_section_properties_IPE([x[index], x[index+1]])
                #member2.custom_section_properties_IPE([x[index], x[index+1]])                
                # Symmetrical joint rigidities
                member2.alpha2 = member1.alpha1
                member2.alpha1 = member1.alpha2                
                member1.alpha_to_Sj()
                member2.alpha_to_Sj()
                index += 2
            value += self.storeys           
        # If there're beams on the symmetry axis
        if self.bays % 2 != 0:
            for i in range(self.storeys):
                key = "Beam" + str(i + value)
                member = self.members[key]
                #member.profile = "Custom"   
                member.h = x[index]
                #member.custom_section_properties_IPE([x[index], x[index+1]])                
                # Symmetrical joint rigidities
                member.alpha1 = member.alpha2
                member.alpha_to_Sj()
                index += 2                   
        # Change column's profile
        # Symmetrical columns
        value = 0
        for j in range(math.ceil(self.bays/2)):
            for i in range(self.storeys):
                key1 = "Column" + str(i + value)
                key2 = "Column" + str(i + self.storeys * (self.bays - j))                
                member1 = self.members[key1]
                member2 = self.members[key2]                
                #member1.profile = "Custom"
                #member2.profile = "Custom"  
                member1.h = x[index]
                member2.h = x[index]
                #member1.custom_section_properties_HEA([x[index], x[index+1]])
                #member2.custom_section_properties_HEA([x[index], x[index+1]])
                index += 2                
            value += self.storeys            
        # If there're columns on the symmetry axis
        # uses i from previous loop
        if self.bays % 2 == 0:
                for _ in range(self.storeys):
                    i += 1
                    key = "Column" + str(i)
                    member = self.members[key]
                    member.profile = "Custom"
                    member.h=x[index]
                    #member.custom_section_properties_HEA([x[index], x[index+1]])
                    index += 2            
        #self.generate_frame()
        self.calculate()


    def constraints(self, x=[]):
        """ Constraints for optimization
            Satisfactory when all values in the array >= 0
            
            Parameters
            ----------
            x : 1-D-array of float
                Array containing optimized members values
            x = [] in the function call allows to check these constraints
                without x-array
                
            Returns
            -------
            cons : 2-D-array of float
                Array of utilization ratio's
        """
        if len(x) != 0:
            self.func(x)
            #self.func2(x)
        # Sway constraint 
        top_node = self.num_elements * self.storeys
        height = self.storey_height * self.storeys
        max_sway = height / 400 # SWAY LIMIT
        sway = self.nodal_displacements[top_node][0]
        cons1 = max_sway - sway        
        # Beam deflection
        cons2 =[]
        num_beams = self.bays * self.storeys
        for i in range(num_beams):
            member = self.members["Beam" + str(i)]
            max_v = 0
            for node in member.nodal_displacements.keys():
                v_val = member.nodal_displacements[node][1]
                if abs(v_val) > max_v:
                    max_v = abs(v_val)
            cons2.append(max_v)
        cons2 = np.asarray(cons2)
        max_deflection = self.bay_length / 300 # DEFLECTION LIMIT
        cons2 = max_deflection - cons2             
        cons = [cons1] + list(cons2)
        # Cross-section, buckling and lateral-buckling strength
        for key in self.members.keys():
            member = self.members[key]
            for r in member.r:
                cons.append(1-r)                           
        return cons
    
    
    
    def log_optimization(self, parameters, result, time):
        """ Function for logging optimization results
        
            Parameters
            ----------
            parameters : dict
                Parameters used to get optimization result
            result : float
                Optimized frame's weight
            time : float
                Time elapsed for optimization
        """
        file_name = self.optimizer + "_optimization_log.txt"
        print("Optimization results logged to {}".format(file_name))
        with open(file_name, 'a') as f:          
            f.write("Frame properties: \n")
            f.write("    Bays: {} \n".format(self.bays))
            f.write("    Bay length: {} \n".format(self.bay_length))
            f.write("    Storeys: {} \n".format(self.storeys))
            f.write("    Storey height: {} \n".format(self.storey_height))
            f.write("    Elements: {} \n".format(self.num_elements))
            if self.line_loads != {}:
                f.write("Loads: \n")
                for loadid in self.line_loads.keys():
                    f.write("    LoadID: {} \n".format(loadid))
                    f.write("    Line loads: \n")
                    for loads in self.line_loads[loadid]:
                        member = loads[0]
                        q_value = loads[1]
                        direction = loads[2]
                        f.write("        Member: {}\n".format(member))
                        f.write("        q value [kN/m]: {}\n".format(q_value))
                        f.write("        Direction: {}\n".format(direction))
                    if self.point_loads[loadid]:
                        f.write("    Point loads: \n")
                        for point_loads in self.point_loads[loadid]:
                            node = point_loads[0]
                            value = point_loads[1]
                            factor = point_loads[2]
                            f.write("        Node: {}\n".format(node))
                            f.write("        Value [Fx,Fy,Mz] [kN]: {}\n"
                                    .format(value))
                            f.write("        Factor: {}\n".format(factor))
            f.write("Optimizer : {}\n".format(self.optimizer))
            if not self.optimize_joints:
                f.write("Joint rigidities NOT optimized")
            f.write("Parameters: \n")
            for key in parameters.keys():
                f.write("    {} : {} \n".format(key, parameters[key]))
            f.write("Optimization result: {0:.2f} kg \n".format(result))
            f.write("Optimization result is feasible: {}\n"
                    .format(((min(self.constraints()) > 0))))
            f.write("    MAX r: {:.3f}% \n"
                    .format((1-min(self.constraints()))*100))
            f.write("Optimized profiles: \n")
            for key in self.members.keys():
                member = self.members[key]
                f.write("    {} : {} r:{:.2f} %\n"
                        .format(key,member.profile, (max(member.r)*100)))
                if member.mtype == 'beam':
                    f.write("    Sj1 : {:.2f}  Sj2 : {:.2f} \n"
                            .format(member.Sj1, member.Sj2))
            f.write("Time elapsed: {0:.2f} s \n".format(time))
            f.write("{}\n".format('-'*81))
    
    
    def initial_x0(self):
        """ Creates initial guess x0 from current frame's member's properties
        """        
        num_beams = self.bays * self.storeys
        sym_beams = int(self.storeys*(self.bays -self.bays%2)/2)
        x0 = []
        index = 0
        # Append beam's joint rigidities
        # Symmetrical beams
        i = 0
        for i in range(sym_beams):
            key = "Beam" + str(i)
            member = self.members[key]
            x0.append(member.alpha1)
            x0.append(member.alpha2)
            index += 2
        # Beams on symmetry axis  
        for j in range(num_beams - 2*sym_beams):
            key = "Beam" + str(i+j)
            member = self.members[key]
            x0.append(member.alpha1)
            index += 1
        # Append beam's h and b
        # Symmetrical beams
        value = 0
        for _ in range(math.floor(self.bays/2)):
            for i in range(self.storeys):
                key = "Beam" + str(i + value)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2
            value += self.storeys            
        # If there're beams on the symmetry axis
        if self.bays % 2 != 0:
            for i in range(self.storeys):
                key = "Beam" + str(i + value)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2
        # Append column's h and b
        # Symmetrical columns
        value = 0
        for _ in range(math.ceil(self.bays/2)):
            for i in range(self.storeys):
                key = "Column" + str(i + value)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2               
            value += self.storeys           
        # If there're columns on the symmetry axis
        if self.bays % 2 == 0:
            for i in range(self.storeys):
                key = "Column" + str(i + self.bays)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2   
        return x0
    
 
      
        
    def optimize_frame(self, debug=True, print_data=False, optimizer="",
                       swarmsize=50, maxiter=20,omega=0.5, phip=0.5,
                       phig=0.5, penalty_val=1000, joints=False, log_data=False,
                       draw=True, lb=[80, 46], ub=[500,246]):
        """ Optimizes frame's members and beam's joints rotational stifness
            
            Parameters
            ----------
            debug : bool (Default: True)
                boolean value, if True prints optimization data 
            print_data : bool (Default: False)
                boolean value, if True prints members' info
            optimizer : string
                Optimization algorithm used to optimize the frame
                    "PSO" -- Particle Swarm optimization
                    "SLSQP" -- Sequential least squares
                    "DIF-EVOL" -- Differential evolution
            swarmsize : int (Default: 50)
                PSO number of particles used for optimization
            maxiter : int (Default: 20)
                number of maximum iterations
            omega : float (Default: 0.5)
                PSO Particle velocity scaling factor 
            phip : float (Default: 0.5)
                PSO Scaling factor to search away from the particle's best
                known position 
            phig : float (Default: 0.5)
                PSO Scaling factor to search away from the swarm's
                best known position
            joints : bool (Default: True)
                boolean value, if True optimizes joint rigidities
            penalty_val : float (Default: 1000)
                DIF-EVOL Penalty value that's added to cost function's value 
                for every non-satisfactory constraint
            log_data : bool (Default: False)
                Logs optimization data
            draw : bool (Default: True)
                Draws frame after optimization
            lb : 1-D-array of float (Default: [80, 46])
                array of lower boundary values [h, b] [mm, mm]
            ub : 1-D-array of float (Default: [500, 246])
                array of upper boundary values [h, b] [mm, mm]               
        """
        self.optimize_joints = joints
        # If no optimizer is given in the fuction call, use last used optimizer
        # Also checks if given optimizer is in the optimizer list and raises
        # error if not, otherwise changes frame's optimizer to that algorithm
        if optimizer == "":
            optimizer = self.optimizer
         
        else:
            optimizer = optimizer.upper()
            if optimizer not in OPTIMIZERS:
                raise ValueError("There's no optimizer called {}"
                                 .format(optimizer))
            else:
                self.optimizer = optimizer
           
        # Set penalty value
        if self.penalty_val != penalty_val:
            self.penalty_val = penalty_val
        # Set boundary arrays
        #if len(ub) != 2 or len(lb) != 2:
        #   raise ValueError("Boundary arrays must have exactly two values.")           
        # Set true for printing optimized values
        self.is_optimized = True
        # Start timer
        start = timer()
        # Number of beams and symmetrical beams
        num_beams = self.bays * self.storeys
        sym_beams = self.storeys*(self.bays -self.bays%2)/2       
        # Objective function
        obj_func = self.costf_weight
        # Constraint function
        cons = self.constraints        
        # Symmetrical beams
        # Joint rigidity boundaries according to SFS-EN 1993-1-8
        lb_j = [0.5]*2*int(sym_beams)
        ub_j = [25]*2*int(sym_beams)        
        # Beams on symmetry axis
        lb_j += [0.5]*1*(num_beams - 2*int(sym_beams))
        ub_j += [25]*1*(num_beams - 2*int(sym_beams))
        if not joints:
            lb_j = []
            ub_j = []
        # Symmetrical members
        # [h, b]
        num_sym_members = self.bays*self.storeys + self.storeys
        lb = lb_j + lb*num_sym_members
        ub = ub_j + ub*num_sym_members  
        # Particle Swarm Optimization
        if optimizer == "PSO": 
            """
            # Designs frame and uses those values as initial guess
            # also scales upper and lower boundary values according to
            # designed members
            x0 = np.asarray(self.initial_x0())
            lb_x0 = x0[len(lb_j):]*0.5
            ub_x0 = x0[len(ub_j):]*1.5
            lb = lb_j + list(lb_x0)
            ub = ub_j + list(ub_x0)
            """           
            xopt, fopt = pso(obj_func,
                             lb,
                             ub,
                             swarmsize=swarmsize,
                             maxiter=maxiter,
                             f_ieqcons=cons,
                             debug=debug,
                             omega=omega,
                             phip=phip,
                             phig=phig)
            if xopt == []:
                return
            self.func(xopt)
            # Parameters saved to dict for logger
            parameters = {"swarmsize": swarmsize,
                          "maxiter": maxiter,
                          "omega": omega,
                          "phip": phip,
                          "phig": phig,
                          "result": xopt}
        # Sequential Least SQuares 
        if optimizer == "SLSQP":     
            """
            # Designs frame and uses those values as initial guess
            # also scales upper and lower boundary values according to
            # designed members
            x0 = np.asarray(self.initial_x0())
            lb_x0 = x0[len(lb_j):]*0.5
            ub_x0 = x0[len(ub_j):]*1.5
            lb = lb_j + list(lb_x0)
            ub = ub_j + list(ub_x0)
            """
            bounds = []
            for i in range(len(lb)):
                bounds.append((lb[i], ub[i]))
            # initial guess
            x0 = np.random.rand(len(ub)) * np.asarray(ub)
            #x0 = self.initial_x0()           
            if debug:
                print("x0: ",x0)
                iprint = 2
            else:
                iprint = 0
            out, fx, its, imode, smode = fmin_slsqp(obj_func,
                                                    x0,
                                                    f_ieqcons=cons,
                                                    bounds=bounds,
                                                    iprint=iprint,
                                                    iter=maxiter,
                                                    full_output=True)
            self.constraints(out)
            parameters = {
                    "x0": x0,
                    "iterations": its,
                    "maxiter": maxiter,
                    "exit mode": smode,
                    "result": out
                    }
        # Differential evolution
        if optimizer == "DIF_EVOL":
            """
            x0 = np.asarray(self.initial_x0())
            lb_x0 = x0[len(lb_j):]-100
            ub_x0 = x0[len(ub_j):]+100
            lb = lb_j + list(lb_x0)
            ub = ub_j + list(ub_x0)
            """
            x0 = np.random.rand(len(ub)) * np.asarray(ub)
            bounds = []
            
            if debug:
                print("x0: ",x0)
                disp = True
            else:
                disp = False
            
            for i in range(len(lb)):
                bounds.append((lb[i], ub[i]))
            result = differential_evolution(obj_func,
                                            bounds,
                                            disp=disp,
                                            maxiter=maxiter)
            self.constraints(result.x)
            parameters = {
                    "x0": x0,
                    "maxiter": maxiter,
                    "penalty": self.penalty_val
                    }            
        # Change optimized continuous members to discrete
        self.continuous_to_discrete() 
        # End timer
        end = timer()
        print("Time elapsed: ", end - start, "s")
        if log_data:
            self.log_optimization(parameters, self.weight, end - start)       
        if print_data:
            self.print_member_info()
        if draw:
            self.draw_frame()
             
    def print_member_info(self):
        """ Prints frame's members' info 
        """
        for key in self.members.keys():
            member = self.members[key]
            print()
            print(key)
            if member.mtype == "beam":
                if self.is_optimized:
                    print("Optimal joint rigidities: ")
                print("Sj1 = {:0.2f}, Sj2 = {:0.2f}"
                      .format(member.Sj1, member.Sj2))
            print("MAX Utilization ratio: {:0.2f} %"
                  .format((max(member.r)*100)))
            
            
    def rigid_joints(self):
        """ Set all beam-column joints to rigid
        """
        for key in self.members.keys():
            member = self.members[key]
            if member.mtype == "beam":
                member.alpha1 = np.inf
                member.alpha2 = np.inf
                member.Sj1 = np.inf
                member.Sj2 = np.inf
        self.generate_frame()
        
    def hinge_joints(self):
        """ Set all beam-column joints to hinges
        """
        for key in self.members.keys():
            member = self.members[key]
            if member.mtype == "beam":
                member.alpha1 = 0
                member.alpha2 = 0
                member.Sj1 = 0
                member.Sj2 = 0
        self.generate_frame()

#-----------------------------------------------------------------------------   
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
        SemiRigidFrame
        FrameFEM
        I_sections
        SteelSection
        SteelMember
            
    """
    
    def __init__(self, coordinates, mtype, length, mem_id, f):
        
        self.f = f
        self.elements = {}
        self.nodes = []
        self.nodal_coordinates = []
        self.loc = []
        self.coordinates = coordinates
        self.material = "S355"
        self.profile = "IPE 100"
       
        if mtype == 'beam': 
            self.cross_section = IPE(100, 355)
        elif mtype == 'column':
            self.cross_section= HEA(100, 355)
        
        
        self.profile_idx = PROFILES.index(self.profile)
        
        self.length = length            
        self.mtype = mtype
        self.mem_id = mem_id
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.loads = {}
        self.alpha1 = 25
        self.alpha2 = 25
        self.Sj1 = np.inf
        self.Sj2 = np.inf
        self.ned = 0                 
        self.ved = 0                    
        self.med = 0               
        self.has_load = False
        self.is_designed = False
        self.is_strong_enough = False
        self.r = 1
        self.self_weight = False
    
    @property
    def weight (self):
        weight = self.cross_section.A * self.length * self.rho
        return weight
    
    
    @property
    def h(self):
        return self.__h
    
    @h.setter
    def h(self, val):
        self.__h = val
        if self.mtype == 'beam': 
            self.cross_section = IPE(val, self.fy, False)
        elif self.mtype == 'column':
            self.cross_section = HEA(val, self.fy, False)
            
        if len(self.elements) > 0:
            for element in self.elements.values():
                element.section.A = self.cross_section.A * 1e-6
                element.section.Iy = self.cross_section.I[0] * 1e-12
    @property
    def material(self):
        return self.__material

    @material.setter
    def material(self, val):
        self.__material = val
        self.fy = tables_and_tuples.mat[self.material]['f_y']
        self.E = tables_and_tuples.mat[self.material]['E']
        self.fu = tables_and_tuples.mat[self.material]['f_u']
        self.nu = tables_and_tuples.mat[self.material]['v']
        self.rho = tables_and_tuples.mat[self.material]['rho']        

    @property
    def profile(self):
        return self.__profile

    @profile.setter
    def profile(self, val):
        self.__profile = val.upper()
        profiles = self.profile.split(" ")
        profile_type = profiles[0]
        kwargs = tables_and_tuples.profile[self.profile]

        if profile_type == 'IPE':
            self.cross_section = IPE(int(profiles[1]), self.fy)
        elif profile_type == 'RHS':
            pass
        elif profile_type == 'CHS':
            pass
        elif profile_type == 'HE':
            if profiles[2] == 'A':
                self.cross_section = HEA(int(profiles[1]), self.fy)
            pass
        elif profile_type == 'SHS':
            pass
        else:
            raise ValueError('{} is not valid profile type!'.format(profile_type))

        if len(self.elements) > 0:
            for element in self.elements.values():
                element.section.A = self.cross_section.A * 1e-6
                element.section.Iy = self.cross_section.I[0] * 1e-12
    
    def alpha_to_Sj(self):
        """ Change alpha value to Sj
        """
        # SFS-EN 1993-1-8 page 60, hinge 0.5 < alpha < 25 rigid
        self.Sj1 = self.alpha1 * self.E * self.cross_section.I[0]*1e-12 / self.length
        self.Sj2 = self.alpha2 * self.E * self.cross_section.I[0]*1e-12 / self.length
        
        """fixity factor alpha from 96semi_rigid.pdf, page 2  
        if self.alpha1 == 1:
            self.Sj1 = fem.kRigid
        elif self.alpha1 == 0:
            self.Sj1 = fem.kHinge
        else:    
            self.Sj1 = - self.alpha1*3*self.E*self.I_y*1e-12/
                        (self.length*(self.alpha1-1))
            
        if self.alpha2 == 1:
            self.Sj2 = fem.kRigid
        elif self.alpha2 == 0:
            self.Sj2 = fem.kHinge
        else:   
            self.Sj2 = - self.alpha2*3*self.E*self.I_y*1e-12/
                        (self.length*(self.alpha2-1))
        """   
    def Sj_to_alpha(self):
        """ Change Sj value to alpha
        """
        # SFS-EN 1993-1-8 alpha
        self.alpha1 = self.Sj1 * self.length / (self.E * self.I_y*1e-12)
        self.alpha2 = self.Sj2 * self.length / (self.E * self.I_y*1e-12)
        """ fixity factor alpha
        self.alpha1 = 1 / (1 + 3 * self.E * self.I_y*1e-12 / 
                           (self.Sj1*self.length))
        self.alpha2 = 1 / (1 + 3 * self.E * self.I_y*1e-12 / 
                           (self.Sj2*self.length))
        """
    def calc_nodal_coordinates(self, num_elements):
        """
            Calculates node locations along member
            Coordinates used are global coordinates
            Adds locations to a list where they can be accessed later
        """
        start_node = self.coordinates[0]
        end_node = self.coordinates[1]

        x = start_node[0]
        y = start_node[1]
        
        for i in range(num_elements):

            self.nodal_coordinates.append([x, y])
        
            x += (end_node[0] - start_node[0]) / num_elements
            y += (end_node[1] - start_node[1]) / num_elements
            self.loc.append(i/num_elements)

        self.nodal_coordinates.append([end_node[0], end_node[1]])
        self.loc.append(1)
        
        
    def add_nodes(self, frame_nodal_coordinates):
        """ Creates nodes to previously calculated locations
        """
        for coordinate in self.nodal_coordinates:
            self.nodes.append(frame_nodal_coordinates.index(coordinate))
        
    def generate_elements(self, index):
        """ Generates elements between nodes
            For beam member's first and last element, creates semi-rigid end
            elements. Elements are added to a list
        """ 
        # EBBeam -elements
        if self.mtype == "column":
            for i in range(len(self.nodes)-1):
                n1 = self.nodes[i]
                n2 = self.nodes[i+1]
                self.elements[index] =\
                fem.EBBeam(self.f.nodes[n1], self.f.nodes[n2],\
                                    self.f.sections[self.mem_id],\
                                    self.f.materials[self.mem_id])
                
                self.f.add_element(self.elements[index])
                index += 1                    
        # EBBeam -elements
        if self.mtype == "beam":
            for i in range(len(self.nodes)-1):
                n1 = self.nodes[i]
                n2 = self.nodes[i+1]
                # EBSemiRigid -element
                if i == 0:
                    self.elements[index] =\
                fem.EBSemiRigidBeam(self.f.nodes[n1], self.f.nodes[n2],\
                                    self.f.sections[self.mem_id],\
                                    self.f.materials[self.mem_id],\
                                    rot_stiff=[self.Sj1, np.inf])
                # EBSemiRigid -element
                elif i == (len(self.nodes)-2):
                    self.elements[index] =\
                fem.EBSemiRigidBeam(self.f.nodes[n1], self.f.nodes[n2],\
                                    self.f.sections[self.mem_id],\
                                    self.f.materials[self.mem_id], \
                                    rot_stiff=[np.inf, self.Sj2])
                
                else:
                    self.elements[index] =\
                fem.EBBeam(self.f.nodes[n1], self.f.nodes[n2],\
                                    self.f.sections[self.mem_id],\
                                    self.f.materials[self.mem_id])
                
                self.f.add_element(self.elements[index])
                index += 1
        

    def calc_nodal_forces(self):
        """ Calculates nodal forces on member's nodes          
        """
        
        i = 0
        for key in self.elements.keys():
            element = self.elements[key]
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
               
    def calc_nodal_displacements(self):
        """ Calculates nodal displacements and saves them to a dict
        """
        for node in self.nodes:
            self.nodal_displacements[node] = self.f.nodes[node].u


    def bmd(self):
        """ Plots member's bending moment diagram
            Plots currently only horizontal diagrams
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
                i += self.length/len(self.elements)
            x2 = [0]*len(x1)
            max_val = max(y1)
            min_val = min(y1)
            if abs(max_val) > abs(min_val) and max_val != y1[0] \
                                            and max_val != y1[-1]:
                val = min_val
            else:
                val = max_val
            max_loc = int(y1.index(val))
            plt.gca().invert_yaxis()
            plt.fill_between(x1,y1, color='lightgray')
            plt.plot(x1,y1)
            plt.plot(x1,x2, color='black')
            plt.text(x1[0],y1[0], str(y1[0])+ "kNm")
            plt.text(x1[-1],y1[-1], str(y1[-1]) + "kNm")
            plt.text(x1[max_loc],y1[max_loc], str(val) + "kNm")
        except ValueError:
            print("Error! Calculate results first.")
       
    def add_line_load(self, load_id, value, direction):
        """
            Adds line load to a member, but not to the calculation model
            Add_load -function adds load to calculation model
            load_id -- load ID
            value -- kN/m
            direction -- 'x' or 'y'
        """
        if not isinstance(value, list):
            value = [value, value]
        self.loads[load_id] = [value, direction]   
        self.has_load = True
        
    def add_load(self, load_id):
        """ Adds load to the calculation model
        """
        if load_id in self.loads.keys():
            value = self.loads[load_id][0]
            direction = self.loads[load_id][1]
            for eid in self.elements.keys():
                load = fem.LineLoad(load_id,                # Load id
                                    self.f.elements[eid],   # element
                                    [0.0,1.0],              # local coordinates
                                    [value[0], value[1]],   
                                    direction)
                self.f.add_load(load)
    
    def add_self_weight(self):
        """ Adds self-weight to the member's loads
        """
        if not self.self_weight:
            self.self_weight = True
            load_id = "self_weight"
            # self.weight is kg's, multiplier changes it to kN's
            multiplier = 1e-2
            value = -1 * multiplier *self.weight / self.length
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



    """
    def custom_section_properties_IPE(self, x):
    
        Makes a custom section and calculates it's properties
            
            Parameters
            ----------
            x : 1-D-array of float
                Array containing h and b values
    """
        
        # IPE values
        #K = [0.0256, 3.2831, 0.0155, 2.4921, 0.3486, 30.178, 0.0155, 2.4921]
       
        # HEA values
        #Q = [0.0294, 5.7651,0.014, 4.2949,  1.0323, 2.5368, 0.0781, 2,9206]
       
    """  
        self.h = x[0]
       
        if len(x) == 2:
            
            self.t_f = K[0]*self.h + K[1]
            self.t_w = K[2]*self.h + K[3]
            self.b = K[4] * self.h + K[5]
            self.round = K[6] * self.h + K[7]
            
        elif len(x) == 5:
            self.t_f = x[2]
            self.t_w = x[3]
            self.round = x[4]
                
        beam = I_sections.Beam(h=self.h,
                               b=self.b,
                               t_f=self.t_f,
                               t_w=self.t_w,
                               r=self.round, # r-value of smallest IPE profile
                               material=self.material,
                               L=self.length)
        
        self.weight = beam.g * self.length * 1000    # kg
        self.A = beam.A                              # mm^2
        self.A_v = beam.A_v                          # mm^2
        self.I_y = beam.I_y                          # mm^4
        self.I_z = beam.I_z                          # mm^4
        self.I_t = beam.I_t                          # mm^4
        self.I_w = beam.I_w                          # mm^6
        self.W_pl_y = beam.W_pl_y                    # mm^3
        self.W_pl_z = beam.W_pl_z                    # mm^3
        self.W_el_y = beam.W_el_y                    # mm^3
        self.W_el_z = beam.W_el_z                    # mm^3
        self.E = beam.E*1e3                          # kN/m^2
        self.nu = beam.v                             # dimensionless
        self.rho = beam.rho*1e9                      # kg/m^3
        self.f_y = beam.f_y
            
        beam.capacity()
        self.N_Rd = beam.N_Rd                        # N
        self.V_y_Rd = beam.V_y_Rd                    # N
        self.M_Rd_y = beam.M_Rd_y / 1e3              # Nm
        self.M_Rd_z = beam.M_Rd_z / 1e3              # Nm
     
    def custom_section_properties_HEA(self, x):
     
    Makes a custom section and calculates it's properties
            
            Parameters
            ----------
            x : 1-D-array of float
                Array containing h and b values
    """
        # IPE values
       # K = [0.0155, 2.4921, 0.0256, 3.2831, 0.3486, 30.178, 0.0155, 2.4921]
       
        # HEA values
        #Q = [0.0294, 5.7651, 0.014, 4.2949,  1.0323, 2.5368, 0.0781, 2,9206]
       
    """  
        self.h = x[0]
       
        if len(x) == 2:
            #jo posaria igual a 1 perque ara la b es troba a partir de h
            self.t_f = Q[0]*self.h + Q[1]
            self.t_w = Q[2]*self.h + Q[3]
            self.b = Q[4] * self.h + Q[5]
            self.round = Q[6] * self.h + Q[7]
            
        elif len(x) == 5:
            self.t_f = x[2]
            self.t_w = x[3]
            self.round = x[4]
                
        column = anna_HEA_sections.Beam(h=self.h,
                               b=self.b,
                               t_f=self.t_f,
                               t_w=self.t_w,
                               r=self.round, # r-value of smallest IPE profile
                               material=self.material,
                               L=self.length)
        
        self.weight = column.g * self.length * 1000    # kg
        self.A = column.A                              # mm^2
        self.A_v = column.A_v                          # mm^2
        self.I_y = column.I_y                          # mm^4
        self.I_z = column.I_z                          # mm^4
        self.I_t = column.I_t                          # mm^4
        self.I_w = column.I_w                          # mm^6
        self.W_pl_y = column.W_pl_y                    # mm^3
        self.W_pl_z = column.W_pl_z                    # mm^3
        self.W_el_y = column.W_el_y                    # mm^3
        self.W_el_z = column.W_el_z                    # mm^3
        self.E = column.E*1e3                          # kN/m^2
        self.nu = column.v                             # dimensionless
        self.rho = column.rho*1e9                      # kg/m^3
        self.f_y = column.f_y
            
        column.capacity()
        self.N_Rd = column.N_Rd                        # N
        self.V_y_Rd = column.V_y_Rd                    # N
        self.M_Rd_y = column.M_Rd_y / 1e3              # Nm
        self.M_Rd_z = column.M_Rd_z / 1e3              # Nm
    """



    """
    def steel_section_properties_IPE(self):
        Calculates cross-section properties 
        
        beam = I_sections.Beam(size=self.profile,
                               material=self.material,
                               L=self.length)
        self.weight = beam.g * self.length * 1000    # kg
        self.h = beam.h                              # mm
        self.b = beam.b                              # mm
        self.t_f = beam.t_f                          # mm
        self.A = beam.A                              # mm^2
        self.A_v = beam.A_v                          # mm^2
        self.I_y = beam.I_y                          # mm^4
        self.I_z = beam.I_z                          # mm^4
        self.I_t = beam.I_t                          # mm^4
        self.I_w = beam.I_w                          # mm^6
        self.W_pl_y = beam.W_pl_y                    # mm^3
        self.W_pl_z = beam.W_pl_z                    # mm^3
        self.W_el_y = beam.W_el_y                    # mm^3
        self.W_el_z = beam.W_el_z                    # mm^3
        self.E = beam.E*1e3                          # kN/m^2
        self.nu = beam.v                             # dimensionless
        self.rho = beam.rho*1e9                      # kg/m^3
        self.f_y = beam.f_y
            
        beam.capacity()
        self.N_Rd = beam.N_Rd                        # N
        self.V_y_Rd = beam.V_y_Rd                    # N
        self.M_Rd_y = beam.M_Rd_y / 1e3              # Nm
        self.M_Rd_z = beam.M_Rd_z / 1e3              # Nm

    
    
    
    
    
    def steel_section_properties_HEA(self):
        Calculates cross-section properties 
        
        column = anna_HEA_sections.Beam(size=self.profile,
                               material=self.material,
                               L=self.length)
        self.weight = column.g * self.length * 1000    # kg
        self.h = column.h                              # mm
        self.b = column.b                              # mm
        self.t_f = column.t_f                          # mm
        self.A = column.A                              # mm^2
        self.A_v = column.A_v                          # mm^2
        self.I_y = column.I_y                          # mm^4
        self.I_z = column.I_z                          # mm^4
        self.I_t = column.I_t                          # mm^4
        self.I_w = column.I_w                          # mm^6
        self.W_pl_y = column.W_pl_y                    # mm^3
        self.W_pl_z = column.W_pl_z                    # mm^3
        self.W_el_y = column.W_el_y                    # mm^3
        self.W_el_z = column.W_el_z                    # mm^3
        self.E = column.E*1e3                          # kN/m^2
        self.nu = column.v                             # dimensionless
        self.rho = column.rho*1e9                      # kg/m^3
        self.f_y = column.f_y
            
        column.capacity()
        self.N_Rd = column.N_Rd                        # N
        self.V_y_Rd = column.V_y_Rd                    # N
        self.M_Rd_y = column.M_Rd_y / 1e3              # Nm
        self.M_Rd_z = column.M_Rd_z / 1e3              # Nm
    """
    
    def add_material(self):
        """ Adds member's material information to calculation model
            Young's modulus kN/m^2
            nu is dimensionless
            Density kg/m^3
        """
        self.f.add_material(self.E, self.nu, self.rho)
        
    def add_section(self):
        """ Units: m**2 A and I_y are in mm**2 and mm**4, respectively"""
        s = fem.BeamSection(self.cross_section.A*1e-6, self.cross_section.I[0]*1e-12)
        self.f.add_section(s)
        
    def check_steel_section_IPE(self):
        """ Checks if cross-section is strong enough.
            Gives boolean value for self.is_strong_enough
            Saves list of stress ratios to self.r
        """
        #if self.profile != "Custom":
        #    self.steel_section_properties_IPE()
        self.r = []
          
        #steel_section = IPE(100,250,False)
        """ The old way
        steel_section = IPE(self.h,355,False)
        steel_member = SteelMember(steel_section,self.length,
                                   Lcr= [1.0,1.0],mtype=self.mtype)
        """
        # The new way
        steel_member = SteelMember(self.cross_section,self.length,
                                   Lcr= [1.0,1.0],mtype=self.mtype)
        i = 0
        for node in self.nodal_forces.keys():
            forces = self.nodal_forces[node]
            ned = forces[0] *1e3        # kN to N
            ved = forces[1] *1e3        # kN to N
            med = forces[2] *1e6        # kNm to Nmm
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
       
        
    def check_steel_section_HEA(self):
        """ Checks if cross-section is strong enough.
            Gives boolean value for self.is_strong_enough
            Saves list of stress ratios to self.r
        """
        #if self.profile != "Custom":
        #    self.steel_section_properties_HEA()
        self.r = []
        #steel_section = HEA(100,250,False)
        """ The old way
        steel_section = HEA(self.h,355,False)
        steel_member = SteelMember(steel_section,self.length,
                                   Lcr= [1.0,1.0],mtype=self.mtype)
        """
        # The new way
        steel_member = SteelMember(self.cross_section,self.length,
                                   Lcr= [1.0,1.0],mtype=self.mtype)
        i = 0
        for node in self.nodal_forces.keys():
            forces = self.nodal_forces[node]
            ned = forces[0] *1e3        # kN to N
            ved = forces[1] *1e3        # kN to N
            med = forces[2] *1e6        # kNm to Nmm
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
        i = max(0, self.profile_idx-5)
        self.is_strong_enough = False
        while(self.is_strong_enough == False):
            i += 1
            self.profile = PROFILES[i]
            self.check_steel_section()
            self.profile_idx = i

#-----------------------------------------------------------------------------
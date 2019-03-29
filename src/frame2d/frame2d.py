"""
@author: Jaakko Huusko
"""
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import framefem.framefem  as fem

from framefem.elements import EBBeam, EBSemiRigidBeam
from sections.steel.ISection import IPE, HEA, HEAA, HEB, HEC, HEM, CustomISection
from sections.steel.RHS import RHS, SHS
from sections.steel.CHS import CHS
from sections.steel.catalogue import mat as MATERIALS
from sections.steel.catalogue import ipe_profiles, h_profiles, rhs_profiles, shs_profiles, chs_profiles
from structures.steel.steel_member import SteelMember

# Eounding precision
PREC = 3



# -----------------------------------------------------------------------------
class Frame2D:
    """ Class for 2D frame

        Parameters:
        ------------
            :param simple: creates a simple frame with given list [storeys, bays, storey height, bay length]
            :param num_elements: number of elements per member
            :param supports: Generates supports to simple frame 'FIXED', 'XHINGED', 'YHINGED', 'XYHINGED'
            :param fem: FrameFEM -instance, used to share fem model with other Frame2D -instances
                        default: creates a new  FrameFEM -instance

            :type simple: list
            :type num_elements: int
            :type supports: string
            :type fem: FrameFEM


            Variables:
            ----------
            :ivar nodes:
            :ivar nodal_coordinates:
            :ivar num_elements:
            :ivar f:
            :ivar alpha_cr:
            :ivar num_members:
            :ivar support_nodes:
            :ivar point_loads:
            :ivar line_loads:
            :ivar nodal_foces:
            :ivar nodal_displacements:
            :ivar joints:
            :ivar r:
            :ivar is_generated:
            :ivar is_calculated:
            :ivar self_weight:
            :ivar truss:
            :ivar simple:

    
    """

    def __init__(self, simple=None, num_elements=None, supports=None, fem_model=None):

        if fem_model:                
            self.f = fem_model
        else:
            self.f = fem.FrameFEM()
        self.alpha_cr = []
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
        self.joints = {}
        self.r = []
        self.is_generated = False
        self.is_calculated = False
        self.self_weight = False
        self.truss = None
        self.simple = simple

        if simple:
            self.storeys = simple[0]
            self.bays = simple[1]
            self.storey_height = simple[2]
            self.bay_length = simple[3]
            self.generate_members()

        if supports:
            self.generate_supports(supports)

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
                    
            if this.nodal_coordinates not in self.nodal_coordinates:
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

        # TRUSS
        elif isinstance(this, Frame2D):
            self.truss = this
            
            # Bottom chord coordinates
            for bchord in this.bottom_chords:
                c0, c1 = bchord.coordinates
                for member in self.members.values():
                    coord = member.line_intersection(bchord.coordinates)
                    if isinstance(coord, list):
                        member.add_node_coord(coord)
                        # Create eccentricity to connection
                        if member.mtype == "column":
                            joint0 = [j for j in bchord.joints if j.loc <= 0.02]                          
                            joint1 = [j for j in bchord.joints if j.loc >= 0.98]
                            if coord == c0 and joint0:
                                joint0 = joint0[0]
                                web = list(joint0.webs.values())[0]
                                theta = abs(joint0.chord.angle - web.angle)
                                ecc_x = member.h / 2 + abs(joint0.g1*math.cos(joint0.chord.angle) +
                                        web.h/2 * math.sin(theta))
                                # m to mm
                                ecc_x = ecc_x/1000
                                x0, y0 = joint0.coordinate
                                # Change joint's coordinate
                                joint0.coordinate = [x0 + ecc_x, y0]
                            if coord == c0:
                                # Calculate chord's new start coordinate
                                c0 = [c0[0] + member.h/2000, c0[1]]
                                # Add eccentricity element's coordinates to list
                                member.ecc_coordinates.append([coord, c0])  
                                bchord.columns[c0[0]] = member
                            if coord == c1 and joint1:
                                joint1 = joint1[0]
                                web = list(joint1.webs.values())[0]
                                theta = abs(joint1.chord.angle - web.angle)
                                ecc_x = member.h / 2 + abs(joint1.g1*math.cos(joint1.chord.angle) +
                                        web.h/2 * math.sin(theta))
                                # m to mm
                                ecc_x = ecc_x/1000
                                x1, y1 = joint1.coordinate
                                joint1.coordinate = [x1 - ecc_x, y1]
                            if coord == c1:
                                c1 = [c1[0] - member.h/2000, c1[1]]
                                member.ecc_coordinates.append([coord, c1])
                                bchord.columns[c1[0]] = member
                # Change chord's end coordinates
                
                bchord.coordinates = [c0, c1]
                bchord.calc_nodal_coordinates()
                      
            # Top chord coordinates
            for tchord in this.top_chords:               
                c0, c1 = tchord.coordinates
                for member in self.members.values():
                    coord = member.line_intersection(tchord.coordinates)
                    if isinstance(coord, list):
                        member.add_node_coord(coord)
                        # Create eccentricity to connection
                        if member.mtype == "column":
                            
                            joint0 = [j for j in tchord.joints if j.loc <= 0.02]                      
                            joint1 = [j for j in tchord.joints if j.loc >= 0.98]

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
                                #ecc_x = member.h / 2 + abs(joint0.g1*math.cos(joint0.chord.angle) +
                                #        web.h/2 * math.sin(theta))
                                ecc_x = member.h / 2000 + joint0.g1
   
                                coord0 = np.asarray(joint0.coordinate)
                                # Change joint's coordinate
                                joint0.coordinate = list(c0 + ecc_x * v)
                            if coord == c0:
                                # Calculate chord's new start coordinate
                                #c0 = [c0[0] + tchord.h / 2000 * u[0], c0[1] + tchord.h / 2000 * u[1]]
                                ecc_y = tchord.h / 2000 * u[1]
                                #c0 = [c0[0] + member.h/2000, c0[1]]
                                # Add eccentricity elements coordinates to list
                                coord[1] -= ecc_y
                                member.ecc_coordinates.append([coord, c0])
                                tchord.columns[c0[0]] = member
                                temp_coords = sorted(member.coordinates, key=lambda x: x[1])
                                temp_coords[1] = coord
                                member.coordinates = temp_coords
                                member.calc_nodal_coordinates()
                           
                            # End coordinate
                            if coord == c1 and joint1:
                                joint1 = joint1[0]
                                web = list(joint1.webs.values())[0]
                                theta = abs(joint1.chord.angle - web.angle)
                                #ecc_x = member.h / 2 + abs(joint1.g1*1000*math.cos(joint1.chord.angle) +
                                #        web.h/2 * math.sin(theta))  
                                ecc_x = member.h / 2000 + joint1.g1
                                coord1 = np.asarray(joint1.coordinate)
                                joint1.coordinate = list(c1 - ecc_x*v)
                            if coord == c1:
                                # Calculate chord's new end point
                                #c1 = [c1[0] + tchord.h/2000 * u[0], c1[1] + tchord.h/2000 * u[1]]
                                #c1 = [c1[0] - member.h/2000, c1[1]]
                                #member.ecc_coordinates.append([coord, c1])
                                #tchord.columns[c1[0]] = member
                                ecc_y = tchord.h / 2000 * u[1]
                                #c0 = [c0[0] + member.h/2000, c0[1]]
                                # Add eccentricity elements coordinates to list
                                coord[1] -= ecc_y
                                member.ecc_coordinates.append([coord, c1])
                                tchord.columns[c1[0]] = member
                                temp_coords = sorted(member.coordinates, key=lambda x: x[1])
                                temp_coords[1] = coord
                                member.coordinates = temp_coords
                                member.calc_nodal_coordinates()
                                
                
                # Change chord's ends' coordinates
                tchord.coordinates = [c0, c1]
                tchord.calc_nodal_coordinates()

            # Add truss's members to frame's members dict
            for key in this.members.keys():
                this.members[key].mem_id = len(self.members)
                self.members[key] = this.members[key]
                
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
            :param load_id: Id of the loads to be added in the calculation model
            
            :type load_id: int
        """
        if self.is_calculated == False:
            self.is_calculated = True
            self.f.add_loadcase(supp_id=1, load_id=load_id)
            self.f.nodal_dofs()

        self.f.linear_statics()
        self.calc_nodal_forces()
        self.calc_nodal_displacements()
        self.check_members_strength()
        #self.alpha_cr, _ = self.f.linear_buckling(k=4)

    def design_members(self, prof_type="IPE"):
        """ Desgins frame's members
        """
        self.is_designed = True
        for member in self.members.values():
            member.design_member(prof_type)

        self.check_members_strength()

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

    def generate(self):
        """ Generates the frame and truss
        """
        # If the frame hasn't been created, create members and calculate nodal
        # coordinates, otherwise initialize frame.
        if not self.is_generated:
            self.is_generated = True
        for member in self.members.values():
            member.round_coordinates()
            member.generate(self.f)
            for coord in member.nodal_coordinates:
                if coord not in self.nodal_coordinates:
                    self.nodal_coordinates.append(coord)
            self.nodes.extend(list(member.nodes.values()))
        # Generate TrussJoints
        if self.truss:
            for joint in self.truss.joints.values():
                joint.generate(self.f)

        # Generate eccentricity elements
        for member in self.members.values():
            member.generate_eccentricity_elements(self.f)
                   
        # Remove duplicate nodes
        self.nodes = list(set(self.nodes))
        for support in self.supports.values():
            support.add_support(self.f)

        for pointLoad in self.point_loads.values():
            pointLoad.add_load(self.f)

        for lineLoad in self.line_loads.values():
            member = lineLoad.member
            lineLoad.element_ids = member.lineload_elements(lineLoad.coordinates)
            lineLoad.add_load(self.f)
            
    
    
    def to_robot(self, filename, num_frames=1, s=1, brace_profile="SHS 50x50x2"):
        """  Creates an Aurodesk Robot Structural Analysis .str -file
        
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
                    n1 = [n1[0], i*s, n1[1]]
                    n2 = [n2[0], i*s, n2[1]]
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
                    if nodes.index(n1) +1 not in vert_braces:
                        vert_braces.append(nodes.index(n1) +1)
                    
                    if nodes.index(n2) +1 not in vert_braces:
                        vert_braces.append(nodes.index(n2) +1)
                
                if profile_type == 'RHS':
                    profile = 'RRHS ' + splitted_val[1]
                elif profile_type == 'HE':
                    profile = splitted_val[0] + splitted_val[2] + splitted_val[1]
                elif profile_type == "SHS":
                    dims = splitted_val[1].split("X")
                    profile = 'RRHS ' + str(dims[0]) + 'X' + str(dims[0]) + 'X' + str(dims[1])
                else:
                    profile = member.profile
                # Webs
                if member.mtype == "web" or member.mtype == "beam" or\
                    member.mtype=="top_chord" or member.mtype == "bottom_chord":
                    Sj1 = min(1e10, member.Sj1)
                    Sj2 = min(1e10, member.Sj2)
                    releases[elements.index(elem) +1] = f'ORIgin RY Hy={Sj1}  END RY Hy={Sj2}'
    
                profiles.append(profile)
                material.append(member.material)
                
                # Eccentricity elements
                if len(member.ecc_coordinates):
                    for coords in member.ecc_coordinates:
                        n1, n2 = coords
                        if num_frames != 1:
                            n1 = [n1[0], i*s, n1[1]]
                            n2 = [n2[0], i*s, n2[1]]
                        if n1 not in nodes:
                            nodes.append(n1)
                        if n2 not in nodes:
                            nodes.append(n2)
                        elem = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                        elements.append(elem)
                        profiles.append("HEA 1000")
                        material.append("S355")
                        #if n1 < n2:
                        #    releases[elements.index(elem) +1] = 'END RY'
                        #else:
                        #    releases[elements.index(elem) +1] = 'ORIgin RY'
                        
                        
    
            try:
                if self.truss:
                    truss = self.truss
                else:
                    truss = self
                for joint in truss.joints.values():
                    # Y joint
                    if len(joint.nodal_coordinates) == 2:
                        n1, n2 = sorted(joint.nodal_coordinates)
                        if num_frames != 1:
                            n1 = [n1[0], i*s, n1[1]]
                            n2 = [n2[0], i*s, n2[1]]
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
                    elif len(joint.nodal_coordinates) == 5:
                        coords = sorted(joint.nodal_coordinates)
                        n1, n2, n3, n4, n5 = coords
                        if num_frames != 1:
                            n1 = [n1[0], i*s, n1[1]]
                            n2 = [n2[0], i*s, n2[1]]
                            n3 = [n3[0], i*s, n3[1]]
                            n4 = [n4[0], i*s, n4[1]]
                            n5 = [n5[0], i*s, n5[1]]
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
                        #releases[elements.index(elem1)+1] = "END RY"
                        #releases[elements.index(elem2)+1] = "END RY"
                    else:
                        pass
                    
            except:
                pass
                    
            for pointload in self.point_loads.values():
                try:
                    n1 = pointload.coordinate
                    if num_frames != 1:
                        n1 = [n1[0], i*s, n1[1]]
                    idx = nodes.index(n1) + 1
                except ValueError:
                    n1 = pointload.coordinate
                    if num_frames != 1:
                        n1 = [n1[0], i*s, n1[1]]
                    nodes.append(n1)
                    idx = nodes.index(n1) +1
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
        for i in range(len(vert_braces) -3):
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
                profile = 'RRHS ' + str(dims[0]) + 'X' + str(dims[0]) + 'X' + str(dims[2])
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
            f.write("LENgth=m	Force=kN \n")
            f.write("NODes \n")
            for i, node in enumerate(nodes):
                if num_frames != 1:
                    f.write(f'{i+1}   {node[0]}    {node[1]}    {node[2]} \n')
                else:
                    f.write(f'{i+1}   {node[0]}    {node[1]} \n')

            f.write('ELEments \n')
            for i, element in enumerate(elements):
                f.write(f' {i+1}  {element[0]}    {element[1]} \n')

            f.write('PROperties \n')
            for i in range(len(elements)):
                f.write(f' "{material[i]}" \n')
                f.write(f' {i+1}  ')
                f.write(f' {profiles[i]}    \n')

            f.write("SUPports \n")
            if self.is_calculated:
                for sup in self.supports.values():
                    for i in range(num_frames):
                        n1 = sup.coordinate
                        if num_frames != 1:
                            n1 = [n1[0], i*s, n1[1]]
                        idx = nodes.index(n1) + 1
                        if sup.dofs[0] == [-1]:
                            f.write(f'{idx}  UX  \n')
                        elif sup.dofs[0:1] == [-1, -1]:
                            f.write(f'{idx}  UX UZ \n')
                        elif sup.dofs == [-1, -1, -1]:
                            f.write(f'{idx}  \n')
                        elif sup.dofs[1] == [-1]:
                            f.write(f'{idx}  UZ  \n')
                        elif sup.dofs[1:2] == [-1, -1]:
                            f.write(f'{idx}  UZ RY  \n')
                        else:
                            f.write(f'{idx}  RY \n')
            else:
                for sup in self.supports.values():
                    for i in range(num_frames):
                        n1 = sup.coordinate
                        if num_frames != 1:
                            n1 = [n1[0], i*s, n1[1]]
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
                            n1 = [n1[0], i*s, n1[1]]
                            n2 = [n2[0], i*s, n2[1]]
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
             loads=True, color=False):
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
        if self.is_calculated and color:
            color = True
              
        # Plot members
        for member in self.members.values():
            member.plot(print_text, color)
        
        # Plot joints
        for joint in self.joints.values():
            joint.plot(color=color)

        # Plot supports
        for support in self.supports.values():
            node_coord = support.coordinate
            if support.dofs == [-1, -1, -1] or support.dofs == [1,1,1]:
                marker = 's'
            elif support.dofs == [1]:
                if node_coord[0] == 0:
                    marker = '3'
                else:
                    marker = '4'
            else:
                marker = '2'
            plt.scatter(node_coord[0], node_coord[1], s=50, c='k', marker=marker)
       

        #if self.truss:
        #    self.truss.plot(show=False, print_text=print_text, color=color)
        if loads:
            self.plot_loads()
        plt.axis('equal')
        if show:
            plt.axis('equal')
            plt.show()

    def plot_loads(self):
        """ Plots loads
        """

        for load in self.point_loads.values():
            x, y = load.coordinate
            plt.scatter(x, y, c='b', marker='*')

        for lineload in self.line_loads.values():
            c0, c1 = lineload.coordinates
            plt.plot([c0[0], c1[0]], [c0[1], c1[1]], c='b')
            
            #x0, y0 = self.f.elements[lineload.element_ids[0]].nodes[0].x
            #x1, y1 = self.f.elements[lineload.element_ids[-1]].nodes[1].x
            #plt.plot([x0, x1], [y0, y1], c='b')

    def plot_deflection(self, scale=1, prec=4, show=True):
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
        #if self.truss:
        #    self.truss.plot_deflection(scale, show=False)
        
        self.plot(print_text=False, show=False)
        self.calc_nodal_displacements()
        for member in self.members.values():
            X = []
            Y = []

            member.calc_nodal_displacements(self.f)
            max_x = 0
            max_y = 0
            for i, node in enumerate(member.nodes.values()):
                x0, y0 = node.x
                x1 = node.u[0]
                y1 = node.u[1]
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

            elif member.mtype == "column":
                plt.plot(loc_max_x, loc_max_y, 'ro')
                plt.text(loc_max_x, loc_max_y,
                         (str(max_x)[0:prec + 1] + " mm"))

        if show:
            plt.show()
        
        
    def plot_buckling(self, scale=1, k=4, show=True):
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
        w, v = self.f.linear_buckling(k=k)
        self.alpha_cr = w
        for j in range(v.shape[1]):
            if show:
                self.plot(print_text=False, show=False, loads=False, color=False)
            for member in self.members.values():
                X = []
                Y = []
                for i, node in enumerate(member.nodes.values()):
                    x0, y0 = node.x
                    if not isinstance(node.v[0], int):
                        x1 = node.v[0][j]
                        y1 = node.v[1][j]
                    else:
                        x1 = node.v[0]
                        y1 = node.v[1]
                    x = x0 + x1 * (scale)
                    y = y0 + y1 * (scale)
                    X.append(x)
                    Y.append(y)
                    
                plt.plot(X, Y, color='m')
            plt.title(f'Buckling shape {j+1},  ' r"$\alpha_{cr}$" f' = {w[j]*1000:.2f}')
            #if self.truss:
            #    self.truss.plot_buckling(show=False)
            
            plt.show()

    def bmd(self, scale=1):
        """ Draws bending moment diagram
            
            Parameters
            ----------
            :param scale: Scaling factor
            :type scale : int
        """
        self.plot(print_text=False, show=False, color=False)
        for member in self.members.values():
            member.bmd_test(scale)
        if self.truss:
            self.truss.bmd(scale)
        plt.show()
        
    def plot_normal_force(self, show=True):
        """ Plots normal force and utilization ratio
            
            Parameters
            ------------
            :param show: Shows the plotted diagram, allows to plot multiple
                        diagrams to be plotted on same picture
            :type show: bool
        
        
        """
        for member in self.members.values():
            member.plot_normal_force()
        #if self.truss:
        #    for member in self.truss.members.values():
        #        member.plot_normal_force()
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
        MIN_VAL = 0
        for member in self.members.values():
            if member.mtype == "beam":
                member.alpha1 = MIN_VAL
                member.alpha2 = MIN_VAL
                member.Sj1 = MIN_VAL
                member.Sj2 = MIN_VAL

    
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
        SemiRigidFrame
        FrameFEM
        I_sections
        CrossSection
        SteelMember
            
    """

    def __init__(self, coordinates, mem_id="", profile="IPE 100",
                 material="S355",num_elements=5,
                 Sj1=np.inf, Sj2=np.inf, mtype=""):

        self.element_ids = []
        self.elements = {}
        self.ecc_elements = {}
        self.ecc_coordinates = []
        self.ecc_element_nodes = {}
        self.steel_member = None
        # strat node, FEMNode object
        self.n1 = None
        # end node, FEMNode object
        self.n2 = None
        self.nodes = {}
        self.added_coordinates = []
        self.nodal_coordinates = []
        self.loc = []
        self.__coordinates = None
        self.coordinates = sorted(coordinates, key=lambda x: x[0])
        self.material = material
        self.cross_section = None
        self.__profile = profile
        self.profile = profile
        #self.profile_idx = PROFILES.index(self.profile)
        #self.length = self.calc_length()
        self.mem_id = mem_id
        self.mtype = mtype
        self.nodal_forces = {}
        self.nodal_displacements = {}
        self.num_elements = num_elements
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
        self.loads = []
        self.is_designed = False
        self.is_strong_enough = False
        self.r = np.ones(7)
        self.self_weight = False
        self.is_generated = False
        # Lineload values
        self.q0 = 0
        self.q1 = 0
        # Create SteelMember object
        self.steel_member = SteelMember(self.cross_section, self.length,
                                   Lcr=[1.0, 1.0], mtype=self.mtype)
    @property
    def angle(self):
        """
        Returns member's angle in radians
        """
        start_node, end_node = self.coordinates
        x0, y0 = start_node
        x1, y1 = end_node
        if (x1-x0) == 0:
            angle = math.radians(90)
        else:
            angle = math.atan((y1 - y0) / (x1 - x0))
        return angle


    @property
    def coordinates(self):
        if self.n1 and self.n2:
            if self.__coordinates != [list(self.n1.x), list(self.n2.x)]:
                self.__coordinates = [list(self.n1.x), list(self.n2.x)]
                #self.calc_nodal_coordinates(self.num_elements)
        return self.__coordinates
    
    @coordinates.setter
    def coordinates(self, val):
        if not self.n1:
            self.__coordinates = val
        else:
            self.__coordinates = val
            x1, y1 = val[0]
            x2, y2 = val[1]
            self.n1.x = np.array([x1, y1])
            self.n2.x = np.array([x2, y2])
            for mem in self.n1.parents:
                mem.calc_nodal_coordinates()
            for mem in self.n2.parents:
                mem.calc_nodal_coordinates()
            
    @property
    def MRd(self):
        """ Returns bending resistance
            Units: Nmm
        """
        return self.cross_section.MRd
    
    @property
    def NRd(self):
        """ Returns axial force resistance
            Units: N
        """
        return self.cross_section.NRd
    
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
        
        weight = self.cross_section.A * self.length * self.rho * 1000
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
        Sets cross-section's height to givel value.
        Changes profile and sets new cross-sectional properties to member's elements
        """
        
        splitted = self.profile.split()
        if len(splitted) == 2:
            profile_type = self.profile.split()[0].upper()
            if profile_type == "IPE":
                self.profile = profile_type + " " + str(val)
            elif profile_type == "SHS":
                val = max(val, 25)
                self.profile = profile_type + " " + str(val) + 'X'+ str(val) + 'X' + str(max(round(val / 20), 2))
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
        self.fy = MATERIALS[self.material]['f_y']
        self.E = MATERIALS[self.material]['E']
        self.fu = MATERIALS[self.material]['f_u']
        self.nu = MATERIALS[self.material]['v']
        self.rho = MATERIALS[self.material]['rho']

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
                raise ValueError('{} is not valid profile type!'.format(profile_type))
        # Change member's elements' properties
        if len(self.elements) > 0:
            for element in self.elements.values():
                element.section.A = self.cross_section.A * 1e-6
                element.section.Iy = self.cross_section.I[0] * 1e-12
        # Change steel_member objects properties
        if self.steel_member:
            self.steel_member.profile = self.cross_section

                      
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

    @property
    def length(self):
        start_node, end_node = self.coordinates
        x0, y0 = start_node
        x1, y1 = end_node
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
            Calculates node locations along member
            Coordinates used are global coordinates
            Adds locations to a list where they can be accessed later
            :param num_elements: int, number of elements to be created
        """
        if num_elements != 0:
            self.num_elements = num_elements
        self.nodal_coordinates = []
        self.nodal_coordinates.extend(self.added_coordinates)
        start_node, end_node = self.coordinates
        x0, y0 = start_node
        Lx = end_node[0] - start_node[0]
        Ly = end_node[1] - start_node[1]
        dx = Lx / self.num_elements
        dy = Ly / self.num_elements
        for i in range(self.num_elements):
            x = i * dx + x0
            y = i * dy + y0
            loc = i / self.num_elements
            node = [x, y]
            if node not in self.nodal_coordinates:
                self.nodal_coordinates.append(node)
            if loc not in self.loc:
                self.loc.append(loc)
        if [end_node[0], end_node[1]] not in self.nodal_coordinates:
            self.nodal_coordinates.append([end_node[0], end_node[1]])
        if 1 not in self.loc:
            self.loc.append(1)

        if self.is_generated:
            for j, node in enumerate(self.nodes.values()):
                #print(self.nodal_coordinates)
                #print(j, len(self.nodal_coordinates))
                node.x = np.array(self.nodal_coordinates[j])
                
        self.nodal_coordinates.sort()
    def add_node_coord(self, coord):
        """
        Adds coordinate to memeber's coordinates and
        calculates new coordinates local location
        If coordinate is not between member's coordinates, reject it
        :param coord: array of two float values, node's coordinates
        """      
        if coord not in self.added_coordinates and self.point_intersection(coord):                  
            self.added_coordinates.append(coord)
            self.nodal_coordinates.append(coord)
            self.nodal_coordinates.sort()
            start_node, end_node = self.coordinates
            x, y = coord
            x1, y1 = start_node
            dx = x - x1
            dy = y - y1
            dz = math.sqrt((dx ** 2 + dy ** 2))
            z = dz / self.length
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
                n2 = node_ids[i+1]
            
                self.elements[index] = \
                    EBBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                           fem_model.sections[self.mem_id], \
                           fem_model.materials[self.mem_id])

                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
                index += 1

        if self.mtype == "beam":
            if len(self.nodes) == 2:
                n1 = node_ids[0]
                n2 = node_ids[1]
                
                self.elements[index] = \
                        EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                                        fem_model.sections[self.mem_id], \
                                        fem_model.materials[self.mem_id], \
                                        rot_stiff=[self.Sj1, self.Sj2])
                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
            else:
                for i in range(len(self.nodes) - 1):
                    n1 = node_ids[i]
                    n2 = node_ids[i+1]
                 
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
        sect = fem.BeamSection(38880*1e-6, 6299657109*1e-12)
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
                                     #rot_stiff=[np.inf, 0])
            fem_model.add_element(self.ecc_elements[idx])


    def calc_nodal_forces(self):
        """ Gets calculated nodal forces on member's elements and
            adds them to dict where node is the key and froces are value
            self.nodal_forces[node_id] = [Fx, Fy, Mz]
        """
        self.med = 0
        self.ved = 0
        self.ned = 0

        i = 0
        node_ids = list(self.nodes.keys())
        for element in self.elements.values():
            node_id = node_ids[i]
            axial_force = element.axial_force[0]
            if abs(axial_force) > abs(self.ned):
                self.ned = axial_force

            shear_force = element.shear_force[0]
            if abs(shear_force) > abs(self.ved):
                self.ved = shear_force

            bending_moment = element.bending_moment[0]
            if abs(bending_moment) > abs(self.med):
                self.med = bending_moment

            self.nodal_forces[node_id] = [axial_force,
                                                shear_force,
                                                bending_moment]
            i += 1
            
            if i == len(self.elements):
                
                node_id = node_ids[i]
                axial_force = element.axial_force[1]
                if abs(axial_force) > abs(self.ned):
                    self.ned = axial_force

                shear_force = element.shear_force[1]
                if abs(shear_force) > abs(self.ved):
                    self.ved = shear_force

                bending_moment = element.bending_moment[1]
                if abs(bending_moment) > abs(self.med):
                    self.med = bending_moment

                self.nodal_forces[node_id] = [axial_force,
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
        
        if not self.steel_member.nsect():
            for i, node in enumerate(self.nodal_forces.keys()):
                forces = self.nodal_forces[node]
                ned = forces[0] * 1e3  # kN to N
                ved = forces[1] * 1e3  # kN to N
                med = forces[2] * 1e6  # kNm to Nmm
                loc = self.loc[i]
                self.steel_member.add_section(ned=ned, vzed=ved,
                                              myed=med, loc=loc)
        else:
            for i, node in enumerate(self.nodal_forces.keys()):
                forces = self.nodal_forces[node]
                self.steel_member.ned[i] = forces[0] * 1e3  # kN to N
                self.steel_member.vzed[i] = forces[1] * 1e3  # kN to N
                self.steel_member.myed[i] = forces[2] * 1e6  # kNm to Nmm
        
        # Cross-sectional stress ratios in member's nodes' locations
        self.r[:4] = self.steel_member.check_sections()
        # Buckling about y - and z - axis 
        buckling_r = self.steel_member.check_buckling()
        self.r[4] = buckling_r[0]
        self.r[5] = buckling_r[1]
        # Lateral-Torsional buckling
        self.r[6] = self.steel_member.check_LT_buckling()
        # Set cross-section's maximum forces for getting results later       
        self.cross_section.Ned = self.ned * 1e3  # kN to N
        self.cross_section.Ved = self.ved * 1e3  # kN to N
        self.cross_section.Med = self.med * 1e6  # kNm to Nmm

        if max(self.r) > 1.0:
            self.is_strong_enough = False
        else:
            self.is_strong_enough = True

    def design_member(self, prof_type):
        """ Goes through all profiles in list and stops iterating when 
            cross-section can bear it's loads.
            If member is previously designed, iterating starts from 3 profiles
            before currrent profile.
        """
        prof_type = prof_type.upper()
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

        for profile in profiles:
            self.profile = profile
            self.check_cross_section()
            if max(self.r) <= 1.0:
                break



    def line_intersection(self, coordinates):
        """
        Calculates coordinate where two members intersect
        :param coordinates: array of two arrays of two float values, [[x1,y1],[x2,y2]]
        :return [px, py]: array of two float values, coordinates for intersection
        source: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
        """
        start_node, end_node = self.coordinates
        x1, y1 = start_node
        x2, y2 = end_node
        x3, y3 = coordinates[0]
        x4, y4 = coordinates[1]
        
        
        if ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)) == 0:
            return None
        else:
            px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / (
                (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
            py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / (
                (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))

            return [round(px,5), round(py,5)]


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
                return math.isclose(y, k*x + b, rel_tol=1e-3) #y == k * x + b
        return False
        

    def plot(self, print_text=True, c='k'):

        X = self.coordinates
        if c:
            if self.is_strong_enough:
                color = 'green'
            else:
                color = 'red'
        else:
            color = 'k'
        # Plot members
        plt.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color)
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
        if self.mtype == 'beam':
            x = (X[0][0] + X[1][0]) / 2
            y = (X[1][1] + X[0][1]) / 2
            horzalign = 'center'
            vertalign = 'bottom'
        elif self.mtype == 'column':
            x = (X[0][0] + X[1][0]) / 2
            y = (X[1][1] + X[0][1]) / 2
            horzalign = 'right'
            vertalign = 'center'
            
        else:
            x = (X[0][0] + X[1][0]) / 2
            y = (X[1][1] + X[0][1]) / 2
            horzalign = 'right'
            vertalign = 'center'
            
        if print_text:
            plt.text(x, y, str(self.mem_id) + ": " + self.profile,
                     rotation=rot, horizontalalignment=horzalign,
                     verticalalignment=vertalign)
        else:
            plt.text(x, y, str(self.mem_id),
                     rotation=rot, horizontalalignment=horzalign,
                     verticalalignment=vertalign)


    def plot_results(self):
        
        UN,UV,UM,UMN = self.cross_section.section_resistance(return_list=True)
        B = max(self.steel_member.check_buckling())
        LT = self.steel_member.check_LT_buckling()
        
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
        plt.xticks(x,X)
        plt.title(f'Member {self.mem_id} stress ratios')
        line_x = np.arange(-0.5,len(Y))
        line_y = [1]*len(line_x)
        plt.plot(line_x, line_y,'--', c='k')
        plt.show()
        

    def bmd(self, scale):
        """
        Plots bending moment diagram
        :param scale: float, scaling factor
        """
        X = []
        Y = []
        self.calc_nodal_forces()

        moment_values = [x[2] for x in self.nodal_forces.values()]
        node_ids = list(self.nodes.keys())
        x = self.nodal_coordinates[0][0]
        y = self.nodal_coordinates[0][1]
        X.append(x)
        Y.append(y)
        for i in range(len(self.nodes)):
            node = node_ids[i]
            if self.mtype == "beam":
                x0, y0 = self.nodal_coordinates[i+1].x
                bending_moment = self.nodal_forces[node][2]
                y1 = bending_moment / (1000 / scale)
                x = x0
                y = y0 - y1
                X.append(x)
                Y.append(y)
                if i == 0 or i == len(self.nodes)-1 or bending_moment == max(moment_values) or\
                                bending_moment == min(moment_values):
                    if bending_moment > 0:
                        vert = 'top'
                    else:
                        vert = 'bottom'
                    plt.text(x, y, f'{bending_moment:.2f} kNm', verticalalignment=vert)

            elif self.mtype == "column":
                x0 = self.nodal_coordinates[i+1][0]
                y0 = self.nodal_coordinates[i+1][1]
                bending_moment = self.nodal_forces[node][2]
                x1 = bending_moment / (1000 / scale)
                x = x0 + x1
                y = y0
                X.append(x)
                Y.append(y)
                if i == 0 or i == len(self.nodes)-1 or bending_moment == max(moment_values) or\
                                bending_moment == min(moment_values):
                    if bending_moment > 0:
                        horz = 'left'
                    else:
                        horz = 'right'
                    plt.text(x, y, f'{bending_moment:.2f} kNm', horizontalalignment=horz)
        X.append(self.nodal_coordinates[i+1][0])
        Y.append(self.nodal_coordinates[i+1][1])
        plt.plot(X, Y, color='gray')

    def bmd_test(self, scale=1):

        # Member's position vector
        v = np.array([math.cos(self.angle), math.sin(self.angle)])
        # Vector perpendicular to v
        u = np.array([-math.sin(self.angle), math.cos(self.angle)])
        X = []
        Y = []
        start_coord, end_coord = self.coordinates
        x01, y01 = start_coord
        x02, y02 = end_coord
        X.append(x01)
        Y.append(y01)
        self.calc_nodal_forces()

        moment_values = [x[2] for x in self.nodal_forces.values()]
        
        for elem in self.elements.values():
            x0, y0 = elem.nodes[0].x
            bending_moment = elem.bending_moment[0]
            val = bending_moment / (1000 / scale)
            x, y = np.array([x0, y0]) - val * u
            X.append(x)
            Y.append(y)
            horzalign = 'center'
            vertalign = 'center'
            if bending_moment == max(moment_values) or bending_moment == min(moment_values):
                plt.text(x, y, f'{bending_moment:.2f} kNm', horizontalalignment=horzalign)
        x0, y0 = elem.nodes[1].x
        bending_moment = elem.bending_moment[1]
        val = bending_moment / (1000 / scale)
        x, y = np.array([x0, y0]) - val * u
        X.append(x)
        Y.append(y)
        X.append(x02)
        Y.append(y02)
        plt.text(x, y, f'{bending_moment:.2f} kNm', horizontalalignment=horzalign)
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
        X = self.coordinates
        if self.ned < 0:          
            r = min(1, max(self.r[4:5]))
            alpha = max(0.1, r)
            color = (1,0.5-r/2,0, alpha)
        else:
            r = min(1, self.r[0])
            alpha = max(0.1, r)
            color = (0,0.5-r/2,1, alpha)
        # Plot members
        plt.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color=color, linewidth=2)
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
        plt.text(x, y, f'{self.ned:.2f}  kN,\nr: {r*100:.2f} %',
                 rotation=rot, horizontalalignment=horzalign,
                 verticalalignment=vertalign)

    def lineload_elements(self, coords):
    
    #Returns the elements that are under line load
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
        
        for i, coord in enumerate(self.nodal_coordinates):
            self.nodal_coordinates[i] = [round(c, prec) for c in coord]

class SteelBeam(FrameMember):
    def __init__(self, coordinates, alpha1=25, alpha2=25, mem_id="",
                 profile="IPE 100", material="S355",num_elements=4):
        super().__init__(coordinates, mem_id, profile, material, num_elements)

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


class SteelColumn(FrameMember):
    def __init__(self, coordinates, mem_id="", profile="IPE 100",
                 material="S355", num_elements=4):
        super().__init__(coordinates, mem_id, profile, material,num_elements)

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
    def __init__(self, coordinates, values, direction, load_id=2, f=1.0,
                 ltype='live', name='LineLoad'):
        super().__init__(load_id, ltype, name, f)
        if isinstance(coordinates, FrameMember):
            self.member = coordinates
            self.member.has_load = True
            self.member.q0 = values[0]
            self.member.q1 = values[1]
            coordinates = self.member.coordinates
            self.member.loads.append(self)
        else:
            self.mem_id = None
            self.member = None
        self.coordinates = coordinates
        self.values = values
        self.direction = direction
        self.f = f
        self.element_ids = None

    def calc_k(self):
        v0, v1 = self.values
        k = (v1 - v0) / len(self.member.elements)
        return k
    
    def add_load(self, fem_model):
        k = self.calc_k()
        v0, v1 = self.values
        for i, elem_id in enumerate(self.member.elements.keys()):
            v1 = (i*k) + v0
            load = fem.LineLoad(self.load_id,
                                fem_model.elements[elem_id],
                                [0.0, 1.0],
                                [v0, v1],
                                self.direction)
            fem_model.add_load(load)
            v0 = v1

# --------------------------- SUPPORT CLASSES ----------------------------
class Support:
    def __init__(self, coordinate, dofs, supp_id=1):
        
        self.node = None
        self.node_id = None
        self.coordinate = coordinate
        self.supp_id = supp_id
        self.dofs = dofs
        
    
    @property
    def coordinate(self):
        if self.node and list(self.node.x) != self.__coordinate:
            self.__coordinate = list(self.node.x)          
        return self.__coordinate
    
    @coordinate.setter
    def coordinate(self, val):
        self.__coordinate = val
        if self.node:
            self.node.x = np.array(val)
 
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
        super().__init__(coordinate, [1, 1, 1], supp_id)


class XHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [1, 0, 0], supp_id)


class YHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0, 1, 0], supp_id)


class XYHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [1, 1, 0], supp_id)


class Hinge(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0, 0, 0], supp_id)

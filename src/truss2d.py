# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:11:35 2018

@author: huuskoj
"""
from frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, Support
from fem.frame.elements.ebbeam import EBBeam
from fem.frame.elements.eb_semi_rigid_beam import EBSemiRigidBeam
from fem.frame.frame_fem import FrameFEM, Material, BeamSection

import math
import numpy as np
import matplotlib.pyplot as plt


class Truss2D(Frame2D):
    def __init__(self, simple=None, num_elements=2, fem=FrameFEM()):
        super().__init__(num_elements=num_elements, fem=fem)
        self.top_chord = None
        self.bottom_chord = None
        self.webs = {}
        self.joints = {}

        if simple:
            coords1, coords2, num = simple
            self.top_chord = TopChord(coords1)
            self.bottom_chord = BottomChord(coords2)


    def add(self, this):
        # Top chord
        if isinstance(this, TopChord):
            self.top_chord = this
            this.mem_id = len(self.members)
            self.members["T" + str(this.mem_id)] = this
            this.calc_nodal_coordinates(self.num_elements)
        # Bottom chord    
        elif isinstance(this, BottomChord):
            self.bottom_chord = this
            this.mem_id = len(self.members)
            self.members["B" + str(this.mem_id)] = this
            this.calc_nodal_coordinates(self.num_elements)
        # Web and truss joint
        elif isinstance(this, TrussWeb):
            # Calculate web's coordinates
            c01 = self.bottom_chord.local(this.bot_loc)
            c02 = self.top_chord.local(this.top_loc)
            # Change web's coordinates to global coordinates
            this.bot_coord = c01
            this.top_coord =  c02
            # Check if there's already a joint at given location
            j1, j2 = None, None
            for joint in self.joints.values():
                if joint.coordinate == c01:
                    j1 = joint
                elif joint.coordinate == c02:
                    j2 = joint    
            # Set joints to chord
            if not j1:
                j1 = TrussJoint(self.bottom_chord, this.bot_loc)
            if not j2:
                j2 = TrussJoint(self.top_chord, this.top_loc)
            # Assign member id
            this.mem_id = len(self.webs)
            # Assign joint id's   
            if j1.jid == None:
                j1.jid = len(self.joints)
                self.joints[j1.jid] = j1
            if j2.jid == None:
                j2.jid = len(self.joints)
                self.joints[j2.jid] = j2

            j1.add_web(this, self.num_elements)
            j2.add_web(this, self.num_elements)
     
            #this.mem_id = len(self.members)
            self.webs[this.mem_id] = this
            self.members["W" + str(this.mem_id)] = this
            
            # Make web's joints hinged
            this.Sj1 = 0
            this.Sj2 = 0

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

    def generate(self):
        """ Generates the truss
        """
        # Generate members' nodes and elements
        for member in self.members.values():
            member.generate(self.f)
            for coord in member.nodal_coordinates:
                if coord not in self.nodal_coordinates:
                    self.nodal_coordinates.append(coord)
            self.nodes.extend(list(member.nodes.values()))
        # Remove duplicate nodes
        self.nodes = list(set(self.nodes))
        # Generate joints' nodes and elements
        for joint in self.joints.values():
            joint.generate(self.f)
        for support in self.supports.values():
            support.add_support(self.f)

        for pointLoad in self.point_loads.values():
            pointLoad.add_load(self.f)

        for lineLoad in self.line_loads.values():
            member = lineLoad.member
            lineLoad.element_ids = member.lineload_elements(lineLoad.coordinates)
            lineLoad.add_load(self.f)


class TrussMember(FrameMember):
    def __init__(self, coordinates, mem_id="",
                 profile="SHS 100x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

    def local(self, value):
        """
        Returns chord's local coordinate in global coordinates
        :param value: local coordinate, 0 <= value <= 1
        """
       
        if 0 > value > 1:
            raise ValueError('Value must be between 0 and 1')
        start_node, end_node = self.coordinates

        x0, y0 = start_node
        Lx = end_node[0] - start_node[0]
        Ly = end_node[1] - start_node[1]
        try:
            k = Ly / Lx
        except ZeroDivisionError:
            k = 0
        return [Lx * value + x0, k * value * Lx + y0]
    
    def shape(self, x):
        start_node, end_node = self.coordinates
        x0, y0 = start_node
        Lx = end_node[0] - start_node[0]
        Ly = end_node[1] - start_node[1]
        try:
            k = Ly / Lx
        except ZeroDivisionError:
            k = 0
        return k*x + y0
        
        
        

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

    def generate_elements(self, fem_model):
        """ Generates elements between nodes
            For beam member's first and last element, creates semi-rigid end
            elements. Elements are added to a list
            :param fem_model: FrameFEM -object
        """
        index = fem_model.nels()
        node_ids = list(self.nodes.keys())
        # EBBeam -elements
        if self.mtype == "top_chord" or self.mtype == 'bottom_chord':
            for i in range(len(self.nodes) - 1):
                n1 = node_ids[i]
                n2 = node_ids[i + 1]

                self.elements[index] = \
                    EBBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                           fem_model.sections[self.mem_id], \
                           fem_model.materials[self.mem_id])

                fem_model.add_element(self.elements[index])
                self.element_ids.append(index)
                index += 1

        if self.mtype == "web":
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
                    n2 = node_ids[i + 1]

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

    def plot(self, print_text, c):

        X = self.coordinates
        if c:
            if self.is_strong_enough:
                color = 'green'
            else:
                color = 'red'
        color = 'k'
        # Plot members
        plt.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color)
        # Plot text
        if print_text:
            # Calculate text location
            delta_y = X[1][1] - X[0][1]
            delta_x = X[1][0] - X[0][0]
            if delta_x == 0:
                rot = 90
            elif delta_y == 0:
                rot = 0
            else:
                rot = math.degrees(math.atan(delta_y / delta_x))

            x = (
                    X[0][0] + X[1][0]) / 2
            y = (X[1][1] + X[0][1]) / 2
            horzalign = 'center'
            vertalign = 'bottom'

            plt.text(x, y, str(self.mem_id) + ": " + self.profile,
                     rotation=rot, horizontalalignment=horzalign,
                     verticalalignment=vertalign)


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
            if self.mtype == "top_chord" or self.mtype == "bottom_chord":
                x0 = self.nodal_coordinates[i][0]
                y0 = self.nodal_coordinates[i][1]
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

            elif self.mtype == "web":
                x0 = self.nodal_coordinates[i][0]
                y0 = self.nodal_coordinates[i][1]
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
        X.append(self.nodal_coordinates[i][0])
        Y.append(self.nodal_coordinates[i][1])
        plt.plot(X, Y, color='gray')


class TopChord(TrussMember):
    def __init__(self, coordinates, mem_id="",
                 profile="SHS 100x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'top_chord'


class BottomChord(TrussMember):
    def __init__(self, coordinates, mem_id="",
                 profile="SHS 100x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)
        self.mtype = 'bottom_chord'


class TrussWeb(TrussMember):
    def __init__(self, bot_loc, top_loc, mem_id="",
                 profile="SHS 50x5", material="S355", Sj1=0, Sj2=0):
        super().__init__([[0,0],[2,2]], profile=profile) # Member's coordinates are changed in add-function
        # Chord's local coordinates
        self.bot_loc = bot_loc
        self.top_loc = top_loc
        # Global coordinates
        self.__top_coord = None
        self.__bot_coord = None
        self.top_coord = self.coordinates[1]
        self.bot_coord = self.coordinates[0]
        
        self.__Sj1 = Sj1
        self.__Sj2 = Sj2
        self.mtype = 'web'
        
     
    @property
    def bot_coord(self):
        return self.__bot_coord
    
    @bot_coord.setter
    def bot_coord(self, val):
        self.__bot_coord = val
        if self.n1:
            self.n1.x = val
        self.coordinates = [val, self.top_coord]
        
        
    @property
    def top_coord(self):
        return self.__top_coord
    
    @top_coord.setter
    def top_coord(self, val):
        self.__top_coord = val
        if self.n2:
            self.n2.x = val
        self.coordinates = [self.bot_coord, val]
        

class TrussJoint():
    def __init__(self, chord, loc, joint_type="N", g1=0.1, g2=0.1):

        self.chord = chord
        self.jid = None
        self.__coordinate = chord.local(loc)
        self.cnode = None # central node
        self.joint_type = joint_type
        self.g1 = g1
        self.g2 = g2
        # nodes in order [center, web, chord, web, chord]
        self.nodal_coordinates = []
        self.nodes = {}
        self.elements = {}
        self.element_ids = []
        self.webs = {}
        self.is_generated = False
        self.num_elements = None
        # Index for eccentricity elements material and section
        self.idx = None

        
    def add_node_coord(self, coord):
        if coord not in self.nodal_coordinates:
            self.nodal_coordinates.append(coord)
    
    def add_web(self, web, num_elements):
        self.webs[web.mem_id] = web
        self.calc_nodal_coordinates(num_elements)

    @property
    def coordinate(self):
        return self.__coordinate

    @coordinate.setter
    def coordinate(self, val):
        self.__coordinate = self.chord.local(val)
        if self.num_elements:
            self.calc_nodal_coordinates(self.num_elements)
            
    
    def calc_nodal_coordinates(self, num_elements):
        
        
        v = np.array([math.cos(self.chord.angle), math.sin(self.chord.angle)])
        u = np.array([-math.sin(self.chord.angle), math.cos(self.chord.angle)])
        
        self.num_elements = num_elements
        # Initialize nodal coordinates as an empty array
        self.nodal_coordinates = []
        central_coord = np.asarray(self.coordinate)
        # Eccentricity along y-axis
        ecc_y = self.chord.h / 2000 # div by 2 and 1000 (mm to m)
        # Add joint's central coordinate x0,y0 to joint and chord
        self.add_node_coord(list(central_coord))
        self.chord.add_node_coord(list(central_coord))
        # If joint type is Y
        if len(self.webs) == 1:
            web = list(self.webs.values())[0]
            if self.chord.mtype == 'top_chord':
                coord = list(central_coord - [0,ecc_y])
                
                self.add_node_coord(coord)
                web.top_coord = coord
            else:
                coord = list(central_coord + [0,ecc_y])
                self.add_node_coord(coord)
                web.bot_coord = coord
            web.calc_nodal_coordinates(num_elements)
            
        # If joint type is K or KT
        elif len(self.webs) >= 2:
            # Iterate trhough every web member in the joint
            for web in self.webs.values():
                # angle between chord and web
                theta = abs(self.chord.angle - web.angle)
                # Eccentricity along x-axis
                ecc_x = self.g1/2 *math.cos(self.chord.angle) + web.h/2000 * math.cos(theta)             
                if len(self.webs) >=3:
                    if 0 < web.angle < math.radians(90):
                        ecc_x = self.g1/2 *math.cos(self.chord.angle) + web.h/2000 * math.cos(theta)
                    elif -math.radians(90) < web.angle < 0:
                        ecc_x = self.g2/2 *math.cos(self.chord.angle) + web.h/2000 * math.cos(theta)
                    else:
                        ecc_x = 0
                # TOP CHORD
                if self.chord.mtype == 'top_chord':                   
                    if 0 < web.angle < math.radians(90):
                        ecc_x = -ecc_x                    
                    # Web's eccentricity node on chord's center line
                    coord1 = central_coord + ecc_x*v
                    if list(coord1)[0] < self.chord.coordinates[0][0] or\
                        list(coord1)[0] > self.chord.coordinates[1][0]:
                        coord1 = central_coord
                    # Web's eccentricity node on chord's edge
                    coord2 = coord1 - ecc_y*u
                    if web.angle == math.radians(90):
                        coord2 = central_coord - [0, ecc_y]
                    # Move chord's coordinate to newly calculated location
                    web.top_coord = list(coord2)
                # BOTTOM CHORD
                else:                   
                    if -math.radians(90) < web.angle < 0 or web.bot_loc == 0 and web.top_loc == 0:
                        ecc_x = -ecc_x
                    # Web's eccentricity node on chord's center line
                    coord1 = central_coord + ecc_x*v
                    if list(coord1)[0] < self.chord.coordinates[0][0] or\
                        list(coord1)[0] > self.chord.coordinates[1][0]:
                        coord1 = central_coord
                    # Web's eccentricity node on chord's edge
                    coord2 = coord1 + ecc_y*u
                    if web.angle == math.radians(90):
                        coord2 = central_coord + [0, ecc_y]
                    # Move chord's coordinate to newly calculated location
                    web.bot_coord = list(coord2)
                                       
                # Calculate web's new nodal locations 
                web.calc_nodal_coordinates(num_elements)
                # Add eccentricity nodes' coordinates to joint and chord
                # Node on the chord's center line needs to be exactly on the
                # center line -> calculate coordinates y-coord with
                # chord's shape function
                coord1 = list([coord1[0], self.chord.shape(coord1[0])])
                # Node coordinate on chord's center line
                self.add_node_coord(coord1)
                self.chord.add_node_coord(coord1)
                # Node coordinate on web's end
                self.add_node_coord(list(coord2))                
                      
        # Move existing nodes, nodes in order [centre, web, chord, web, chord]
        if self.nodes:
            for i, node in enumerate(self.nodes.values()):
                # Move only central and chrod's nodes
                # webs nodes are moved earlier
                if i%2 == 0:
                    node.x = self.nodal_coordinates[i]
        
        # Set joint type
        if len(self.nodal_coordinates) == 0:
            pass
        elif len(self.nodal_coordinates) == 2:
            self.joint_type = 'Y'
        elif len(self.nodal_coordinates) == 4:
            self.joint_type = 'N'
        elif len(self.nodal_coordinates) == 5:
            self.joint_type = 'K'
        elif len(self.nodal_coordinates) == 6:
            self.joint_type = 'KT'
        else:
            raise ValueError('Too many nodes for one joint')
            
            
                
    def generate(self, fem_model):
        """
        """
        self.generate_nodes(fem_model)
        self.add_section_and_material(fem_model)
        self.generate_eccentricity_elements(fem_model)
        
                
    def generate_nodes(self, fem_model):
        """
        Generates eccentricity nodes to chord
        :param fem_model: FrameFEM -instance, model where nodes are added
        """
        if not self.is_generated:
            self.is_generated = True
            for coord in self.nodal_coordinates:
                try:
                    idx = fem_model.nodal_coords.index(coord)
                    self.nodes[idx] = fem_model.nodes[idx]
                    if coord == self.coordinate:
                        self.cnode = fem_model.nodes[idx]
                except ValueError:
                    pass
    
    def add_section_and_material(self, fem_model):
        """
        Add's eccentricity elements section and material properties
        to the fem model.
        Eccentricity elements are modelled as infinitely rigid elements
        """
        
        
        self.idx = len(fem_model.materials)
        self.idx = self.chord.mem_id
        # Material(E, nu, rho)
        fem_model.add_material(210e3, 0.3, 7850e-9)
        # BeamSection(A, Iy)
        sect = BeamSection(1e20, 1e20)
        fem_model.add_section(sect)
        
    
    
    def generate_eccentricity_elements(self, fem_model):
        # Create a list of node instances sorted by nodes' x-coordinate
        nodes = [x for _,x in sorted(zip(self.nodal_coordinates,
                                         list(self.nodes.values())))]
        index = fem_model.nels()
        if self.joint_type == 'Y':
            cn, wn = nodes
            self.elements[index] = \
                EBBeam(cn,
                        wn,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            
        elif self.joint_type == 'N':
            cn1, wn1, cn2, wn2 = nodes
            self.elements[index] = \
                EBBeam(cn1,
                        wn1,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            # Chrod and web
            self.elements[index+1] = \
                EBBeam(cn2,
                        wn2,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+1])
            self.element_ids.append(index+1)
            
            
        elif self.joint_type == 'K':
            print(self.nodal_coordinates)
            cn1, wn1, cn2, wn2, cn3 = nodes
            # Chrod and web
            self.elements[index] = \
                EBBeam(cn1,
                        wn1,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            # Chrod and web
            self.elements[index+1] = \
                EBBeam(cn3,
                        wn2,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+1])
            self.element_ids.append(index+1)
            
        elif self.joint_type == 'KT':
            cn1, wn1, cn2, wn2, cn3, wn3 = nodes
            # Chrod and web
            self.elements[index] = \
                EBBeam(cn1,
                        wn1,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            # Central and web
            self.elements[index+1] = \
                EBBeam(cn2,
                        wn2,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+1])
            self.element_ids.append(index+1)
            # Chrod and web
            self.elements[index+2] = \
                EBBeam(cn3,
                        wn3,
                        fem_model.sections[self.idx], 
                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+2])
            self.element_ids.append(index+2)
            
        
    
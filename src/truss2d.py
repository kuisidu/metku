# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:11:35 2018

@author: huuskoj
"""
from frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, Support
from fem.frame.elements.ebbeam import EBBeam
from fem.frame.elements.eb_semi_rigid_beam import EBSemiRigidBeam
from fem.frame.frame_fem import FrameFEM

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
            if not j1.jid:
                j1.jid = len(self.joints)
                self.joints[j1.jid] = j1
            if not j2.jid:
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



class TrussMember(FrameMember):
    def __init__(self, coordinates, mem_id="",
                 profile="SHS 50x5", material="S355"):
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

    @property
    def angle(self):
        """
        Returns member's angle in radians
        """
        x0, y0 = self.coordinates[0]
        x1, y1 = self.coordinates[1]
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
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'top_chord'


class BottomChord(TrussMember):
    def __init__(self, coordinates, mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)
        self.mtype = 'bottom_chord'


class TrussWeb(TrussMember):
    def __init__(self, bot_loc, top_loc, mem_id="",
                 profile="SHS 50x5", material="S355", Sj1=0, Sj2=0):
        super().__init__([[0,0],[2,2]]) # Member's coordinates are changed in add-function
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
        self.joint_type = joint_type
        self.g1 = g1
        self.g2 = g2
        self.nodal_coordinates = []
        self.nodes = {}
        self.elements = {}
        self.webs = {}
        
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
        
    def calc_coord(self, web_id):
        webs = self.webs.values()
        x0, y0 = self.coordinate
        x = 0
        for web in webs:
            if len(self.webs) == 1:
                x = 0
            
            gamma = web.angle - self.chord.angle
            coords1, coords2 = web.coordinates
            if len(self.webs) > 1:
                if self.chord.mtype == 'top_chord':
                    coords = min(coords1, coords2)
                    y = 0.5 * self.chord.h/1000 # mm to m
                    if math.degrees(gamma) < 0:
                        x = -(0.5 * web.h/1000 * math.sin(gamma) + self.g1/2)
                    else:
                        x = 0.5 * web.h/1000 * math.sin(gamma) + self.g1/2
                else:
                    coords = max(coords1, coords2)
                    y = -0.5 * self.chord.h/1000 # mm to m
                    
                    if math.degrees(gamma) < 0:
                        x = 0.5 * web.h/1000 * math.sin(gamma) + self.g1/2
                    else:    
                        x = -1*(0.5 * web.h/1000 * math.sin(gamma) + self.g1/2)
                
        
        """
        if len(self.webs) == 1:
            x = 0
        elif len(self.webs) == 2:
            if self.chord.mtype == 'top_chord':
                x = 0.5 * web.h/1000 * math.sin(gamma) + self.g1
            else:
                x = 0.5 * web.h/1000 * math.sin(gamma) + self.g1
            
        elif len(self.webs) == 3:
            if self.chord.mtype == 'top_chord':
                x = 0.5 * web.h/1000 * math.sin(gamma) + self.g2
            else:
                x = -1*(0.5 * web.h/1000 * math.sin(gamma) + self.g2)
        """
        
        #return [x0, y0]
        if x0-x > self.chord.coordinates[1][0] or x0-x < self.chord.coordinates[0][0]:
            return [x0,y0]
    
        return [x0-x, y0] 
    
    def calc_nodal_coordinates(self, num_elements):
        # Initialize nodal coordinates as an empty array
        self.nodal_coordinates = []
        x0, y0 = self.coordinate
        ecc_y = self.chord.h / 2000
        self.add_node_coord([x0, y0])
        self.chord.add_node_coord([x0, y0])
        if len(self.webs) == 1:
            web = list(self.webs.values())[0]
            if self.chord.mtype == 'top_chord':
                self.add_node_coord([x0, y0 - ecc_y])
                web.top_coord = [x0, y0 - ecc_y]
            else:
                self.add_node_coord([x0, y0 + ecc_y])
                web.bot_coord = [x0, y0 + ecc_y]
            web.calc_nodal_coordinates(num_elements)
                            
        elif len(self.webs) == 2:
            for web in self.webs.values():
                gamma = web.angle - self.chord.angle
                print("GAMMA: ",web.mem_id,  math.degrees(gamma))
                ecc_x = self.g1/2 + web.h/2000 * math.cos(gamma)
                if self.chord.mtype == 'top_chord':
                    if web.angle == math.radians(90):
                        x, y = x0, y0 - ecc_y
                    
                    elif web.angle > 0:
                        x, y = x0 - ecc_x, y0 - ecc_y        
                        
                    else:
                        x, y = x0 + ecc_x, y0 - ecc_y
                    """
                    # Coordinates need to be rotated to align with the chord
                    dx = x * (1 - math.cos(web.angle))
                    dy = x * math.sin(web.angle)
                    if web.angle < 0:
                        dx = -dx
                    else:
                        dy = -dy
                    """
                    
                    web.top_coord = [x, y]
                    
                else:
                    if web.angle == math.radians(90):
                        x, y = x0, y0 + ecc_y
                    elif math.degrees(gamma) < 0:
                        x, y = x0 - ecc_x, y0 + ecc_y
                    else:
                        x, y = x0 + ecc_x, y0 + ecc_y
                    web.bot_coord = [x, y]
                    
                web.calc_nodal_coordinates(num_elements)               
                self.add_node_coord([x, y])
                self.add_node_coord([x, y0])
                self.chord.add_node_coord([x, y0])
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:11:35 2018

@author: huuskoj
"""
from src.frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, Support
from src.framefem.elements import EBBeam, EBSemiRigidBeam
#from fem.elements.eb_semi_rigid_beam import EBSemiRigidBeam
from src.sections.steel import HEA

from src.framefem import FrameFEM, BeamSection
from src.structures.steel.steel_member import SteelMember
from src.eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint, RHSYJoint


import math
import numpy as np
import matplotlib.pyplot as plt

PREC = 3
L_sys = 0.8

class Truss2D(Frame2D):

    def __init__(self, origo=[0,0], simple={}, num_elements=2, fem_model=None):
        super().__init__(num_elements=num_elements, fem_model=fem_model)
        self.truss = self
        self.origo = origo
        self.top_chords = []
        self.bottom_chords = []
        self.webs = {}
        self.joints = {}
        self.simple = {'H0': 0,
                      'H1': 0,
                      'H2': 0,
                      'H3': 0,
                      'L1': 0,
                      'L2': 0,
                      'n': 0,
                       'dx': 0}
        for key in simple:
            self.simple[key] = simple[key]

        if len(simple):
            self.H0 = self.simple['H0']
            self.H1 = self.simple['H1']
            self.H2 = self.simple['H2']
            if not self.simple['H3']:
                self.H3 = self.H1
            else:
                self.H3 = self.simple['H3']
            self.L1 = self.simple['L1']
            if not self.simple['L2']:
                self.L2 = self.L1
            else:
                self.L2 = self.simple['L2']
            if not self.simple['n']:
                self.n = 10
            else:
                self.n = self.simple['n']
            self.dx = self.simple['dx']
            self.generate_chords()
            self.generate_webs()
            
    def generate_chords(self):
        """
        Generates chords
        :return:
        """
        x, y = self.origo

        if self.H2:
            self.add(TopChord([[x,y+self.H0 + self.H1],
                               [x+self.L1,y+self.H0 + self.H2]]))
            self.add(TopChord([[x+self.L1,y+self.H0 + self.H2],
                               [x+self.L1 + self.L2,y+ self.H0 + self.H3]]))
    
        else:
            self.add(TopChord([[x,y+self.H0 + self.H1],
                               [x+self.L1*2,y+self.H0 + self.H1]]))

        self.add(BottomChord([[x+ self.dx, y+self.H0],
                              [x+self.L1 + self.L2 - self.dx, y+self.H0]]))

                 
    def generate_webs(self):
        """
        Generates webs symmetrically
        :return:
        """
        c1 = 0.0
        c2 = 0.0

        # If true first web starts from top chord else from bottom
        if self.n // 2 % 2:
            flip = True
        else:
            flip = False


        for i in range(1, int(self.n/2)+1):
            if i%2 == 0:
                c1 = round(i/self.n, 4)
                if c1 > 1:
                    c1 = 1
                if self.H2 and not flip:
                    c1 = min(1, c1*2)

            elif i == 1 and self.dx:
                c2 = 0

            elif i:
                c2 = round(i/self.n, 4)
                if c2 > 1:
                    c2 = 1
                if flip:
                    c2 = min(1, c2*2)

            if flip:
                self.add(TrussWeb(c1, c2))
                self.add(TrussWeb(1-c1, 1-c2))
            else:
                self.add(TrussWeb(c2, c1))
                self.add(TrussWeb(1-c2, 1-c1))




    def add(self, this):
        """
        Adds object to truss
        :param this: object to be added
        :return:
        """
        # TOP CHORD
        if isinstance(this, TopChord):
            self.top_chords.append(this)
            this.mem_id = len(self.members)
            self.members["T" + str(this.mem_id)] = this
            this.calc_nodal_coordinates(self.num_elements)
            
        # BOTTOM CHORD   
        elif isinstance(this, BottomChord):
            self.bottom_chords.append(this)
            this.mem_id = len(self.members)
            self.members["B" + str(this.mem_id)] = this
            this.calc_nodal_coordinates(self.num_elements)
            
        # WEB
        elif isinstance(this, TrussWeb):
            if this.coord_type == "local":
                # Calculate web's global coordinates
                bottom_chord = self.bottom_chords[0]
                c01 = bottom_chord.to_global(this.bot_loc)
                for top_chord in self.top_chords:
                    if c01[0] <= top_chord.coordinates[1][0]:
                        c02 = top_chord.to_global(this.top_loc)
                        break
                # Change web's coordinates to global coordinates
                this.bot_coord = c01
                this.top_coord = c02
            else:
                # Sort by y-coordinate
                c01, c02 = sorted(this.coordinates, key=lambda x: x[1])
                # Find corresponding bottom chord
                for chord in self.bottom_chords:
                    if chord.point_intersection(c01):
                        bottom_chord = chord
                        this.bot_loc = bottom_chord.global_coord(c01)
                        break 
                else:
                #if UnboundLocalError and not bottom_chord:
                    raise ValueError("Web's coordinates doesn't match bottom chord's coordinates")                     
                # Find corresponding top chord
                for chord in self.top_chords:   
                    if chord.global_coord(c02) <= 1 and c01[0] <= chord.coordinates[1][0]:
                        top_chord = chord
                        this.top_loc = top_chord.global_coord(c02)
                        break   
                else:
                #if UnboundLocalError and not top_chord:
                    raise ValueError("Web's coordinates doesn't match top chord's coordinates")
            # Check if there's already a joint at given location
            # if not, create a new joint
            j1, j2 = None, None 
            R = 0.03 
            for joint in bottom_chord.joints:    
                if np.sqrt((joint.loc - this.bot_loc)**2)  <= R:
                    j1 = joint
                    this.j1 = j1
                    break
            else:
                j1 = TrussJoint(bottom_chord, this.bot_loc)
                bottom_chord.joints.append(j1)
                this.j1 = j1
                    
                    
            for joint in top_chord.joints:
                if np.sqrt((joint.loc - this.top_loc)**2) <= R:
                    j2 = joint
                    this.j2 = j2
                    break
            else:
                j2 = TrussJoint(top_chord, this.top_loc)
                top_chord.joints.append(j2)
                this.j2 = j2
                
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
     
            self.webs[this.mem_id] = this
            self.members["W" + str(this.mem_id)] = this

        else:
            super().add(this)



    def calculate(self, load_id=2):
        super().calculate(load_id)
        #for joint in self.joints.values():
        #    joint.calculate()

    def generate(self):
        """ Generates the truss
        """
        if not self.is_generated:
            self.is_generated = True
        # Generate members' nodes and elements

        for tc in self.top_chords:
            tc.calc_nodal_coordinates()
            tc.round_coordinates()
            tc.generate(self.f)

        for bc in self.bottom_chords:
            bc.calc_nodal_coordinates()
            bc.round_coordinates()
            bc.generate(self.f)

        for w in self.webs.values():
            w.calc_nodal_coordinates()
            w.round_coordinates()
            w.generate(self.f)


        for member in self.members.values():
            # if member.mtype != "web":
            #     member.calc_nodal_coordinates()
            # member.round_coordinates()
            # member.generate(self.f)
            # # Save coordinates to truss' list
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

    def mirror(self):
        top_coords = self.top_chord.coordinates
        top_coords = [[top_coords[0][0], top_coords[1][1]],
                      [top_coords[1][0], top_coords[0][1]]]
        
        self.top_chord.coordinates = top_coords
        bottom_coords = self.bottom_chord.coordinates
        bottom_coords = [[bottom_coords[0][0], bottom_coords[1][1]],
                      [bottom_coords[1][0], bottom_coords[0][1]]]
        self.bottom_chord.coordinates = bottom_coords
        
class TrussMember(FrameMember):
    def __init__(self,
                 coordinates,
                 mem_id="",
                 profile="SHS 100x5", 
                 material="S355",
                 num_elements=2,
                 Sj1=np.inf,
                 Sj2=np.inf):

        super().__init__(coordinates,
                         mem_id,
                         profile,
                         material,
                         num_elements,
                         Sj1,
                         Sj2)
        self.joints = []
        # FrameMember -objects that chord is connected to
        self.columns = {}
        # List of chords' members between joints
        self.steel_members = []


    # def add_node_coord(self, coord):
    #     """
    #     Adds coordinate to memeber's coordinates and
    #     calculates new coordinates local location
    #     If coordinate is not between member's coordinates, reject it
    #     :param coord: array of two float values, node's coordinates
    #     """
    #
    #     #if self.shape(coord[0]) == coord[1]:
    #     if list(coord) not in self.added_coordinates and self.point_intersection(coord):
    #
    #         self.added_coordinates.append(coord)
    #         self.nodal_coordinates.append(coord)
    #         c0, c1 = self.coordinates
    #         self.nodal_coordinates.sort(
    #             key=lambda c: np.linalg.norm(np.asarray(c0) - np.asarray(c)))
    #         x, y = coord
    #         x1, y1 = self.coordinates[0]
    #         dx = x - x1
    #         dy = y - y1
    #         dz = math.sqrt((dx ** 2 + dy ** 2))
    #         z = dz / self.length
    #         self.loc.append(z)
    #         self.loc.sort()
    
    def calc_nodal_coordinates(self, num_elements=0):
        """
            Calculates node locations along member
            Coordinates are global coordinates
            Adds locations to a list where they can be accessed later
            :param num_elements: int, number of elements to be created
        """
        # if num_elements:
        #     self.num_elements = num_elements
        #
        # self.nodal_coordinates = []
        # coords = [list(c) for c in self.coordinates]
        # start_node, end_node = sorted(coords)
        # x0, y0 = start_node
        # v = np.array([math.cos(self.angle), math.sin(self.angle)])
        #
        # for joint in self.joints:
        #         joint.coordinate = joint.loc
        #         #joint.calc_nodal_coordinates()
        #
        # joint_coords = [joint.coordinate for joint in self.joints]
        # joint_coords.sort()
        #
        # # calculate nodal coordinates between joints
        # if self.mtype != 'web' and len(self.joints) >= 1:
        #     for coord in joint_coords:
        #         x, y = coord
        #         # Calculate distance between joints
        #         L = np.sqrt((x - x0)**2 + (y - y0)**2)
        #         # Calculate nodal coordinates between joints
        #         for i in range(self.num_elements):
        #             node = np.asarray([x0, y0]) + (i*L) / self.num_elements * v
        #             node = list(node)
        #             loc = self.global_coord(node)
        #             if node not in self.nodal_coordinates:
        #                 self.nodal_coordinates.append(node)
        #             if loc not in self.loc:
        #                 self.loc.append(loc)
        #         x0, y0 = x, y
        #     x, y = end_node
        #     # Distance between joints
        #     L = np.sqrt((x - x0)**2 + (y - y0)**2)
        #
        #     for i in range(self.num_elements):
        #             node = np.asarray([x0, y0]) + (i*L) / self.num_elements * v
        #             node = list(node)
        #
        #             loc = self.global_coord(node)
        #             if node not in self.nodal_coordinates:
        #                 self.nodal_coordinates.append(node)
        #             if loc not in self.loc:
        #                 self.loc.append(loc)
        #     if [end_node[0], end_node[1]] not in self.nodal_coordinates:
        #         self.nodal_coordinates.append([end_node[0], end_node[1]])
        #
        #     if 1 not in self.loc:
        #         self.loc.append(1)
        #
        # # Calculate nodal coordinates along member
        # else:
        #     for i in range(self.num_elements):
        #
        #         x, y = np.asarray([x0, y0]) + (i*self.length) / self.num_elements * v
        #         loc = i / self.num_elements
        #         node = [x, y]
        #         if node not in self.nodal_coordinates:
        #             self.nodal_coordinates.append(node)
        #         if loc not in self.loc:
        #             self.loc.append(loc)
        #     if [end_node[0], end_node[1]] not in self.nodal_coordinates:
        #         self.nodal_coordinates.append(end_node)
        #     if 1 not in self.loc:
        #         self.loc.append(1)
        # # move previously created nodes
        # if self.is_generated:
        #     #for joint in self.joints:
        #     #     joint.coordinate = joint.loc
        #     for j, node in enumerate(self.nodes.values()):
        #         node.coord = np.array(self.nodal_coordinates[j])
        #
        # else:
        #     for joint in self.joints:
        #         joint.coordinate = joint.loc
        #         #joint.calc_nodal_coordinates()
        #
        # #self.nodal_coordinates.extend(self.added_coordinates)
        # #self.nodal_coordinates.sort()
        #
        self.nodal_coordinates = []
        self.nodal_coordinates.extend(self.added_coordinates)
        if num_elements:
            self.num_elements = num_elements

        joint_locs = [j.loc for j in self.joints]

        if len(joint_locs):
            if 0 not in joint_locs:
                joint_locs.append(0)
            if 1 not in joint_locs:
                joint_locs.append(1)
            joint_locs.sort()
            for i in range(len(joint_locs) -1):

                j_loc_0 = joint_locs[i]
                j_loc_1 = joint_locs[i + 1]

                d_loc = j_loc_1 - j_loc_0
                d_loc /= self.num_elements

                for j in range(self.num_elements):
                    loc = j_loc_0 + d_loc * j
                    coord = list(self.to_global(loc))
                    if loc <= 1 and coord not in self.nodal_coordinates:
                        self.nodal_coordinates.append(coord)
            last_joint_coord = list(self.to_global(joint_locs[-1]))
            if last_joint_coord not in self.nodal_coordinates:
                self.nodal_coordinates.append(last_joint_coord)
            if list(self.to_global(1)) not in self.nodal_coordinates:
                self.nodal_coordinates.append(list(self.to_global(1)))

            c0, c1 = self.coordinates
            self.nodal_coordinates.sort(key=lambda c: np.linalg.norm(np.asarray(c0) - np.asarray(c)))
        else:
            super().calc_nodal_coordinates(num_elements)

        if self.is_generated:
            nodes = list(self.nodes.values())
            c0, c1 = self.coordinates
            nodes = sorted(nodes, key= lambda n: np.linalg.norm(np.asarray(c0) - np.asarray(n.coord)))
            for j, node in enumerate(nodes):
                node.coord = np.array(self.nodal_coordinates[j])





    def check_cross_section(self):
        super().check_cross_section()
        self.steel_members = []
        if self.mtype != 'web':  # If member is chord
            # Sort joints using joints' local coordinate
            self.joints.sort(key=lambda x: x.loc)
            for i in range(len(self.joints) -1):
                j0 = self.joints[i]
                j1 = self.joints[i + 1]
                x0, y0 = j0.coordinate
                x1, y1 = j1.coordinate
                L = np.sqrt((x0 - x1)**2 + (y0 - y1)**2)
                steel_mem = SteelMember(self.cross_section, L, [L_sys, L_sys])
                for id in range(j0.cnode.nid, j1.cnode.nid):
                    N, V, M = self.nodal_forces[id]
                    steel_mem.add_section(ned=N, myed=M, vyed=V)
                self.steel_members.append(steel_mem)

            for steel_member in self.steel_members:
                for id in range(j0.cnode.nid, j1.cnode.nid):
                    N, V, M = self.nodal_forces[id]
                    steel_member.add_section(ned=N, myed=M, vyed=V)
                r = []
                r_sect = steel_member.check_sections()
                r_buckling = steel_member.check_buckling()
                r_lt_buckling = steel_member.check_LT_buckling()
                r.extend(r_sect)
                r.extend(r_buckling)
                r.append(r_lt_buckling)


        
    def local(self, value):
        """
        Returns chord's local coordinate in global coordinates
        :param value: local coordinate, 0 <= value <= 1
        """      
        if 0 >= value >= 1:
            raise ValueError('Value must be between 0 and 1')
        start_node, end_node = self.coordinates
        
        v = np.array([math.cos(self.angle), math.sin(self.angle)])
        
        x0, y0 = start_node
        Lx = self.length * value

        return list(np.asarray(start_node) + Lx * v)
    
    
    def global_coord(self, value):
        """
        Return's global coordinate in chord's local coordinates
        """
        x, y = value
        x0, y0 = self.coordinates[0]
        
        L = np.sqrt((x0-x)**2 + (y0-y)**2)
        local = L / self.length
        
        return local
        
    
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
        except ZeroDivisionError:
            k = 0

        return k*(x-x0) + y0



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
       
        for i in range(len(self.nodes) - 1):
            
            n1 = node_ids[i]
            n2 = node_ids[i + 1]       
            # EBSemiRigid -element, first element
            if i == 0:
                self.elements[index] = \
                    EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                                    self.cross_section, \
                                    self.material, \
                                    rot_stiff=[self.Sj1, np.inf])
            # EBSemiRigid -element, last element
            elif i == (len(self.nodes) - 2):
                self.elements[index] = \
                    EBSemiRigidBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                                    self.cross_section, \
                                    self.material, \
                                    rot_stiff=[np.inf, self.Sj2])
            # EBBeam -elements
            else:              
                self.elements[index] = \
                    EBBeam(fem_model.nodes[n1], fem_model.nodes[n2], \
                           self.cross_section, \
                           self.material)
        
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            index += 1

    def plot(self, print_text=True, color=True):

        X = self.coordinates
        if color:
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
            #rot = self.angle

        x = (X[0][0] + X[1][0]) / 2
        y = (X[1][1] + X[0][1]) / 2
        horzalign = 'center'
        vertalign = 'center'
        if print_text:
            plt.text(x, y, str(self.mem_id) + ": " + self.profile,
                     rotation=rot, horizontalalignment=horzalign, rotation_mode='anchor')
        else:
            plt.text(x, y, str(self.mem_id),
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




class Chord(FrameMember):

    def __init__(self, coordinates, **kwargs):
        super().__init__(coordinates, **kwargs)

        self.members = []
        self.joints = []


    def generate(self, fem_model):

        locs = [0, 1]
        for joint in self.joints:
            if joint.loc not in locs:
                locs.append(joint.loc)

        locs.sort()
        for i in range(len(locs) -1):
            c0 = self.to_global(locs[i])
            c1 = self.to_global(locs[i+1])
            mem = FrameMember([c0, c1])
            self.members.append(mem)
            mem.generate(fem_model)





    def check_cross_section(self):
        for mem in self.members:
            mem.check_cross_section()

class TopChord(TrussMember):
    def __init__(self,
                 coordinates,
                 **kwargs):
        super().__init__(coordinates,
                         **kwargs)

        self.mtype = 'top_chord'




class BottomChord(TrussMember):

    def __init__(self,
                 coordinates,
                 **kwargs):
        super().__init__(coordinates,
                         **kwargs)
        self.mtype = 'bottom_chord'

    @property
    def perpendicular(self):
        return -1 * super().perpendicular

class TrussWeb(TrussMember):
    def __init__(self,
                 bot_loc,
                 top_loc,
                 coord_type="local",
                 mem_id="",
                 profile="SHS 50x2",
                 material="S355",
                 num_elements=2,
                 Sj1=0,
                 Sj2=0):

        if coord_type == "global" or coord_type == "g":
            super().__init__([bot_loc, top_loc], profile=profile,
                 material=material, num_elements=num_elements)  
                        # Global coordinates
            self.__top_coord = top_loc
            self.__bot_coord = bot_loc
            # Chord's local coordinates
            self.bot_loc = None
            self.top_loc = None
        else:
            # Member's coordinates are changed in add-function
            super().__init__([[0,0],[1,1]], profile=profile,
                 material=material, num_elements=num_elements)
            # Chord's local coordinates
            self.bot_loc = bot_loc
            self.top_loc = top_loc
            # Global coordinates
            self.__top_coord = [1,1]
            self.__bot_coord = [0,0]


        self.coord_type = coord_type   
        self.top_coord = self.coordinates[1]
        self.bot_coord = self.coordinates[0]
        
        self.j1 = None
        self.j2 = None
        
        self.Sj1 = Sj1
        self.Sj2 = Sj2
        self.mtype = 'web'
        
     
    @property
    def bot_coord(self):
        return self.__bot_coord
    
    @bot_coord.setter
    def bot_coord(self, val):
        self.__bot_coord = val
        if self.n1:
            self.n1.coord = val
        self.coordinates = [val, self.top_coord]
        
        
    @property
    def top_coord(self):
        return self.__top_coord
    
    @top_coord.setter
    def top_coord(self, val):
        self.__top_coord = val
        if self.n2:
            self.n2.coord = val
        self.coordinates = [self.bot_coord, val]


class TrussJoint:
    def __init__(self, chord, loc, joint_type="N", g1=30, g2=30):

        self.chord = chord
        self.chord_elements = {}
        self.jid = None
        self.__loc = loc
        self.__coordinate = chord.to_global(loc)
        self.cnode = None # central node
        self.chord_nodes = []
        self.joint_type = joint_type
        self.__g1 = g1
        self.__g2 = g2
        self.nodal_coordinates = {}
        self.nodes = {}
        self.elements = {}
        self.element_ids = []
        self.webs = {}
        self.is_generated = False
        # Index for eccentricity elements' material and section
        self.idx = None
        self.rhs_joint = None
        self.r = []
                
    def add_node_coord(self, coord):
        coord = [round(c, PREC) for c in coord]
        if coord not in self.nodal_coordinates:
            self.nodal_coordinates.append(coord)
    
    def add_web(self, web, num_elements):
        self.webs[web.mem_id] = web
        #if self.is_generated:
        self.calc_nodal_coordinates()

    @property
    def loc(self):
        return round(self.__loc, PREC)
    
    @loc.setter
    def loc(self, val):
        val = max(min(val, 1), 0)
        self.__loc = round(val, PREC)
        if self.chord.mtype == 'top_chord':
            for web in self.webs.values():
                web.top_loc = self.loc
        else:
            for web in self.webs.values():
                web.bot_loc = self.loc

        self.calc_nodal_coordinates()



    @property
    def coordinate(self):     
        return [round(c, PREC) for c in self.chord.to_global(self.loc)]

    @coordinate.setter
    def coordinate(self, val):
        # (x, y) -coordinate
        if isinstance(val, list):            
            local = self.chord.global_coord(val)
            self.coordinate = local
        # local coordinate
        elif isinstance(val, float):
            if 0 <= val <= 1:
                self.loc = val
                self.calc_nodal_coordinates()


    @property
    def e(self):
        if self.joint_type == 'K':
            w1, w2 = list(self.webs.values())
            h1 = w1.h
            h2 = w2.h
            h0 = self.chord.h

            theta1 = min(abs(w1.angle - self.chord.angle),
                         abs(np.pi - w1.angle - self.chord.angle))
            theta2 = min(abs(w2.angle - self.chord.angle),
                         abs(np.pi - w2.angle - self.chord.angle))

            e = h1/(2*np.sin(theta1)) + h2/(2*np.sin(theta2)) + self.g1
            e *= (np.sin(theta1) * np.sin(theta2)) / np.sin(theta1 + theta2)
            e -= h0/2

            return e

        return 0


    @e.setter
    def e(self, val):

        if self.joint_type == 'K':
            w1, w2 = self.webs.values()
            h1 = w1.h
            h2 = w2.h
            h0 = self.chord.h

            theta1 = min(abs(w1.angle - self.chord.angle),
                         abs(np.pi - w1.angle - self.chord.angle))
            theta2 = min(abs(w2.angle - self.chord.angle),
                         abs(np.pi - w2.angle - self.chord.angle))

            g = val + h0/2
            g *= np.sin(theta1 + theta2) / (np.sin(theta1) * np.sin(theta2))
            g -= h1 / (2*np.sin(theta1))
            g -= h2 / (2 * np.sin(theta2))
            self.g1 = g



    @property
    def g1(self):
        return self.__g1
    
    @g1.setter
    def g1(self, val):
        self.__g1 = val
        self.calc_nodal_coordinates()
        
    @property
    def g2(self):
        return self.__g2
    
    @g2.setter
    def g2(self, val):
        self.__g2 = val
        self.calc_nodal_coordinates()
        
    def calculate(self):
        """
        """
        NEd0 = [elem.axial_force for elem in self.chord.elements.values() if
                self.cnode in elem.nodes]
        # item for sublist in list for item in sublist
        NEd0 = [N for forces in NEd0 for N in forces]
        maxN = max(NEd0)
        minN = min(NEd0)
        if abs(minN) > maxN:
            NEd0 = minN
        else:
            NEd0 = maxN
        
        if self.joint_type == "K" or self.joint_type == 'N':
            NEd1, NEd2 = [w.ned for w in self.webs.values()]
            wNEd = np.array([NEd1, NEd2])
            # Resistances
            r1 = self.chord_face_failure()
            r2 = self.punching_shear()
            r3 = self.brace_failure()
            r4 = self.chord_shear()
            # Stress ratios
            r1 = abs(wNEd / r1)
            r2 = abs(wNEd / r2)
            r3 = abs(wNEd / r3)
            r5 = abs(NEd0 / r4[1]) # chord
            r4 = abs(wNEd / r4[0])

            #self.r = [r1, r2, r3, r4]
            maxR_for_webs = np.max((r1, r2, r3, r4), axis=0)
            maxR = np.append(maxR_for_webs, r5)
            self.r = maxR
            #return maxR

        elif self.joint_type == 'Y':
            NEd1 = [w.ned for w in self.webs.values()][0]
            
            r1 = self.chord_face_failure()
            r2 = self.chord_web_buckling()
            r3 = self.brace_failure()
            
            r1 = abs(NEd0 / r1)
            r2 = abs(NEd0 / r2)
            r3 = abs(NEd1 / r3)
            #self.r = [r1, r2, r3]

            maxR = [r3, max(r1, r2)]
            self.r = maxR
            #return maxR
        elif self.joint_type == 'KT':
            pass
        
    def chord_face_failure(self):
        if self.rhs_joint:
            results = self.rhs_joint.chord_face_failure() 
            return results / 1000
    
    def punching_shear(self):
        if self.rhs_joint and self.joint_type == 'K':
            results = self.rhs_joint.punching_shear() 
            return results / 1000
        
    def brace_failure(self):
        if self.rhs_joint:           
            results = self.rhs_joint.brace_failure()
            return results / 1000 # N to kN
    
    def chord_shear(self):
        if self.rhs_joint and self.joint_type == 'K':
            V0 = [elem.shear_force for elem in self.chord.elements.values() if
                  self.cnode in elem.nodes]
            # item for sublist in list for item in sublist
            V0 = [abs(V) for forces in V0 for V in forces]
            V0 = max(V0)
            self.rhs_joint.V0 = V0 * 1e3 # kN to N
            results = [r / 1000 for r in self.rhs_joint.chord_shear()]
            return results
        
    def chord_web_buckling(self):
        if self.rhs_joint and self.joint_type == 'Y':
            return self.rhs_joint.chord_web_buckling() / 1000
        
    def calc_nodal_coordinates(self):
        """
        Calculates eccentricity nodes coordinates and adds them to chord's
        node-dict, changes webs end coordinates and saves all calculated
        nodal coordinates to own nodal_coordinates -list
        
        coord1 : eccentricity coordinate on chord
        coord2 : web's new end coordinate
        v : chord's position vector
        u : vector perpendicular to v
        
        1: Calc ecc_y, distance from chord's center line to edge
        2: Calc ecc_x, distance from joint to web's edge
        3: Calc coord1, ecc_x * v
        4: Calc coord2, coord1 + ecc_y * u
        5: Move web to coord2
        6: Save coordinates 
        
        """
        
        # Chord's position vector
        v = np.array([math.cos(self.chord.angle), math.sin(self.chord.angle)])
        # Vector perpendicular to v
        u = np.array([-math.sin(self.chord.angle), math.cos(self.chord.angle)])
        # Remove previously calculated coordinates from chord
        for coord in self.nodal_coordinates:
            if coord in self.chord.added_coordinates:
                self.chord.added_coordinates.remove(coord)
            #if coord in self.chord.nodal_coordinates and \
            #    coord not in self.chord.coordinates:
            #    self.chord.nodal_coordinates.remove(coord)
        # Initialize nodal coordinates as an empty array
        self.nodal_coordinates = []
        central_coord = np.asarray(self.coordinate)
        # Eccentricity along y-axis
        ecc_y = self.chord.h / 2

        # Add joint's central coordinate to joint's and chord's list
        self.add_node_coord(list(central_coord))
        # If joint type is Y or T
        if len(self.webs) == 1:
            web = list(self.webs.values())[0]
            # angle between chord and web
            theta = abs(self.chord.angle - web.angle)
            # Eccentricity along x-axis
            ecc_x = abs(self.g1 * math.cos(self.chord.angle)
                    + web.h/2 / math.sin(theta))

            ecc_x = round(ecc_x, PREC)
            # TOP CHORD
            if self.chord.mtype == 'top_chord':
                """
                if self.loc == 0:
                    coord1 = list(central_coord + v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) - u *ecc_y)
                elif self.loc == 1:
                    coord1 = list(central_coord - v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) - u *ecc_y)
                else:
                """
                coord1 = list(central_coord)
                coord2 = list(np.asarray(coord1) - u *ecc_y)
                self.add_node_coord(coord2)
                self.chord.add_node_coord(coord1)
                web.top_coord = coord2
                
            # BOTTOM CHORD
            else:
                """
                if self.loc == 0:
                    coord1 = list(central_coord + v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) + u *ecc_y)
                elif self.loc == 1:
                    coord1 = list(central_coord - v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) + u *ecc_y)
                else:
                """
                coord1 = list(central_coord)
                coord2 = list(np.asarray(coord1) + u *ecc_y)

                self.add_node_coord(coord2)
                self.chord.add_node_coord(coord1)
                web.bot_coord = coord2
            web.calc_nodal_coordinates()    
            
        # If joint type is K or KT
        elif len(self.webs) >= 2:
            self.chord.add_node_coord(list([central_coord[0], self.chord.shape(central_coord[0])]))
            # Calculate coordinates n-times because web's angle changes
            n = 10
            for i in range(n):                
                # Iterate through every web member in the joint
                for web in self.webs.values():
                    # angle between chord and web
                    theta = abs(self.chord.angle - web.angle)
                    # Eccentricity along x-axis
                    ecc_x = abs(self.g1/2 * math.cos(self.chord.angle)
                            + web.h/2 / math.sin(theta))
                    # If joint type is KT, there's two gap values g1 and g2
                    if len(self.webs) == 3:
                        # Find the leftmost, rightmost and middle web with x-values
                        # leftmost is the one with smallest x (miv_val)
                        # rightmost is the one with highest x (max_val)
                        # middle is the one between those two
                        min_val = min([x.bot_coord[0] for x in self.webs.values()])
                        max_val = max ([x.bot_coord[0] for x in self.webs.values()])
                        middle_web = [x for x in self.webs.values() if x.bot_coord[0] != min_val and x.bot_coord[0] != max_val][0]
                        # eccentricity for the leftmost web
                        if web.bot_coord[0] == min_val:
                            ecc_x = self.g1 *math.cos(self.chord.angle) + web.h/2 / math.sin(theta) + middle_web.h/2
                        # eccentricity for the rightmost web
                        elif web.bot_coord[0] == max_val:
                            ecc_x = self.g2 *math.cos(self.chord.angle) + web.h/2 / math.sin(theta) + middle_web.h/2
                        # eccentricity for the middle web
                        else:
                            ecc_x = 0
                    # TOP CHORD
                    if self.chord.mtype == 'top_chord':
                        # if the web is the leftmost web
                        #if web.bot_coord[0] == min([x.bot_coord[0] for x in self.webs.values()]):   
                        if web.bot_loc == min([x.bot_loc for x in self.webs.values()]):
                            ecc_x = -ecc_x                    
                        # Web's eccentricity node on chord's center line
                        coord1 = central_coord + ecc_x*v
                        coord1 = [round(c, PREC) for c in coord1]
                        # Check that coordinate is between chord's coordinates
                        if list(coord1)[0] < self.chord.coordinates[0][0] or\
                            list(coord1)[0] > self.chord.coordinates[1][0]:
                            coord1 = central_coord
                        # Web's eccentricity node on chord's edge
                        coord2 = coord1 - ecc_y*u
                        # Move chord's coordinate to newly calculated location
                        web.top_coord = list(coord2)
                    # BOTTOM CHORD
                    else:                  
                        # if the web is the leftmost web
                        #if web.top_coord[0] == min([x.top_coord[0] for x in self.webs.values()]):
                        if web.top_loc == min([x.top_loc for x in self.webs.values()]):
                            ecc_x = -ecc_x
                        # Web's eccentricity node on chord's center line
                        coord1 = central_coord + ecc_x*v
                        # Check that coordinate is between chord's coordinates
                        if list(coord1)[0] < self.chord.coordinates[0][0] or\
                            list(coord1)[0] > self.chord.coordinates[1][0]:
                            coord1 = central_coord
                        # Web's eccentricity node on chord's edge
                        coord2 = coord1 + ecc_y*u
                        # Move chord's coordinate to newly calculated location
                        coord2 = [round(c, 3) for c in coord2]
                        web.bot_coord = list(coord2)                    
                    # Calculate web's new nodal locations 
                    web.calc_nodal_coordinates()
                    # Add eccentricity nodes' coordinates to joint and chord
                    # Node on the chord's center line needs to be exactly on the
                    # center line, so we need to calculate the y-coord with
                    # chord's shape function
                    if i == n-1:
                        coord1 = list([coord1[0], self.chord.shape(coord1[0])])
                        coord1 = [round(c,PREC) for c in coord1]
                        coord2 = [round(c,PREC) for c in coord2]
                        # Node coordinate on chord's center line
                        self.add_node_coord(coord1)
                        self.chord.add_node_coord(coord1)
                        # Node coordinate on web's end
                        self.add_node_coord(list(coord2))
                        
        # Move existing nodes
        if len(self.nodes):
            for i, node in enumerate(self.nodes.values()):
                node.coord = self.nodal_coordinates[i]
          
        # Set joint type
        if len(self.nodal_coordinates) == 0:
            pass
        elif len(self.nodal_coordinates) == 2:
            self.joint_type = 'Y'
            self.rhs_joint = RHSYJoint(self.chord.cross_section,
                                       [w.cross_section for w in self.webs.values()][0],
                                       theta)
        elif len(self.nodal_coordinates) == 3:
            self.joint_type = 'Y'
            self.rhs_joint = RHSYJoint(self.chord.cross_section,
                                       [w.cross_section for w in self.webs.values()][0],
                                       theta)
        elif len(self.nodal_coordinates) == 4:
            self.joint_type = 'N'
            self.rhs_joint = RHSKGapJoint(self.chord.cross_section,
                             [w.cross_section for w in self.webs.values()],
                             [math.degrees(self.chord.angle - w.angle) for w in self.webs.values()],
                             self.g1) # m to mm
        elif len(self.nodal_coordinates) == 5:
            self.joint_type = 'K'
            self.rhs_joint = RHSKGapJoint(self.chord.cross_section,
                                         [w.cross_section for w in self.webs.values()],
                                         [math.degrees(self.chord.angle - w.angle) for w in self.webs.values()],
                                         self.g1) # m to mm
        elif len(self.nodal_coordinates) == 6:
            self.joint_type = 'KT'
        else:
            raise ValueError('Too many nodes for one joint')
            


    def calc_nodal_coordinates(self):
        """

        :return:
        """
        coords = []

        if len(self.webs) == 1:
            self.joint_type = "Y"
            w1 = list(self.webs.values())[0]
            c0 = self.chord.to_global(self.loc)
            c1 = c0 + self.chord.h / 2 * self.chord.perpendicular
            if self.chord.mtype == 'top_chord':
                w1.top_coord = c1
            else:
                w1.bot_coord = c1

            coords = [c0, c1]

            self.nodal_coordinates['c0'] = c0
            self.nodal_coordinates['c1'] = c1

        elif len(self.webs) == 2:
            self.joint_type = 'K'
            c0 = self.chord.to_global(self.loc)
            ec = c0 - self.e * self.chord.perpendicular

            if self.chord.mtype == 'top_chord':
                w1, w2 = sorted(self.webs.values(), key= lambda w: w.bot_loc)
                c01 = self.chord.line_intersection([w1.bot_coord, ec])
                c02 = self.chord.line_intersection([w2.bot_coord, ec])

            else:
                w1, w2 = sorted(self.webs.values(), key=lambda w: w.top_loc)
                c01 = self.chord.line_intersection([w1.top_coord, ec])
                c02 = self.chord.line_intersection([w2.top_coord, ec])

            theta1 = w1.angle - self.chord.angle
            d1 = abs(self.chord.h / 2   / np.tan(theta1))
            c01 -= self.chord.unit * d1

            theta2 = w2.angle - self.chord.angle
            d2 = abs(self.chord.h / 2 / np.tan(theta2))
            c02 += self.chord.unit * d2

            if self.chord.mtype == 'top_chord':
                c1 = np.asarray(c01) + self.chord.h / 2 * self.chord.perpendicular
                c2 = np.asarray(c02) + self.chord.h / 2 * self.chord.perpendicular
                w1.top_coord = c1
                w2.top_coord = c2
            else:
                c1 = np.asarray(c01) + self.chord.h / 2 * self.chord.perpendicular
                c2 = np.asarray(c02) + self.chord.h / 2 * self.chord.perpendicular
                w1.bot_coord = c1
                w2.bot_coord = c2

            coords = [c0, c01, c02, c1, c2]

            self.nodal_coordinates['c0'] = c0
            self.nodal_coordinates['c01'] = c01
            self.nodal_coordinates['c02'] = c02
            self.nodal_coordinates['c1'] = c1
            self.nodal_coordinates['c2'] = c2

            #
            #
            # plt.plot([c0[0], ec[0]], [c0[1], ec[1]])
            # w1.plot()
            # w2.plot()
            # plt.scatter(*w1.line_intersection(w2), marker='x')
            # plt.scatter(*ec)
            # plt.scatter(*c0, marker='s')
            # plt.scatter(*c01)
            # plt.scatter(*c02)
            # plt.scatter(*(np.asarray(c01) + self.chord.h / 2 * self.chord.perpendicular))
            # plt.scatter(*(np.asarray(
            #     c02) + self.chord.h / 2 * self.chord.perpendicular))
            # plt.axis('equal')
            # plt.show()

        elif len(self.webs) == 3:
            self.joint_type = 'KT'
            c0 = self.chord.to_global(self.loc)
            ec = c0 - self.e * self.chord.perpendicular
            if self.chord.mtype == 'top_chord':
                w1, w2, w3 = sorted(self.webs.values(), key= lambda w: w.bot_loc)
            else:
                w1, w2, w3 = sorted(self.webs.values(), key=lambda w: w.top_loc)
            c01 = self.chord.line_intersection(w1)
            theta1 = w1.angle - self.chord.angle
            d1 = abs(self.chord.h / 2   / np.tan(theta1))
            c01 -= self.chord.unit * d1

            c02 = self.chord.line_intersection(w2)
            theta2 = w2.angle - self.chord.angle
            d2 = abs(self.chord.h / 2 / np.tan(theta2))
            c02 += self.chord.unit * d2

        for c in coords:
            self.chord.add_node_coord(c)

        if self.is_generated:
            for key, coord in self.nodal_coordinates.items():
                self.nodes[key].coord = coord


    def round_coordinates(self, prec=PREC):
        for name, coord in self.nodal_coordinates.items():
            self.nodal_coordinates[name] = [round(c, prec) for c in coord]
            
    def generate(self, fem_model):
        """
        """
        self.round_coordinates()
        self.generate_nodes(fem_model)
        self.add_section_and_material(fem_model)
        self.generate_eccentricity_elements(fem_model)
        self.get_chord_elements()
        
                
    def generate_nodes(self, fem_model):
        """
        Generates eccentricity nodes to chord
        :param fem_model: FrameFEM -instance, model where nodes are added
        """
        if not self.is_generated:
            self.is_generated = True
            for name, coord in self.nodal_coordinates.items():
                idx = fem_model.nodal_coords.index(coord)
                self.nodes[name] = fem_model.nodes[idx]
                if coord == self.coordinate:
                    self.cnode = fem_model.nodes[idx]
                if coord in self.chord.nodal_coordinates:
                    self.chord_nodes.append(fem_model.nodes[idx])
        # sort chord nodes from leftmost to rightmost
        self.chord_nodes = sorted(self.chord_nodes, key= lambda node: node.x)
                    
    def get_chord_elements(self):
        """
        Gets the joint's elements which are on chord
        """
        if len(self.chord_nodes) > 2:
            for i in range(2):
                nodes = self.chord_nodes[i:i+2]
                for key in self.chord.elements.keys():
                    element = self.chord.elements[key]
                    if nodes[0] in element.nodes and nodes[1] in element.nodes:
                        self.chord_elements[key] = element
        
        
    def add_section_and_material(self, fem_model):
        """
        Add's eccentricity elements section and material properties
        to the fem model.
        Eccentricity elements are modelled as infinitely rigid elements
        """
        self.idx = len(fem_model.sections)
        
        # Material(E, nu, rho)
        fem_model.add_material(210e3, 0.3, 7850e-9)
        # BeamSection(A, Iy) [m^2, m^4]
        sect = BeamSection(38880*1e-6, 6299657109*1e-12)
        fem_model.add_section(sect)
        
        # USE CHORD'S PROPERTIES
        #self.idx = self.chord.mem_id
  
    
    def generate_eccentricity_elements(self, fem_model):
        """
        Generates eccentricity elements
        """


        section = HEA(1000)

        index = len(fem_model.elements)
        if self.joint_type == 'Y':
            self.elements[index] = EBBeam(self.nodes['c0'],
                                        self.nodes['c1'],
                                        section,
                                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)


        elif self.joint_type == 'K':
            # Chrod and web
            self.elements[index] = EBBeam(self.nodes['c01'],
                                        self.nodes['c1'],
                                        section,
                                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            # Chrod and web
            self.elements[index+1] = EBBeam(self.nodes['c02'],
                                            self.nodes['c2'],
                                            section,
                                            fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+1])
            self.element_ids.append(index+1)

        # TODO: KT
        # elif self.joint_type == 'KT':
        #     cn1, cn2, cn3 = chord_nodes
        #     wn1, wn2, wn3 = web_nodes
        #     # Chrod and web
        #     self.elements[index] = EBBeam(cn1,
        #                                 wn1,
        #                                 section,
        #                                 fem_model.materials[self.idx])
        #     # Add it to fem model
        #     fem_model.add_element(self.elements[index])
        #     self.element_ids.append(index)
        #     # Central and web
        #     self.elements[index+1] = EBBeam(cn2,
        #                                     wn2,
        #                                     section,
        #                                     fem_model.materials[self.idx])
        #     # Add it to fem model
        #     fem_model.add_element(self.elements[index+1])
        #     self.element_ids.append(index+1)
        #     # Chrod and web
        #     self.elements[index+2] = EBBeam(cn3,
        #                                     wn3,
        #                                     section,
        #                                     fem_model.materials[self.idx])
        #     # Add it to fem model
        #     fem_model.add_element(self.elements[index+2])
        #     self.element_ids.append(index+2)
            
        
    def plot(self, print_text=True, color=True):

        x, y = self.coordinate
        if color:
            if self.is_strong_enough:
                color = 'green'
            else:
                color = 'red'
        else:
            color = 'cyan'

        horzalign = 'center'
        vertalign = 'center'
        plt.scatter(x, y,s=100, c=color)
        plt.text(x, y, str(self.jid), horizontalalignment=horzalign,
                 verticalalignment=vertalign, fontsize=12)
        
    def plot_joint(self, length=0.3, show=True):
        """
        Plots joint
        """
        X, Y = self.coordinate  
        end_line = None
        col = None
        # CHORD
        # Chord's position vector
        v = np.array([math.cos(self.chord.angle), math.sin(self.chord.angle)])
        # Vector perpendicular to v
        u = np.array([-math.sin(self.chord.angle), math.cos(self.chord.angle)])
        if X-length*2 >= self.chord.local(0)[0]:
            start = X - length*2
        else:
            start = self.chord.coordinates[0][0]
            end_line = self.chord.coordinates[0][0]
        
        if X + length*2 <= self.chord.local(1)[0]:
            end = X + length*2
        else:
            end = self.chord.coordinates[1][0]  
            end_line = self.chord.coordinates[1][0]
            
        X0 = [start, X, end]
        Y0 = [Y-(X-start)*math.tan(self.chord.angle),
              Y,
              Y+(end-X)*math.tan(self.chord.angle)]
        X1 = [start, X, end]
        Y1 = [Y-(X-start)*math.tan(self.chord.angle) - self.chord.h/2000,
              Y-self.chord.h/2000,
              Y+(end-X)*math.tan(self.chord.angle)-self.chord.h/2000]
        X2 = [start, X,end]
        Y2 = [Y-(X-start)*math.tan(self.chord.angle)+self.chord.h/2000,
              Y+self.chord.h/2000,
              Y+(end-X)*math.tan(self.chord.angle)+self.chord.h/2000]
        
        plt.plot(X0, Y0, 'k--')
        plt.plot(X1, Y1, 'k')
        plt.plot(X2, Y2, 'k')
        if end_line:
            plt.plot([end_line, end_line],[Y-self.chord.h/2000, Y+self.chord.h/2000],'k')
            
        # COLUMN
        if start in self.chord.columns.keys():
            col = self.chord.columns[start]
            end = start - col.h/2000
            inner = start + col.h/2000
            Y1 = Y- length
            Y2 = Y + length
            if Y2 >= col.coordinates[1][1]:
                Y2 = col.coordinates[1][1]
            # outer
            plt.plot([end, end], [Y1, Y2], 'k')
            # mid
            plt.plot([start, start],[Y1, Y2], 'k--')
            # inner
            plt.plot([inner, inner], [Y1, Y2], 'k')
            
        elif end in self.chord.columns.keys():
            col = self.chord.columns[end]
            out = end + col.h/2000
            inner = end - col.h/2000
            Y1 = Y- length
            Y2 = Y + length
            if Y2 >= col.coordinates[1][1]:
                Y2 = col.coordinates[1][1]
            # outer
            plt.plot([out, out], [Y1, Y2], 'k')
            # mid
            plt.plot([end, end],[Y1, Y2], 'k--')
            # inner
            plt.plot([inner, inner], [Y1, Y2], 'k')
            
        
        # WEBS
        for web in self.webs.values():
            theta = abs(self.chord.angle - web.angle)
            if self.chord.mtype == "top_chord":
                coord = web.top_coord
                k = -1
            else:
                coord = web.bot_coord
                k = 1
            if web.angle < 0:
                X0 = [coord[0], coord[0] - k*length]
                Y0 = [coord[1], coord[1] - k*(math.tan(web.angle)*length)]
                X1 = [coord[0]+web.h/2000 / math.sin(theta) , coord[0]+web.h/2000 / math.sin(theta) - k*length]
                Y1 = [coord[1], coord[1] - k*(math.tan(web.angle)*length)]
                X2 = [coord[0]-web.h/2000 / math.sin(theta) , coord[0]-web.h/2000 / math.sin(theta) - k*length]
                Y2 = [coord[1], coord[1] - k*(math.tan(web.angle)*length)]
            
            elif math.degrees(web.angle) == 90:
                X0 = [coord[0], coord[0]]
                Y0 = [coord[1], coord[1] + k*2*length]
                X1 = [coord[0]+web.h/2000 , coord[0]+web.h/2000]
                Y1 = [coord[1], coord[1] + k*2*length]
                X2 = [coord[0]-web.h/2000 , coord[0]-web.h/2000]
                Y2 = [coord[1], coord[1] + k*2*length]
            else:
                X0 = [coord[0], coord[0] + k*length]
                Y0 = [coord[1], coord[1] + k*(math.tan(web.angle)*length)]
                X1 = [coord[0] + web.h/2000 / math.sin(theta) , coord[0]+web.h/2000 / math.sin(theta) + k*length]
                Y1 = [coord[1], coord[1] + k*(math.tan(web.angle)*length)]
                X2 = [coord[0] - web.h/2000 / math.sin(theta) , coord[0]-web.h/2000 / math.sin(theta) + k*length]
                Y2 = [coord[1], coord[1] + k*(math.tan(web.angle)*length)]

            plt.plot(X0, Y0, 'k--')
            plt.plot(X1, Y1, 'k')
            plt.plot(X2, Y2, 'k')
            
        # NODES
        for key in self.nodes.keys():
            node = self.nodes[key]
            x, y = node.coord
            plt.scatter(x, y, s=100, c='pink')
            plt.text(x, y, str(key), horizontalalignment='center',
                     verticalalignment='center', fontsize=15)
            if x < X:
                x -= length
            elif x > X:
                x += length/1.5
            if y < Y:
                y -= length
            elif y > Y:
                y += length

            #plt.text(x, y+0.08, f'Fx {node.Fx:.0f} kN', fontsize=11)
            #plt.text(x, y+0.05, f'Fy {node.Fy:.0f} kN', fontsize=11)
            #plt.text(x, y+0.02, f'Mz {node.Mz:.0f} kNm', fontsize=11)
            
        # ECCENTRICITY ELEMENTS
        for elem in self.elements.values():
            n1, n2 = elem.nodes
            x1, y1 = n1.coord
            x2, y2 = n2.coord
            plt.plot([x1,x2], [y1, y2], c='b')
            
        plt.axis('equal')
        plt.title("Joint " + str(self.jid))
        if show:
            plt.show()



class TrussColumn(Truss2D):
    
    def __init__(self, coordinate, H1, H2, L1, L2, n, flip=False):
        super().__init__()
        self.flip = flip
        self.coordinate = coordinate
        self.L1 = L1
        self.L2 = L2
        self.H1 = H1
        self.H2 = H2
        self.n = n
        
        self.generate_chords()
        self.generate_webs()
        
    def add(self, this):
        # TOP CHORD
        if isinstance(this, TopChord):
            self.top_chords.append(this)
            this.mem_id = len(self.members)
            self.members["T" + str(this.mem_id)] = this
            this.calc_nodal_coordinates(self.num_elements)
            
        # BOTTOM CHORD   
        elif isinstance(this, BottomChord):
            self.bottom_chords.append(this)
            this.mem_id = len(self.members)
            self.members["B" + str(this.mem_id)] = this
            this.calc_nodal_coordinates(self.num_elements)
            
        # WEB
        elif isinstance(this, TrussWeb):
            if this.coord_type == "local":
                # Calculate web's global coordinates
                bottom_chord = self.bottom_chords[0]
                top_chord = self.top_chords[0]
                c01 = self.bottom_chords[0].local(this.bot_loc)
                c02 = self.top_chords[0].local(this.top_loc)
                # Change web's coordinates to global coordinates
                this.bot_coord = c01
                this.top_coord =  c02
            else:
                # Sort by y-coordinate
                c01, c02 = sorted(this.coordinates, key=lambda x: x[1])
                # Find corresponding bottom chord
                for chord in self.bottom_chords:
                    if chord.point_intersection(c01):
                        bottom_chord = chord
                        this.bot_loc = bottom_chord.global_coord(c01)
                        break 
                else:
                #if UnboundLocalError and not bottom_chord:
                    raise ValueError("Web's coordinates doesn't match bottom chord's coordinates")                     
                # Find corresponding top chord
                for chord in self.top_chords:   
                    if chord.global_coord(c02) <= 1 and c01[0] <= chord.coordinates[1][0]:
                        top_chord = chord
                        this.top_loc = top_chord.global_coord(c02)
                        break   
                else:
                #if UnboundLocalError and not top_chord:
                    raise ValueError("Web's coordinates doesn't match top chord's coordinates")
            # Check if there's already a joint at given location
            # if not, create a new joint
            j1, j2 = None, None 
            R = 0.03 
            for joint in bottom_chord.joints:    
                if np.sqrt((joint.loc - this.bot_loc)**2)  <= R:
                    j1 = joint
                    this.j1 = j1
                    break
            else:
                j1 = TrussJoint(bottom_chord, this.bot_loc)
                bottom_chord.joints.append(j1)
                this.j1 = j1
                    
                    
            for joint in top_chord.joints:
                if np.sqrt((joint.loc - this.top_loc)**2) <= R:
                    j2 = joint
                    this.j2 = j2
                    break
            else:
                j2 = TrussJoint(top_chord, this.top_loc)
                top_chord.joints.append(j2)
                this.j2 = j2
                
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
     
            self.webs[this.mem_id] = this
            self.members["W" + str(this.mem_id)] = this

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
        elif isinstance(this, Truss2D):
            self.truss = this


        else:
            print(type(this), " is not supported.")
            raise TypeError
        
        
    def generate_chords(self):
        
        x0, y0 = self.coordinate
        
        top_chord = TopChord([[x0, y0], [x0, y0 + self.H1]])
        bottom_chord = BottomChord([[x0 + self.L1, y0], [x0 + self.L1 + self.L2, y0 + self.H2]])
        
        self.add(top_chord)
        self.add(bottom_chord)
        
    def generate_webs(self):
        
        c1 = 0.0
        c2 = 0.0
        
        #self.add(TrussWeb(1,1))
        
        for i in range(1, self.n +1):
            if i%2 == 0:
                c1 = round(i/self.n, 4)
                if c1 > 1:
                    c1 = 1
            elif i!=0:
                c2 = round(i/self.n, 4)
                if c2 > 1:
                    c2 = 1
            
            if self.flip:
                self.add(TrussWeb(c1, c2))
            else:
                self.add(TrussWeb(c2, c1))
        
        


if __name__ == '__main__':

    from src.frame2d.frame2d import *
    frame = Frame2D(simple=[1,1, 5000, 5000], supports='fixed', beams=False)

    simple = dict(
        H0=5000,
        H1=1500,
        H2=2000,
        L1=12500,
        H3=1500,
        dx=1000,
        n=16
    )

    truss = Truss2D(simple=simple)
    for tc in truss.top_chords:
        truss.add(LineLoad(tc, [-30, -30], 'y'))

    truss.add(XYHingedSupport([0,simple['H0'] + simple['H1']]))
    truss.add(YHingedSupport([simple['L1'] * 2, simple['H0'] + simple['H3']]))
    frame.add(truss)
    truss.generate()
    # truss.calculate()
    truss.f.draw()
    for j in truss.joints.values():
        j.e = -100

    truss.f.draw()



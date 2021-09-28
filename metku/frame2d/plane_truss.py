# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 15:17:02 2021

Class for plane (steel) trusses

@author: kmela
"""
try:
    from metku.frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, \
        Support, XYHingedSupport, YHingedSupport
    from metku.framefem.elements import EBBeam, EBSemiRigidBeam
    # from fem.elements.eb_semi_rigid_beam import EBSemiRigidBeam
    from metku.sections.steel import HEA

    from metku.framefem import FrameFEM, BeamSection
    from metku.structures.steel.steel_member import SteelMember
    from metku.eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint, \
        RHSYJoint
except:
    from frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, \
        Support, XYHingedSupport, YHingedSupport
    from framefem.elements import EBBeam, EBSemiRigidBeam
    # from fem.elements.eb_semi_rigid_beam import EBSemiRigidBeam
    from sections.steel import HEA

    from framefem import FrameFEM, BeamSection
    from structures.steel.steel_member import SteelMember
    from eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint, RHSYJoint

from loadIDs import LoadIDs

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Arc

PREC = 3
L_sys = 0.9


class PlaneTruss(Frame2D):

    def __init__(self, origin=(0, 0), span=10000, left_height=1000,
                 mid_height=1500, right_height=1000, bot_chord_dx = 500,
                 num_elements=2,
                 topology='K',
                 ndiv=4,
                 fem_model=None, q=None):
        """ Constructor 
            Parameters:
            ------------
            :param origin: location of the origin (tuple)
            :param num_elements: number of elements per member
            :param fem_model: FrameFEM object
            :param q: load?
    
            :type origin: tuple
            :type num_elements: int
            :type fem_model: FrameFEM (optional)
    
    
            Variables:
            ----------
            :ivar truss:
            :ivar origin:
            :ivar top_chords: list of members constituting the top chord
            :ivar bottom_chords: list of bottom chord members
            :ivar braces: dict of brace members
            :ivar joints: dict of truss joints
            
            :ivar _H0: y-coordinate of the bottom chord
            :ivar _H1: height at
            :ivar _H2: height at
            :ivar _H3: height at
            :ivar L1: span
            :ivar L2: span
            
        
        """
        super().__init__(num_elements=num_elements, fem_model=fem_model)
        
        self.origin = origin
        self.top_chords = []
        self.bottom_chords = []
        self.braces = {}
        self.joints = {}
               
        self.span = span
        
        # Top and bottom chord nodes
        self.top_coords = []
        self.bottom_coords = []
        
        if mid_height < left_height:
            raise ValueError('Truss height at mid-span cannot be smaller than height at left support.')
            
        if mid_height < right_height:
            raise ValueError('Truss height at mid-span cannot be smaller than height at right support.')
        
        self._H0 = left_height
        self._H1 = mid_height
        self._H2 = right_height
        
        self.n = 1
        self._dx = bot_chord_dx
        
        self.top = topology
        self.ndiv = ndiv
        
        self.generate_topology(topology,ndiv)
        self.generate_joints()
        #self.generate_chords()
        #self.generate_webs()

        #if q is not None:
        #    for tc in self.top_chords:
        #        self.add(LineLoad(tc, [q, q], 'y'))

    @property
    def H0(self):
        return self._H0

    @property
    def H1(self):
        return self._H1

    @H1.setter
    def H1(self, val):
        """
        Sets top chord's height of the leftmost eave
        :param val: distance between top and bottom chord
        """
        
        if len(self.top_chords) > 0:
            tc = self.top_chords[0]
            locs = tc.locs
            tc.n1.coord = np.asarray(self.origin) + np.array([0, self.H0 + val])
            for loc, node in zip(locs, tc.nodes.values()):
                node.coord = tc.to_global(loc)
                for j in tc.joints:
                    j.calc_nodal_coordinates()

        if self.H1 == self.H3:
            self.H3 = val
        self._H1 = val

    @property
    def H2(self):
        return self._H2

    @H2.setter
    def H2(self, val):
        """
        Sets top chord's height on ridge
        :param val: distance between top and bottom chord
        """
        if len(self.top_chords) > 0:
            tc = self.top_chords[0]
            tc1 = self.top_chords[1]
            locs = tc.locs
            locs1 = tc1.locs
            # Set the ridge's height
            tc.n2.coord = np.asarray(self.origin) + np.array(
                [self.L1, self.H0 + val])
    
            # Left side
            for loc, node in zip(locs, tc.nodes.values()):
                node.coord = tc.to_global(loc)
            for j in tc.joints:
                j.calc_nodal_coordinates()
    
            # Right side
            for loc, node in zip(locs1, tc1.nodes.values()):
                node.coord = tc1.to_global(loc)
            for j in tc1.joints:
                j.calc_nodal_coordinates()

        self._H2 = val

    @property
    def H3(self):
        return self._H3

    @H3.setter
    def H3(self, val):
        """
        Sets top chord's height of the rightmost eave
        :param val: distance between top and bottom chord
        """
        if len(self.top_chords) > 0:
            tc = self.top_chords[1]
            locs = tc.locs
            tc.n2.coord = np.asarray(self.origin) + np.array(
                [self.L1 + self.L2, self.H0 + val])
            for loc, node in zip(locs, tc.nodes.values()):
                node.coord = tc.to_global(loc)
            for j in tc.joints:
                j.calc_nodal_coordinates()
        self._H3 = val

    @property
    def dx(self):
        return self._dx

    @dx.setter
    def dx(self, val):

        if len(self.bottom_chords) > 0:
            bc = self.bottom_chords[0]
            locs = bc.locs
            bc.n1.coord = np.asarray(self.origin) + np.array([val, self.H0])
            bc.n2.coord = np.asarray(self.origin) + np.array(
                [self.L1 + self.L2 - val, self.H0])
            for loc, node in zip(locs, bc.nodes.values()):
                node.coord = bc.to_global(loc)
            for j in bc.joints:
                j.calc_nodal_coordinates()
        self._dx = val

    @property
    def top_chord_angle(self):
        """ Angle of the top chord with respect to horizontal axis """
        
        a = (self.H1-self.H0)/(0.5*self.span)
        
        return np.arctan(a)
    
    def top_chord_direction_vector(self,half='left'):
        """ Unit vector in the direction of the top chord """
        
        a = self.top_chord_angle
        
        if half == 'left':
            v = np.array([np.cos(a),np.sin(a)])
        else:
            v = np.array([np.cos(a),-np.sin(a)])
        
        return v

    def generate_topology(self,topology='K',ndiv=4):
        """ Generates truss members for a given topology 
        
            :param: topology .. 'K', 'N', or 'KT'
            :param: ndiv .. number of divisions on one half of the truss
            
        """
        
        if topology == 'K':
            # Generate nodes first
            
            # Top chord nodes:
            S = np.linspace(0,1,ndiv+1)
                        
            X0 = np.array([self.origin[0],self.origin[0]+self.H0])
            X1 = np.array([self.origin[0]+0.5*self.span,self.origin[0]+self.H1])            
            X2 = np.array([self.origin[0]+self.span,self.origin[0]+self.H0])
            
            Xtop = [X0+s*(X1-X0) for s in S]
            Xtop += [X1+s*(X2-X1) for s in S[1:]]
            
            # Bottom chord nodes:
            S = np.linspace(0,1,ndiv)
            dx = 0.5*(Xtop[1][0]-Xtop[0][0])            
            X0 = np.array([self.origin[0]+dx,self.origin[0]])
            X1 = np.array([self.origin[0]+0.5*self.span-dx,self.origin[0]])  
            X2 = np.array([self.origin[0]+self.span-dx,self.origin[0]])
            
            Xbot = [X0+s*(X1-X0) for s in S]
            X1 = np.array([self.origin[0]+0.5*self.span+dx,self.origin[0]])  
            Xbot += [X1+s*(X2-X1) for s in S]
        
                
        
            # Generate top chord members
            for i in range(len(Xtop)-1):                
                self.add(TopChord([Xtop[i],Xtop[i+1]]))
        
            # Generate bottom chord members
            for i in range(len(Xbot)-1):
                self.add(BottomChord([Xbot[i],Xbot[i+1]]))
        
            # Generate braces
            
            # First half
            for i in range(ndiv):
                self.add(TrussBrace([Xtop[i],Xbot[i]]))
                self.add(TrussBrace([Xbot[i],Xtop[i+1]]))
            
            # Second half            
            for i in range(ndiv):                
                self.add(TrussBrace([Xtop[ndiv+i],Xbot[ndiv+i]]))
                self.add(TrussBrace([Xbot[ndiv+i],Xtop[ndiv+i+1]]))
            
            
        """    
        fig, ax = plt.subplots(1)
        
        ax.plot([X[0] for X in Xtop], [X[1] for X in Xtop],'ob')
        ax.plot([X[0] for X in Xbot], [X[1] for X in Xbot],'ob')
        ax.set_aspect('equal')
        """  
        self.top_coords = Xtop
        self.bottom_coords = Xbot    
            
    def generate_joints(self):
        """ Generates joints for the truss """

        # Top chord joints        
        for node in self.top_coords:            
            joint_braces = []
            for brace in self.braces.values():            
                if list(node) in brace.coordinates:
                    joint_braces.append(brace)
            
            joint_chord = []
            for chord in self.top_chords:                                
                if all(np.isin(node,chord.coordinates)):                    
                    joint_chord.append(chord)
            
            self.add(TrussJoint(node,joint_chord,joint_braces))
        
        # Bottom chord joints        
        for node in self.bottom_coords:            
            joint_braces = []
            for brace in self.braces.values():            
                if list(node) in brace.coordinates:
                    joint_braces.append(brace)
            
            joint_chord = []
            for chord in self.bottom_chords:                                
                if all(np.isin(node,chord.coordinates)):                    
                    joint_chord.append(chord)
            
            self.add(TrussJoint(node,joint_chord,joint_braces))

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
            self.members["TC" + str(this.mem_id)] = this
            """ Add nodal coordinates for the member """
            this.calc_nodal_coordinates(self.num_elements)

        # BOTTOM CHORD   
        elif isinstance(this, BottomChord):
            self.bottom_chords.append(this)
            this.mem_id = len(self.members)
            self.members["BC" + str(this.mem_id)] = this
            """ Add nodal coordinates for the member """
            this.calc_nodal_coordinates(self.num_elements)

        # BRACE
        elif isinstance(this, TrussBrace):            
            
            # Check if there's already a joint at given location
            # if not, create a new joint
            
            """
            j1, j2 = None, None
            R = 0.03
            for joint in bottom_chord.joints:
                if np.sqrt((joint.loc - this.bot_loc) ** 2) <= R:
                    j1 = joint
                    this.j1 = j1
                    break
            else:
                j1 = TrussJoint(bottom_chord, this.bot_loc)
                bottom_chord.joints.append(j1)
                this.j1 = j1

            for joint in top_chord.joints:
                if np.sqrt((joint.loc - this.top_loc) ** 2) <= R:
                    j2 = joint
                    this.j2 = j2
                    break
            else:
                j2 = TrussJoint(top_chord, this.top_loc)
                top_chord.joints.append(j2)
                this.j2 = j2
            
            """
            
            # Assign member id
            this.mem_id = len(self.braces)
            # Assign joint id's   
            """
            if j1.jid == None:
                j1.jid = len(self.joints)
                self.joints[j1.jid] = j1
            if j2.jid == None:
                j2.jid = len(self.joints)
                self.joints[j2.jid] = j2            

            j1.add_web(this, self.num_elements)
            j2.add_web(this, self.num_elements)
            """
            self.braces[this.mem_id] = this
            self.members["B" + str(this.mem_id)] = this

        elif isinstance(this,TrussJoint):
            
            this.jid = len(self.joints)
            self.joints[this.jid] = this

        else:
            super().add(this)

    
    def calculate(self, load_id=LoadIDs.ULS):
        super().calculate(load_id)
        for joint in self.joints.values():
            joint.update_forces()
    
    def generate0(self):
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

        for w in self.braces.values():
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
            
        # Add point loads (if any)
        for pointLoad in self.point_loads.values():
            pointLoad.add_load(self.f)
            # Creates new loadcase if one with same load_id does not exist
            lcase_ids = [lc.load for lc in self.f.loadcases.values()]
            if pointLoad.load_id not in lcase_ids:
                self.f.add_loadcase(supp_id=1,
                                    load_id=pointLoad.load_id)

        """
        if self_weight:
            self.add_self_weight()
        """
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

        """
        for pointLoad in self.point_loads.values():
            pointLoad.add_load(self.f)

        for lineLoad in self.line_loads.values():
            member = lineLoad.member
            lineLoad.element_ids = member.lineload_elements(
                lineLoad.coordinates)
            lineLoad.add_load(self.f)
        """

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
                 Sj2=np.inf, **kwargs):

        self.joints = []
        # FrameMember -objects that chord is connected to
        self.columns = {}
        # List of chords' members between joints
        self.steel_members = []
        super().__init__(coordinates,
                         mem_id,
                         profile,
                         material,
                         num_elements,
                         Sj1,
                         Sj2, **kwargs)

    @property
    def profile(self):
        return super().profile

    @profile.setter
    def profile(self, val):
        super(TrussMember, self.__class__).profile.fset(self, val)
        for joint in self.joints:
            joint.calc_nodal_coordinates()

    def calc_nodal_coordinates(self, num_elements=0):
        """
            Calculates node locations along member
            Coordinates are global coordinates
            Adds locations to a list where they can be accessed later
            :param num_elements: int, number of elements to be created
        """
        self.nodal_coordinates.clear()
        self.loc.clear()
        self.nodal_coordinates.extend(self.added_coordinates)
        if num_elements > 0:
            self.num_elements = num_elements
        c00, c01 = self.coordinates
        
        # Create calculation nodes for joints (if joints are present)
        if len(self.joints) > 0:
            joint_locs = [j.loc for j in self.joints]
            if 0 not in joint_locs:
                joint_locs.append(0)
            if 1 not in joint_locs:
                joint_locs.append(1)
            joint_locs.sort()
            for i in range(len(joint_locs) - 1):
                j_loc_0 = joint_locs[i]
                j_loc_1 = joint_locs[i + 1]

                d_loc = j_loc_1 - j_loc_0
                d_loc /= self.num_elements
                for j in range(self.num_elements + 1):
                    loc = j_loc_0 + d_loc * j
                    loc = round(loc, 5)
                    coord = list(self.to_global(loc))
                    loc = min(loc, 1)
                    if loc not in self.loc:
                        self.loc.append(loc)
                    if loc <= 1 and coord not in self.nodal_coordinates:
                        self.nodal_coordinates.append(coord)

            last_joint_coord = list(self.to_global(joint_locs[-1]))
            if last_joint_coord not in self.nodal_coordinates:
                self.nodal_coordinates.append(last_joint_coord)
            if list(self.to_global(1)) not in self.nodal_coordinates:
                self.nodal_coordinates.append(list(self.to_global(1)))
            c0, c1 = self.coordinates
            self.nodal_coordinates.sort(
                key=lambda c: np.linalg.norm(np.asarray(c0) - np.asarray(c)))
            joint_locs.sort()
            coords = [self.to_global(loc) for loc in self.loc if
                      loc not in joint_locs]

        else:
            # If no joints are present, use the method of the FrameMember class
            #print("use FrameMember calc method")
            super().calc_nodal_coordinates(num_elements)

        if self.is_generated:
            nodes = [n for n in self.nodes.values()]
            nodes = sorted(nodes, key=lambda n: np.linalg.norm(
                np.asarray(c00) - np.asarray(n.coord)))
            self.loc.sort()

            for loc, node in zip(self.locs, nodes):
                node.coord = self.to_global(loc)

    def check_cross_section(self):

        for smem in self.steel_members:
            smem.profile = self.cross_section

        if not isinstance(self.TrussBrace): # self.mtype != 'web':  # If member is chord
            # Sort joints using joints' local coordinate
            self.joints.sort(key=lambda x: x.loc)
            if not len(self.steel_members):
                for i in range(len(self.joints) - 1):
                    j0 = self.joints[i]
                    j1 = self.joints[i + 1]
                    x0, y0 = j0.coordinate
                    x1, y1 = j1.coordinate
                    L = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
                    steel_mem = SteelMember(self.cross_section, L,
                                            [L_sys, L_sys])
                    for id in range(j0.nodes['c0'].nid, j1.nodes['c0'].nid):
                        N, V, M = self.nodal_forces[id]
                        steel_mem.add_section(ned=N, myed=M, vyed=V)
                    self.steel_members.append(steel_mem)
            else:
                for i in range(len(self.joints) - 1):
                    j0 = self.joints[i]
                    j1 = self.joints[i + 1]
                    x0, y0 = j0.coordinate
                    x1, y1 = j1.coordinate
                    L = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
                    steel_mem = self.steel_members[i]
                    steel_mem.length = L
                    steel_mem.clear_sections()
                    for id in range(j0.nodes['c0'].nid,
                                    j1.nodes['c0'].nid):
                        N, V, M = self.nodal_forces[id]
                        steel_mem.add_section(ned=N, myed=M, vyed=V)
                # r = []
                # r_sect = steel_member.check_sections()
                # r_buckling = steel_member.check_buckling()
                # r_lt_buckling = steel_member.check_LT_buckling()
                # r.extend(r_sect)
                # r.extend(r_buckling)
                # r.append(r_lt_buckling)
        else:
            super().check_cross_section()

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

        L = np.sqrt((x0 - x) ** 2 + (y0 - y) ** 2)
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

        return k * (x - x0) + y0

    # @property
    # def angle(self):
    #     """
    #     Returns member's angle in radians
    #     """
    #     start_node, end_node = self.coordinates
    #     x0, y0 = start_node
    #     x1, y1 = end_node
    #     if (x1 - x0) == 0:
    #         angle = math.radians(90)
    #     else:
    #         angle = math.atan((y1 - y0) / (x1 - x0))
    #     return angle


    def plot(self, print_text=True, color=True, axes=None):

        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes        

        X = self.coordinates
        if color:
            if self.is_strong_enough:
                color = 'green'
            else:
                color = 'red'
        else:
            color = 'k'
        # Plot members
        ax.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color)
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
            # rot = self.angle

        x = (X[0][0] + X[1][0]) / 2
        y = (X[1][1] + X[0][1]) / 2
        horzalign = 'center'
        vertalign = 'center'
        if print_text:
            ax.text(x, y, str(self.mem_id) + ": " + str(self.cross_section),
                     rotation=rot, horizontalalignment=horzalign,
                     rotation_mode='anchor')
        else:
            ax.text(x, y, str(self.mem_id),
                     rotation=rot, horizontalalignment=horzalign,
                     verticalalignment=vertalign)

    def bmd0(self, scale):
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
                if i == 0 or i == len(self.nodes) - 1 or bending_moment == max(
                        moment_values) or \
                        bending_moment == min(moment_values):
                    if bending_moment > 0:
                        vert = 'top'
                    else:
                        vert = 'bottom'
                    plt.text(x, y, f'{bending_moment:.2f} kNm',
                             verticalalignment=vert)

            elif self.mtype == "web":
                x0 = self.nodal_coordinates[i][0]
                y0 = self.nodal_coordinates[i][1]
                bending_moment = self.nodal_forces[node][2]
                x1 = bending_moment / (1000 / scale)
                x = x0 + x1
                y = y0
                X.append(x)
                Y.append(y)
                if i == 0 or i == len(self.nodes) - 1 or bending_moment == max(
                        moment_values) or \
                        bending_moment == min(moment_values):
                    if bending_moment > 0:
                        horz = 'left'
                    else:
                        horz = 'right'
                    plt.text(x, y, f'{bending_moment:.2f} kNm',
                             horizontalalignment=horz)
        X.append(self.nodal_coordinates[i][0])
        Y.append(self.nodal_coordinates[i][1])
        plt.plot(X, Y, color='gray')


class Chord(FrameMember):

    def __init__(self, coordinates, **kwargs):
        super().__init__(coordinates, **kwargs)

        self.members = []
        self.joints = []

    def generate0(self, fem_model):

        locs = [0, 1]
        for joint in self.joints:
            if joint.loc not in locs:
                locs.append(joint.loc)

        locs.sort()
        for i in range(len(locs) - 1):
            c0 = self.to_global(locs[i])
            c1 = self.to_global(locs[i + 1])
            mem = FrameMember([c0, c1])
            self.members.append(mem)
            mem.generate(fem_model)

    def check_cross_section(self):
        for mem in self.members:
            mem.check_cross_section()
    
    def calc_nodal_coordinates(self, num_elements=2):
        """
            Calculates node locations along a chord member
            Coordinates are global coordinates
            Adds locations to a list where they can be accessed later
            :param num_elements: int, number of elements to be created
        """
        self.nodal_coordinates.clear()
        self.loc.clear()
        print(self.added_coordinates)
        self.nodal_coordinates.extend(self.added_coordinates)
        if num_elements:
            self.num_elements = num_elements
        c00, c01 = self.coordinates
        joint_locs = [j.loc for j in self.joints]
        if len(joint_locs):
            if 0 not in joint_locs:
                joint_locs.append(0)
            if 1 not in joint_locs:
                joint_locs.append(1)
            joint_locs.sort()
            for i in range(len(joint_locs) - 1):
                j_loc_0 = joint_locs[i]
                j_loc_1 = joint_locs[i + 1]

                d_loc = j_loc_1 - j_loc_0
                d_loc /= self.num_elements
                for j in range(self.num_elements + 1):
                    loc = j_loc_0 + d_loc * j
                    loc = round(loc, 5)
                    coord = list(self.to_global(loc))
                    loc = min(loc, 1)
                    if loc not in self.loc:
                        self.loc.append(loc)
                    if loc <= 1 and coord not in self.nodal_coordinates:
                        self.nodal_coordinates.append(coord)

            last_joint_coord = list(self.to_global(joint_locs[-1]))
            if last_joint_coord not in self.nodal_coordinates:
                self.nodal_coordinates.append(last_joint_coord)
            if list(self.to_global(1)) not in self.nodal_coordinates:
                self.nodal_coordinates.append(list(self.to_global(1)))
            c0, c1 = self.coordinates
            self.nodal_coordinates.sort(
                key=lambda c: np.linalg.norm(np.asarray(c0) - np.asarray(c)))
            joint_locs.sort()
            coords = [self.to_global(loc) for loc in self.loc if
                      loc not in joint_locs]

        else:
            super().calc_nodal_coordinates(num_elements)

        if self.is_generated:
            nodes = [n for n in self.nodes.values()]
            nodes = sorted(nodes, key=lambda n: np.linalg.norm(
                np.asarray(c00) - np.asarray(n.coord)))
            self.loc.sort()

            for loc, node in zip(self.locs, nodes):
                node.coord = self.to_global(loc)


class TopChord(TrussMember):
    def __init__(self,
                 coordinates,
                 **kwargs):
        super().__init__(coordinates,
                         **kwargs)

        self.mtype = 'top_chord'
        self.steel_members = []


class BottomChord(TrussMember):

    def __init__(self,
                 coordinates,
                 **kwargs):
        super().__init__(coordinates,
                         **kwargs)
        self.mtype = 'bottom_chord'
        self.steel_members = []

    @property
    def perpendicular(self):
        return -1 * super().perpendicular


class TrussBrace(TrussMember):
    """ Class for truss web members, or braces """
    def __init__(self,
                 coordinates,
                 mem_id="",
                 profile="SHS 50x2",
                 material="S355",
                 num_elements=1,
                 Sj1=0,
                 Sj2=0):
        """ Constructor
        
            :param: bot_loc .. location of the member end at the bottom chord
            :param: top_loc .. location of the member end at the top chord
            :param: coord_type .. 'local' (default) or 'global' coordinates
        """

        self.j1 = None
        self.j2 = None
        
        self.top_joint = None
        self.bottom_joint = None

        super().__init__(coordinates, profile=profile,material=material, num_elements=num_elements)
   
        self.top_coord = self.coordinates[1]
        self.bot_coord = self.coordinates[0]
        self.Sj1 = Sj1
        self.Sj2 = Sj2
        self.mtype = 'brace'
        self.steel_members = [self.member] 


class TrussJoint:
    """ Class for welded tubular joints """

    def __init__(self, node, chord, brace, reverse=False):
    
        for b in brace:
            b.joints.append(self)
            
            if isinstance(chord[0],TopChord):
                b.j1 = self
            else:
                b.j2 = self
            
        for c in chord:
            c.joints.append(self)

        self.chord = chord
        self.helper_chord = chord.copy()
        self.chord_elements = {}
        self.jid = None
 
        self.cnode = node
        self.chord_nodes = []
        self.joint_type = None
        self.nodal_coordinates = {}
        self.nodes = {}
        self.elements = {}
        self.element_ids = []
        self.braces = brace
        self.is_generated = False
        # Index for eccentricity elements' material and section
        self.idx = None
        self.rhs_joint = None
        self.r = []
        self.__e = 0
        self.g_set = False
        self.reverse = reverse

    def add_node_coord(self, coord):
        coord = [round(c, PREC) for c in coord]
        if coord not in self.nodal_coordinates:
            self.nodal_coordinates.append(coord)

    def add_brace(self, brace):
        self.braces[brace.mem_id] = brace

        # if self.is_generated:
        self.calc_nodal_coordinates()

    @property
    def loc(self):
        return round(self.__loc, PREC)

    @loc.setter
    def loc(self, val):
        # if self.chord.mtype == 'bottom_chord' or 0.01 <= self.loc <= 0.99:
        min_val = self.g1 / self.chord.length
        max_val = 1 - self.g1 / self.chord.length

        min_val = 0
        max_val = 1

        val = min(max(min_val, val), max_val)
        val = round(val, PREC)

        # locs = [j.loc for j in self.chord.joints]
        # idx = locs.index(self.loc)
        # if idx == 0:
        #     min_val = 1e-2
        # else:
        #     min_val = locs[idx - 1] + 0.05
        # if idx == len(locs) -1:
        #     max_val = 1-1e-2
        # else:
        #     max_val = locs[idx + 1] - 0.05
        #
        # val = max(min(val, max_val), min_val)

        if self.reverse:
            val = 1 - val

        if val != self.loc:
            self.__loc = val
            # if self.chord.mtype == 'top_chord':
            #     for web in self.webs.values():
            #         web.top_loc = self.loc
            # else:
            #     for web in self.webs.values():
            #         web.bot_loc = self.loc

            self.calc_nodal_coordinates()

    @property
    def coordinate(self):
        if self.chord.mtype == 'bottom_chord':
            return [round(c, PREC) for c in
                    self.helper_chord.to_global(self.loc)]

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
    def theta1(self):
        if self.chord.mtype == 'top_chord' and not self.chord.reverse:
            webs = sorted(self.webs.values(), key=lambda w: w.bot_loc)
        else:
            webs = sorted(self.webs.values(), key=lambda w: w.top_loc)
        w1 = webs[0]
        theta1 = min(abs(w1.angle - self.chord.angle),
                     abs(np.pi - w1.angle - self.chord.angle))

        return theta1

    @property
    def theta2(self):
        if len(self.webs) >= 2:
            if self.chord.mtype == 'top_chord' and not self.chord.reverse:
                webs = sorted(self.webs.values(), key=lambda w: w.bot_loc)
            else:
                webs = sorted(self.webs.values(), key=lambda w: w.top_loc)
            w2 = webs[1]
            theta2 = min(abs(w2.angle - self.chord.angle),
                         abs(np.pi - w2.angle - self.chord.angle))
            return theta2
        return 0

    @property
    def e(self):
        """ Eccentrictity of the connection """

        # if self.joint_type == 'K':
        #     w1, w2 = self.webs.values()
        #     ec = w1.line_intersection(w2)
        #     # TODO: distance from ec to chord's axis

        if self.g_set:
            if self.joint_type == 'K':
                w1, w2 = list(self.webs.values())
                h1 = w1.h
                h2 = w2.h
                h0 = self.chord.h
                theta1 = min(abs(w1.angle - self.chord.angle),
                             abs(np.pi - w1.angle - self.chord.angle))
                theta2 = min(abs(w2.angle - self.chord.angle),
                             abs(np.pi - w2.angle - self.chord.angle))
                e = h1 / (2 * np.sin(theta1)) + h2 / (
                        2 * np.sin(theta2)) + self.__g1
                e *= (np.sin(theta1) * np.sin(theta2)) / np.sin(
                    theta1 + theta2)
                e -= h0 / 2

                return e

            return 0

        else:
            return self.__e

    @e.setter
    def e(self, val):
        self.g_set = False
        delta = val - self.e
        # self.calc_nodal_coordinates()
        if self.joint_type == 'K':
            w1, w2 = self.webs.values()
            h1 = w1.h
            h2 = w2.h
            h0 = self.chord.h

            X1 = np.cos(self.theta1) * w1.length
            X2 = np.cos(self.theta2) * w2.length

            theta1 = np.arctan(self.theta1 + delta / X1)
            theta2 = np.arctan(self.theta2 + delta / X2)

            # theta1 = min(abs(w1.angle - self.chord.angle),
            #              abs(np.pi - w1.angle - self.chord.angle))
            # theta2 = min(abs(w2.angle - self.chord.angle),
            #              abs(np.pi - w2.angle - self.chord.angle))

            g = val + h0 / 2
            g *= np.sin(theta1 + theta2) / (np.sin(theta1) * np.sin(theta2))
            g -= h1 / (2 * np.sin(theta1))
            g -= h2 / (2 * np.sin(theta2))
            self.__g1 = g

            if self.chord.mtype == 'top_chord':
                w1.j1.calc_nodal_coordinates()
                w2.j1.calc_nodal_coordinates()
            else:
                w1.j2.calc_nodal_coordinates()
                w2.j2.calc_nodal_coordinates()

        self.__e = val

    @property
    def g1(self):

        if self.joint_type == 'K':
            w1, w2 = self.webs.values()
            h1 = w1.h
            h2 = w2.h
            h0 = self.chord.h

            theta1 = min(abs(w1.angle - self.chord.angle),
                         abs(np.pi - w1.angle - self.chord.angle))
            theta2 = min(abs(w2.angle - self.chord.angle),
                         abs(np.pi - w2.angle - self.chord.angle))

            g = self.e + h0 / 2
            g *= np.sin(theta1 + theta2) / (np.sin(theta1) * np.sin(theta2))
            g -= h1 / (2 * np.sin(theta1))
            g -= h2 / (2 * np.sin(theta2))
            return g

        else:
            return self.__g1

    @g1.setter
    def g1(self, val):
        self.g_set = True
        self.__g1 = val
        # if self.joint_type == 'K':
        #     w1, w2 = list(self.webs.values())
        #     h1 = w1.h
        #     h2 = w2.h
        #     h0 = self.chord.h
        #
        #     e = h1 / (2 * np.sin(self.theta1)) + h2 / (2 * np.sin(self.theta2)) + self.__g1
        #     e *= (np.sin(self.theta1) * np.sin(self.theta2)) / np.sin(self.theta1 + self.theta2)
        #     e -= h0 / 2
        #
        #     self.e = e

        self.calc_nodal_coordinates()

    @property
    def g2(self):
        return self.__g2

    @g2.setter
    def g2(self, val):
        self.__g2 = val
        self.calc_nodal_coordinates()

    @property
    def N0(self):
        """
        Returns maximum normal force on chord on gap

        :return: max N
        """
        if self.joint_type == 'Y':
            N0 = self.chord.nodal_forces[self.nodes['c0'].nid][0]


        elif self.joint_type == 'K':
            N01 = self.chord.nodal_forces[self.nodes['c01'].nid][0]
            N02 = self.chord.nodal_forces[self.nodes['c02'].nid][0]
            N0 = max(abs(N01), abs(N02))
            # N0 *= np.sign(N01)
        else:
            # TODO: KT
            N0 = 0
        return N0
        # N = [abs(elem.shear_force[0]) for elem in self.elements.values()]
        # return max(N)

    @property
    def V0(self):
        """
        Returns maximum shear force on chord on gap

        :return: max N
        """
        V0 = self.chord.nodal_forces[self.nodes['c0'].nid][1]
        return abs(V0)

        # V = [abs(elem.axial_force[0]) for elem in self.elements.values()]
        # return max(V)

    @property
    def M0(self):
        if self.joint_type == 'Y':
            M0 = self.chord.nodal_forces[self.nodes['c0'].nid][2]

        elif self.joint_type == 'K':
            M01 = self.chord.nodal_forces[self.nodes['c01'].nid][2]
            M02 = self.chord.nodal_forces[self.nodes['c02'].nid][2]
            M0 = max(abs(M01), abs(M02))
        else:
            # TODO: KT
            pass
        return abs(M0)

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
            r5 = abs(NEd0 / r4[1])  # chord
            r4 = abs(wNEd / r4[0])

            # self.r = [r1, r2, r3, r4]
            maxR_for_webs = np.max((r1, r2, r3, r4), axis=0)
            maxR = np.append(maxR_for_webs, r5)
            self.r = maxR
            # return maxR

        elif self.joint_type == 'Y':
            NEd1 = [w.ned for w in self.webs.values()][0]

            r1 = self.chord_face_failure()
            r2 = self.chord_web_buckling()
            r3 = self.brace_failure()

            r1 = abs(NEd0 / r1)
            r2 = abs(NEd0 / r2)
            r3 = abs(NEd1 / r3)
            # self.r = [r1, r2, r3]

            maxR = [r3, max(r1, r2)]
            self.r = maxR
            # return maxR
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
            return results / 1000  # N to kN

    def chord_shear(self):
        if self.rhs_joint and self.joint_type == 'K':
            V0 = [elem.shear_force for elem in self.chord.elements.values() if
                  self.cnode in elem.nodes]
            # item for sublist in list for item in sublist
            V0 = [abs(V) for forces in V0 for V in forces]
            V0 = max(V0)
            self.rhs_joint.V0 = V0 * 1e3  # kN to N
            results = [r / 1000 for r in self.rhs_joint.chord_shear()]
            return results

    def chord_web_buckling(self):
        if self.rhs_joint and self.joint_type == 'Y':
            return self.rhs_joint.chord_web_buckling() / 1000

    def calc_nodal_coordinates(self):
        """
        Calculates nodal coordinates

        Y:
        -----
        TODO

        K :
        -----
        c0 -- joint's center coordinate
        ec -- coordinate of web's axis' intersection
        w1 -- leftmost web
        c01 -- web's end coordinate's projection to chord's axis
        c1 -- web's end coordinate on chord's flange
        w2 -- rightmost web
        c02 -- web's end coordinate's projection to chord's axis
        c2 -- web's end coordinate on chord's flange

        KT:
        -----
        TODO

        :return:
        """
        coords = []

        # Remove coordinates from chord
        for coord in self.nodal_coordinates.values():
            if list(coord) in self.chord.added_coordinates:
                self.chord.added_coordinates.remove(list(coord))

        if list(self.coordinate) in self.chord.added_coordinates:
            self.chord.added_coordinates.remove(list(self.coordinate))

        # Y
        if len(self.webs) == 1:
            if self.chord.mtype == "top_chord":
                if self.loc == 0:
                    self.loc = self.g1 / self.chord.length
                if self.loc == 1:
                    self.loc = 1 - self.g1 / self.chord.length
            self.joint_type = "Y"
            w1 = list(self.webs.values())[0]
            c0 = self.coordinate
            c1 = c0 + self.chord.h / 2 * self.chord.perpendicular
            if self.chord.mtype == 'top_chord':
                w1.top_coord = c1
            else:
                w1.bot_coord = c1
            coords = [c0, c1]

            self.nodal_coordinates['c0'] = c0
            self.nodal_coordinates['c1'] = c1

            theta1 = min(abs(w1.angle - self.chord.angle),
                         abs(np.pi - w1.angle - self.chord.angle))
            # Set design class
            self.rhs_joint = RHSYJoint(self.chord.cross_section,
                                       w1.cross_section,
                                       np.degrees(theta1))

        # K
        elif len(self.webs) == 2:
            # Set joint type
            self.joint_type = 'K'
            # Center coordinate
            c0 = self.coordinate
            # Eccentricity coordinate
            ec = c0 - self.e * self.chord.perpendicular
            # Calculate webs' and chord's axis' intersections
            if self.chord.mtype == 'top_chord' and not self.chord.reverse:
                w1, w2 = sorted(self.webs.values(), key=lambda w: w.bot_loc)
                c01 = self.chord.line_intersection([w1.bot_coord, ec])
                c02 = self.chord.line_intersection([w2.bot_coord, ec])
            else:
                w1, w2 = sorted(self.webs.values(), key=lambda w: w.top_loc)
                c01 = self.chord.line_intersection([w1.top_coord, ec])
                c02 = self.chord.line_intersection([w2.top_coord, ec])
            # Calculate distance between axis' intersection
            # and projection from web to chord's axis
            theta1 = min(abs(w1.angle - self.chord.angle),
                         abs(np.pi - w1.angle - self.chord.angle))
            theta2 = min(abs(w2.angle - self.chord.angle),
                         abs(np.pi - w2.angle - self.chord.angle))

            d1 = abs(self.chord.h / 2 / np.tan(theta1))
            c01 -= self.chord.unit * d1

            d2 = abs(self.chord.h / 2 / np.tan(theta2))
            c02 += self.chord.unit * d2

            if np.linalg.norm(
                    self.chord.to_global(1) - c01) > self.chord.length:
                d01 = np.linalg.norm(self.chord.to_global(0) - c01)
                c01 += self.chord.unit * d01
                c0 += self.chord.unit * d01
                c02 += self.chord.unit * d01

            elif np.linalg.norm(
                    self.chord.to_global(0) - c02) > self.chord.length:
                d02 = np.linalg.norm(self.chord.to_global(1) - c02)
                c01 -= self.chord.unit * d02
                c0 -= self.chord.unit * d02
                c02 -= self.chord.unit * d02

            if self.chord.mtype == 'top_chord':
                c1 = np.asarray(
                    c01) + self.chord.h / 2 * self.chord.perpendicular
                c2 = np.asarray(
                    c02) + self.chord.h / 2 * self.chord.perpendicular
                w1.top_coord = c1
                w2.top_coord = c2
            else:
                c1 = np.asarray(
                    c01) + self.chord.h / 2 * self.chord.perpendicular
                c2 = np.asarray(
                    c02) + self.chord.h / 2 * self.chord.perpendicular
                w1.bot_coord = c1
                w2.bot_coord = c2

            coords = [c0, c01, c02, c1, c2]

            self.nodal_coordinates['c0'] = c0
            self.nodal_coordinates['c01'] = c01
            self.nodal_coordinates['c02'] = c02
            self.nodal_coordinates['c1'] = c1
            self.nodal_coordinates['c2'] = c2

            # Set design class
            self.rhs_joint = RHSKGapJoint(self.chord.cross_section,
                                          [w1.cross_section,
                                           w2.cross_section],
                                          [np.degrees(theta1),
                                           np.degrees(theta2)],
                                          self.g1)

        # KT
        elif len(self.webs) == 3:
            self.joint_type = 'KT'
            c0 = self.chord.to_global(self.loc)
            ec = c0 - self.e * self.chord.perpendicular
            if self.chord.mtype == 'top_chord':
                w1, w2, w3 = sorted(self.webs.values(),
                                    key=lambda w: w.bot_loc)
            else:
                w1, w2, w3 = sorted(self.webs.values(),
                                    key=lambda w: w.top_loc)
            c01 = self.chord.line_intersection(w1)
            theta1 = w1.angle - self.chord.angle
            d1 = abs(self.chord.h / 2 / np.tan(theta1))
            c01 -= self.chord.unit * d1

            c02 = self.chord.line_intersection(w2)
            theta2 = w2.angle - self.chord.angle
            d2 = abs(self.chord.h / 2 / np.tan(theta2))
            c02 += self.chord.unit * d2

            # TODO

        # Move existing nodes
        if self.is_generated:
            for key, coord in self.nodal_coordinates.items():
                self.nodes[key].coord = coord
        # Add nodes to chord
        else:
            for c in coords:
                self.chord.add_node_coord(c)

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
        self.chord_nodes = sorted(self.chord_nodes, key=lambda node: node.x)

    def get_chord_elements(self):
        """
        Gets the joint's elements which are on chord
        """
        if len(self.chord_nodes) > 2:
            for i in range(2):
                nodes = self.chord_nodes[i:i + 2]
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
        sect = BeamSection(38880 * 1e-6, 6299657109 * 1e-12)
        fem_model.add_section(sect)

        # USE CHORD'S PROPERTIES
        # self.idx = self.chord.mem_id

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
            self.elements[index + 1] = EBBeam(self.nodes['c02'],
                                              self.nodes['c2'],
                                              section,
                                              fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index + 1])
            self.element_ids.append(index + 1)

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
        plt.scatter(x, y, s=100, c=color)
        plt.text(x, y, str(self.jid), horizontalalignment=horzalign,
                 verticalalignment=vertalign, fontsize=12)

    def plot_joint(self, length=0.02, show=True):
        """
        Plots joint
        :param length: local coordinate, defines how much members are show
        :param show: boolean, whether to plot or not
        """
        X, Y = self.coordinate
        end_line = None
        col = None
        # CHORD
        x, y = self.chord.to_global(self.loc)
        h = self.chord.h
        C0 = self.chord.to_global(max(self.loc - length, 0))
        C01 = self.chord.to_global(min(self.loc + length, 1))

        d = np.linalg.norm(C01 - C0)

        C1 = C0 + h / 2 * self.chord.perpendicular
        C2 = C0 - h / 2 * self.chord.perpendicular

        C11 = C1 + d * self.chord.unit
        C21 = C2 + d * self.chord.unit

        plt.plot([C0[0], C01[0]], [C0[1], C01[1]], '--', color='lightblue')
        plt.plot([C1[0], C11[0]], [C1[1], C11[1]], 'k')
        plt.plot([C2[0], C21[0]], [C2[1], C21[1]], 'k')

        # WEBS
        for w in self.webs.values():
            d = abs(d)
            h = w.h
            if self.chord.mtype == 'top_chord':
                W0 = w.top_coord
            else:
                W0 = w.bot_coord

            # if np.degrees(w.angle) >= 0:
            #     d += 500
            # else:
            #     d -= 500

            W1 = W0 + h / 2 * w.perpendicular
            W2 = W0 - h / 2 * w.perpendicular
            if self.chord.mtype == 'top_chord':
                W01 = W0 - d * w.unit
                W11 = W1 - d * w.unit
                W21 = W2 - d * w.unit
                W1 = line_intersection([C1, C11], [W1, W11])
                W2 = line_intersection([C1, C11], [W2, W21])
            else:
                W01 = W0 + d * w.unit
                W11 = W1 + d * w.unit
                W21 = W2 + d * w.unit
                W1 = line_intersection([C1, C11], [W1, W11])
                W2 = line_intersection([C1, C11], [W2, W21])

            plt.plot([W0[0], W01[0]], [W0[1], W01[1]], '--', color='lightblue')
            plt.plot([W1[0], W11[0]], [W1[1], W11[1]], 'k')
            plt.plot([W2[0], W21[0]], [W2[1], W21[1]], 'k')

        if self.joint_type == 'K':
            if self.chord.mtype == 'top_chord' and not self.chord.reverse:
                w1, w2 = sorted(self.webs.values(), key=lambda w: w.bot_loc)
            else:
                w1, w2 = sorted(self.webs.values(), key=lambda w: w.top_loc)

            theta0 = np.degrees(self.chord.angle)
            theta1 = np.degrees(self.theta1)
            theta2 = np.degrees(self.theta2)

            if self.chord.mtype == "bottom_chord":
                xy1 = w1.bot_coord  # - w1.cross_section.H / 2 * self.chord.unit
                xy2 = w2.bot_coord  # + w2.cross_section.H / 2 * self.chord.unit
                arc1 = Arc(xy1, 100, 100, theta1=theta0 + 180 - theta1,
                           theta2=theta0 + 180)
                arc2 = Arc(xy2, 100, 100, theta1=theta0,
                           theta2=theta0 + theta2)

                x1, y1 = (xy1 + np.array([-50, -10])) + 50 * w1.unit
                plt.text(x1, y1, f"{theta1:.1f}\N{DEGREE SIGN}")

                x2, y2 = (xy2 + np.array([50, -10])) + 50 * w1.unit
                plt.text(x2, y2, f"{theta2:.1f}\N{DEGREE SIGN}")
            else:
                xy1 = w1.top_coord  # - w1.cross_section.H / 2 * self.chord.unit
                xy2 = w2.top_coord  # + w2.cross_section.H / 2 * self.chord.unit

                if theta0 < 90:
                    arc1 = Arc(xy1, 100, 100, theta1=theta0 + 180,
                               theta2=theta0 + 180 + theta1)
                    arc2 = Arc(xy2, 100, 100, theta1=360 - theta0 - theta2,
                               theta2=theta0)
                else:
                    arc1 = Arc(xy1, 100, 100, theta1=theta0,
                               theta2=theta0 + theta1)
                    arc2 = Arc(xy2, 100, 100, theta1=180 + theta0 - theta2,
                               theta2=180 + theta0)

                x1, y1 = (xy1 + np.array([-50, -10])) - 50 * w1.unit
                plt.text(x1, y1, f"{theta1:.1f}\N{DEGREE SIGN}")

                x2, y2 = (xy2 + np.array([50, -10])) - 50 * w1.unit
                plt.text(x2, y2, f"{theta2:.1f}\N{DEGREE SIGN}")

            ax = plt.subplot(111)
            ax.add_patch(arc1)
            ax.add_patch(arc2)

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
                x += length / 1.5
            if y < Y:
                y -= length
            elif y > Y:
                y += length

            # plt.text(x, y+0.08, f'Fx {node.Fx:.0f} kN', fontsize=11)
            # plt.text(x, y+0.05, f'Fy {node.Fy:.0f} kN', fontsize=11)
            # plt.text(x, y+0.02, f'Mz {node.Mz:.0f} kNm', fontsize=11)

        x, y = self.nodal_coordinates['c0'] - self.e * self.chord.perpendicular
        plt.scatter(x, y, s=100, c='pink')
        plt.text(x, y, "ec", horizontalalignment='center',
                 verticalalignment='center', fontsize=15)

        # ECCENTRICITY ELEMENTS
        for elem in self.elements.values():
            n1, n2 = elem.nodes
            x1, y1 = n1.coord
            x2, y2 = n2.coord
            plt.plot([x1, x2], [y1, y2], c='b')

        plt.axis('equal')
        plt.title("Joint " + str(self.jid))
        if show:
            plt.show()

    def update_forces(self):
        """
        Updates RHS_joint's forces
        """
        self.rhs_joint.M0 = self.M0
        self.rhs_joint.N0 = self.N0
        self.rhs_joint.V0 = self.V0

        if self.joint_type == 'Y':
            self.angle = np.degrees(self.theta1)

        elif self.joint_type == 'K':
            self.rhs_joint.angles = [np.degrees(self.theta1),
                                     np.degrees(self.theta2)]


        elif self.joint_type == 'KT':
            # TODO
            pass

def line_intersection(coordinates1, coordinates2):
    """
    Calculates coordinate where two members intersect
    :param coordinates: array of two arrays of two float values, [[x1,y1],[x2,y2]]
    :return [px, py]: array of two float values, coordinates for intersection
    source: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    """
    start_node, end_node = coordinates1
    x1, y1 = start_node
    x2, y2 = end_node
    x3, y3 = coordinates2[0]
    x4, y4 = coordinates2[1]

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


if __name__ == '__main__':

    """
    from metku.frame2d.frame2d import *

    n_bays = 1
    n_storeys = 1
    L_bay = 10000
    H_storey = 5000

    simple_frame = [n_storeys, n_bays, H_storey, L_bay]

    frame = Frame2D(simple=simple_frame, create_beams=False, supports='fixed')

    simple_truss = dict(
        H0=H_storey,
        H1=1500,
        H2=2000,
        L1=L_bay / 2,
        H3=1500,
        dx=1000,
        n=16
    )
    """
    truss = PlaneTruss(span=10500, left_height=1000, mid_height=1150, right_height=1000,
                       num_elements=4)

    
    for tc in truss.top_chords:
        truss.add(LineLoad(tc, [-30, -30], 'y',load_id=LoadIDs.ULS))
    
    
    truss.add(XYHingedSupport(truss.top_coords[0]))
    truss.add(YHingedSupport(truss.top_coords[-1]))
    
    #truss.plot()
    truss.generate()
    truss.calculate()
    #truss.bmd(scale=100,loads=False)

    # truss.H1 = 1000
    # truss.H2 = 1500
    # truss.H3 = 1300
    # truss.dx = 100
    #
    # truss.f.draw()

    """
    joint = truss.joints[4]
    w1, w2 = joint.webs.values()

    print("Before: ")
    print('e0', joint.e)
    print('g0', joint.g1)
    print('e1', w1.j1.e)
    print('g1', w1.j1.g1)
    print('e2', w2.j1.e)
    print('g2', w2.j1.g1)

    print("theta1", np.degrees(joint.theta1))
    print("theta2", np.degrees(joint.theta2))
    print("-" * 20)
    val = 50
    print("Changed g0 to:", val)
    print("-" * 20)
    joint.g1 = val
    print("After: ")
    print('e0', joint.e)
    print('g0', joint.g1)
    print('e1', w1.j1.e)
    print('g1', w1.j1.g1)
    print('e2', w2.j1.e)
    print('g2', w2.j1.g1)
    print("theta1", np.degrees(joint.theta1))
    print("theta2", np.degrees(joint.theta2))
    """
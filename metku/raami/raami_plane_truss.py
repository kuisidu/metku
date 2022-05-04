# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Wed Dec 29 20:45:07 2021

Class for steel plane trusses, primarily tubular trusses

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt
from copy import copy

from raami.raami import Raami
from raami.frame_node import FrameNode
from raami.frame_member import FrameMember, SteelFrameMember, MultiSpanSteelMember, MemberGroup
from raami.frame_loads import PointLoad, LineLoad, LoadIDs
from raami.frame_supports import XYHingedSupport, YHingedSupport
from raami.exports import AbaqusOptions

from framefem.elements.ebbeam import EBBeam, EBBeam3D
from framefem.elements.rod import Rod
from sections.steel.RHS import RHS, SHS
from sections.steel.ISection import IPE
from eurocodes.en1993.en1993_1_8.rhs_joints import Line, RHSKGapJoint, RHSYJoint, RHSKTGapJoint
from eurocodes.en1993.en1993_1_8.en1993_1_8 import full_strength_weld_tube

from cost.workshop import Workshop


#from loadIDs import LoadIDs

def unit_vector(x1,x2):
    """ Creates a unit vector between points x1 and x2 """
    
    return (x2-x1)/np.linalg.norm(x2-x1)

class PlaneTruss(Raami):

    def __init__(self, origin=(0, 0), span=10000, left_height=1000,
                 mid_height=1500, right_height=1000, 
                 bot_chord_dx_left = 500, bot_chord_dx_right = 500,
                 rotation = 0,
                 num_elements=6,
                 name="Raami Plane Truss",
                 workshop=Workshop()):
        """ Constructor 
            Parameters:
            ------------
            :param origin: location of the origin (tuple)
            :param span: span of the truss [mm]
            :param left_height: hei
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
        super().__init__(name)
        
        self.origin = origin
        self.top_chord = []
        self.bottom_chord = []
        self.braces = {}
        self.joints = {}
        
        self.member_groups['top_chord'] = MemberGroup(members=[],name="Top Chord")
        self.member_groups['bottom_chord'] = MemberGroup(members=[],name="Bottom Chord")
        
        self.span = span
        
        # Top and bottom chord nodes
        self.top_nodes = []
        self.bottom_nodes = []
        
        if mid_height < left_height:
            raise ValueError('Truss height at mid-span cannot be smaller than height at left support.')
            
        if mid_height < right_height:
            raise ValueError('Truss height at mid-span cannot be smaller than height at right support.')
        
        self._H0 = left_height
        self._H1 = mid_height
        self._H2 = right_height
        
        self.n = 1
        self._dx_left = bot_chord_dx_left
        self._dx_right = bot_chord_dx_right
        
        # Workshop for cost calculations
        
        self.workshop = workshop
        self.costs = {}
        
        # Flag for stating whether or not the braces will be modelled as beams (True)
        # or rod elements (False)
        self.braces_as_beams = False
        #self.generate_outline()
        #self.top = topology
        #self.ndiv = ndiv
        
        #self.generate_topology(topology,ndiv)
        #self.generate_joints()
        #self.generate_chords()
        #self.generate_webs()

        #if q is not None:
        #    for tc in self.top_chords:
        #        self.add(LineLoad(tc, [q, q], 'y'))
    
    @property
    def Hleft(self):
        """ Height of the truss at left end """
        return self._H0
    
    @property
    def Hmid(self):
        """ Height of the truss at midspan """
        return self._H1
    
    @property
    def Hright(self):
        """ Height of the truss at right end """
        return self._H2

    @property
    def top_chord_profile(self):
        """ Returns top chord profile """
        return self.top_chord[0].cross_section
    
    @top_chord_profile.setter
    def top_chord_profile(self,val):
        """ Sets top chord profile """
        
        for chord in self.top_chord:
            chord.cross_section = val
    
    @property
    def bottom_chord_profile(self):
        """ Returns bottom chord profile """
        return self.bottom_chord[0].cross_section
    
    @bottom_chord_profile.setter
    def bottom_chord_profile(self,val):
        """ Sets bottom chord profile """
        
        for chord in self.bottom_chord:
            chord.cross_section = val
    
    def bottom_chord_dy(self,dy):
        """ Move the bottom chord nodes vertically by the distance dy """
        
        for node in self.bottom_nodes:
            node.y += dy
    
    def xmid(self):
        """ coordinates of the top chord at midspan """        
        return np.array([self.origin[0]+0.5*self.span,self.origin[1]+self.Hmid-self.Hleft])
    
    def xend(self):
        """ coordinates of the right end of top chord """
        return np.array([self.origin[0]+self.span,self.origin[1]+self.Hright-self.Hleft])
    
    def xbottom_left(self):
        """ coordinates of the bottom chord left node """
        return np.array([self.origin[0]+self._dx_left,self.origin[1]-self.Hleft])
    
    def xbottom_right(self):
        """ coordinates of the bottom chord left node """
        return np.array([self.origin[0]+self.span-self._dx_right,self.origin[1]-self.Hright])
    
    def dir_vec_origin_to_midspan(self):
        """ Unit vector from origin to the point at midspan """
        
        x0 = np.array(self.origin)
        
        return unit_vector(x0,self.xmid())
    
    def dir_vec_midspan_to_end(self):
        """ Unit vector from origin to the point at midspan """
                        
        return unit_vector(self.xmid(),self.xend())
    
    def generate_outline(self):
        """ Creates the nodes outlining the truss """
        
        self.add(FrameNode(self.origin))        
        self.add(FrameNode(self.xmid()))
        self.add(FrameNode(self.xend()))
        self.add(FrameNode(self.xbottom_left()))
        self.add(FrameNode(self.xbottom_right()))
    
    def add(self, this):
        """
        Adds object to truss
        :param this: object to be added
        :return:
        """
        # TOP CHORD
        if isinstance(this, TopChord):             
            this.mem_id = len(self.members)
            self.members["TC" + str(this.mem_id)] = this
            this.frame = self
            
            self.top_chord.append(this)
            self.member_groups['top_chord'].add_member(this)
            """ Add nodal coordinates for the member """
            #this.calc_nodal_coordinates(self.num_elements)

        # BOTTOM CHORD   
        elif isinstance(this, BottomChord):            
            this.mem_id = len(self.members)
            self.members["BC" + str(this.mem_id)] = this
            this.frame = self
            
            self.bottom_chord.append(this)
            self.member_groups['bottom_chord'].add_member(this)
                        
            """ Add nodal coordinates for the member """
            #this.calc_nodal_coordinates(self.num_elements)

        # BRACE
        elif isinstance(this, TrussBrace):            
            
            # Assign member id
            this.mem_id = len(self.members)
            # Assign joint id's   
            self.braces[this.mem_id] = this
            self.members["B" + str(this.mem_id)] = this
            this.frame = self
            
            if self.braces_as_beams:
                # If braces are modelled as beam elements, add hinges to
                # both ends
                this.mtype = 'beam'
                this.add_hinge(loc=0)
                this.add_hinge(loc=1)
            
        # JOINT
        elif isinstance(this, TubularJoint):
            this.joint_id = len(self.joints)
            self.joints[this.joint_id] = this
            
        else:
            super().add(this)
            
    
    def cost(self, verb=False):
        """ Calculates the cost of the truss """
        Ctot = 0
        Cbraces = 0.0
        Ctop = 0.0
        Cbottom = 0.0
        
        costs = {'material':0.0,'blasting':0.0,'painting':0.0,'welding':0.0, 'sawing':0.0}
        
        # Cost of braces
        for mem in self.braces.values():
            C = mem.cost()
            Cbraces += sum([val for val in C.values()])
            for key, value in mem.costs.items():
                costs[key] += value
        
        # Cost of top chord
        for mem in self.top_chord:
            C = mem.cost()
            Ctop += sum([val for val in C.values()])
            for key, value in mem.costs.items():
                costs[key] += value
        
        # Cost of bottom chord
        for mem in self.bottom_chord:
            C = mem.cost()
            Cbottom += sum([val for val in C.values()])    
            for key, value in mem.costs.items():
                costs[key] += value
                
        # Cost of welded joints
        Cjoints = 0.0
        for joint in self.joints.values():
            Cjoints += joint.cost()
            for key, value in joint.costs.items():
                costs[key] += value
        
        Ctot = Cbraces + Ctop + Cbottom + Cjoints
                
        
        if verb:
            print("Truss costs: ")
            print(f"Total cost: {Ctot:.2f} €")
            print(f"Top chord: {Ctop:.2f} € ({Ctop/Ctot*100:.1f} %)")
            print(f"Bottom chord: {Cbottom:.2f} € ({Cbottom/Ctot*100:.1f} %)")
            print(f"Braces: {Cbraces:.2f} € ({Cbraces/Ctot*100:.1f} %) ")
            print(f"Welded joints: {Cjoints:.2f} € ({Cjoints/Ctot*100:.1f} %)")
            
            print('\nCost Distribution:')
            for key, value in costs.items():
                print(f"{key.capitalize()}: {value:.2f} € ({value/Ctot*100:.1f}%)")
            
        
        self.costs = costs
        
        return Ctot
    
    def determine_symmetry_joints(self):
        
        for joint in self.joints.values():
            joint.symmetry_joint(self.joints)
    
    def design_joints(self,load_id,verb=False):
        """ Performs design of joints in the truss """
        
        for joint in self.joints.values():
            r = joint.design(load_id)
            
            if verb:
                if isinstance(r,float):
                    print(f"Joint {joint.__repr__()}: {r:.2f}")
                else:
                    s = f"Joint {joint.__repr__()}:"
                    for R in r:
                        s += f" {R:.2f} "
                    print(s)
    
    def to_abaqus(self,target_dir='C:/Users/kmela/Data/',filename="Raami",partname='Truss',options=AbaqusOptions()):
        """ Exports truss to abaqus 
            The export method of 'raami' is mainly used. The main point here is to
            write the MPCs corresponding to eccentricity elements
        """
        
        # options
        #opts = {'x_monitor':0.5*self.span, 'n_monitored':2, 'mpc':[],'ecc_elements':[]}
        
        # Create element sets for top chord and bottom chord
        # top chord is halved
        L2 = 0.5*self.span
        options.elsets['Top_chord_left'] = []
        options.elsets['Top_chord_right'] = []
        options.elsets['Bottom_chord'] = []
        
        for mem in self.top_chord:
            for el in mem.fem_elements:
                if all([n.x <= L2 for n in el.nodes]):
                    options.elsets['Top_chord_left'].append(el)
                else:
                    options.elsets['Top_chord_right'].append(el)
        
        for mem in self.bottom_chord:
            for el in mem.fem_elements:
                options.elsets['Bottom_chord'].append(el)
        
        for joint in self.joints.values():
            if isinstance(joint,TubularKGapJoint):
                if not(joint.fem_elements['ecc'] is None):
                    # This applies if the model 'en1993' was used
                    pass
                elif not(joint.fem_elements['left_ecc'] is None):
                    # Determine which node of the left eccentricity element
                    # is in the brace member. That is the slave node
                    n = joint.fem_elements['left_ecc'].nodes[0]
                    if n in joint.left_brace.fem_nodes:
                        options.mpc.append([n.nid+1,joint.fem_elements['left_ecc'].nodes[1].nid+1])
                    else:
                        options.mpc.append([joint.fem_elements['left_ecc'].nodes[1].nid+1,n.nid+1])
                    
                    n = joint.fem_elements['right_ecc'].nodes[0]
                    if n in joint.right_brace.fem_nodes:
                        options.mpc.append([n.nid+1,joint.fem_elements['right_ecc'].nodes[1].nid+1])
                    else:
                        options.mpc.append([joint.fem_elements['right_ecc'].nodes[1].nid+1,n.nid+1])
                    
                    options.ecc_elements.append(joint.fem_elements['left_ecc'])
                    options.ecc_elements.append(joint.fem_elements['right_ecc'])
                    
                    if joint.chord == "top":
                        if joint.ridge:
                            options.top_gap_elements += joint.fem_elements['gap']
                        else:
                            options.top_gap_elements.append(joint.fem_elements['gap'])
                    elif joint.chord == 'bottom':
                        options.bottom_gap_elements.append(joint.fem_elements['gap'])
        
        for el in options.top_gap_elements:
            if all([n.x <= L2 for n in el.nodes]):
                options.elsets['Top_chord_left'].append(el)
            else:
                options.elsets['Top_chord_right'].append(el)
        
        for el in options.bottom_gap_elements:
            options.elsets['Bottom_chord'].append(el)
            
        
        super().to_abaqus(target_dir,filename,partname,options)
                    
    

class SlopedTruss(PlaneTruss):
    """ Class for roof truss with a slope """
    
    def __init__(self,L1=10000,L2=10000,dx1=500,dx2=500,**kwargs):
        """
        Constructor

        Parameters
        ----------
        L1 : TYPE, optional
            DESCRIPTION. The default is 10000.
        L2 : TYPE, optional
            DESCRIPTION. The default is 10000.
        dx1 : TYPE, optional
            DESCRIPTION. The default is 500.
        dx2 : TYPE, optional
            DESCRIPTION. The default is 500.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.L1 = L1
        self.L2 = L2
        
        span = L1 + L2
        
        
        # The shape of the truss can be given by providing two of the 
        # following three dimensions. The third will be 
        h1 = 0
        h2 = 0
        self.slope = 0
        
        
        for key, value in kwargs.items():
            if key == "h1":
                h1 = value
            elif key == "h2":
                h2 = value
            elif key == "slope":
                self.slope = value
            #else:
            #    raise ValueError(f"Erroneous key {key}.")
        
        if h1 == 0:
            # User has given the height at ridge and slope
            dy = self.height_of_slope_left()
            h1 = h2-dy
        elif h2 == 0:
            # User has given the height from bottom chord to left
            # end of top chord and slope
            dy = self.height_of_slope_left()
            h2 = h1+dy
        else:
            # User has given the two heights, so slope must be
            # calculated
            dy = h2-h1
            self.slope = np.degrees(np.arctan(dy/self.L1))
            
        
        PlaneTruss.__init__(self,origin=(0, 0), span=span, left_height=h1,
                     mid_height=h2, right_height=h1, 
                     bot_chord_dx_left = dx1, bot_chord_dx_right = dx2,
                     rotation = 0,
                     num_elements=6,
                     name="Raami Sloped Roof Truss")
        
        # Joint at the ridge
        self.ridge_joint = None
    
    def height_of_slope_left(self):
        """ Height of the truss above the origin """
        return self.L1*np.tan(np.radians(self.slope))
    
    def weight_detailed(self,verb=False):
        """ Returns detailed weight 
            The weight of chords and braces are separated.
        """
        Wtop = 0
        Wbottom = 0
        Wbraces = 0
        
        for mem in self.members.values():
            W = mem.weight()
            
            if isinstance(mem,TopChord):
                Wtop += W
            elif isinstance(mem,BottomChord):
                Wbottom += W
            else:
                Wbraces += W
        
        Wtot = Wtop + Wbottom + Wbraces
        
        if verb:
            print(f"Truss weight: {Wtot:.2f} kg")
            print(f"Top Chord: {Wtop:.2f} kg ({Wtop/Wtot*100:.1f}%)")
            print(f"Bottom Chord: {Wbottom:.2f} kg ({Wbottom/Wtot*100:.1f}%)")
            print(f"Braces: {Wbraces:.2f} kg ({Wbraces/Wtot*100:.1f}%)")
        
        return Wtop, Wbottom, Wbraces
            
        
    def generate_topology(self,topology='K',ndiv=4, 
                          top_chord_profile = SHS(160,5),
                          bottom_chord_profile = SHS(140,5),
                          brace_profile = SHS(90,5),
                          first_diagonal_up=False,
                          edge_verticals=True,
                          nel_chord=4,
                          nel_brace=1):
        """ Generates truss members for a given topology 
        
            :param: topology .. 'K', 'N', or 'KT'
            :param: ndiv .. number of divisions on one half of the truss
            
        """
        
        
        if topology == 'K':
            # Generate nodes first
            # Bottom chord nodes:            
            dx1 = self._dx_left            
            dx2 = self._dx_right
            ycoord = self.origin[1]-self.Hleft
            
            if dx1 > 0:
                # If the first bottom chord node is not at x = 0, do the
                # following:
                # 1) divide half of the span into ndiv trios of two diagonals
                #   and a vertical.
                # 2) each trio is symmterical with respect to the vertical.
                # 3) for the first trio, the distance in x direction between
                #    nodes is dx1 and for the others it is Dx
                S = np.linspace(0,1,ndiv)
                # Distance between nodes in x direction
                #Dx = (self.L1-2*dx1)/(ndiv-1)
                Dx = (self.L1-dx1)/(ndiv-0.5)
            else:
                S = np.linspace(0,1,ndiv+1)
                if not first_diagonal_up:
                    S = np.linspace(0,1,ndiv)
                    Xbot0 = FrameNode(np.array([self.origin[0],ycoord]))
                    Dx = self.L1/(ndiv-0.5)
                else:
                    # Distance between nodes in x direction
                    Dx = self.L1/ndiv
            
            # First half
            if dx1 == 0 and not first_diagonal_up:
                X0 = np.array([self.origin[0]+0.5*Dx,ycoord])
                X2 = np.array([self.origin[0]+self.span-0.5*Dx,ycoord])
            else:
                X0 = np.array([self.origin[0]+dx1,ycoord])
                X2 = np.array([self.origin[0]+self.span-dx2,ycoord])
            X1 = np.array([self.origin[0]+self.L1-0.5*Dx,ycoord])  
            
            Xbot = [FrameNode(X0+s*(X1-X0)) for s in S]
        
            # Second half
            X1[0] += Dx
            
            Xbot += [FrameNode(X1+s*(X2-X1)) for s in S]
        
            if dx1 == 0 and not first_diagonal_up:
                Xbot.insert(0,Xbot0)
                Xbot.append(FrameNode(np.array([self.origin[0]+self.span,ycoord])))
        
            # Top chord nodes:
            # Local coordinates of the top chord
            #            
            X0 = np.array([self.origin[0],self.origin[1]])
            X2 = np.array([self.origin[0]+self.span,self.origin[1]])
                                 
            X1 = self.xmid()
            #X2 = self.xend()
            
            v1 = X1-X0
            v2 = X2-X1
            
            #Xtop = []
            
            if dx1 == 0:
                if first_diagonal_up:
                    Xtop = [FrameNode(X0)]
                    v0 = v1/np.linalg.norm(v1)
                    top_line = Line(v=v0,p1=X0)
                    
                    for Xb2, Xb1 in zip(Xbot[:ndiv+1],Xbot[1:ndiv+2]):
                        X = np.array([0.5*(Xb2.x+Xb1.x),Xb2.y])
                        bottom_line = Line(v=np.array([0.0,1.0]),p1=X)                
                        t, Xtop1 = top_line.intersect(bottom_line)
                        
                        #print(top_line, bottom_line)                    
                        Xtop.append(FrameNode(Xtop1))
                    
                    v0 = v2/np.linalg.norm(v2)
                    top_line = Line(v=v0,p1=Xtop[-1].coords)
                    
                    for Xb1, Xb2 in zip(Xbot[ndiv+1:],Xbot[ndiv+2:]):
                        X = np.array([0.5*(Xb2.x+Xb1.x),Xb2.y])
                        bottom_line = Line(v=np.array([0.0,1.0]),p1=X)                
                        t, Xtop1 = top_line.intersect(bottom_line)
                        
                        #print(top_line, bottom_line)                    
                        Xtop.append(FrameNode(Xtop1))
                    
                    Xtop.append(FrameNode(X2))
                    
                else:
                    S = np.linspace(0,1,ndiv+1)
                
                    Xtop = [FrameNode(X0+s*v1) for s in S]
                    Xtop += [FrameNode(X1+s*v2) for s in S[1:]]
                    
            else:
                Xtop = [FrameNode(X0)]
                v0 = v1/np.linalg.norm(v1)
                top_line = Line(v=v0,p1=X0)
                
                for Xb2, Xb1 in zip(Xbot[:ndiv],Xbot[1:ndiv+1]):
                    X = np.array([0.5*(Xb2.x+Xb1.x),Xb2.y])
                    bottom_line = Line(v=np.array([0.0,1.0]),p1=X)                
                    t, Xtop1 = top_line.intersect(bottom_line)
                    
                    #print(top_line, bottom_line)                    
                    Xtop.append(FrameNode(Xtop1))
                
                v0 = v2/np.linalg.norm(v2)
                top_line = Line(v=v0,p1=Xtop[-1].coords)
                
                for Xb1, Xb2 in zip(Xbot[ndiv:],Xbot[ndiv+1:]):
                    X = np.array([0.5*(Xb2.x+Xb1.x),Xb2.y])
                    bottom_line = Line(v=np.array([0.0,1.0]),p1=X)                
                    t, Xtop1 = top_line.intersect(bottom_line)
                    
                    #print(top_line, bottom_line)                    
                    Xtop.append(FrameNode(Xtop1))
                
                Xtop.append(FrameNode(X2))
            
            for Xt in Xtop:
                self.add(Xt) 
                        
            for Xb in Xbot:
                self.add(Xb) 
                                    
        
            # Generate top chord members
            for i in range(len(Xtop)-1):                
                self.add(TopChord([Xtop[i],Xtop[i+1]],copy(top_chord_profile),nel=nel_chord))
        
            # Generate bottom chord members
            for i in range(len(Xbot)-1):
                self.add(BottomChord([Xbot[i],Xbot[i+1]],copy(bottom_chord_profile),nel=nel_chord))
        
            
            # Generate braces
            
            if dx1 == 0:
                # Diagonals
                if first_diagonal_up:
                    # Diagonals run from bottom chord to top chord
                    for i in range(1,2*ndiv+2):
                        self.add(TrussBrace([Xtop[i],Xbot[i-1]],copy(brace_profile),nel=nel_brace))
                        self.add(TrussBrace([Xtop[i],Xbot[i]],copy(brace_profile),nel=nel_brace))
                else:
                    # Diagonals run from top chord to bottom chord                    
                    for i in range(1,2*ndiv+1):
                        self.add(TrussBrace([Xbot[i],Xtop[i-1]],copy(brace_profile),nel=nel_brace))
                        self.add(TrussBrace([Xbot[i],Xtop[i]],copy(brace_profile),nel=nel_brace))
                
                # Close with verticals
                if edge_verticals:
                    self.add(TrussBrace([Xbot[0],Xtop[0]],copy(brace_profile),nel=nel_brace))
                    self.add(TrussBrace([Xbot[-1],Xtop[-1]],copy(brace_profile),nel=nel_brace))
            
            else:
                # If the first node of the bottom chord is not at the same
                # x location as the first node of the top chord,
                # the first diagonal goes from top to bottom chord.
                for i in range(0,2*ndiv):
                    self.add(TrussBrace([Xbot[i],Xtop[i]],copy(brace_profile),nel=nel_brace))
                    self.add(TrussBrace([Xbot[i],Xtop[i+1]],copy(brace_profile),nel=nel_brace))
                
                
           
        elif topology == 'N':
            # Generate nodes first
            
            # Bottom chord nodes:            
            dx1 = self._dx_left            
            dx2 = self._dx_right
            ycoord = self.origin[1]-self.Hleft
            
            if dx1 > 0:
                S = np.linspace(0,1,ndiv)
            else:
                S = np.linspace(0,1,ndiv+1)
            
            # First half
            X0 = np.array([self.origin[0]+dx1,ycoord])
            X1 = np.array([self.origin[0]+0.5*self.span,ycoord])  
            X2 = np.array([self.origin[0]+self.span-dx2,ycoord])
                        
            Xbot = [FrameNode(X0+s*(X1-X0)) for s in S]
            
            # Second half
            #X1 = np.array([self.origin[0]+0.5*self.span+dx1,ycoord])  
            Xbot += [FrameNode(X1+s*(X2-X1)) for s in S[1:]]
        
                    
            # Top chord nodes:
            # Local coordinates of the top chord
            
            
            #            
            X0 = np.array([self.origin[0],self.origin[0]])
            X2 = np.array([self.origin[0]+self.span,self.origin[0]])
                                 
            X1 = self.xmid()
            X2 = self.xend()
            
            v1 = X1-X0
            v2 = X2-X1
            
            if dx1 == 0:
                S = np.linspace(0,1,ndiv+1)
                Xtop = [FrameNode(X0+s*v1) for s in S]
                Xtop += [FrameNode(X1+s*v2) for s in S[1:]]
            else:
                S = np.linspace(0,1,ndiv)
                Xtop = [FrameNode(X0)]
                v0 = v1/np.linalg.norm(v1)
                top_line = Line(v=v0,p1=X0)
                
                for X in Xbot[:ndiv]:                    
                    bottom_line = Line(v=np.array([0.0,1.0]),p1=X.coords)                
                    t, Xtop1 = top_line.intersect(bottom_line)
                    
                    #print(top_line, bottom_line)                    
                    Xtop.append(FrameNode(Xtop1))
            
                v0 = v2/np.linalg.norm(v2)
                top_line = Line(v=v0,p1=Xtop[-1].coords)
                
                for X in Xbot[ndiv:]:
                    bottom_line = Line(v=np.array([0,1]),p1=X.coords)                
                    t, Xtop1 = top_line.intersect(bottom_line)                
                    Xtop.append(FrameNode(Xtop1))
                
                Xtop.append(FrameNode(X2))
                        
            for Xt in Xtop:
                self.add(Xt) 
                        
            for Xb in Xbot:
                self.add(Xb)    
            
            
            # Generate top chord members
            for i in range(len(Xtop)-1):                
                self.add(TopChord([Xtop[i],Xtop[i+1]],copy(top_chord_profile),nel=nel_chord))
        
            # Generate bottom chord members
            for i in range(len(Xbot)-1):
                self.add(BottomChord([Xbot[i],Xbot[i+1]],copy(bottom_chord_profile),nel=nel_chord))
            
            # Generate braces
            if dx1 == 0:
                # Verticals
                if edge_verticals:
                    for X, Y in zip(Xtop,Xbot):
                        self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                else:
                    for X, Y in zip(Xtop[1:-1],Xbot[1:-1]):
                        self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                # Diagonals
                if first_diagonal_up:
                    # Diagonals run from bottom chord to top chord
                    for X, Y in zip(Xtop[1:ndiv+1],Xbot[:ndiv]):
                        self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                        
                    for X, Y in zip(Xtop[ndiv:-1],Xbot[ndiv+1:]):
                        self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                else:
                    # Diagonals run from top chord to bottom chord
                    for X, Y in zip(Xtop[:ndiv],Xbot[1:ndiv+1]):
                        self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                        
                    for X, Y in zip(Xtop[ndiv+1:],Xbot[ndiv:-1]):
                        self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
            else:
                # Lower chord is shorter than top chord. Now the first
                # diagonal is downward anyway.
                
                # Verticals
                for X, Y in zip(Xtop[1:-1],Xbot):
                    self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                    
                # Diagonals run from top chord to bottom chord
                for X, Y in zip(Xtop[:ndiv],Xbot[:ndiv]):
                    self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                    
                for X, Y in zip(Xtop[ndiv+1:],Xbot[ndiv-1:]):
                    self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
        
        elif topology == "KT":
            # KT truss:
            # Generate nodes first
            
            # Bottom chord nodes:            
            dx1 = self._dx_left            
            dx2 = self._dx_right
            ycoord = self.origin[1]-self.Hleft
            
            if dx1 > 0:
                # If the first bottom chord node is not at x = 0, do the
                # following:
                # 1) divide half of the span into ndiv trios of two diagonals
                #   and a vertical.
                # 2) each trio is symmterical with respect to the vertical.
                # 3) for the first trio, the distance in x direction between
                #    nodes is dx1 and for the others it is Dx
                S = np.linspace(0,1,ndiv)
                # Distance between nodes in x direction
                #Dx = (self.L1-2*dx1)/(ndiv-1)
                Dx = (self.L1-dx1)/(ndiv-0.5)
            else:
                # If the first bottom chord node is not at x = 0, do the
                # following:
                # 1) divide half of the span into ndiv trios of two diagonals
                #   and a vertical.
                # 2) each trio is symmterical with respect to the vertical.
                # 3) if the first diagonal is up, the distance of it in x direction
                #    is 0.5*Dx. If the first diagonal is down, all, trios have
                #    the distance in x direction between nodes 0.5*Dx
                S = np.linspace(0,1,ndiv+1)
                if not first_diagonal_up:
                    S = np.linspace(0,1,ndiv)
                    Xbot0 = FrameNode(np.array([self.origin[0],ycoord]))
                    Dx = self.L1/ndiv
                else:
                    # Distance between nodes in x direction
                    Dx = 2*self.L1/(2*ndiv+1)
            
            # First half
            if dx1 == 0 and not first_diagonal_up:
                X0 = np.array([self.origin[0]+0.5*Dx,ycoord])
                X2 = np.array([self.origin[0]+self.span-0.5*Dx,ycoord])
            else:
                X0 = np.array([self.origin[0]+dx1,ycoord])
                X2 = np.array([self.origin[0]+self.span-dx2,ycoord])
            X1 = np.array([self.origin[0]+self.L1-0.5*Dx,ycoord])  
            
            
            Xbot = [FrameNode(X0+s*(X1-X0)) for s in S]
            
            
            
            # Second half
            X1[0] += Dx
            
            Xbot += [FrameNode(X1+s*(X2-X1)) for s in S]
        
            if dx1 == 0 and not first_diagonal_up:
                Xbot.insert(0,Xbot0)
                Xbot.append(FrameNode(np.array([self.origin[0]+self.span,ycoord])))
        
            # Top chord nodes:
            # Local coordinates of the top chord
            #            
            X0 = np.array([self.origin[0],self.origin[0]])
            X2 = np.array([self.origin[0]+self.span,self.origin[0]])
                                 
            X1 = self.xmid()
            X2 = self.xend()
            
            v1 = X1-X0
            v2 = X2-X1
            
            #Xtop = []
            
            if dx1 == 0:
                if first_diagonal_up:
                    S = np.linspace(0,1,2*ndiv+2)                    
                else:
                    S = np.linspace(0,1,2*ndiv+1)
                
                Xtop = [FrameNode(X0+s*v1) for s in S]
                Xtop += [FrameNode(X1+s*v2) for s in S[1:]]
                    
            else:
                S = np.linspace(0,1,ndiv)
                Xtop = [FrameNode(X0)]
                v0 = v1/np.linalg.norm(v1)
                top_line = Line(v=v0,p1=X0)
                
                for X in Xbot[:ndiv]:                    
                    bottom_line = Line(v=np.array([0.0,1.0]),p1=X.coords)                
                    t, Xtop1 = top_line.intersect(bottom_line)
                    
                    #print(top_line, bottom_line)                    
                    Xtop.append(FrameNode(Xtop1))
                    
                    p1New = X.coords + np.array([0.5*Dx,0])
                    bottom_line = Line(v=np.array([0.0,1.0]),p1=p1New)
                    t, Xtop2 = top_line.intersect(bottom_line)
                    
                    #print(top_line, bottom_line)                    
                    Xtop.append(FrameNode(Xtop2))
                    
                v0 = v2/np.linalg.norm(v2)
                top_line = Line(v=v0,p1=Xtop[-1].coords)
                
                for X in Xbot[ndiv:-1]:
                    bottom_line = Line(v=np.array([0,1]),p1=X.coords)                
                    t, Xtop1 = top_line.intersect(bottom_line)                
                    Xtop.append(FrameNode(Xtop1))
                    
                    p1New = X.coords + np.array([0.5*Dx,0])
                    bottom_line = Line(v=np.array([0.0,1.0]),p1=p1New)
                    t, Xtop2 = top_line.intersect(bottom_line)
                    
                    #print(top_line, bottom_line)                    
                    Xtop.append(FrameNode(Xtop2))
                
                bottom_line = Line(v=np.array([0,1]),p1=Xbot[-1].coords)                
                t, Xtop1 = top_line.intersect(bottom_line)                
                Xtop.append(FrameNode(Xtop1))
                
                Xtop.append(FrameNode(X2))
            
            for Xt in Xtop:
                self.add(Xt) 
                        
            for Xb in Xbot:
                self.add(Xb) 
                
            # Generate top chord members
            for i in range(len(Xtop)-1):                
                self.add(TopChord([Xtop[i],Xtop[i+1]],copy(top_chord_profile),nel=nel_chord))
        
            # Generate bottom chord members
            for i in range(len(Xbot)-1):
                self.add(BottomChord([Xbot[i],Xbot[i+1]],copy(bottom_chord_profile),nel=nel_chord))
            
            # Generate braces
            if dx1 == 0:
                if first_diagonal_up:
                    # Verticals
                    if edge_verticals:
                        for X, Y in zip(Xtop[::2],Xbot):
                            self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                    else:
                        for X, Y in zip(Xtop[2::2],Xbot[1:-1]):
                            self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
                
                    # Diagonals
                    self.add(TrussBrace([Xbot[0],Xtop[1]],copy(brace_profile),nel=nel_brace))                    
                    
                    for i in range(1,len(Xbot)-1):
                        self.add(TrussBrace([Xtop[2*i-1],Xbot[i]],copy(brace_profile),nel=nel_brace))
                        self.add(TrussBrace([Xtop[2*i+1],Xbot[i]],copy(brace_profile),nel=nel_brace))
                        
                    
                    self.add(TrussBrace([Xbot[-1],Xtop[-2]],copy(brace_profile),nel=nel_brace))
                else:   
                    if edge_verticals:
                        self.add(TrussBrace([Xbot[0],Xtop[0]],copy(brace_profile),nel=nel_brace))
                        self.add(TrussBrace([Xbot[-1],Xtop[-1]],copy(brace_profile),nel=nel_brace))
                   
                    for i in range(1,len(Xbot)-1):                           
                        self.add(TrussBrace([Xtop[2*i-2],Xbot[i]],copy(brace_profile),nel=nel_brace))
                        self.add(TrussBrace([Xtop[2*i-1],Xbot[i]],copy(brace_profile),nel=nel_brace))
                        self.add(TrussBrace([Xtop[2*i],Xbot[i]],copy(brace_profile),nel=nel_brace))
                                       
                    
                   
                    """
                    if first_diagonal_up:
                        for i in range(1,len(Xbot)-2):                           
                            self.add(TrussBrace([Xtop[2*i-2],Xbot[i]],brace_profile))
                            self.add(TrussBrace([Xtop[2*i-1],Xbot[i]],brace_profile))
                            self.add(TrussBrace([Xtop[2*i+1],Xbot[i]],brace_profile))
                    else:
                        for i in range(1,len(Xbot)-1):                           
                            self.add(TrussBrace([Xtop[2*i-2],Xbot[i]],brace_profile))
                            self.add(TrussBrace([Xtop[2*i-1],Xbot[i]],brace_profile))
                            self.add(TrussBrace([Xtop[2*i],Xbot[i]],brace_profile))
                       
                        self.add(TrussBrace([Xbot[-1],Xtop[-1]],brace_profile))
                    """
            else:
                # Lower chord is shorter than top chord. Now the first
                # diagonal is downward anyway.
                
                # Verticals
                #for X, Y in zip(Xtop[1::2],Xbot):
                #    self.add(TrussBrace([X,Y],brace_profile))
                    
                for i in range(len(Xbot)):
                    #print(i)                       
                    self.add(TrussBrace([Xtop[2*i],Xbot[i]],copy(brace_profile),nel=nel_brace))
                    self.add(TrussBrace([Xtop[2*i+1],Xbot[i]],copy(brace_profile),nel=nel_brace))
                    self.add(TrussBrace([Xtop[2*i+2],Xbot[i]],copy(brace_profile),nel=nel_brace))
        """    
        fig, ax = plt.subplots(1)
        
        ax.plot([X[0] for X in Xtop], [X[1] for X in Xtop],'ob')
        ax.plot([X[0] for X in Xbot], [X[1] for X in Xbot],'ob')
        ax.set_aspect('equal')
        """  
        self.top_nodes = Xtop
        self.bottom_nodes = Xbot
    def generate_supports(self):
        """ Creates supports to the truss """
        
        if self._dx_left > 0:
            # If the first node of the bottom chord is not at the left edge,
            # the first node of the top chord is supported
            self.add(XYHingedSupport(self.top_nodes[0]))
        else:
            # Support the first node of the bottom chords
            self.add(XYHingedSupport(self.bottom_nodes[0]))
        
        if self._dx_right > 0:
            self.add(YHingedSupport(self.top_nodes[-1]))
        else:
            self.add(YHingedSupport(self.bottom_nodes[-1]))
    
    def generate_uniform_load(self,q,chord="top",load_id=LoadIDs['ULS'],ltype='live'):
        """ Creates a uniform line load to one of the chords """

        vals = [q,q]

        if chord == "top":
            for mem in self.top_chord:
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, 'LineLoad', "global"))
        elif chord == "bottom":            
            for mem in self.bottom_chord:
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, 'LineLoad', "global"))
    
    def generate_joints(self):
        """ Creates joints for the truss. These are welded tubular joints by default. """
        
        joints = {}
        
        # First, determine which members are connected to each node
        # This information is stored in 'joints' dict.
        for node in self.nodes:
            new_joint = {'chord':[],'braces':[]}
            for mem in self.top_chord:
                if node in mem.nodes:
                    new_joint['chord'].append(mem)
                    if len(new_joint['chord'])==2:
                        break                    
            
            for mem in self.bottom_chord:
                if node in mem.nodes:
                    new_joint['chord'].append(mem)
                    if len(new_joint['chord'])==2:
                        break
            
            for mem in self.braces.values():
                if node in mem.nodes:
                    new_joint['braces'].append(mem)
                    # The brace knows to which joints it belongs.
                    if node in self.top_nodes:
                        mem.top_joint = new_joint
                    else:
                        mem.bottom_joint = new_joint
            
            
            joints[node] = new_joint
        
        # Create joints. As a first guess, the joint type is
        # based on number of braces. This may not be accurate for
        # determining X joints, because both X and K joints have two braces.
        # The user must probably provide some additional information here to
        # distinguish the proper joint type, or the type will be based on
        # the results of structural analysis.
        njoint = 0
        for node, joint in joints.items():
            
            nbraces = len(joint['braces'])
            
            if nbraces == 3:
                newJoint = TubularKTGapJoint(node, joint['chord'], joint['braces'])                 
                self.add(newJoint)                
                njoint += 1
                
                if node in self.top_nodes and node.x == self.L1:
                    self.ridge_joint = newJoint
                
            elif nbraces == 2:
                # K joint
                newJoint = TubularKGapJoint(node, joint['chord'], joint['braces'])                 
                self.add(newJoint)                
                njoint += 1
                
                if node in self.top_nodes and node.x == self.L1:
                    self.ridge_joint = newJoint
                    
            elif nbraces == 1:                
                newJoint = TubularYJoint(node, joint['chord'], joint['braces'])
                self.add(newJoint)
                njoint += 1
            elif nbraces == 0:
                # Number of braces can be 0, if top or bottom chord is supported
                # at the node
                newJoint = None
            
            # Attach the joint to the corresponding braces
            for brace in joint['braces']:
                if node in self.top_nodes:
                    brace.top_joint = newJoint
                else:
                    brace.bottom_joint = newJoint
            
            # Attach the joint to the corresponding chord members
            if not newJoint is None:
                for chord in joint['chord']:                
                    chord.joints.append(newJoint)
        
        # Set the actual geometry for all joints
        for joint in self.joints.values():                      
            if isinstance(joint,TubularYJoint):
                angle = joint.brace_angle(joint.braces[0],geometry="actual")
                joint.joint.angle = angle
            elif isinstance(joint,TubularKGapJoint):
                # The gap calculation routine updates the
                # brace angles and gap
                #print(joint.gap)
                joint.calc_gap()                
                #print(joint.gap)
    
    def clear_fem(self):
        """ Clears FEM data """
        super().clear_fem()
        
        for joint in self.joints.values():
            joint.clear_fem()
        
    def generate_fem(self,model="no_eccentricity"):
        """
        Creates FEM model

        Parameters
        ----------
        model : string, optional
            Type of FEM model. Allowed values are:
                "simple": all members are pin-jointed bars
                "no_eccentricity": braces meet at originally created nodes
                "en1993": for K joints, the braces meet at the intersection point of their
                        centerlines, and an eccentricity element is added
                "braces2chord_faces": braces end at chord faces.
                The default is "no_eccentricity".

        Returns
        -------
        None.

        """
                
        
        if model == "no_eccentricity":
            # Here, we can use the fem model generation of Raami class
            super().generate_fem()
            
            for joint in self.joints.values():
                # Generate fem nodes for the joints
                x = joint.node.coords
                for node in self.fem.nodes:
                    if np.linalg.norm(node.coord-x) < 1e-6:
                        joint.fem_nodes["xc"] = node
                        break
                    
        elif model == "en1993":
            #print("FEM Model according to EN 1993")
            # This is more complicated:
            # 1) braces are bars that meet at the intersection points of centerlines
            for joint in self.joints.values():
                if isinstance(joint,TubularYJoint):
                    # For Y joints, a single FEM node is created and this coincides
                    # with the point of the FrameNode object, which is located at the
                    # center line
                    newNode = self.fem.add_node(joint.node.x,joint.node.y)                    
                    joint.fem_nodes['xc'] = newNode
                    joint.node.fem_node = newNode
                elif isinstance(joint,TubularKGapJoint):
                    # For K gap joints, we introduce two FEM nodes:
                    # 1) for point xp (intersection of brace centerlines)
                    # 2) for point xc on centerline of the chord that is closest
                    #    to xp. If xc and xp coincide, they will be the one and the
                    #    same node.
                    xp = joint.xp_coord
                    xpNode = self.fem.add_node(xp[0],xp[1])
                    if joint is self.ridge_joint:
                        xc = joint.node.coords
                    else:
                        xc = joint.xchord()
                    
                    if np.linalg.norm(xp-xc) < 1e-6:
                        xcNode = xp
                    else:
                        xcNode = self.fem.add_node(xc[0],xc[1])
                        
                    joint.fem_nodes['xc'] = xcNode
                    joint.fem_nodes['xp'] = xpNode
                    joint.node.fem_node = newNode
                elif isinstance(joint,TubularKTGapJoint):
                    # For KT gap joints, 
                    newNode = self.fem.add_node(joint.node.x,joint.node.y)                    
                    joint.fem_nodes['xc'] = newNode
                    joint.node.fem_node = newNode
            
            for mem in self.braces.values():
                #Generate member elements                
                mem.generate_elements(self.fem,model,self.braces_as_beams)
            
            for mem in self.top_chord:
                mem.generate_elements(self.fem,model)
        
            for mem in self.bottom_chord:
                mem.generate_elements(self.fem,model)
            
            # TODO! Change eccentricity elements to rigid links!
            # (Is it needed, or do IPE 500 elements suffice?)
            for joint in self.joints.values():
                if isinstance(joint,TubularKGapJoint):
                    n1 = joint.fem_nodes['xp']
                    n2 = joint.fem_nodes['xc']
                    cs = IPE(500)
                                 
                    newElement = EBBeam(n1,n2,cs,cs.material)                    
                    
                    self.fem.add_element(newElement)
                    joint.fem_elements['ecc'] = newElement
            
            # Generate supports
            for supp in self.supports.values():
                supp.add_support(self.fem)
                            
            # Generate loads
            for load in self.loads.values():
                load.add_load(self.fem)
                # Creates new loadcase if one with same load_id does not exist
                lcase_ids = [lc.load for lc in self.fem.loadcases.values()]
                if load.load_id not in lcase_ids:
                    self.fem.add_loadcase(supp_id=1,load_id=load.load_id)
            
            # Set nodal degrees of freedom
            self.fem.nodal_dofs()
        elif model == "ecc_elements":
            # Eccentricity elements are added for K joints
            # from the intersection of brace centerline and chord
            # face to the chord centerline.
            for joint in self.joints.values():
                if isinstance(joint,TubularYJoint):      
                    # For Y joints, no eccentricity elements are used.
                    newNode = self.fem.add_node(joint.node.x,joint.node.y)                    
                    joint.fem_nodes['xc'] = newNode
                    joint.node.fem_node = newNode
                elif isinstance(joint,TubularKGapJoint):
                    # For K gap joints, we introduce the following FEM nodes:
                    # 1) intersection of each brace and chord face
                    # 2) chord centerline points closest to the chord face points.
                    
                    
                    """
                    if joint is self.ridge_joint:
                        xc = joint.node.coords
                        xcNode = self.fem.add_node(xc[0],xc[1])
                        joint.fem_nodes['xc'] = xcNode
                        
                        xface, xcenter = joint.brace_chord_face_x()
                        print(xface,xcenter)
                    else:
                        xface, xcenter = joint.brace_chord_face_x()
                        
                        for brace, xf in xface.items():
                            newNode = self.fem.add_node(xf[0],xf[1])
                            joint.fem_nodes['xcf'][brace] = newNode
                        
                        for brace, xc in xcenter.items():
                            newNode = self.fem.add_node(xc[0],xc[1])
                            joint.fem_nodes['xch'][brace] = newNode
                    """
                    xface, xcenter = joint.brace_chord_face_x()
                    for brace, xf in xface.items():
                        newNode = self.fem.add_node(xf[0],xf[1])
                        joint.fem_nodes['xcf'][brace] = newNode
                    
                    for brace, xc in xcenter.items():
                        newNode = self.fem.add_node(xc[0],xc[1])
                        joint.fem_nodes['xch'][brace] = newNode
                    
                    # For the apex, create a FEM node to the intersection
                    # of chord centerlines
                    if joint is self.ridge_joint:
                        xc = joint.node.coords
                        xcNode = self.fem.add_node(xc[0],xc[1])
                        joint.fem_nodes['xc'] = xcNode
                        joint.node.fem_node = xcNode
                    
                    """
                    xpNode = self.fem.add_node(xp[0],xp[1])
                    if joint is self.ridge_joint:
                        xc = joint.node.coords
                    else:
                        xc = joint.xchord()
                    """
                    
                    """
                    if np.linalg.norm(xp-xc) < 1e-6:
                        xcNode = xp
                    else:
                        xcNode = self.fem.add_node(xc[0],xc[1])
                        
                    joint.fem_nodes['xc'] = xcNode
                    joint.fem_nodes['xp'] = xpNode
                    joint.node.fem_node = newNode
                    """
                elif isinstance(joint,TubularKTGapJoint):
                    # For KT gap joints, currently only a node in the
                    # initial FrameNode location is created.
                    newNode = self.fem.add_node(joint.node.x,joint.node.y)                    
                    joint.fem_nodes['xc'] = newNode
                    joint.node.fem_node = newNode
            
            for mem in self.braces.values():
                #Generate member elements                
                mem.generate_elements(self.fem,model,self.braces_as_beams)
            
            for mem in self.top_chord:
                mem.generate_elements(self.fem,model)
        
            for mem in self.bottom_chord:
                mem.generate_elements(self.fem,model)
            
            # TODO! Change eccentricity elements to rigid links!
            # (Is it needed, or do IPE 500 elements suffice?)
            for joint in self.joints.values():
                if isinstance(joint,TubularKGapJoint):
                    # Eccentricity elements
                    n1 = joint.fem_nodes['xcf'][joint.left_brace]
                    nc1 = joint.fem_nodes['xch'][joint.left_chord]
                    cs = IPE(500)
                                 
                    newElement = EBBeam(n1,nc1,cs,cs.material)                    
                    
                    self.fem.add_element(newElement)
                    joint.fem_elements['left_ecc'] = newElement
                    
                    n1 = joint.fem_nodes['xcf'][joint.right_brace]
                    nc2 = joint.fem_nodes['xch'][joint.right_chord]
                                 
                    newElement = EBBeam(n1,nc2,cs,cs.material)
                    
                    joint.fem_elements['right_ecc'] = newElement
                    
                    self.fem.add_element(newElement)
                    
                    # Add gap elements
                    if joint.ridge:
                        nApex = joint.fem_nodes['xc']
                        cs_chord = joint.chords[0].cross_section
                        newElement1 = EBBeam(nc1,nApex,cs_chord,cs_chord.material)
                        self.fem.add_element(newElement1)
                        newElement2 = EBBeam(nApex,nc2,cs_chord,cs_chord.material)
                        self.fem.add_element(newElement2)
                        joint.fem_elements['gap'] = [newElement1,newElement2]
                    else:
                        cs_chord = joint.chords[0].cross_section
                        newElement = EBBeam(nc1,nc2,cs_chord,cs_chord.material)
                        self.fem.add_element(newElement)
                        joint.fem_elements['gap'] = newElement
            
            # Generate supports
            for supp in self.supports.values():
                supp.add_support(self.fem)
                            
            # Generate loads
            for load in self.loads.values():
                load.add_load(self.fem)
                # Creates new loadcase if one with same load_id does not exist
                lcase_ids = [lc.load for lc in self.fem.loadcases.values()]
                if load.load_id not in lcase_ids:
                    self.fem.add_loadcase(supp_id=1,load_id=load.load_id)
            
            # Set nodal degrees of freedom
            self.fem.nodal_dofs()
            
        else:
            print("Weird FEM model.")
    
    def plot(self, geometry=False, print_text=True, show=True,
             loads=True, color=False, axes=None, save=False, mem_dim=False):
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
        if not geometry:
            super().plot(print_text, show, loads, color, axes, save, mem_dim)
        else:
            
            if axes is None:
                fig, ax = plt.subplots(1)
            else:
                ax = axes
    
            # Plot members
            for mem in self.top_chord:
                mem.plot(print_text, color, ax, mem_dim=True)
            
            for mem in self.bottom_chord:
                mem.plot(print_text, color, ax, mem_dim=True)
            
            for mem in self.braces.values():                            
                mem.plot(print_text, color, ax, mem_dim=True, geometry=geometry)
    
            # Plot joints
            #for joint in self.joints.values():
            #    joint.plot(color=color)
    
            # Plot supports
            
            for support in self.supports.values():
                node_coord = support.node.coords
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
                
            ax.axis('equal')
            if save:
                plt.savefig('default.svg', format='svg')
            if show:
                plt.axis('equal')
                plt.show()
                
    def optimize_members(self, prof_type="CURRENT",verb=False,**kwargs):    
        """ Finds minimum weight profiles for members """
        
        # TODO!
        # Optimoinnissa voidaan haluta:
        # Sauvoille tietty materiaali (yläpaarre, alapaarre, uumasauvat)
        # Sauvoille tietty poikkileikkausluokka (1 tai 2)
        #
        top = {'material': 'S355', 'class': 2}
        bottom = {'material': 'S355', 'class': 2}
        braces = {'material': 'S355', 'class': 2}
            
        for key, value in kwargs.items():
            if key == 'top':
                top = value
            elif key == 'bottom':
                bottom = value
            elif key == 'braces':
                braces = value
        
        kmax = 10
        k = 1
                        
        while k < kmax:            
            explored = []
            MEM_PROFILES_CHANGED = []
            self.structural_analysis('all','REM')
            if verb:
                print(f"Optimize profiles, iteration {k}")
                
            # Go through groups first
            for name, group in self.member_groups.items():
                if name == "top_chord":
                    material = top['material']
                    sec_class = top['class']
                elif name == 'bottom_chord':
                    material = bottom['material']
                    sec_class = bottom['class']
                    
                group.optimum_design(prof_type,verb,material,sec_class)
                
                for mem in group.members:
                    explored.append(mem)
                
            
            for member in self.members.values():
                #print(member)                
                if not member in explored:
                    if verb:
                        print(member)
                    explored.append(member)
                    if isinstance(member,TrussBrace):
                        material = braces['material']
                        sec_class = braces['class']
                    else:
                        material = 'S355'
                        sec_class = 2
                    MEM_PROFILES_CHANGED.append(member.optimum_design(prof_type,verb,material,sec_class))
                                        
                    if not member.symmetry_pair is None:
                        explored.append(member.symmetry_pair)
                    
                    #print(member.cross_section)
            
            # If none of the members changed profile, exit the loop
            if not any(MEM_PROFILES_CHANGED):
                break
                
            k += 1
        
        if verb:
            print(f"Number of iterations: {k}.")
        
        # First run the optimization method of Raami class
        # super().optimize_members(prof_type,verb)
        
        
        
    
class TrussMember(SteelFrameMember):
    
    def __init__(self,nodes,section,mem_type='bar',mem_id="",nel=4,hinges=None):

        # List of joints to which the member is attached
        self.joints = []
        
        if hinges == None:
            hinges = [False,False]

        SteelFrameMember.__init__(self,nodes,section,mem_type,mem_id,nel,hinges)

class TopChord(SteelFrameMember):
    def __init__(self,nodes,section,mem_id="",nel=4,hinges=None,lcr=[0.9,1.0]):

        mem_type = 'top_chord'
        
        if hinges == None:
            hinges = [False,False]
                
        SteelFrameMember.__init__(self,nodes,section,mem_type,mem_id,nel,hinges=hinges,lcr=lcr)

        self.mtype = 'top_chord'
        self.steel_members = []

        self.joints = []
    
    @property
    def perpendicular(self):
        """
        Returns vector perpendicular to member.
        NOTE! This needs to be modified for 3D!        
        :return:
        """
        unit = self.dir_vector            

        return np.array([-unit[1], unit[0]])
    
    def cost(self):
        """ Cost of top chord member. This includes:
            i) material cost
            ii) blasting cost
            iii) painting cost
        """
        Ctot = {'material':self.material_cost(),
                'blasting':self.blasting_cost(),
                'painting':self.painting_cost()}
        
        self.costs = Ctot
        
        return Ctot

    def generate_elements(self,fem,model="no_eccentricity"):
        """ Generate finite elements """
        
        if model == "no_eccentricity":
            super().generate_elements(fem)
        elif model == "en1993":
            # First end node
            self.fem_nodes.append(self.joints[0].fem_nodes['xc'])
            
            # Generate internal nodes:
            # global_node_coord include also the end node coordinates,
            # so they are not used in the iteration
            for x in self.global_node_coords[1:-1]:                
                if self.frame.dim == 2:
                    newNode = fem.add_node(x[0],x[1])    
                else:
                    newNode = fem.add_node(x[0],x[1],x[2])
                    
                self.fem_nodes.append(newNode)
            
            # Last end node
            self.fem_nodes.append(self.joints[1].fem_nodes['xc'])
            
            # Create members
            # The zip command allows to create two iterables, where the
            # items are actually consecutive items of the list fem_nodes.
            for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):
                if self.frame.dim == 2:                   
                    newElement = EBBeam(n1,n2,self.cross_section,self.material)     
                else:
                    newElement = EBBeam3D(n1,n2,self.cross_section,self.material)     
                
                fem.add_element(newElement)
                self.fem_elements.append(newElement)
        
        elif model == "ecc_elements":
            if isinstance(self.joints[0],TubularYJoint):
                self.fem_nodes.append(self.joints[0].fem_nodes['xc'])                
            elif isinstance(self.joints[0],TubularKGapJoint):
                # in case of K gap joint, the node is created at the
                # intersection between brace center line and chord face            
                self.fem_nodes.append(self.joints[0].fem_nodes['xch'][self])
            
            # Generate internal nodes:
            # global_node_coord include also the end node coordinates,
            # so they are not used in the iteration
            for x in self.global_node_coords[1:-1]:                
                if self.frame.dim == 2:
                    newNode = fem.add_node(x[0],x[1])    
                else:
                    newNode = fem.add_node(x[0],x[1],x[2])
                    
                self.fem_nodes.append(newNode)
        
            
            if isinstance(self.joints[1],TubularYJoint):
                self.fem_nodes.append(self.joints[1].fem_nodes['xc'])                
            elif isinstance(self.joints[1],TubularKGapJoint):
                # in case of K gap joint, the node is created at the
                # intersection between brace center line and chord face            
                self.fem_nodes.append(self.joints[1].fem_nodes['xch'][self])
            
            # Add elements between first and last node.
            for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):
                if self.frame.dim == 2:                   
                    newElement = EBBeam(n1,n2,self.cross_section,self.material)                    
                else:
                    newElement = EBBeam3D(n1,n2,self.cross_section,self.material)
                    
                fem.add_element(newElement)
                self.fem_elements.append(newElement)
        
        if self.hinges[0]:
            """ There is a hinge in the first end """
            self.fem_elements[0].releases = [2]
        
        if self.hinges[1]:
            """ There is a hinge in the last end """
            self.fem_elements[-1].releases = [5]
    
    
        
    
class MultiSpanTopChord(MultiSpanSteelMember):
    """ Class for multi-span top chord """
    
    def __init__(self,nodes,section,mem_id="",nel=2,hinges=None):
        
        mtype='top_chord'
        
        if hinges == None:
            hinges = [False,False]
        
        super().__init__(nodes,section,mtype,mem_id,nel,hinges)

class BottomChord(SteelFrameMember):

    def __init__(self,nodes,section,mem_id="",nel=4,hinges=None,lcr=[0.9,1.0]):
        
        mem_type = 'bottom_chord'
        
        if hinges == None:
            hinges = [False,False]
        
        SteelFrameMember.__init__(self,nodes,section,mem_type,mem_id,nel,hinges,lcr) 
        
        self.steel_members = []

        self.joints = []
    
    @property
    def perpendicular(self):
        """
        Returns vector perpendicular to member.
        NOTE! This needs to be modified for 3D!        
        :return:
        """
        unit = self.dir_vector    

        return np.array([unit[1], -unit[0]])
    
    def cost(self):
        """ Cost of chord member. This includes:
            i) material cost
            ii) blasting cost
            iii) painting cost

        """
        Ctot = {'material':self.material_cost(),
                'blasting':self.blasting_cost(),
                'painting':self.painting_cost()}
        
        self.costs = Ctot
        
        return Ctot

    #@property
    #def perpendicular(self):
    #    return -1 * super().perpendicular

    def generate_elements(self,fem,model="no_eccentricity"):
        """ Generate finite elements """
        
        if model == "no_eccentricity":
            super().generate_elements(fem)
        elif model == "en1993":
            # First end node
            self.fem_nodes.append(self.joints[0].fem_nodes['xc'])
            
            # Generate internal nodes:
            # global_node_coord include also the end node coordinates,
            # so they are not used in the iteration
            for x in self.global_node_coords[1:-1]:                
                if self.frame.dim == 2:
                    newNode = fem.add_node(x[0],x[1])    
                else:
                    newNode = fem.add_node(x[0],x[1],x[2])
                    
                self.fem_nodes.append(newNode)
            
            # Last end node
            self.fem_nodes.append(self.joints[1].fem_nodes['xc'])
            
            # Create members
            # The zip command allows to create two iterables, where the
            # items are actually consecutive items of the list fem_nodes.
            for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):
                if self.frame.dim == 2:                   
                    newElement = EBBeam(n1,n2,self.cross_section,self.material)                    
                else:
                    newElement = EBBeam3D(n1,n2,self.cross_section,self.material)
                    
                fem.add_element(newElement)
                self.fem_elements.append(newElement)
        elif model == "ecc_elements":
            if isinstance(self.joints[0],TubularYJoint):
                self.fem_nodes.append(self.joints[0].fem_nodes['xc'])                
            elif isinstance(self.joints[0],TubularKGapJoint):
                # in case of K gap joint, the node is created at the
                # intersection between brace center line and chord face            
                self.fem_nodes.append(self.joints[0].fem_nodes['xch'][self])
            elif isinstance(self.joints[0],TubularKTGapJoint):
                # in case of K gap joint, the node is created at the
                # intersection between brace center line and chord face            
                self.fem_nodes.append(self.joints[0].fem_nodes['xc'])
            
            # Generate internal nodes:
            # global_node_coord include also the end node coordinates,
            # so they are not used in the iteration
            for x in self.global_node_coords[1:-1]:                
                if self.frame.dim == 2:
                    newNode = fem.add_node(x[0],x[1])    
                else:
                    newNode = fem.add_node(x[0],x[1],x[2])
                    
                self.fem_nodes.append(newNode)
        
            if isinstance(self.joints[1],TubularYJoint):
                self.fem_nodes.append(self.joints[1].fem_nodes['xc'])                
            elif isinstance(self.joints[1],TubularKGapJoint):
                # in case of K gap joint, the node is created at the
                # intersection between brace center line and chord face            
                self.fem_nodes.append(self.joints[1].fem_nodes['xch'][self])
            elif isinstance(self.joints[1],TubularKTGapJoint):
                # in case of K gap joint, the node is created at the
                # intersection between brace center line and chord face            
                self.fem_nodes.append(self.joints[1].fem_nodes['xc'])
            
            for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):
                if self.frame.dim == 2:                   
                    newElement = EBBeam(n1,n2,self.cross_section,self.material)                    
                else:
                    newElement = EBBeam3D(n1,n2,self.cross_section,self.material)
                    
                fem.add_element(newElement)
                self.fem_elements.append(newElement)
        
        if self.hinges[0]:
            """ There is a hinge in the first end """
            self.fem_elements[0].releases = [2]
        
        if self.hinges[1]:
            """ There is a hinge in the last end """
            self.fem_elements[-1].releases = [5]
            

class TrussBrace(SteelFrameMember):
    """ Class for truss web members, or braces """
    
    def __init__(self,nodes,section,mem_id="",nel=1,hinges=[False,False],lcr=[1.0,1.0]):
        """ Constructor
        
        """
        mem_type='brace'
                
        SteelFrameMember.__init__(self,nodes,section,mem_type,mem_id,nel,hinges,lcr)
        
        self.top_joint = None
        self.bottom_joint = None
        
        self.top_chord_face_point = None
        self.bottom_chord_face_point = None
        
        # See if the member is vertical
        # This used with joints.
        if abs(nodes[0].x-nodes[1].x) < 1e-6:
            self.vertical = True
        else:
            self.vertical = False
    
    def cost(self,verb=False):
        """ Cost of brace. This includes:
            i) material cost
            ii) blasting cost
            iii) painting cost
            iv) sawing cost
            v) welding cost
        """
        Ctot = {'material':self.material_cost(),
                'blasting':self.blasting_cost(),
                'painting':self.painting_cost(),
                'sawing':self.sawing_cost()}
        
        self.costs = Ctot
        
        if verb:
            print(f"Brace {self.mem_id} cost:")
            print(f"Material cost: {Ctot['material']:.2f}")
            print(f"Blasting cost: {Ctot['blasting']:.2f}")
            print(f"Painting cost: {Ctot['painting']:.2f}")
            print(f"Sawing cost: {Ctot['sawing']:.2f}")
        
        return Ctot
    
    def sawing_cost(self):
        """ Cost of sawing """
        
        h = self.cross_section.B
        steel_grade = self.material.__repr__()
        bevels = 0
        
        T = [self.cross_section.T,self.cross_section.T]
        Af = 2*self.cross_section.H*self.cross_section.T
        AH = [Af,Af]
        AT = [self.cross_section.A,self.cross_section.A]
        
        L = self.length()
        
        return self.frame.workshop.cost_centres['sawing'].cost(L,bevels, steel_grade, [h,h], T, AH, AT)
    
    
    def robot_releases(self):
        """ Returns a string indicating releases that Robot Structural Analysis
            can use in creating a str file.
            
            For truss braces, create releases to both ends
        """
        return f'ELEments {self.mem_id+1} ORIgin RY END RY\n'
    
    def generate_elements(self,fem,model="no_eccentricity",beams=False):
        """ Generate finite elements """
        
        if model == "no_eccentricity":
            super().generate_elements(fem)
        elif model == "en1993":
            # Attach FEM nodes to the current member
            if isinstance(self.top_joint,TubularYJoint):
                n1 = self.top_joint.fem_nodes['xc']
                #self.fem_nodes.append(self.top_joint.fem_nodes['xc'])                
            elif isinstance(self.top_joint,TubularKGapJoint):
                #self.fem_nodes.append(self.top_joint.fem_nodes['xp'])
                n1 = self.top_joint.fem_nodes['xp']
            elif isinstance(self.top_joint,TubularKTGapJoint):
                #self.fem_nodes.append(self.top_joint.fem_nodes['xc'])
                n1 = self.top_joint.fem_nodes['xc']
            
            self.fem_nodes.append(n1)
            
            if isinstance(self.bottom_joint,TubularYJoint):
                n2 = self.bottom_joint.fem_nodes['xc']
                #self.fem_nodes.append(self.bottom_joint.fem_nodes['xc'])                
            elif isinstance(self.bottom_joint,TubularKGapJoint):
                n2 = self.bottom_joint.fem_nodes['xp']
                #self.fem_nodes.append(self.bottom_joint.fem_nodes['xp'])
            elif isinstance(self.bottom_joint,TubularKTGapJoint):
                n2 = self.bottom_joint.fem_nodes['xc']
                #self.fem_nodes.append(self.bottom_joint.fem_nodes['xc'])
            
            if beams:
                # Generate internal nodes:
                # global_node_coord include also the end node coordinates,
                # so they are not used in the iteration
                if self.frame.dim == 2:
                    r1 = np.array([n1.x,n1.y])
                    r2 = np.array([n2.x,n2.y])
                else:
                    r1 = np.array([n1.x,n1.y,n1.z])
                    r2 = np.array([n2.x,n2.y,n2.z])
                #r0 = self.top_chord_face_point
                #r1 = self.bottom_chord_face_point
                X = [r1 + t*(r2-r1) for t in self.local_node_coords[1:-1]]
                
                for x in X:
                    if self.frame.dim == 2:
                        newNode = fem.add_node(x[0],x[1])    
                    else:
                        newNode = fem.add_node(x[0],x[1],x[2])
                        
                    self.fem_nodes.append(newNode)
            
            self.fem_nodes.append(n2)
            
            """
            if beams:
                # Generate internal nodes:
                # global_node_coord include also the end node coordinates,
                # so they are not used in the iteration
                for x in self.global_node_coords[1:-1]:                
                    if self.frame.dim == 2:
                        newNode = fem.add_node(x[0],x[1])    
                    else:
                        newNode = fem.add_node(x[0],x[1],x[2])
                        
                    self.fem_nodes.append(newNode)
            
            
            if isinstance(self.bottom_joint,TubularYJoint):
                self.fem_nodes.append(self.bottom_joint.fem_nodes['xc'])                
            elif isinstance(self.bottom_joint,TubularKGapJoint):
                self.fem_nodes.append(self.bottom_joint.fem_nodes['xp'])
            elif isinstance(self.bottom_joint,TubularKTGapJoint):
                self.fem_nodes.append(self.bottom_joint.fem_nodes['xc'])
            """
            # Create elements
            # It can be desired, that also braces are modelled with beam elements.
            # In this case, 'beams' is True.
            if beams:
                for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):
                    if self.frame.dim == 2:        
                        newElement = EBBeam(n1,n2,self.cross_section,self.material)
                    else:
                        newElement = EBBeam3D(n1,n2,self.cross_section,self.material)
                    
                    fem.add_element(newElement)
                    self.fem_elements.append(newElement)
            else:
                n1 = self.fem_nodes[0]
                n2 = self.fem_nodes[1]
                newElement = Rod(n1,n2,self.cross_section,self.material)
            
                fem.add_element(newElement)
                self.fem_elements.append(newElement)
            
            """
            if beams:
                for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):
                    if self.frame.dim == 2:        
                        newElement = EBBeam(n1,n2,self.cross_section,self.material)
                    else:
                        newElement = EBBeam3D(n1,n2,self.cross_section,self.material)
                    
                    fem.add_element(newElement)
                    self.fem_elements.append(newElement)
            else:
                n1 = self.fem_nodes[0]
                n2 = self.fem_nodes[1]
                newElement = Rod(n1,n2,self.cross_section,self.material)
            
                fem.add_element(newElement)
                self.fem_elements.append(newElement)
            """
        elif model == "ecc_elements":
            # Attach FEM nodes to the current member
            
            # First node is the node connected to the top chord. This depends
            # on the joint type.
            if isinstance(self.top_joint,TubularYJoint):
                n1 = self.top_joint.fem_nodes['xc']                
            elif isinstance(self.top_joint,TubularKGapJoint):
                # in case of K gap joint, the node is created at the
                # intersection between brace center line and chord face 
                n1 = self.top_joint.fem_nodes['xcf'][self]                
            elif isinstance(self.top_joint,TubularKTGapJoint):
                n1 = self.top_joint.fem_nodes['xc']
                            
            self.fem_nodes.append(n1)
            
            # Last node is the node connected to the bottom chord.
            # This will be appended after the internal nodes in case
            # of using beam elements.
            if isinstance(self.bottom_joint,TubularYJoint):
                n2 = self.bottom_joint.fem_nodes['xc']                
            elif isinstance(self.bottom_joint,TubularKGapJoint):
                n2 = self.bottom_joint.fem_nodes['xcf'][self]                
            elif isinstance(self.bottom_joint,TubularKTGapJoint):
                n2 = self.bottom_joint.fem_nodes['xc']                
            
            if beams:
                # Generate internal nodes:
                # global_node_coord include also the end node coordinates,
                # so they are not used in the iteration
                if self.frame.dim == 2:
                    r1 = np.array([n1.x,n1.y])
                    r2 = np.array([n2.x,n2.y])
                else:
                    r1 = np.array([n1.x,n1.y,n1.z])
                    r2 = np.array([n2.x,n2.y,n2.z])
                #r0 = self.top_chord_face_point
                #r1 = self.bottom_chord_face_point
                X = [r1 + t*(r2-r1) for t in self.local_node_coords[1:-1]]
                
                for x in X:
                    if self.frame.dim == 2:
                        newNode = fem.add_node(x[0],x[1])    
                    else:
                        newNode = fem.add_node(x[0],x[1],x[2])
                        
                    self.fem_nodes.append(newNode)
            
            self.fem_nodes.append(n2)
            
            """
            if isinstance(self.bottom_joint,TubularYJoint):
                self.fem_nodes.append(self.bottom_joint.fem_nodes['xc'])                
            elif isinstance(self.bottom_joint,TubularKGapJoint):
                #self.fem_nodes.append(self.bottom_joint.fem_nodes['xp'])
                self.fem_nodes.append(self.bottom_joint.fem_nodes['xcf'][self])
            elif isinstance(self.bottom_joint,TubularKTGapJoint):
                self.fem_nodes.append(self.bottom_joint.fem_nodes['xc'])
            """
            # Create element
            
            if beams:
                for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):
                    if self.frame.dim == 2:        
                        newElement = EBBeam(n1,n2,self.cross_section,self.material)
                    else:
                        newElement = EBBeam3D(n1,n2,self.cross_section,self.material)
                    
                    fem.add_element(newElement)
                    self.fem_elements.append(newElement)
            else:
                n1 = self.fem_nodes[0]
                n2 = self.fem_nodes[1]
                newElement = Rod(n1,n2,self.cross_section,self.material)
            
                fem.add_element(newElement)
                self.fem_elements.append(newElement)
        
        if self.hinges[0]:
            """ There is a hinge in the first end """
            self.fem_elements[0].releases = [2]
        
        if self.hinges[1]:
            """ There is a hinge in the last end """
            self.fem_elements[-1].releases = [5]
    
    def angle_top(self):
        """ Angle between the brace and top chord """
        if not self.top_joint is None:
            # If the top chord joint has been defined,
            # use it to calculate the angle
            angle = self.top_joint.brace_angle(self)
        else:
            # If the top chord joint has not been defined,
            # use the centerlines of the members
            if not self.frame is None:
                # Take top chord direction from the parent frame
                # It is assume that the first top chord member
                # is on the left-hand side of the ridge and the last
                # top chord member is on the right-hand side.
                if self.node.x < self.frame.L1:                    
                    top_dir = self.frame_top_chord[0].dir_vector
                else:
                    top_dir = self.frame_top_chord[-1].dir_vector
                
                angle = np.degrees(np.arccos(np.dot(top_dir,self.dir_vector)))
                if angle > 90:
                    angle = 180-angle
            else:
                raise ValueError("Brace does not know about top chord.")
        
        return angle
    
    def angle_bottom(self):
        """ Angle between the brace and bottom chord """
        if not self.bottom_joint is None:
            # If the bottom chord joint has been defined,
            # use it to calculate the angle
            angle = self.bottom_joint.brace_angle(self)
        else:
            # If the bottom chord joint has not been defined,
            # use the centerlines of the members.
            # It is assumed that the bottom chord has always the same direction.
            bottom_dir = self.frame.bottom_chord[0].dir_vector
                
            angle = np.degrees(np.arccos(np.dot(bottom_dir,self.dir_vector)))
            if angle > 90:
                angle = 180-angle
            
        return angle
    
    def coords_actual_geometry(self):
        """ Returns the coordinates according to the actual truss geometry,
            including eccentricities, where relevant.
        """
        if isinstance(self.top_joint,TubularKGapJoint):
            Xt = self.top_joint.xp.coords
        elif isinstance(self.top_joint,TubularYJoint):
            Xt = self.top_joint.node.coords
        else:
            Xt = self.top_joint.node.coords
        
        if isinstance(self.bottom_joint,TubularKGapJoint):
            Xb = self.bottom_joint.xp.coords
        elif isinstance(self.bottom_joint,TubularYJoint):
            Xb = self.bottom_joint.node.coords
        else:
            Xb = self.bottom_joint.node.coords
            
        return Xt, Xb
    
    def dir_vector_xp(self):
        """ Direction vector of the brace using the eccentricity nodes
            where relevant.
        """
        Xt, Xb = self.coords_actual_geometry()
        
        if np.linalg.norm(Xt-Xb) < 1e-6:
            print(self)
            print(Xt,Xb)
        
        return (Xt-Xb)/np.linalg.norm(Xt-Xb)

    def plot(self, print_text=True, c='k', axes=None, mem_dim=False, geometry=False):
        
        if not geometry:
            super().plot(print_text,c,axes,mem_dim)
        else:
            if self.active:
                if axes is None:
                    fig, ax = plt.subplots(1)
                else:
                    ax = axes
        
                if geometry:
                    X0, X1 = self.coords_actual_geometry()
                else:
                    X0, X1 = self.coords()
                
                if c:
                    if self.resistance_check:
                        color = 'green'
                    else:
                        color = 'red'
                else:
                    color = 'k'
                
                # Plot members
                if mem_dim:
                    ax.plot([X0[0], X1[0]], [X0[1], X1[1]], color, linestyle='dashdot')
                    v = self.perpendicular
                    h = 0.5*self.cross_section.H*v
                    ax.plot([X0[0]+h[0], X1[0]+h[0]], [X0[1]+h[1], X1[1]+h[1]], color, linestyle='solid')
                    ax.plot([X0[0]-h[0], X1[0]-h[0]], [X0[1]-h[1], X1[1]-h[1]], color, linestyle='solid')
                else:
                    ax.plot([X0[0], X1[0]], [X0[1], X1[1]], color)
                                
                horzalign = 'center'
                vertalign = 'center'
                
                x, y = self.local_to_global(0.5) - self.perpendicular * 50
                rot = np.degrees(self.angle)
        
                if rot > 90:
                    rot = -(180-rot)
        
                if print_text:
                    ax.text(x, y, str(self.mem_id) + ": " + str(self.cross_section),
                             rotation=rot, horizontalalignment=horzalign,
                             verticalalignment=vertalign)
                else:
                    ax.text(x, y, str(self.mem_id),
                             rotation=rot, horizontalalignment=horzalign,
                             verticalalignment=vertalign)

class TubularJoint:
    """ Class for welded tubular joints """
    
    def __init__(self,node,chords,braces):
        """
        Constructor

        Parameters
        ----------
        node : FrameNode
            Node at which the joint is created.
        chords : List of chord members
            chord members appearing at joint. One or two members
        braces : List of TrussBrace objects
            Braces meeting at the joint. from one to three members are allowed

        Returns
        -------
        None.

        """
                
        self.joint_id = None
        self.node = node
        self.braces = braces
        self.chords = chords
    
        self.joint = None
        self.sym_joint = None
        
        self.fem_nodes = {}
        
        self.costs = {'welding':0.0}
        
        self.utilization = {}
        
        if isinstance(chords[0],TopChord):
            self.chord = "top"
        elif isinstance(chords[0],BottomChord):
            self.chord = "bottom"
        else:
            self.chord = None
        
    @property
    def h0(self):
        """ Height of the chord """
        return self.chords[0].cross_section.H
    
    @property
    def b0(self):
        """ Height of the chord """
        return self.chords[0].cross_section.B
    
    @property
    def h(self):
        """ Height of the braces """
        return np.array([brace.cross_section.H for brace in self.braces])
    
    @property
    def t(self):
        """ Wall thickness of the braces """
        return np.array([brace.cross_section.T for brace in self.braces])
    
    def weld_size(self,brace):
        """ Size of full strength fillet weld """
        a = full_strength_weld_tube(brace.material)
        return a*brace.cross_section.T
    
    def brace_angle(self,brace=None,geometry="init"):
        """ Angle between braces and chord """
        
        if brace is None:
            brace = self.braces[0]
        
        if geometry == "init":
            # Angle is calculated according to the initial geometry
            angle = np.degrees(np.arccos(np.dot(self.chords[0].dir_vector,brace.dir_vector)))
        elif geometry == "actual":
            # Angle is calculated according to the actual geometry, taking into account
            # eccentricities when relevant.
            
            # Chord direction vector
            vc = self.chords[0].dir_vector
            
            # Brace direction vector 
            vb = brace.dir_vector_xp()
            
            angle = np.degrees(np.arccos(np.dot(vc,vb)))
        
        # Take only sharp angles
        if angle > 90:
            angle = 180-angle
        
        return angle    
    
    def cost(self):
        """ Cost of welding """
        C = 0.0
        self.costs['welding'] = 0.0
        for brace in self.braces:            
            Hproj = brace.cross_section.H/np.sin(np.radians(self.brace_angle(brace)))            
            B = brace.cross_section.B
            R = brace.cross_section.R
            weld_length = 2*(Hproj+B) - 8*R + 2*np.pi*R
            weld_size = self.weld_size(brace)
            C += brace.frame.workshop.cost_centres['assembly_welding'].cost(weld_size,weld_length,weld_type='fillet')
        
        self.costs['welding'] = C
        return C
    
    def symmetry_joint(self,joints):
        """ Determine the joint lying symmetrically with respect to the 'node'.
            Requires symmetry_node attribute of the 'node' to have been set.
        """
        
        sym_joint = None
        
        if not self.node.symmetry_node is None:                  
            for joint in joints.values():
                if joint != self and joint.node == self.node.symmetry_node:
                    sym_joint = joint
                    self.sym_joint = sym_joint
                    break
                        
        return sym_joint
            
            
    
    def design(self,load_id):
        return 0.0
    
    def clear_fem(self):
        pass
        
        
class TubularYJoint(TubularJoint):
    """ Class for Y joints """
    
    def __init__(self,node,chords,braces,gap=0):
        """
        Constructor

        Parameters
        ----------
        node : FrameNode
            Node at which the joint is created.
        chords : List of chord members
            chord members appearing at joint. One or two members
        braces : List of TrussBrace objects
            Braces meeting at the joint. from one to three members are allowed

        Returns
        -------
        None.

        """
        super().__init__(node,chords,braces)
        
        # Initial angle is set according to the initial geometry
        angle = self.brace_angle(braces[0])
        
        # Could the chord and braces be members instead of cross-sections?
        self.joint = RHSYJoint(chords[0].cross_section,braces[0].cross_section,
                                  angle)
        
        self.fem_nodes['xc'] = None
    
    def __repr__(self):
        
        s = f"Y joint: Node {self.node.node_id}. {self.chords[0].cross_section.__repr__()}/{self.braces[0].cross_section.__repr__()}"
        
        return s
    
    def design(self,load_id):
        """ Design the joint in load case 'load_id' """
        
        """ The the load in the brace """
        self.joint.brace.Ned = self.braces[0].NEd[load_id]
        
        """ Load in the chord """
        N0 = 0.0
        for chord in self.chords:
            # Here, the FEM node can be something else besides 'xc',
            # depending on the FE model.
            n = chord.fem_nodes.index(self.fem_nodes['xc'])
            N = chord.nodal_forces[load_id]['N'][n]
            if abs(N) > abs(N0):
                N0 = N
        
        self.joint.N0 = N0

        r = self.joint.design()
        
        self.utilization[load_id] = r
        
        return r

    def clear_fem(self):
        
        self.fem_nodes['xc'] = None

class TubularKGapJoint(TubularJoint):
    """ Class for K gap joints """
    
    def __init__(self,node,chords,braces,gap=0):
        """
        Constructor

        Parameters
        ----------
        node : FrameNode
            Node at which the joint is created.
        chords : List of chord members
            chord members appearing at joint. One or two members
        braces : List of TrussBrace objects
            Braces meeting at the joint. Two braces are required.
            

        Returns
        -------
        None.
        
        The tricky bit is handling the gap and eccentricity, because the
        geometry of the joint depends on the joints on the opposite chord
        as well.
        
        One option is to consider the location of 'node' such that the
        gap is divided equally on each side. However, this might cause some difficulties
        when spread across the truss.
        
        Perhaps the easiest approach to implement is to take the
        intersection point of the brace centerlines as the variable and calculate
        everything from them. For a joint, the intersection points define eccentricity
        and the brace angles, from which the gap can be calculated. In order
        to find a suitable layout that satisifes gap and eccentricity conditions,
        the intersection points need to be altered, perhaps through an optimization
        problem. 
        

        """
        super().__init__(node,chords,braces)
        
        self.ridge = False
        if len(self.chords) > 1 and abs(abs(self.chords[0].dir_vector.dot(self.chords[1].dir_vector)) -1) > 1e-8:
            self.ridge = True
        
        # yloc and xloc are the local coordinates of the point of
        # intersection of the braces. 
        # Local coordinate axis of the joint is such that origin
        # is at 'node',
        # xloc is the coordinate along the chord and
        # yloc is the coordinate perpendicular to the chord. 'yloc' is
        # positive when moving away from the joint side.
        self.__yloc = 0
        self.__xloc = 0
        
        # How to determine gap and angles?
        if gap == 0:
            gap = 20
        
        # Take as initial angles the centerline angles
        angles = [self.brace_angle(brace) for brace in braces]
           
        # Could the chord and braces be members instead of cross-sections?
        self.joint = RHSKGapJoint(chords[0].cross_section,
                                  [braces[0].cross_section,braces[1].cross_section],
                                  angles,gap)
        
        min_gap = max(sum(self.joint.t),0.5*(1-self.joint.beta())*self.joint.b0)
        
        if self.joint.gap < min_gap:            
            self.joint.gap = min_gap
        
        # Determine which brace and chord member are on the left (right) side
        # of the joint node.
        if self.braces[0].opposite_node(self.node).x < self.braces[1].opposite_node(self.node).x:
            self.left_brace = self.braces[0]
            self.right_brace = self.braces[1]        
        else:
            self.left_brace = self.braces[1]
            self.right_brace = self.braces[0]
        
        if len(self.chords) == 1:
            # This can happen for bottom chord joints where the chord
            # end at the joint.
            if self.chords[0].opposite_node(self.node).x < self.node.x:
                self.left_chord = self.chords[0]
                self.right_chord = None
            else:
                self.left_chord = None
                self.right_chord = self.chords[0]
        else:        
            if self.chords[0].opposite_node(self.node).x < self.chords[1].opposite_node(self.node).x:
                self.left_chord = self.chords[0]
                self.right_chord = self.chords[1]        
            else:
                self.left_chord = self.chords[1]
                self.right_chord = self.chords[0]
        
        # Point of intersection of brace centerlines
        # This point is for the ideal model where the brace centerlines
        # meet at the centerline of the chord
        self.xp = FrameNode(node.coords)
        
        # Set initial eccentricity
        self.yloc = 0.25*self.chords[0].cross_section.H
        """
        if len(self.chords) == 2 and \
            np.linalg.norm(self.chords[0].dir_vector-self.chords[1].dir_vector) > 1e-6:
            # In this case, there are two chord members meeting at the joint and
            # the chord members are not parallel, i.e. the joint is located at
            # the ridge joint. Then, the point of intersection is
            xp = node.coords + np.array([0,self.yloc])
            #xp[1] += 0.25*self.chords[0].cross_section.H
            #print(xp)
        else:
            xp = node.coords - self.yloc*self.chords[0].perpendicular
        
        self.xp = FrameNode(xp)
        """
        
        # FEM nodes for structural analysis
        # Other fem nodes may need to be created.
        self.fem_nodes['xp'] = None
        self.fem_nodes['xc'] = None
        
        # FEM nodes for the model, where nodes are
        # added at chord face and on the chord center line
        # perpendicular to the chord face nodes
        self.fem_nodes['xcf'] = {braces[0]: None, braces[1]: None}
        #self.fem_nodes['xcf2'] = None
        
        if len(chords) > 1:
            self.fem_nodes['xch'] = {chords[0]: None, chords[1]: None}
        else:
            self.fem_nodes['xch'] = {chords[0]: None}
        #self.fem_nodes['xch2'] = None
        # Line passing in the direction of chord through its centerline
        #self.chord_line = Line(v=self.chords[0].dir_vector,p1=self.node.coords)
        
        
        self.fem_elements = {'ecc':None, 'left_ecc':None, 'right_ecc': None, 'gap': None}
    
    def __repr__(self):
        
        s = f"K gap joint: Node {self.node.node_id}."
        s += f" {self.chords[0].cross_section.__repr__()}/{self.braces[0].cross_section.__repr__()}/{self.braces[1].cross_section.__repr__()}"
        
        return s
    
    @property
    def xp_coord(self):
        """ Coordinates of the intersection point of brace centerlines """
        return self.xp.coords
    
    @property
    def gap(self):
        """ Returns the gap """
        return self.joint.gap
    
    @property
    def angles(self):
        """ Returns brace angles """
        return self.joint.angles
    
    @property
    def yloc(self):
        """ Location of the intersection point in the local
            coordinate system of the joint, y coordinate. This is equal
            to eccentricity of the joint
        """
        return self.__yloc
    
    @yloc.setter
    def yloc(self,val):
        """ Set the local y coordinate of the intersection point of the
            brace center lines.
        """
        self.__yloc = val
        
        #if len(self.chords) == 2 and \
        #    np.linalg.norm(self.chords[0].dir_vector-self.chords[1].dir_vector) > 1e-6:
        if self.ridge:
            # In this case, the joint is at the apex of the truss.
            # In this case, there are two chord members meeting at the joint and
            # the chord members are not parallel, i.e. the joint is located at
            # the ridge joint. Then, the point of intersection is
            self.xp.coords = self.node.coords + np.array([0,self.__yloc])            
        else:            
            self.xp.coords = self.node.coords + self.__xloc*self.chords[0].dir_vector + self.__yloc*self.chords[0].perpendicular
            if not self.sym_joint is None:
                self.sym_joint.__yloc = val
                self.sym_joint.xp.coords = self.sym_joint.node.coords + self.sym_joint.__xloc*self.sym_joint.chords[0].dir_vector + self.sym_joint.__yloc*self.sym_joint.chords[0].perpendicular
        
    @property
    def xloc(self):
        """ Location of the intersection point in the local
            coordinate system of the joint, x coordinate. This is the distance
            of the intersection point along the chord line.
        """
        return self.__xloc
    
    @xloc.setter
    def xloc(self,val):
        """ Set the local x coordinate of the intersection point of the
            brace center lines.
        """
        self.__xloc = val
        
        self.xp.coords = self.node.coords + self.__xloc*self.chords[0].dir_vector + self.__yloc*self.chords[0].perpendicular
    
        if not self.sym_joint is None:
            self.sym_joint.__xloc = -val
            self.sym_joint.xp.coords = self.sym_joint.node.coords + self.sym_joint.__xloc*self.sym_joint.chords[0].dir_vector + self.sym_joint.__yloc*self.sym_joint.chords[0].perpendicular
            
    
    
    def xc_face(self):
        """ Coordinates of the chord face point closest to the coordinates
            of 'node' on the side of the joint
        """
        if self.ridge:
            xc = self.node.coords - np.array([0,0.5*self.joint.h0/np.cos(self.chords[0].angle)])
        else:
            xc = self.node.coords - 0.5*self.joint.h0*self.chords[0].perpendicular
        
        return xc
    
    def brace_chord_face_x(self):
        """ Determine the points of interesection between braces and chord face 
            Additionally, determine the points on the chord center line closest to
            the surface points.
        """
        
        xint = {}
        xchord = {}
        
        if self.ridge:
            # Direction vector of chord on the left hand side of the ridge
            vleft = self.left_chord.dir_vector
            xcf = self.xc_face()
            
            xp = self.xp_coord
            
            chord_face_line = Line(v=vleft,p1=xcf)
            chord_center_line = Line(v=vleft,p1=self.node.coords)
            
            b_v = self.left_brace.dir_vector_xp()     
            b_line = Line(v=b_v,p1=xp)
            tb, xb = b_line.intersect(chord_face_line)
            xint[self.left_brace] = xb
            
            self.left_brace.top_chord_face_point = xb
            
            d, xc = chord_center_line.distance(xb)
            xchord[self.left_chord] = xc
            
            # Right side brace:
            vright = self.right_chord.dir_vector
            
            chord_face_line = Line(v=vright,p1=xcf)
            chord_center_line = Line(v=vright,p1=self.node.coords)
            
            b_v = self.right_brace.dir_vector_xp()     
            b_line = Line(v=b_v,p1=xp)
            tb, xb = b_line.intersect(chord_face_line)
            xint[self.right_brace] = xb
            
            self.right_brace.top_chord_face_point = xb
            
            d, xc = chord_center_line.distance(xb)
            xchord[self.right_chord] = xc
            
        else:
            # Make a line in the direction of the chord 
            # passing through the chord face point
            v = self.chords[0].dir_vector
            xcf = self.xc_face()        
            chord_face_line = Line(v=v,p1=xcf)
            chord_center_line = Line(v=v,p1=self.node.coords)
            
            xp = self.xp_coord
            
            for brace in self.braces:
                b_v = brace.dir_vector_xp()        
                b_line = Line(v=b_v,p1=xp)
            
                # calculate intersection point between brace
                # center line and chord face line
                tb, xb = b_line.intersect(chord_face_line)
                #t.append(tb)
                #xint.append(xb)
                # keys of xint are the brace objects. These will be
                # used later in generating the fem model.
                xint[brace] = xb
                
                if self.chord == "top":
                    brace.top_chord_face_point = xb
                elif self.chord == "bottom":
                    brace.bottom_chord_face_point = xb
                
                d, xc = chord_center_line.distance(xb)
                #D.append(d)
                #xchord.append(xc)
                # TODO! PAARTEEN KESKILINJAN SOLMU PITÄISI LIITTÄÄ
                # KO. PAARRESAUVAAN, KOSKA SOLMUSTA TEHDÄÄN PAARRESAUVALLE
                # FEM-SOLMU.
                # MITEN KYTKETÄÄN LIITOKSEN UUMASAUVA brace OIKEAAN PAARRESAUVAAN?
                if brace == self.left_brace:
                    xchord[self.left_chord] = xc
                else:
                    xchord[self.right_chord] = xc
        #print(t,xint)
        
        return xint, xchord
        
    
    def xchord(self):
        """ Finds the point on the centerline of the chord closest to the
            intersection point of the braces
        """
        
        v = self.chords[0].dir_vector
        
        chord_line = Line(v=v,p1=self.node.coords)
        xp_dist, xc = chord_line.distance(self.xp_coord)
        
        return xc
    
    def eccentricity(self):
        """ Evaluates eccentricity based on the intersection point xp """
        
        """
        v = self.chords[0].dir_vector
        
        chord_line = Line(v=v,p1=self.node.coords)
        xp_dist, _ = chord_line.distance(self.xp_coord)
        
        return xp_dist
        """
        return self.yloc
    
    def calc_gap(self):
        """ Calculates the gap based on eccentricity and brace angles """
        e = self.eccentricity()
        theta = np.array([self.brace_angle(brace,geometry="actual") for brace in self.braces])
        
        # Set angles to the RHS joint object as well
        self.joint.angles = theta
        h0 = self.h0
        h = self.h
        
        sint = np.sin(np.radians(theta))
        
        gap = (e + 0.5*h0)*np.sin(np.radians(sum(theta)))/sint.prod()-0.5*sum(h/sint)
                        
        # Set gap for the RHS joint object
        self.joint.gap = gap
        #self.gap = gap
        
        return gap

    def design(self,load_id):
        """ Design the joint in load case 'load_id' """
        
        """ The load in the brace """
        self.joint.braces[0].Ned = self.braces[0].NEd[load_id]
        self.joint.braces[1].Ned = self.braces[1].NEd[load_id]
        #for joint_brace, brace in zip(self.joint.braces,self.braces):
            
        #    joint_brace.Ned = brace.NEd[load_id]
        
        """ Load in the chord """
        N0 = 0.0
        M0 = 0.0
        for chord in self.chords:
            # Here, the FEM node can be something else besides 'xc',
            # depending on the FE model.
            n = chord.fem_nodes.index(self.fem_nodes['xc'])
            N = chord.nodal_forces[load_id]['N'][n]
            M = chord.nodal_forces[load_id]['My'][n]
            if abs(N) > abs(N0):
                N0 = N
            
            if abs(M) > abs(M0):
                M0 = M
        
        self.joint.N0 = N0

        # The chord stress function in rhs_joints.py is constructed such
        # that positive bending moment causes compression stress in the
        # joint. Then, for top chord joints, by reversing the sign of the
        # bending moment, negative chord moment causes positive compressive
        # stress in the joint.        
        if isinstance(self.chords[0],TopChord):
            M0 = -M0
        self.joint.M0 = M0

        r = self.joint.design()
        
        self.utilization[load_id] = r
        
        return r

    def clear_fem(self):
        
        self.fem_nodes['xc'] = None
        self.fem_nodes['xp'] = None
        
        for fem_node in self.fem_nodes['xcf'].values():
            fem_node = None
        
        for fem_node in self.fem_nodes['xch'].values():
            fem_node = None

        for ele in self.fem_elements.values():
            ele = None

class TubularKTGapJoint(TubularJoint):
    """ Class for KT gap joints """
    
    def __init__(self,node,chords,braces,gaps=[0,0]):
        """
        
        One of the braces is 

        Parameters
        ----------
        node : TYPE
            DESCRIPTION.
        chords : TYPE
            DESCRIPTION.
        braces : TYPE
            DESCRIPTION.
        gaps : TYPE, optional
            DESCRIPTION. The default is [0,0].

        Returns
        -------
        None.

        """
        
        super().__init__(node,chords,braces)
        
        # yloc and xloc are the local coordinates of the point of
        # intersection of the braces. 
        # Local coordinate axis of the joint is such that origin
        # is at 'node',
        # xloc is the coordinate along the chord and
        # yloc is the coordinate perpendicular to the chord. 'yloc' is
        # positive when moving away from the joint side.
        self.__yloc = 0
        self.__xloc = 0
        
        # How to determine gap and angles?
        gaps = [max(20,gap) for gap in gaps]
        
        # Take as initial angles the centerline angles
        angles = [self.brace_angle(brace) for brace in braces]
           
        # Could the chord and braces be members instead of cross-sections?
        self.joint = RHSKTGapJoint(chords[0].cross_section,
                                  [brace.cross_section for brace in self.braces],
                                  angles,gaps)
        
        #min_gap = max(sum(self.joint.t),0.5*(1-self.joint.beta())*self.joint.b0)
        
        #if self.joint.gap < min_gap:
        #    self.joint.gap = min_gap
    
        self.fem_nodes['xp'] = None
        self.fem_nodes['xc'] = None

    def __repr__(self):
        
        s = f"KT gap joint: Node {self.node.node_id}."
        s += f" {self.chords[0].cross_section.__repr__()}/{self.braces[0].cross_section.__repr__()}/{self.braces[1].cross_section.__repr__()}/{self.braces[2].cross_section.__repr__()}"
        
        return s

    def clear_fem(self):
        
        self.fem_nodes['xc'] = None
        self.fem_nodes['xp'] = None
        
        """
        for fem_node in self.fem_nodes['xcf'].values():
            fem_node = None
        
        for fem_node in self.fem_nodes['xch'].values():
            fem_node = None

        for ele in self.fem_elements.values():
            ele = None
        """

def Ktruss_example(h2=2000,h1=1500,dx1=1000,dx2=1000,first=False,edges=True):
        
    t = SlopedTruss(h2=h2,h1=h1,dx1=dx1,dx2=dx2)
    t.generate_topology('K',4,first_diagonal_up=first,edge_verticals=edges)
    #t.plot_nodes()
    #t.plot()
    #print(t.top_nodes)
    #t.generate_supports()
    #t.generate_joints()
    #print(t.top_nodes)
    #t.generate_uniform_load(q=-25)
    #t.generate_fem(model='en1993')
    #t.generate_fem(model='no_eccentricity')
    #t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    #t.bmd(scale=20,load_id=t.load_ids[0],loads=False)
    #t.plot(geometry=True)
    
    return t
    #t.plot(mem_dim=True)

def Ntruss_example(h2=2000,h1=1500,dx1=1000,dx2=1000,first=False,edges=True):
    
    t = SlopedTruss(h2=h2,h1=h1,dx1=dx1,dx2=dx2)
    t.generate_topology('N',4,first_diagonal_up=first,edge_verticals=edges)
    
    t.plot()
    
    return t

def KTtruss_example(h2=2000,h1=1500,dx=1000,first=False,edges=True):
    
    t = SlopedTruss(h2=h2,h1=h1,dx1=dx,dx2=dx)
    t.generate_topology('KT',3,first_diagonal_up=first,edge_verticals=edges)
    
    #t.plot()
    
    return t

if __name__ == "__main__":
    
    
    from timeit import default_timer as timer
    
    #t = Ktruss_example(h1=1500,dx1=1000, dx2=1000,first=True, edges=False)
    #t = Ntruss_example(dx1=1000, dx2=1000, first=True,edges=False)
    #t = Ntruss_example(dx1=1000, dx2=1000, first=False,edges=False)
    #t = Ntruss_example(dx1=1000, dx2=1000, first=False,edges=True)
    #t = Ntruss_example(dx1=1000, dx2=1000, first=True,edges=True)
    #t = Ntruss_example()
    
    dx1 = 1500
    t = KTtruss_example(dx=dx1, first=True,edges=False)
    #t = KTtruss_example(dx=dx1, first=False,edges=False)
    #t = KTtruss_example(dx=dx1, first=False,edges=True)
    #t = KTtruss_example(dx=dx1, first=True,edges=True)
    
    t.symmetry()
    
    """
    #t = PlaneTruss()
    t = SlopedTruss(h2=2000,h1=1500,dx1=1000,dx2=1000)
    t.generate_topology('K',4)    
    #print(t.top_nodes)
    """
    t.generate_supports()
    t.generate_joints()
    #print(t.joints[1].gap)
    #print(t.top_nodes)
    t.generate_uniform_load(q=-25)
    t.generate_fem(model='en1993')
    #t.generate_fem(model='no_eccentricity')
    t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    
    #t.bmd(scale=20,load_id=t.load_ids[0],loads=False)
    #t.cost()
    #start = timer()
    #t.optimize_members(verb=True,top={'material':'S700','class':2},
    #                   bottom={'material':'S700','class':2},
    #                   braces={'material':'S700','class':2})
    #end = timer()
    #print(f"Time elapsed: {end-start:.2f} s")
    """
    t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    
    t.plot(geometry=True)
    #t.plot(mem_dim=True)
    """
    #opts = {'x_monitor':0.5*t.span, 'n_monitored':2}
    #t.to_abaqus(filename='KT-ristikko',partname="KT-ristikko",options=opts)

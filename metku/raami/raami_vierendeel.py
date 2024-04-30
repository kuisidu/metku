# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Tue Mar 26 12:20:46 2024

Class for Vierendeel trusses

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt
from copy import copy

from metku.raami.raami import Raami
from metku.raami.raami_plane_truss import PlaneTruss, TopChord, BottomChord, TrussBrace, TubularYJoint
from metku.raami.raami_functions import divide_line_segment
from metku.raami.frame_node import FrameNode
from metku.raami.frame_member import FrameMember, SteelFrameMember, MultiSpanSteelMember, MemberGroup
from metku.raami.frame_loads import PointLoad, LineLoad, PWLineLoad, LoadIDs
from metku.raami.frame_supports import XYHingedSupport, YHingedSupport
from metku.raami.exports import AbaqusOptions

from metku.framefem import LineLoad as FEMLineLoad
from metku.framefem.elements.ebbeam import EBBeam, EBBeam3D
from metku.framefem.elements.rod import Rod
from metku.sections.steel.RHS import RHS, SHS
from metku.sections.steel.ISection import IPE
from metku.eurocodes.en1993.en1993_1_8.rhs_joints import Line, RHSKGapJoint, RHSYJoint, RHSKTGapJoint
from metku.eurocodes.en1993.en1993_1_8.en1993_1_8 import full_strength_weld_tube

from metku.cost.workshop import Workshop


#from loadIDs import LoadIDs

def unit_vector(x1,x2):
    """ Creates a unit vector between points x1 and x2 """
    
    return (x2-x1)/np.linalg.norm(x2-x1)

class VierendeelTruss(PlaneTruss):

    def __init__(self, origin=(0, 0), span=10000, height=1000, panels=5,
                 num_elements=6,
                 name="Raami Vierendeel Truss",
                 workshop=Workshop()):
        """ Constructor 
            Parameters:
            ------------
            :param origin: location of the origin (tuple)
            :param span: span of the truss [mm]
            :param height: height of the truss [mm]
            :param panels: number of panels, i.e. division of span to vertical braces
            :param num_elements: number of elements per member
            
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
        PlaneTruss.__init__(self,origin=origin, span=span, left_height=height,
                     mid_height=height, right_height=height, 
                     bot_chord_dx_left = 0, bot_chord_dx_right = 0,
                     rotation = 0,
                     num_elements=num_elements,
                     name="Raami Vierendeel Truss")
        
        # Braces are modelled using beam elements
        self.braces_as_beams = True
        
        self._H = height
        self._npanels = panels                
    
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
    
    def generate_topology(self, 
                          top_chord_profile = RHS(160,160,5),
                          bottom_chord_profile = RHS(140,140,5),
                          brace_profile = RHS(90,90,5),
                          nel_chord=4,
                          nel_brace=4):
        """ Generates truss members for a given number of braces 
        
            :param: topology .. 'K', 'N', or 'KT'
            :param: ndiv .. number of divisions on one half of the truss
            
        """

        # Generate nodes first
            
        # Bottom chord nodes:            
        dx = self.span/self._npanels  
        ycoord = self.origin[1]
        
        S = np.linspace(0,1,self._npanels+1)

        X0 = np.array([self.origin[0],ycoord])
        X1 = np.array([self.origin[0]+self.span,ycoord])  
        Xbot = [FrameNode(X0+s*(X1-X0)) for s in S]
            
        
        # Top chord nodes:
        # Local coordinates of the top chord
                        
        X0 = np.array([self.origin[0],ycoord+self._H])
        X1 = np.array([self.origin[0]+self.span,ycoord+self._H])
          
        Xtop = [FrameNode(X0+s*(X1-X0)) for s in S]                      
        
        print(Xtop)
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
        for X, Y in zip(Xtop,Xbot):
            self.add(TrussBrace([X,Y],copy(brace_profile),nel=nel_brace))
        
        for brace in self.braces.values():
            brace.hinges = [False,False]            
        
        self.top_nodes = Xtop
        self.bottom_nodes = Xbot
            
    def generate_supports(self):
        """ Creates supports to the truss """
        
        self.add(XYHingedSupport(self.bottom_nodes[0]))        
        self.add(YHingedSupport(self.bottom_nodes[-1]))
    
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
            
            if nbraces == 1:                
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
            angle = joint.brace_angle(joint.braces[0],geometry="actual")
            joint.joint.angle = angle
    
    def generate_uniform_load(self,q,chord="top",load_id=LoadIDs['ULS'],ltype='live'):
        """ Creates a uniform line load to one of the chords """

        vals = [q,q]

        if chord == "top":
            load_name = 'TC_' + ltype + '_line'
            for mem in self.top_chord:
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, load_name, "global"))
        elif chord == "bottom":          
            load_name = 'BC_' + ltype + '_line'
            for mem in self.bottom_chord:
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, load_name, "global"))
    
    def generate_chord_point_load(self,p,chord="top",load_id=LoadIDs['ULS'],ltype='live'):
        """ Creates point loads to one of the chords """

        #(self, node, v, load_id=LoadIDs['ULS'], ltype='live',
        #name='PointLoad',f=1.0)

        if chord == "top":
            load_name = 'TC_' + ltype + '_point'
            self.add(PointLoad(self.top_nodes[0], [0.0,p[0],0.0], load_id, ltype, load_name + '_1',1.0))
            self.add(PointLoad(self.top_nodes[-1], [0.0,p[0],0.0], load_id,ltype, load_name + '_1',1.0))
            for node in self.top_nodes[1:-1]:                
                self.add(PointLoad(node, [0.0,p[1],0.0], load_id, ltype, load_name + '_2',1.0))
        elif chord == "bottom":          
            load_name = 'BC_' + ltype + '_point'
            # It is assumed that the first and last bottom chord nodes are vertically supported,
            # so no vertical load is inserted there
            for node in self.bottom_nodes[1:-1]:                
                self.add(PointLoad(node, [0.0,p[1],0.0], load_id, ltype, load_name + '_2', 1.0))
    
    def symmetry(self,axis='y',sym_value=None):
        
        # Make general actions for symmetry
        super().symmetry(axis,sym_value)
        
        # Symmetry of joints
        #joint_list = list(self.joints.values())
        for joint in self.joints.values():
            #print(f"Symmetry for joint {joint.joint_id}.")
            sym_node = joint.symmetry_joint(self.joints)
            #print(sym_node)
    
    def clear_fem(self):
        """ Clears FEM data """
        super().clear_fem()
        
        for joint in self.joints.values():
            joint.clear_fem()
            
    def generate_fem(self,model="no_eccentricity",nodes_and_elements_only=False):
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
                            
        elif model == "ecc_elements":
            # Eccentricity elements are added for K joints
            # from the intersection of brace centerline and chord
            # face to the chord centerline.
            for joint in self.joints.values():
                # For Y joints, no eccentricity elements are used.
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
            
       
            if not nodes_and_elements_only:
                # Generate supports
                for supp in self.supports.values():
                    supp.add_support(self.fem)
                                
                # Generate loads
                for load in self.loads.values():
                    # add_load addes the loads on the target member
                    load.add_load(self.fem)
                    # Creates new loadcase if one with same load_id does not exist
                    lcase_ids = [lc.load for lc in self.fem.loadcases.values()]
                    if load.load_id not in lcase_ids:
                        self.fem.add_loadcase(supp_id=1,load_id=load.load_id)
                
                # Set nodal degrees of freedom
                self.fem.nodal_dofs()
                
        else:
            print("Weird FEM model.")
    
    def optimize_members(self, prof_type="CURRENT",verb=False,**kwargs):    
        """ Finds minimum weight profiles for members """
               
        top = {'material': 'S355', 'class': 2, 'utility': 1.0, 'bmin':10, 'bmax':1e5}
        bottom = {'material': 'S355', 'class': 2, 'utility': 1.0,'bmin':10, 'bmax':1e5}
        braces = {'material': 'S355', 'class': 2, 'utility_tens': 1.0, 'utility_comp':1.0, 'bmin':10, 'bmax':1e5}
            
        limit_width = False
        
        for key, value in kwargs.items():
            if key == 'top':
                top = value
            elif key == 'bottom':
                bottom = value
            elif key == 'braces':
                braces = value
            elif key == 'limit_width':
                limit_width = value
        
        kmax = 10
        k = 1
                        
        while k < kmax:            
            explored = []
            MEM_PROFILES_CHANGED = []
            self.structural_analysis('all','REM')
            if verb:
                print(f"Optimize profiles, iteration {k}")
            
            if limit_width:
                b_chord_min = self.max_brace_width()
                #b_brace_max = min(self.top_chord_profile.b,self.bottom_chord_profile.b)
                b_brace_max = 1e5
                b_brace_min = 10
            else:
                b_chord_min = 10
                b_brace_max = 1e5
                b_brace_min = 10
            
            # Go through groups first
            for name, group in self.member_groups.items():
                if name == "top_chord":
                    material = top['material']
                    sec_class = top['class']
                    max_utility = top['utility']                    
                elif name == 'bottom_chord':
                    material = bottom['material']
                    sec_class = bottom['class']
                    max_utility = bottom['utility']                    
                    
                group.optimum_design(prof_type,verb,material,sec_class,max_utility,b_chord_min)
                
                for mem in group.members:
                    explored.append(mem)
                
            # Braces
            for member in self.members.values():
                #print(member)                
                if not member in explored:
                    if verb:
                        print(member)
                    explored.append(member)
                    if isinstance(member,TrussBrace):
                        material = braces['material']
                        sec_class = braces['class']
                        if list(member.NEd.values())[0] >= 0.0:
                            max_utility = braces['utility_tens']
                        else:
                            max_utility = braces['utility_comp']
                        
                        # Minimum brace width is taken from the geometry condition of T-joints:
                        # bi >= 0.25*b0
                        if limit_width:
                            b_brace_min = 0.25*max(member.top_joint.b0,member.bottom_joint.b0)
                            b_brace_max = min(member.top_joint.b0,member.bottom_joint.b0)                        
                        
                    else:
                        material = 'S355'
                        sec_class = 2
                        max_utility = 1.0
                        
                    MEM_PROFILES_CHANGED.append(member.optimum_design(prof_type,verb,material,sec_class,max_utility,bmin=b_brace_min,bmax=b_brace_max))
                                        
                    if not member.symmetry_pair is None:
                        explored.append(member.symmetry_pair)
                    
                    #print(member.cross_section)
            
            # If none of the members changed profile, exit the loop
            if not any(MEM_PROFILES_CHANGED):
                break
                
            k += 1
        
        if verb:
            print(f"Number of iterations: {k}.")
            
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
            
        
        PlaneTruss.__init__(self,origin=(0, 0), span=span, left_height=height,
                     mid_height=height, right_height=height, 
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
                          top_chord_profile = RHS(160,160,5),
                          bottom_chord_profile = RHS(140,140,5),
                          brace_profile = RHS(90,90,5),
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
            X0 = np.array([self.origin[0],self.origin[1]])
            X2 = np.array([self.origin[0]+self.span,self.origin[1]])
                                 
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
            X0 = np.array([self.origin[0],self.origin[1]])
            X2 = np.array([self.origin[0]+self.span,self.origin[1]])
                                 
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
            load_name = 'TC_' + ltype + '_line'
            for mem in self.top_chord:
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, load_name, "global"))
        elif chord == "bottom":          
            load_name = 'BC_' + ltype + '_line'
            for mem in self.bottom_chord:
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, load_name, "global"))
    
    def generate_piecewise_uniform_load(self,q,xzones,load_id=LoadIDs['ULS'],ltype='live'):
        """ Creates piecewise uniform line load which affects the top chord 
            :param q: list of load values
            :param xzones: list of x-coordinates of the wind zones
        """
        load_name0 = 'TC_' + ltype
        
        for mem in self.top_chord:
            if mem.nodes[1].x <= xzones[0]:
                # Member is in the first zone
                vals = [q[0],q[0]]
                load_name = load_name0 + '_line'
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, load_name, "global"))
            elif mem.nodes[0].x > xzones[-1]:
                # Member is in the last zone
                vals = [q[-1],q[-1]]
                load_name = load_name0 + '_line'
                self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, load_name, "global"))
            else:
                
                for i, (x1, x2) in enumerate(zip(xzones[:-1],xzones[1:])):
                    if mem.nodes[0].x < x1 and mem.nodes[1].x > x1:
                        # member is partially outside
                        vals = [q[i],q[i+1]]
                        xdiv = (x1-mem.nodes[0].x)/(mem.nodes[1].x-mem.nodes[0].x)
                        #xdiv = x1
                        load_name = load_name0 + '_pw_line'
                        self.add(PWLineLoad(mem, vals, xdiv, 'y', load_id, 1.0, ltype, load_name,"global"))
                    elif mem.nodes[0].x >= x1 and mem.nodes[1].x > x2:
                        # member is partially outside
                        vals = [q[i+1],q[i+2]]
                        xdiv = (x2-mem.nodes[0].x)/(mem.nodes[1].x-mem.nodes[0].x)
                        #xdiv = x2
                        load_name = load_name0 + '_pw_line'
                        self.add(PWLineLoad(mem, vals, xdiv, 'y', load_id, 1.0, ltype, load_name,"global"))
                    elif mem.nodes[0].x >= x1 and mem.nodes[1].x <= x2:
                        # member is totally inside
                        vals = [q[i+1],q[i+1]]
                        load_name = load_name0 + '_line'
                        self.add(LineLoad(mem, vals, 'y', load_id, 1.0, ltype, load_name,"global"))
    
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
    
    def symmetry(self,axis='y',sym_value=None):
        
        # Make general actions for symmetry
        super().symmetry(axis,sym_value)
        
        # Symmetry of joints
        #joint_list = list(self.joints.values())
        for joint in self.joints.values():
            #print(f"Symmetry for joint {joint.joint_id}.")
            sym_node = joint.symmetry_joint(self.joints)
            #print(sym_node)
    
    def clear_fem(self):
        """ Clears FEM data """
        super().clear_fem()
        
        for joint in self.joints.values():
            joint.clear_fem()
        
    def generate_fem(self,model="no_eccentricity",nodes_and_elements_only=False):
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
            
            if not nodes_and_elements_only:
                # Generate supports
                for supp in self.supports.values():
                    supp.add_support(self.fem)
                                
                # Generate loads
                for load in self.loads.values():
                    # add_load addes the loads on the target member
                    load.add_load(self.fem)
                    # Creates new loadcase if one with same load_id does not exist
                    lcase_ids = [lc.load for lc in self.fem.loadcases.values()]
                    if load.load_id not in lcase_ids:
                        self.fem.add_loadcase(supp_id=1,load_id=load.load_id)
                
                # Add loads to gap elements
                for joint in self.joints.values():
                    if isinstance(joint,TubularKGapJoint):
                        # Add load to gap element
                        for load1 in joint.chords[0].loads:
                            for load2 in joint.chords[1].loads:
                                if load1.load_id == load2.load_id:
                                    if joint == self.ridge_joint:
                                        q = load1.f*load1.values[1]
                                        load = FEMLineLoad(load1.load_id,
                                                           joint.fem_elements['gap'][0],
                                                           [0.0,1.0],
                                                           [q,q],
                                                           load1.direction,
                                                           load1.coord_sys)
                                        self.fem.add_load(load)
                                        q = load1.f*load2.values[0]
                                        load = FEMLineLoad(load1.load_id,
                                                           joint.fem_elements['gap'][1],
                                                           [0.0,1.0],
                                                           [q,q],
                                                           load1.direction,
                                                           load1.coord_sys)
                                        self.fem.add_load(load)
                                    else:
                                        if isinstance(load1,LineLoad) or isinstance(load1,PWLineLoad):
                                            if abs(load1.values[1]) >= abs(load2.values[0]):
                                                q = load1.f*load1.values[1]
                                                
                                            else:
                                                q = load2.f*load2.values[0]
                                                                                    
                                            load = FEMLineLoad(load1.load_id,
                                                               joint.fem_elements['gap'],
                                                               [0.0,1.0],
                                                               [q,q],
                                                               load1.direction,
                                                               load1.coord_sys)                                            
                                            self.fem.add_load(load)
                                    break
                                    
                
                # Set nodal degrees of freedom
                self.fem.nodal_dofs()
                
        else:
            print("Weird FEM model.")
    
    def update_fem_joint_nodes(self):
        """ Updates nodal positions at joints. It is assumed that cross-sections
            have been changed along the way.
        """
        
        for joint in self.joints.values():
            if isinstance(joint,TubularKGapJoint):
                # Calculate brace-to-chord face intersection points and the
                # corresponding projections to the centerline of the chord
                xface, xcenter = joint.brace_chord_face_x()
                for brace, xf in xface.items():                    
                    joint.fem_nodes['xcf'][brace].x = xf[0]
                    joint.fem_nodes['xcf'][brace].y = xf[1]
                
                for brace, xc in xcenter.items():                    
                    joint.fem_nodes['xch'][brace].x = xc[0]
                    joint.fem_nodes['xch'][brace].y = xc[1]
    
    def plot(self, geometry=False, print_text=True, show=True,
             loads=True, color=False, axes=None, save=False, mem_dim=False,
             saveopts={'filename':'default','format':'svg','orientation':'landscape', 'papertype':'a3'}):
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
            plt.rcParams['text.color'] = "red"
            plt.rcParams['font.size'] = 2.5
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
                #plt.savefig('default.svg', format='svg')
                plt.savefig(saveopts['filename'],format=saveopts['format'],\
                            orientation=saveopts['orientation'],papertype=saveopts['papertype'])
            if show:
                plt.axis('equal')
                plt.show()
                
    def optimize_members(self, prof_type="CURRENT",verb=False,**kwargs):    
        """ Finds minimum weight profiles for members """
               
        top = {'material': 'S355', 'class': 2, 'utility': 1.0, 'bmin':10, 'bmax':1e5}
        bottom = {'material': 'S355', 'class': 2, 'utility': 1.0,'bmin':10, 'bmax':1e5}
        braces = {'material': 'S355', 'class': 2, 'utility_tens': 1.0, 'utility_comp':1.0, 'bmin':10, 'bmax':1e5}
            
        limit_width = False
        
        for key, value in kwargs.items():
            if key == 'top':
                top = value
            elif key == 'bottom':
                bottom = value
            elif key == 'braces':
                braces = value
            elif key == 'limit_width':
                limit_width = value
        
        kmax = 10
        k = 1
                        
        while k < kmax:            
            explored = []
            MEM_PROFILES_CHANGED = []
            self.structural_analysis('all','REM')
            if verb:
                print(f"Optimize profiles, iteration {k}")
            
            if limit_width:
                b_chord_min = self.max_brace_width()
                #b_brace_max = min(self.top_chord_profile.b,self.bottom_chord_profile.b)
                b_brace_max = 1e5
            else:
                b_chord_min = 10
                b_brace_max = 1e5
            
            # Go through groups first
            for name, group in self.member_groups.items():
                if name == "top_chord":
                    material = top['material']
                    sec_class = top['class']
                    max_utility = top['utility']                    
                elif name == 'bottom_chord':
                    material = bottom['material']
                    sec_class = bottom['class']
                    max_utility = bottom['utility']                    
                    
                group.optimum_design(prof_type,verb,material,sec_class,max_utility,b_chord_min)
                
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
                        if list(member.NEd.values())[0] >= 0.0:
                            max_utility = braces['utility_tens']
                        else:
                            max_utility = braces['utility_comp']
                    else:
                        material = 'S355'
                        sec_class = 2
                        max_utility = 1.0
                        
                    MEM_PROFILES_CHANGED.append(member.optimum_design(prof_type,verb,material,sec_class,max_utility,bmax=b_brace_max))
                                        
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
        
def vierendeel_example(L=20000,h=2000,panels=8):
        
    t = VierendeelTruss(span=L, height=h, panels=panels)
    t.generate_topology(nel_chord=6,nel_brace=6)
    t.generate_supports()
    t.generate_joints()
    #t.generate_uniform_load(q=-25)
    P = -17e3
    t.generate_chord_point_load(p=[0.5*P,P])
    t.generate_fem(model='no_eccentricity')
    t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    t.bmd(scale=1,load_id=t.load_ids[0],loads=False)
    t.optimize_members(verb=True,top={'material':'S700','class':2,'utility':1.0},
                       bottom={'material':'S700','class':2,'utility':1.0},
                       braces = {'material': 'S700', 'class': 2, 'utility_tens': 1.0, 'utility_comp':1.0},
                       limit_width=True)
    t.bmd(scale=1,load_id=t.load_ids[0],loads=False)
    return t
    #t.plot(mem_dim=True)


if __name__ == "__main__":
    
    
    from timeit import default_timer as timer
    
    t = vierendeel_example(L=6*3000,h=2500,panels=6)

    t.plot(loads=False)
    
    """
    t.symmetry()
    
   
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
    
    t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    
    t.plot(geometry=True)
    #t.plot(mem_dim=True)
    """
    #opts = {'x_monitor':0.5*t.span, 'n_monitored':2}
    #t.to_abaqus(filename='KT-ristikko',partname="KT-ristikko",options=opts)

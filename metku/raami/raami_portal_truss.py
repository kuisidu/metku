# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 11:33:47 2022

Portal truss: two columns and a truss

TODO (17.1.2023)
    
    1. Liitokset pilareihin: epäkeskisyys alapaarteella ja yläpaarteella
    2. FEM-malli EN 1993 -epäkeskisyyksillä.
    
    Nyt create_truss tekee ristikon liitokset. Kun alapaarre tulee kiinni
    pilariin ja ensimmäinen diagonaali lähtee alapaarteelta yläpaarteelle, ei
    yläpaarteen ja pilarin välille tule paarresauvalle liitosta.
    
    Ristikon ja pilarin liitoksille täytyy luoda omat olionsa.
    
    ** Yläpaarteella on kaksi vaihtoehtoa
    a) yläpaarre tulee pilarin päälle. Tällöin ei tarvita epäkeskisyyttä, ja
        FEM-solmut alkavat paarteen ja pilarin keskilinjojen leikkauspisteestä.
    
    b) yläpaarre tulee pilarin kylkeen. Tällöin FEM-solmut tehdään pilarin päähän
        ja kylkeen. Pilarin kylkeen tulevan solmun kohdalle tulee nivel
        
    Jos ensimmäinen diagonaali lähtee yläpaarteelta alapaarteelle, on diagionaalin
    pää tuotava riittävän etäälle pilarin kyljestä.


    ** Jos alapaarre tulee kiinni pilariin, tulee se pilarin kylkeen. Tällöin
    menetellään kuten yläpaarteen kanssa, eli solmut tehdään pilarin keskilinjan
    ja laipan reunan kohdille. Alapaarteelta lähtevä diagonaali on tuotava
    sopivan etäisyyden päähän pilarin reunasta.    
    
    raami_plane_truss-tiedostossa TopChord -luokan generate_elements ei toimi, jos
    eccentricity on 'en1993', koska sauvalla on vain yksi liitos

@author: kmela
"""

from copy import copy
import matplotlib.pyplot as plt
import numpy as np

from metku.raami.raami import Raami
from metku.raami.frame_member import FrameMember, SteelFrameMember, MultiSpanSteelMember, MemberGroup
from metku.raami.frame_node import FrameNode
from metku.raami.frame_loads import PointLoad, LineLoad, LoadIDs
from metku.raami.frame_supports import XYHingedSupport, FixedSupport
from metku.raami.raami_plane_truss import SlopedTruss, TubularYJoint, TubularKGapJoint, TrussBrace

from metku.framefem.elements.ebbeam import EBBeam

from metku.raami.exports import AbaqusOptions

from metku.sections.steel.ISection import HEA, HEB, IPE
from metku.sections.steel.RHS import RHS, SHS

class PortalTruss(Raami):
    """ Class for portal frame trusses """


    def __init__(self,span=20000,height_eaves=10000,column=None,truss=None,name="Portal Truss"):
        
        super().__init__(name)
                
        self.L = span
        self.height = height_eaves
        self.column_profile = column
        self.truss = truss
        self.columns = []
    
        self.bottom_chord_to_column = False
    
        # Create nodes for column bases
        self.add(FrameNode([0,0]))
        self.add(FrameNode([self.L,0]))
    
        self.column_joints = []
    
    def __repr__(self):
        
        return f"{self.__name__}: height = {self.height:.2f} mm, span = {self.L:.2f} mm"
    
    @property
    def col_length(self):
        """ total length of the column """
        return self.height
    
    @property
    def col_profile(self):
        """ cross-section of the profiles """
        if len(self.columns) == 0:
            p = self.column_profile
        else:
            p = self.columns[0].cross_section
        return p
    
    @col_profile.setter
    def col_profile(self,val):
        """ set cross-section of the profiles """
        self.column_profile = val
        
        for col in self.columns:
            col.cross_section = val    
        
    def weight_detailed(self,verb=False):
        """ Returns detailed weight 
            The weight of truss and columns are separated.
        """
        
        Wtop, Wbottom, Wbraces = self.truss.weight_detailed(False)
        
        Wcols = 2*self.columns[0].weight()
        
        Wtruss = Wtop + Wbottom + Wbraces
        Wtot = Wtruss + Wcols
        
        if verb:
            print(f"Portal Truss weight: {Wtot:.2f} kg")
            print(f"Weight of columns (total): {Wcols:.2f} kg ({Wcols/Wtot*100:.1f}%)")            
            print(f"Truss weight: {Wtruss:.2f} kg ({Wtruss/Wtot*100:.1f}%)")            
            print(f"Top Chord: {Wtop:.2f} kg")
            print(f"Bottom Chord: {Wbottom:.2f} kg")
            print(f"Braces: {Wbraces:.2f} kg")
        
        return Wtop, Wbottom, Wbraces
    
    def ridge_height(self):
        """ Height at ridge level """
        return self.height + self.truss.height_of_slope_left()
    
    def add(self, this):
        """
        Adds object to portal truss
        :param this: object to be added
        :return:
        """
        # COLUMN
        if (isinstance(this, SteelFrameMember) or isinstance(this, MultiSpanSteelMember)) and this.mtype == 'column':
            self.columns.append(this)
            this.mem_id = len(self.members)
            self.members["COL" + str(this.mem_id)] = this
            this.frame = self
            
            if isinstance(this, MultiSpanSteelMember):
                for mem in this.members:
                    mem.frame = self
        else:
            super().add(this)
    
    def create_columns(self,profile=HEA(260),hinges=[False,False], nel=[4,2]):
        """ creates columns """
        if not self.truss is None:
            # Truss has already been created
            
            if self.bottom_chord_to_column:
                # Bottom chord is connected to the column
                self.add(MultiSpanSteelMember([self.nodes[0],self.truss.bottom_nodes[0],self.truss.top_nodes[0]],
                                             profile,mem_type='column',nel=nel,hinges=hinges,sway=True))                    
                self.add(MultiSpanSteelMember([self.nodes[1],self.truss.bottom_nodes[-1],self.truss.top_nodes[-1]],profile,
                                                          mem_type='column',nel=nel,hinges=hinges,sway=True))
            else:
                self.add(SteelFrameMember([self.nodes[0],self.truss.top_nodes[0]],profile,
                                                      mem_type='column',nel=nel,hinges=hinges,lcr=[2.0,1.0],sway=True))
            
                self.add(SteelFrameMember([self.nodes[1],self.truss.top_nodes[-1]],profile,
                                                      mem_type='column',nel=nel,hinges=hinges,lcr=[2.0,1.0],sway=True))
        
        col_group = MemberGroup(self.columns,'columns')
        self.member_groups['columns'] = col_group
            
    
    def create_truss(self, topology='K', ndiv=4, truss_type='simple', **kwargs):
        """ Makes the truss based on user data """
        
        if truss_type == 'simple':
            # Simple truss means that the bottom chord is not connected to the
            # column.
            try:
                dx1 = kwargs['dx']
                dx2 = kwargs['dx']
            except:
                dx1 = 1500
                dx2 = 1500
        else:
            # Non-simple truss means that the bottom chord is connected to the column
            dx1 = 0
            dx2 = 0
        
        try:
            first = kwargs['first_diagonal_up']
        except:
            first = True
        
        try:
            nel_chord = kwargs['nel_chord']
        except:
            nel_chord = 4
        
        try:
            nel_brace = kwargs['nel_brace']            
        except:
            nel_brace = 1
        
        if dx1 == 0:
            self.bottom_chord_to_column = True
        
        self.truss = SlopedTruss(L1=0.5*self.L,L2=0.5*self.L,dx1=dx1,dx2=dx2,**kwargs)
        
        if nel_brace > 1:
            self.truss.braces_as_beams = True
        
        self.truss.origin = (0,self.height)
        self.truss.generate_topology(topology,ndiv,first_diagonal_up=first,edge_verticals=False,nel_chord=nel_chord,nel_brace=nel_brace)
        self.truss.generate_joints()
        self.truss.fem = self.fem
        self.truss.symmetry()
        
        self.member_groups = self.truss.member_groups
        
        # Add the truss nodes to the frame:
        n0 = len(self.nodes)        
        
        for node in self.truss.nodes:
            node.node_id += n0
            self.nodes.append(node)
                
        #for mem in self.truss.members.values():
        #self.add(mem)
        self.members.update(self.truss.members)
        
        """
        self.truss.top_chord[0].hinges[0] = True
        self.truss.top_chord[-1].hinges[1] = True
        
        if self.bottom_chord_to_column:
            self.truss.bottom_chord[0].hinges[0] = True
            self.truss.bottom_chord[-1].hinges[1] = True
        """
    def create_truss_to_column_joints(self,y_joint_node=False):
        """ Creates joints between columns and joints 
            input:
                y_joint_node .. True, if the Y-joint at the chord is given a new
                node.
        """
        
        for col, ndx, left in zip(self.columns,[0,-1],[True,False]):
            # Top chord joints
            newJoint = Truss2ColumnJoint(col.nodes[2],self.truss.top_chord[ndx],col,left)
            self.column_joints.append(newJoint)
            
            # Make hinge to the top chord member
            # NOTE: the index '1' works for both ends of the bottom chord, because
            # the element numbering of the left-hand-side top chord member starts from
            # the first K-joint in the top chord and goes to the column.
            # See: raami_plane_truss.py -> TopChord -> generate_elements
            #self.truss.top_chord[ndx].hinges[1] = True
                            
            if y_joint_node:
                # Create new FrameNode which also replaces the node
                # at the Y joint of the chord.
                for joint in newJoint.chord.joints:
                    if isinstance(joint,TubularYJoint):
                        newCoord = copy(joint.node.coords)
                        D = 0.5*newJoint.hc*abs(np.cos(newJoint.chord.angle)) + newJoint.gap
                        if left:
                            newCoord += D*newJoint.chord.dir_vector
                        else:
                            newCoord -= D*newJoint.chord.dir_vector
                                                
                        newNode = FrameNode(newCoord)                    
                        self.add(newNode)
                        joint.node = newNode
            
            if self.bottom_chord_to_column:
                newJoint = Truss2ColumnJoint(col.nodes[1],self.truss.bottom_chord[ndx],col,left)
                self.column_joints.append(newJoint)
                
                #if left:
                #    self.truss.bottom_chord[ndx].hinges[0] = True
                #else:
                #    self.truss.bottom_chord[ndx].hinges[1] = True
                if y_joint_node:
                    # Create new FrameNode which also replaces the node
                    # at the Y joint of the chord.
                    for joint in newJoint.chord.joints:
                        if isinstance(joint,TubularYJoint):
                            newCoord = copy(joint.node.coords)
                            
                            # Determine the index of the node among the
                            # FrameNodes of the brace of the Y-joint
                            y = newCoord[1]
                            if joint.braces[0].nodes[0].y == y:
                                brace_node_nd = 0
                            else:
                                brace_node_nd = 1
                            
                            D = 0.5*newJoint.hc*abs(np.cos(newJoint.chord.angle)) + newJoint.gap
                            if left:
                                newCoord += D*newJoint.chord.dir_vector
                            else:
                                newCoord -= D*newJoint.chord.dir_vector
                            
                            #print(newCoord,joint.node.coords)
                            newNode = FrameNode(newCoord)       
                            self.add(newNode)
                            joint.node = newNode
                            joint.braces[0].nodes[brace_node_nd] = newNode
            
            
    
    def generate_supports(self,col_base="fixed"):
        """ Creates supports to the truss """
        
        if col_base == 'fixed':
            self.add(FixedSupport(self.nodes[0]))
            self.add(FixedSupport(self.nodes[1]))

        elif col_base == 'pinned':
            self.add(XYHingedSupport(self.nodes[0]))
            self.add(XYHingedSupport(self.nodes[1]))

    def generate_fem(self,truss_eccentricity=False,column_eccentricity=False):
        """ Creates the finite element model """
        if not self.fem_generated:
            self.fem_generated = True
        
        if self.bottom_chord_to_column:
            # In this case, the bottom chord is connected to the column
            # First, generate columns                       
            for col in self.columns:
                for node in col.nodes:
                    if self.dim == 2:
                        newNode = self.fem.add_node(node.x,node.y)
                    else:
                        newNode = self.fem.add_node(node.x,node.y,node.z)
                    
                    node.fem_node = newNode
            
            for member in self.columns:            
                member.generate_elements(self.fem)
                
            # Generate eccentricity elements, if any
            if column_eccentricity:
                for joint in self.column_joints:
                    if joint.left:
                        dx = 0.5*joint.hc
                    else:
                        dx = -0.5*joint.hc
                    
                    dy = 0.5*joint.hc*abs(np.tan(joint.chord.angle))
                    newNode = self.fem.add_node(joint.node.x+dx,joint.node.y+dy)
                                        
                    joint.fem_nodes['xc'] = newNode
                    
                    n1 = joint.node.fem_node
                    n2 = joint.fem_nodes['xc']
                    cs = IPE(500)
                                 
                    newElement = EBBeam(n1,n2,cs,cs.material)                    
                    
                    # Add hinge to the end of the member
                    #newElement.releases = [5]
                    
                    self.fem.add_element(newElement)
                    joint.fem_elements = newElement
            
            # Then, generate truss, including loads on trusses
            self.truss.generate_fem("ecc_elements",nodes_and_elements_only=False)
            
            self.truss.top_chord[0].fem_elements[-1].releases = [5]
            self.truss.top_chord[-1].fem_elements[-1].releases = [5]
            
            if self.bottom_chord_to_column:
                # Finally, connect the bottom chords to the eccentricity elements                
                n1 = self.column_joints[1].fem_nodes['xc']
                n2 = self.column_joints[1].chord.joints[0].fem_nodes['xc']                
                cs = self.column_joints[1].chord.cross_section
                newElement = EBBeam(n1,n2,cs,cs.material)                
                self.fem.add_element(newElement)
                self.column_joints[1].chord.fem_elements.insert(0,newElement)
                
                self.column_joints[1].chord.fem_elements[0].releases = [2]                
                
                n2 = self.column_joints[3].fem_nodes['xc']
                n1 = self.column_joints[3].chord.joints[1].fem_nodes['xc']
                cs = self.column_joints[3].chord.cross_section
                newElement = EBBeam(n1,n2,cs,cs.material)                
                self.fem.add_element(newElement)
                self.column_joints[3].chord.fem_elements.append(newElement)
            
                self.column_joints[3].chord.fem_elements[-1].releases = [5]
            
            self.fem.draw()
        else:
            # In this case, the bottom chord is not connected to the column
            pass
            
        
        # Generate loads
        for load in self.loads.values():
            load.add_load(self.fem)
            # Creates new loadcase if one with same load_id does not exist
            lcase_ids = [lc.load for lc in self.fem.loadcases.values()]
            if load.load_id not in lcase_ids:
                self.fem.add_loadcase(supp_id=1,load_id=load.load_id)
        
        # Generate supports
        for supp in self.supports.values():
            supp.add_support(self.fem)

        # Set nodal degrees of freedom        
        self.fem.nodal_dofs()
        
        """
        # Generate finite element model of the truss        
        
        if truss_eccentricity:
            self.truss.generate_fem("en1993")
        else:
            self.truss.generate_fem("no_eccentricity")
        
        
        # Add fem nodes for the column bases
        for node in self.nodes[:2]:
            if self.dim == 2:
                newNode = self.fem.add_node(node.x,node.y)
            else:
                newNode = self.fem.add_node(node.x,node.y,node.z)
            
            node.fem_node = newNode
       
        
        # Generate finite elements for the columns        
        for member in self.columns:            
            member.generate_elements(self.fem)
        
        # Generate supports
        for supp in self.supports.values():
            supp.add_support(self.fem)

        # Set nodal degrees of freedom        
        self.fem.nodal_dofs()
        """
    
    def clear_fem(self):
        
        super().clear_fem()
        self.truss.clear_fem()
        self.truss.fem = self.fem
        
    
    def update_fem_joint_nodes(self):
        """ Updates nodal positions at joints. It is assumed that cross-sections
            have been changed along the way.
        """
        
        # Update column joint FEM-nodes. This means
        # 1) The node at the column face is moved according to the new column section
        # 2) 
        for joint in self.column_joints:
            if joint.left:
                dx = 0.5*joint.hc
            else:
                dx = -0.5*joint.hc
            
            dy = 0.5*joint.hc*abs(np.tan(joint.chord.angle))
            joint.fem_nodes['xc'].x = joint.node.x+dx
            joint.fem_nodes['xc'].y = joint.node.y+dy   

    
        # Update truss joints FEM nodes
        self.truss.update_fem_joint_nodes()    
    
    def generate_uniform_load(self,part='truss',q=-20,member="top",load_id=LoadIDs['ULS'],ltype='live'):
        
        if part == 'truss':
            # Tässä generoidaan ristikolle kuormitustapaus, jonka load_id on "load_id"
            # Jos ristikolla on jo ko. load_id, lisätään tämä kuormitus ko. load_caseen.
            # Ristikolla on "load_cases" attribuutti, johon kuormitus lisätään.
            #
            # Nyt pitäisi sitten lisätä tämä ristikon kuormitus myös PortalTruss-oliolle.
            #
            # a) jos PortalTruss-luokalla on jo load_id:n omaava kuormitustapaus, lisätään
            #    kuorma tälle.
            #    - voidaan kopioida suoraan ristikolta load_id:n load_case myös kehälle.
            #      tässä voi olla se riski, että jos ristikolle halutaan useita kuormia samalla id:llä,
            #      voi tämä kopiointi epäonnistua. Ehkä dict:n update-metodi toimii tässä.
            #    - update toiminee, jos kaikki ristikon kuormat annetaan ennen pilareiden kuormia.
            #
            self.truss.generate_uniform_load(q,member,load_id,ltype)
            
            if not (load_id in self.load_ids):
                self.load_ids.append(load_id)
                self.add(self.truss.load_cases[load_id])
                
            """
            if not (this.load_id in self.load_ids):
                
                NewLoadCase = LoadCase(load_id=this.load_id,loads = [this])
                
                self.add(NewLoadCase)
                self.load_ids.append(this.load_id)
            else:
                # If a load case with this.load_id exists, add the new
                # load to that load case
                self.load_cases[this.load_id].add(this)
                
            if this.name in self.loads.keys():
                this.name += str(len(self.loads))
            self.loads[this.name] = this
            """
            
        elif part == 'column':
            if not (load_id in self.load_ids):
                self.load_ids.append(load_id)
            
            load_name = 'COL'
            if member == 'left':
                mem = self.columns[0]
                load_name += '_left'
            else:
                mem = self.columns[1]
                load_name += '_right'
            
            vals = [q,q]
            
            self.add(LineLoad(mem, vals, 'x', load_id, 1.0, ltype, load_name, "global"))
    
    def generate_piecewise_uniform_load(self,q,xzones,load_id=LoadIDs['ULS'],ltype='live'):
        """ Adds a piecewise uniform load. See the same method for SlopedTruss in 'raami_plane_truss.py' """
        self.truss.generate_piecewise_uniform_load(q,xzones,load_id,ltype)
        
        if not (load_id in self.load_ids):
            self.load_ids.append(load_id)
            self.add(self.truss.load_cases[load_id])
    
    def plot(self,print_text=True, show=True,
             loads=True, color=False, axes=None, save=False, mem_dim=False, geometry=False):
        
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        self.truss.plot(geometry=geometry,print_text=print_text,show=show,
                        loads=loads,color=color,axes=ax,save=False,mem_dim=mem_dim)
        
        for col in self.columns:
            col.plot(print_text=print_text, c=color, axes=ax, mem_dim=mem_dim)
        
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
        

        # if self.truss:
        #    self.truss.plot(show=False, print_text=print_text, color=color)
        if loads:
            self.plot_loads()
    
    def optimize_members(self, prof_type="CURRENT",verb=False,**kwargs):    
        """ Finds minimum weight profiles for members """
               
        top = {'material': 'S355', 'class': 2, 'utility': 1.0, 'bmin':10, 'bmax':1e5}
        bottom = {'material': 'S355', 'class': 2, 'utility': 1.0,'bmin':10, 'bmax':1e5}
        braces = {'material': 'S355', 'class': 2, 'utility_tens': 1.0, 'utility_comp':1.0, 'bmin':10, 'bmax':1e5}
        columns = {'material': 'S355', 'class': 2, 'utility': 1.0, 'bmin':150, 'bmax':1e5} 
            
        limit_width = False
        
        for key, value in kwargs.items():
            if key == 'top':
                top.update(value)
            elif key == 'bottom':
                bottom.update(value)
            elif key == 'braces':
                braces.update(value)
            elif key == 'columns':
                columns.update(value)
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
                b_chord_min = self.truss.max_brace_width()
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
                    b_min = max(top['bmin'],b_chord_min)
                elif name == 'bottom_chord':
                    material = bottom['material']
                    sec_class = bottom['class']
                    max_utility = bottom['utility']
                    b_min = max(bottom['bmin'],b_chord_min)
                elif name == 'columns':
                    material = columns['material']
                    sec_class = columns['class']
                    max_utility = columns['utility']
                    b_min = columns['bmin']
                    
                group.optimum_design(prof_type,verb,material,sec_class,max_utility,b_min)
                
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
            else:
                # Update nodal positions of FEM nodes at joints
                self.update_fem_joint_nodes()
            k += 1
        
        # Update nodal positions of FEM nodes at joints
        self.update_fem_joint_nodes()
        
        if verb:
            print(f"Number of iterations: {k}.")
            
    
    def to_abaqus(self,target_dir='./',filename="Raami",partname='Truss',options=AbaqusOptions()):
        """ Exports truss to abaqus 
            The export method of 'raami' is mainly used. The main point here is to
            write the MPCs corresponding to eccentricity elements
        """
        
        # options
        #opts = {'x_monitor':0.5*self.span, 'n_monitored':2, 'mpc':[],'ecc_elements':[]}
        
        # Create element sets for top chord and bottom chord
        # top chord is halved
        L2 = 0.5*self.L
        options.elsets['Top_chord_left'] = []
        options.elsets['Top_chord_right'] = []
        #options.elsets['Bottom_chord'] = []
        options.elsets['Bottom_chord_left'] = []
        options.elsets['Bottom_chord_right'] = []
        
        options.elsets['Column_left'] = []
        options.elsets['Column_right'] = []
        
        
        options.load_elsets['Top_chord_load'] = []
        options.load_elsets['Column_left_load'] = []
        options.load_elsets['Column_right_load'] = []
        
        options.load_elsets['Wind_GF'] = []
        options.load_elsets['Wind_H'] = []
        options.load_elsets['Wind_I'] = []
        
        options.load_nsets['Column_left_top_node'] = []
        options.load_nsets['Column_right_top_node'] = []
        
        for mem_id, brace in self.truss.braces.items():
            options.included_members.append('B' + str(mem_id))
        
        for mem in self.truss.top_chord:
            for el in mem.fem_elements:
                if all([n.x <= L2 for n in el.nodes]):
                    options.elsets['Top_chord_left'].append(el)
                else:
                    options.elsets['Top_chord_right'].append(el)
        
        for mem in self.truss.bottom_chord:
            for el in mem.fem_elements:
                #options.elsets['Bottom_chord'].append(el)
                if all([n.x <= L2 for n in el.nodes]):
                    options.elsets['Bottom_chord_left'].append(el)
                else:
                    options.elsets['Bottom_chord_right'].append(el)
        
        for mem in self.columns:
            for el in mem.fem_elements:
                #options.elsets['Bottom_chord'].append(el)
                if mem.nodes[0].x < L2:
                    options.elsets['Column_left'].append(el)
                else:
                    options.elsets['Column_right'].append(el)
        
        for joint in self.truss.joints.values():
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
            #options.elsets['Bottom_chord'].append(el)
            if all([n.x <= L2 for n in el.nodes]):
                options.elsets['Bottom_chord_left'].append(el)
            else:
                options.elsets['Bottom_chord_right'].append(el)
        
        
            
        options.load_elsets['Top_chord_load'] = options.elsets['Top_chord_left'] + \
            options.elsets['Top_chord_right'] 
        options.load_elsets['Column_left_load'] = options.elsets['Column_left']
        options.load_elsets['Column_right_load'] = options.elsets['Column_right']
        
        # Define node sets for point loads at column ends
        options.load_nsets['Column_left_top_node'].append(self.columns[0].members[1].fem_nodes[-1])
        options.load_nsets['Column_right_top_node'].append(self.columns[1].members[1].fem_nodes[-1])
        
        # Determine elsets for wind zones F/G, H, and I
        e = 2*self.ridge_height()
        
        xI = 0.5*e
        xG = 0.1*e
        
        for el in options.load_elsets['Top_chord_load']:            
            x = [node.x for node in el.nodes]
            x.sort()
            
            if x[0] >= xI:
                options.load_elsets['Wind_I'].append(el)
            elif x[0] < xG:
                options.load_elsets['Wind_GF'].append(el)
            else:
                options.load_elsets['Wind_H'].append(el)
        
        
        # Add MPCs for eccentricities of truss-to-column connections        
        for joint in self.column_joints:
            # By default, there is a single element in a truss-to-column joint
            if not (joint.fem_elements is None):                
                el = joint.fem_elements
                # The indexing works, because the first FEM node of the
                # connection is always on the column center line
                options.mpc.append([el.nodes[1].nid+1,el.nodes[0].nid+1])
                options.ecc_elements.append(el)
        
        
        
        # Loads
        # Implementation: each load case is 
        # Vaihtoehdot:
        #for load_id, loads in self.load_cases.items():
            
        
        super().to_abaqus(target_dir,filename,partname,options)
            
class Truss2ColumnJoint:
    """ Class for defining joints between chord and column """
    
    def __init__(self,node,chord,column,left=True):
        
        self.gap = 30
        self.node = node
        self.chord = chord
        self.column = column
        self.left = left
        
        # 'col' .. node at the center line of the column
        # 'xc' .. node at the surface of the column face
        self.fem_nodes = {'col': None, 'xc': None}
        self.fem_elements = None
        
        try:
            self.chord.joints.append(self)
        except:
            pass
                
    
    @property
    def hc(self):
        """ Height of the column profile """
        try:
            hc = self.column.cross_section.h
        except:
            hc = self.column.cross_section.H
        
        return hc
 
        
if __name__ == "__main__":
    
    from timeit import default_timer as timer
    
    p = PortalTruss()
    
    p.create_truss(topology='K',ndiv=3,h2=2000,h1=1500,dx=000,first_diagonal_up=True)
    #p.truss.top_chord[0].add_hinge(0)
    #p.truss.top_chord[-1].add_hinge(1)
    #if p.bottom_chord_to_column:
    #    p.truss.bottom_chord[0].add_hinge(0)
    #    p.truss.bottom_chord[-1].add_hinge(1)
        
    p.generate_uniform_load('truss',q=-25,ltype='snow')
    
    #p.truss.generate_supports()
    p.create_columns(profile=SHS(250,8),hinges=[False,False])
    p.create_truss_to_column_joints(y_joint_node=True)
    p.generate_supports()
    p.symmetry()
    
    
    p.generate_uniform_load('column',q=9,member='left',load_id=LoadIDs['ULS'],ltype='wind')
    p.generate_uniform_load('column',q=5,member='right',load_id=LoadIDs['ULS'],ltype='wind')
    
        
    p.generate_uniform_load('truss',q=-18,ltype='snow',load_id=10)
    p.generate_uniform_load('column',q=9,member='left',load_id=10,ltype='wind')
    p.generate_uniform_load('column',q=5,member='right',load_id=10,ltype='wind')
    
    
    #p.generate_uniform_load('column',q=5,member='left')    
    
    #start = timer()
    print(p.weight())
    #end = timer()
    #print(end-start)
    #p.weight_detailed(True)
    #Wtop, Wbottom, Wbraces = p.truss.weight_detailed(True)
    p.generate_fem(truss_eccentricity=True,column_eccentricity=True)
        
    #print(p.weight())
    #p.structural_analysis(load_id=p.load_ids[0],support_method="REM")
    
    p.structural_analysis(load_id='all',support_method="REM")
    
    
    #p.optimize_members(verb=True)
    #p.bmd(scale=20,load_id=p.load_ids[0],loads=False)
    #p.fem.draw()
    #p.plot()
    
    opts = AbaqusOptions(x_monitor = 0.5*p.L, n_monitored = 2)
    p.to_abaqus(target_dir='C:/Users/kmela/Data/',filename='K-portaali',partname="K-portaali",options=opts)
        
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 11:33:47 2022

Portal truss: two columns and a truss

@author: kmela
"""

import matplotlib.pyplot as plt

from raami import Raami
from frame_member import FrameMember, SteelFrameMember, MultiSpanSteelMember
from frame_node import FrameNode
from frame_loads import PointLoad, LineLoad, LoadIDs
from frame_supports import XYHingedSupport, FixedSupport
from raami_plane_truss import SlopedTruss



from sections.steel.ISection import HEA, HEB, IPE

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
    
    
    def add(self, this):
        """
        Adds object to truss
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
    
    def create_columns(self,profile=HEA(260),hinges=[False,False]):
        """ creates columns """
        if not self.truss is None:
            # Truss has already been created
            
            if self.bottom_chord_to_column:
                # Bottom chord is connected to the column
                self.add(MultiSpanSteelMember([self.nodes[0],self.truss.bottom_nodes[0],self.truss.top_nodes[0]],
                                             profile,mem_type='column',nel=4,hinges=hinges,sway=True))                    
                self.add(MultiSpanSteelMember([self.nodes[1],self.truss.bottom_nodes[-1],self.truss.top_nodes[-1]],profile,
                                                          mem_type='column',nel=4,hinges=hinges,sway=True))
            else:
                self.add(SteelFrameMember([self.nodes[0],self.truss.top_nodes[0]],profile,
                                                      mem_type='column',nel=4,hinges=hinges,lcr=[2.0,1.0],sway=True))
            
                self.add(SteelFrameMember([self.nodes[1],self.truss.top_nodes[-1]],profile,
                                                      mem_type='column',nel=4,hinges=hinges,lcr=[2.0,1.0],sway=True))
            
    
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
        
        if dx1 == 0:
            self.bottom_chord_to_column = True
        
        self.truss = SlopedTruss(L1=0.5*self.L,L2=0.5*self.L,dx1=dx1,dx2=dx2,**kwargs)
        
        
        self.truss.origin = (0,self.height)
        self.truss.generate_topology(topology,ndiv,first_diagonal_up=first,edge_verticals=False)
        self.truss.generate_joints()
        self.truss.fem = self.fem
        self.truss.symmetry()
        
        self.member_groups = self.truss.member_groups
        
        # Add the truss nodes to the frame:
        n0 = len(self.nodes)        
        
        for node in self.truss.nodes:
            node.node_id += n0
            self.nodes.append(node)
        
        for mem in self.truss.members.values():
            self.add(mem)
        #self.members = self.truss.members
    
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
        
        
        # Generate finite element model of the truss        
        
        if truss_eccentricity:
            self.truss.generate_fem("en1993")
        else:
            self.truss.generate_fem("no_eccentricity")
        
        
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
        
        
    def generate_uniform_load(self,part='truss',q=-20,chord="top",load_id=LoadIDs['ULS'],ltype='live'):
        
        if part == 'truss':
            if not (load_id in self.load_ids):
                self.load_ids.append(load_id)
                
            self.truss.generate_uniform_load(q,chord,load_id,ltype)
    
    def plot(self,print_text=True, show=True,
             loads=True, color=False, axes=None, save=False, mem_dim=False):
        
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        self.truss.plot(print_text=print_text,show=show,
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

if __name__ == "__main__":
    
    from timeit import default_timer as timer
    
    p = PortalTruss()
    
    p.create_truss(topology='K',ndiv=3,h2=2000,h1=1500,dx=000,first_diagonal_up=True)
    p.truss.top_chord[0].add_hinge(0)
    p.truss.top_chord[-1].add_hinge(1)
    if p.bottom_chord_to_column:
        p.truss.bottom_chord[0].add_hinge(0)
        p.truss.bottom_chord[-1].add_hinge(1)
        
    p.generate_uniform_load('truss',q=-25)
    #p.truss.generate_supports()
    p.create_columns(hinges=[False,False])
    p.generate_supports()
    p.symmetry()
    
        
    
    #start = timer()
    print(p.weight())
    #end = timer()
    #print(end-start)
    #p.weight_detailed(True)
    #Wtop, Wbottom, Wbraces = p.truss.weight_detailed(True)
    p.generate_fem(truss_eccentricity=False)
    
    #p.optimize_members(verb=True)
    #print(p.weight())
    p.structural_analysis(load_id=p.load_ids[0],support_method="REM")
    #p.bmd(scale=20,load_id=p.load_ids[0],loads=False)
    #p.fem.draw()
    #p.plot()
    
        
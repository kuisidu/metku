# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 13:37:34 2021

Raami: for operating frame structures

@author: kmela
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.path as mpath
from datetime import datetime

import framefem.framefem as fem
from frame_node import FrameNode
from frame_member import FrameMember, SteelFrameMember
from frame_loads import PointLoad, LineLoad, LoadIDs, LoadCase, LoadCombination
from frame_supports import Support, FixedSupport, XYHingedSupport, XHingedSupport, YHingedSupport

#from loadIDs import LoadIDs

class Raami:
    """ Class for dealing with frame structures.
        The frame can be 2D or 3D, and the frame members
        can be of any material.
        
        The frame is determined by the following data:
            Members
            Nodes
            Supports
            Loads
            Loading combinations
        
        The main operations are:
            1) Structural analysis
            2) Member design
            3) Plotting

        The class is intended to be used with structural optimization too
        
    """
    
    def __init__(self,name="Raami"):
        
        
        self.__name__ = name
        self.dim = 2
        self.members = {}        
        self.nodes = []
        self.member_groups = {}
        self.node_groups = {}
        self.supports = {}
        self.loads = {}
        self.load_ids = []
        self.load_cases = {}
        self.load_combs = {}
    
        self.fem = fem.FrameFEM()
        self.fem_generated = False # Flag for checking if FEM model has been created
        self.is_analysed = False # Flag for checking if structure has been analysed
        
    def __repr__(self):
        
        s = "Frame " + self.__name__ + f" with {len(self.members)} members"
        
        return s
    
    #@property
    #def load_cases(self):
    #    """ Returns a list of load case labels """
    #    return self.load_ids
    
    def weight(self):
        """ Returns weight of the frame in kg """
        #W = sum([mem.weight() for mem in self.members.values()])
        
        W = 0
        for mem in self.members.values():
            W += mem.weight()
        
        return W
    
    def add(self,this):
        """ Adding different things to the frame """
        
        if isinstance(this,FrameNode):
            this.node_id = int(len(self.nodes))            
            self.nodes.append(this)
    
        elif isinstance(this, FrameMember):
            # Give member a unique id
            # id is used in creating elements with wanted cross-sectional properties
            this.mem_id = int(len(self.members))

            # Create a link between the FrameMember object and the frame
            this.frame = self

            # Generate member coordinates
            #this.calc_nodal_coordinates()

            self.members[this.mem_id] = this
        elif isinstance(this,LoadCase):            
            self.load_cases[this.load_id] = this
            if not (this.load_id in self.load_ids):
                self.load_ids.append(this.load_id)
        elif isinstance(this,LoadCombination):
            # Connect the current frame to load combination
            this.frame = self
            self.load_combs[this.comb_id] = this
        # POINTLOADS
        elif isinstance(this, PointLoad):

            """ If a point load with same 'name' is already included
                in the frame, add a number to the end of the name
                of the new load such that each point load has a unique
                name.
                
                NOTE: using 'len' to provide new id number might fail,
                if point loads are deleted such that the number of point loads
                changes along the way.
            """            
            
            if not (this.load_id in self.load_ids):
                # There is no load case for the new load.
                # Add a new load case for this load
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

            """ If the location of the point load is in between
                member end nodes, add a node to the corresponding member
            
            for member in self.members.values():
                if member.point_intersection(this.coordinate):
                    member.add_node_coord(this.coordinate)
            """
        # LINELOADS
        elif isinstance(this, LineLoad):

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
            
        # SUPPORTS
        elif isinstance(this, Support):
            #print("Adding support")
            # this.supp_id = len(self.supports)
            supp_label = len(self.supports)
            self.supports[supp_label] = this
            # self.supports[this.supp_id] = this
            #self.support_nodes.append(this.coordinate)

            """ If the support is located between end nodes of a member,
                add a node to that member
            
            for member in self.members.values():                
                if member.point_intersection(this.coordinate):                    
                    member.add_node_coord(this.coordinate)
            """
    def generate_fem(self):
        """ Generates FEM model 
            The idea is that the FEM model is generated once and subsequently
            the same model is modified.
        """
        
        if not self.fem_generated:
            self.fem_generated = True
        
        
        # Create nodes at the locations of the
        # FrameNodes
        for node in self.nodes:
            if self.dim == 2:
                newNode = self.fem.add_node(node.x,node.y)
            else:
                newNode = self.fem.add_node(node.x,node.y,node.z)
            
            node.fem_node = newNode
        
        # Generate elements and internal nodes for members
        for member in self.members.values():            
            member.generate_elements(self.fem)
        
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
    
    def structural_analysis(self, load_id=None, support_method='ZERO'):
        """ Calculates forces and displacements
            
            Parameters
            ----------
            :param load_id: Id of the loads to be added in the calculation model
            
            :type load_id: int / str
        """
        
        # if no load_id is given, use the first load_id available.
        if load_id is None:
            load_id = self.load_ids[0]
        
        # print(f'calculation happening with load id {load_id}')
        if not self.is_analysed:
            self.is_analysed = True
            """ Support ID is always 1! """
            #self.fem.nodal_dofs()

        # If load_id == 'ALL' calculates all load cases
        # calls recursively itself for each case
        if str(load_id).upper() == 'ALL':
            #print("Calculate all cases")
            #lcase_ids = self.load_cases.keys()
            for lid in self.load_cases.keys():                
                self.structural_analysis(load_id=lid, support_method=support_method)
        else:
            #print('Calculate case:' + str(load_id))            
            self.fem.linear_statics(support_method=support_method,lcase=load_id)
            # Using the FEM results, calculate nodal forces and displacements
            self.calc_nodal_forces(load_id)
            #self.calc_nodal_displacements(load_id)
            # Assign forces to members
            #self.assign_forces(load_id)            
            # Do member design
            self.design_members(load_id)            
            # self.alpha_cr, _ = self.f.linear_buckling(k=4)
    
    def design_members(self,load_id = LoadIDs['ULS']):
        """ Designs frame members (checks resistance)
        """
        
        for member in self.members.values():            
            member.design(load_id)
    
    def optimize_members(self, prof_type="CURRENT", verb=False):
        """ Find minimum smallest profiles for frame members
            for given loads.
        """
                
        kmax = 10
        k = 1
                        
        while k < kmax:            
            explored = []
            MEM_PROFILES_CHANGED = []
            self.structural_analysis('all','REM')
            if verb:
                print(f"Optimize profiles, iteration {k}")
                
            # Go through groups first
            for group in self.member_groups.values():                
                group.optimum_design(prof_type,verb)
                
                for mem in group.members:
                    explored.append(mem)
                
            
            for member in self.members.values():
                #print(member)                
                if not member in explored:
                    if verb:
                        print(member)
                    explored.append(member)
                    MEM_PROFILES_CHANGED.append(member.optimum_design(prof_type,verb))
                                        
                    if not member.symmetry_pair is None:
                        explored.append(member.symmetry_pair)
                    
                    #print(member.cross_section)
            
            # If none of the members changed profile, exit the loop
            if not any(MEM_PROFILES_CHANGED):
                break
                
            k += 1
        
        if verb:
            print(f"Number of iterations: {k}.")
    
    def calc_nodal_forces(self,lcase=LoadIDs['ULS']):
        """ Calculates nodal forces and saves values to
            self.nodal_forces - dict       
        """
        # Iterate through every frame's member's node and calculate forces
        # acting on that node
        for member in self.members.values():
            member.calc_nodal_forces(lcase)
    
    def calc_nodal_displacements(self,lcase=LoadIDs['ULS']):
        """ Calculates nodal displacements and saves values to
            self.nodal_displacements -dict
        """
        # Iterate through every frame's member's node and calculate 
        # displacements for that node
        for member in self.members.values():
            member.calc_nodal_displacements(self.fem,lcase)
            
    
    def plot(self, print_text=True, show=True,
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
        if axes is None:
            if self.dim == 2:
                fig, ax = plt.subplots(1)
            else:
                fig, ax = plt.subplots(1,subplot_kw=dict(projection='3d'))
        else:
            ax = axes


        #if self.is_calculated and color:
        #    color = True

        # Plot members
        for member in self.members.values():
            if self.dim == 2:
                member.plot(print_text, color, ax, mem_dim)
            else:
                member.plot3d(print_text, color, ax, mem_dim)

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
            
            if self.dim == 2:
                ax.scatter(node_coord[0], node_coord[1], s=50, c='k',
                           marker=marker)
            else:
                ax.scatter(node_coord[0], node_coord[1], node_coord[2], zdir='y', s=50, c='k',
                           marker=marker)

        # if self.truss:
        #    self.truss.plot(show=False, print_text=print_text, color=color)
        if loads:
            self.plot_loads()
        
        if self.dim == 2:
            ax.axis('equal')
        else:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            
        if save:
            plt.savefig('default.svg', format='svg')
        if show:
            if self.dim == 2:
                plt.axis('equal')
            plt.show()
    
    def plot_nodes(self, print_text=True, show=True, color=False, axes=None, save=False):
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
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes


        #if self.is_calculated and color:
        #    color = True

        # Plot members
        for node in self.nodes:            
            node.plot(print_text, 'ok', ax)        
        
        if self.dim == 2:
            ax.axis('equal')
        
        if save:
            plt.savefig('default.svg', format='svg')
        if show:
            if self.dim == 2:
                plt.axis('equal')
            plt.show()
    
    def plot_loads(self):
        """ Plots loads

            Point loads are blue arrows
            Starting point is the given coordinate

            Line loads are red arrows
            Plots loads above given member

            TODO: Load value plotting
        """

        xrange, yrange, _ = self.bounding_box()
        dX = max(xrange)-min(xrange)
        dY = max(yrange)-min(yrange)

        for load in self.loads.values():
            if isinstance(load, PointLoad):
                x, y = load.coordinate                
                scl = max(dX, dY) * 1e-1
                dx, dy, _ = load.v
                plt.arrow(x, y, np.sign(dx) * scl, np.sign(dy) * scl,
                          head_width=scl * 1e-1, ec='b')
            
            elif isinstance(load,LineLoad):
                c0, c1 = load.coordinates
                q1, q2 = load.values
                dq = (q2 - q1) / 10
                x1, y1 = c0
                x2, y2 = c1
                num_arrows = 10
                dx = (x2 - x1) / num_arrows
                dy = (y2 - y1) / num_arrows
                if dx:
                    X = np.arange(x1, x2 + dx, dx)
                else:
                    X = np.ones(11) * x1
                if dy:
                    Y = np.arange(y1, y2 + dy, dy)
                else:
                    Y = np.ones(11) * y1
                scl = max(dX, dY) * 8e-2
                for i, (x, y) in enumerate(zip(X, Y)):
                    q = q1 + dq * i
                    q_scl = q / max(abs(q2), abs(q1))
                    if load.direction == 'y':
                        # Moves arrows above member
                        y += (np.sign(q1) - 2e-1) * scl * q_scl
                        dx = 0
                        dy = -np.sign(q1) * scl * q_scl
                        plt.arrow(x, y, dx, dy,
                                  head_width=scl * 1e-1, ec='r',
                                  head_starts_at_zero=False)
                    else:
                        x -= (np.sign(q1) + 0.25) * scl * q_scl
                        dx = np.sign(q1) * scl * q_scl
                        dy = 0
                        plt.arrow(x, y, dx, dy,
                                  head_width=scl * 1e-1, ec='y')
    def bmd(self, scale=1, save=False, show=True, loads=True, load_id=LoadIDs['ULS']):
        """ Draws bending moment diagram

            Parameters
            ----------
            :param scale: Scaling factor
            :type scale : int
        """
        self.plot(print_text=False, loads=loads, show=False, color=False)
        for member in self.members.values():
            member.bmd(scale, load_id=load_id)
        
        if save:
            plt.savefig('bending moment diagram.svg', format='svg')
        if show:
            plt.show()
    
    def xrange(self):
        """ Range of x coordinate values """
        xrange = [node.coords[0] for node in self.nodes]
        
        return (min(xrange),max(xrange))
    
    def yrange(self):
        """ Range of y coordinate values """
        yrange = [node.coords[1] for node in self.nodes]
        
        return (min(yrange),max(yrange))
    
    def bounding_box(self):
        """ Determines the range of coordinate values for different
            axes.
        """
        xrange = [node.coords[0] for node in self.nodes]
        yrange = [node.coords[1] for node in self.nodes]
        
        if self.dim == 3:
            zrange = [node.coords[2] for node in self.nodes]
        else:
            zrange = None
        
        return xrange, yrange, zrange
    
    def symmetry(self,axis='y',sym_value=None):
        """ Detects nodes and members lying symmetrically with respect
            to a given axis.
        """
        
        if axis == 'y':
            if sym_value is None:
                xvalues = self.xrange()
                sym_value = 0.5*(xvalues[0]+xvalues[1])
            
            
            nodes_explored = 0                 
            for i, node in enumerate(self.nodes[:-1]):
                if node.symmetry_node is None:
                    nodes_explored += 1
                    for node2 in self.nodes[i+1:]:
                        if abs(node2.y-node.y) < 1e-6 and abs(abs(node2.x-sym_value)-abs(node.x-sym_value)) < 1e-6:
                            node.symmetry_node = node2
                            node2.symmetry_node = node
                            break
            
            #print(f"Explored nodes: {nodes_explored}.")
                            
        # Find symmetrically lying members:
        # One of the two conditions must be satisfied:
        # a) The both nodes of the members are symmetry pairs, or
        # b) One pair of nodes is symmetric and the other is a shared node
        
        sym_mem_keys = []
        
        for key, mem in self.members.items():            
            if key in sym_mem_keys:
                continue
            else:
                sym_mem_keys.append(key)
                
                if mem.nodes[0].symmetry_node is None and mem.nodes[1].symmetry_node is None:
                    # In this case neither of the member nodes has a symmetry node.
                    # This happens usually with verticals along the symmetry line.
                    continue
                else:
                    # Determine the nodes that the supposed symmetry member should
                    # be linked with.
                    sym_nodes = []
                    for node in mem.nodes:
                        if node.symmetry_node is None:
                            sym_nodes.append(node)
                        else:
                            sym_nodes.append(node.symmetry_node)
                            
                    # Explore the remaining members
                    for sym_key in self.members.keys():
                        if not sym_key in sym_mem_keys:
                            if all([sym_node in self.members[sym_key].nodes for sym_node in sym_nodes]):
                                mem.symmetry_pair = self.members[sym_key]
                                self.members[sym_key].symmetry_pair = mem
                                sym_mem_keys.append(sym_key)
    
    def to_robot(self,filename="Raami"):
        """ Writes the frame as a Autodesk Robot Structural Analysis
            str file.
        """
        
        fname_str = f' File Name:  {filename}  '
        now = datetime.now()
        date_str = now.strftime(" Date: %d/%m/%Y    %H:%M ")
        robot_str = " ROBOT 97 v.34.0     "
        
        top_row = ";+" + '-'*len(fname_str) + '+' + '-'*len(date_str) + "+" + '-'*len(robot_str) + "+"
        mid_row = ";!" + fname_str + "!" + date_str + "!" + robot_str + "!"
        
        # Add point to the end of the filename, if it does not exist yet.
        if filename[1] != '.':
            filename += '.'
        
        with open(filename + 'str', 'w') as f:
            f.write(top_row + "\n")
            f.write(mid_row + "\n")
            f.write(top_row + "\n")
            f.write("ROBOT97 \n\n")
            if self.dim == 3:
                f.write("FRAme SPAce \n\n")
            else:
                f.write("FRAme PLAne \n\n")
                
            f.write("NUMbering DIScontinuous \n\n")
            
            f.write(f'NODes {len(self.nodes)}  ELEments {len(self.members)} \n\n')
            f.write("UNIts \n")
            f.write("LENgth=mm	Force=N \n\n")
            
            f.write("NODes \n")
            for node in self.nodes:
                # It is assumed here that node_id is the ordinal number of the node.
                f.write(f"{node.node_id+1: ^9} {node.x: ^20} {node.y: ^20}\n")              
        
            f.write('\n')
            f.write('ELEments \n')
            f.write(';+------+-------+-------+\n')
            f.write(';! No.  ! STRT  ! END   !\n')
            f.write(';+------+-------+-------+\n')
            
            for mem in self.members.values():
                f.write(f'  {mem.mem_id+1: <6} {mem.nodes[0].node_id+1: ^7} {mem.nodes[1].node_id+1: ^7} \n')

            f.write('\n')
            
            f.write('PROperties \n')
            
            for mem in self.members.values():
                f.write(f' "{mem.material}" \n')
                f.write(f' {mem.mem_id + 1}  ')
                f.write(f' {mem.cross_section.robot()}    \n')

            f.write('\n')
            
            f.write("SUPports \n\n")
            
            for sup in self.supports.values():
                s = f' {sup.node.node_id +1} '
                # sup_dof lists those degrees of freedom that are NOT fixed.
                # 3D cases need to be checked.
                if isinstance(sup,FixedSupport):
                    sup_dof = ' '                    
                elif isinstance(sup,XYHingedSupport):
                    if self.dim == 2:
                        sup_dof = 'RY'
                    else:
                        sup_dof = 'RY, RZ, RX'    
                elif isinstance(sup,XHingedSupport):
                    if self.dim == 2:
                        sup_dof = 'UZ, RY'
                    else:
                        sup_dof = 'UY, UX, RX, RY, RZ'
                elif isinstance(sup,YHingedSupport):
                    if self.dim == 2:
                        sup_dof = 'UX, RY'
                    else:
                        sup_dof = 'UX, UZ, RX, RY, RZ'
                else:
                    raise ValueError(f"Unidentified support type {type(sup)}.")
                
                f.write(s + sup_dof + '\n')
          
            f.write('\n')
            
            f.write("RELeases \n")
            for mem in self.members.values():
                s = mem.robot_releases()
                if not s is None:
                    f.write(s)
            
            f.write('\n')

            
            """
            f.write("LOAds \n")
            
            loads2robot = []
                        
            for i, load_id in enumerate(self.load_ids):
                
                f.write(f"CASe # {i+1} LC1 \n")
                if len(self.line_loads):
                    f.write("ELEments \n")
                    for i in range(num_frames):
                        for lineload in self.line_loads.values():
                            n1, n2 = lineload.member.coordinates
                            if num_frames != 1:
                                n1 = [n1[0], i * s, n1[1]]
                                n2 = [n2[0], i * s, n2[1]]
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
            
            """
            f.write("END")
                    
                                
                    
            
if __name__ == '__main__':
    
    from sections.steel.RHS import RHS, SHS
    import cProfile
    
    f = Raami()
    
    L = 5000.0
    H = 2000.0
    nodes = [[0.0,0.0],[0.0,H],[L,H],[L,0.0]]
    
    s = RHS(200,200,6)
    
    for node in nodes:
        f.add(FrameNode(node))
    
    f.add(SteelFrameMember(nodes=f.nodes[:2],section=s,mem_type='column',nel=6))
    f.add(SteelFrameMember(nodes=f.nodes[1:3],section=s,mem_type='beam',nel=6,hinges=[True,True]))
    f.add(SteelFrameMember(nodes=f.nodes[2:4],section=s,mem_type='column',nel=6))
    
    f.add(FixedSupport(f.nodes[0]))
    f.add(FixedSupport(f.nodes[-1]))
    
    #f.add(PointLoad(f.nodes[1], [0, -2000, 0],load_id=LoadIDs.ULS))
    
    f.add(LineLoad(f.members[1],[-20,-20],'y',coord_sys='local',load_id=0,ltype='snow',name='Lumi'))
    f.add(LineLoad(f.members[0],[3,3],'x',coord_sys='global',load_id=1,ltype='wind',name='Tuuli x+'))
    
    f.add(LoadCombination(comb_id=2,comb_type='ULS',load_cases=list(f.load_cases.values())))
    
    f.load_combs[2].combine_loads()
    #f.plot(loads=False)
    f.generate_fem()
    #f.to_robot()
    f.structural_analysis(load_id=2,support_method="REM")
    f.bmd(scale=10,load_id=2)
    #f.optimize_members("CURRENT")
    
            
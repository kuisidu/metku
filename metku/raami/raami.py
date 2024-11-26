# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 13:37:34 2021

Raami: for operating frame structures

@author: kmela
"""
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.path as mpath
from datetime import datetime

#import ..framefem.framefem as fem
from metku.framefem import framefem as fem
from metku.raami.frame_node import FrameNode
from metku.raami.frame_member import FrameMember, SteelFrameMember
from metku.raami.frame_loads import PointLoad, LineLoad, PWLineLoad, LoadIDs, LoadCase, LoadCombination
from metku.raami.frame_supports import Support, FixedSupport, XYHingedSupport, XHingedSupport, YHingedSupport
from metku.raami.exports import AbaqusOptions, write_elset, write_nset

#from loadIDs import LoadIDs

# default options for writing to abaqus
# x_monitor: x-coordinate of the fem nodes to follow (or closest)
# n_monitored: number of monitored nodes
def_opts = {'x_monitor':0,'n_monitored':2,'mpc':[],'ecc_elements':[]}

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
        self.self_weight = False
        
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


    def add_node(self, node:FrameNode) -> None:
        """
        Adds node to frame
        Args:
            node: node to be added to frame

        Returns: None

        """
        node.node_id = int(len(self.nodes))
        self.nodes.append(node)

    def add_member(self, member: FrameMember) -> None:
        """
        Adds member to frame
        Args:
            member: member to be added to frame

        Returns: None

        """
        # Give member a unique id
        # id is used in creating elements with wanted cross-sectional properties
        member.mem_id = int(len(self.members))
        # Create a link between the FrameMember object and the frame
        member.frame = self
        self.members[member.mem_id] = member

    def add_loadcase(self, loadcase:LoadCase) -> None:
        """
        Adds loadcase to frame
        Args:
            loadcase: loadcase to be added to frame

        Returns: None

        """

        self.load_cases[loadcase.load_id] = loadcase
        if not (loadcase.load_id in self.load_ids):
            self.load_ids.append(loadcase.load_id)

        for load in loadcase.loads:
            if load not in self.loads:
                self.add(load)

    def add_loadcombination(self, loadcombination: LoadCombination) -> None:
        # Connect the current frame to load combination
        loadcombination.frame = self
        self.load_combs[loadcombination.comb_id] = loadcombination
        loadcombination.combine_loads()

    def add_pointload(self, pointload:PointLoad) -> None:
        """ If a point load with same 'name' is already included
                        in the frame, add a number to the end of the name
                        of the new load such that each point load has a unique
                        name.

                        NOTE: using 'len' to provide new id number might fail,
                        if point loads are deleted such that the number of point loads
                        changes along the way.
                    """

        if not (pointload.load_id in self.load_ids):
            # There is no load case for the new load.
            # Add a new load case for this load
            NewLoadCase = LoadCase(load_id=pointload.load_id, loads=[pointload])
            self.add(NewLoadCase)
            self.load_ids.append(pointload.load_id)
        else:
            # If a load case with this.load_id exists, add the new
            # load to that load case
            self.load_cases[pointload.load_id].add(pointload)

        if pointload.name in self.loads.keys():
            pointload.name += str(len(self.loads))

        if pointload not in self.loads.values():
            self.loads[pointload.name] = pointload

        """ If the location of the point load is in between
            member end nodes, add a node to the corresponding member

        for member in self.members.values():
            if member.point_intersection(this.coordinate):
                member.add_node_coord(this.coordinate)
        """

    def add_lineload(self, lineload: LineLoad | PWLineLoad) -> None:
        if not (lineload.load_id in self.load_ids):

            NewLoadCase = LoadCase(load_id=lineload.load_id, loads=[lineload])

            self.add(NewLoadCase)
            self.load_ids.append(lineload.load_id)
        else:
            # If a load case with this.load_id exists, add the new
            # load to that load case
            self.load_cases[lineload.load_id].add(lineload)

        if lineload.name in self.loads.keys():
            lineload.name += str(len(self.loads))

        self.loads[lineload.name] = lineload

    def add_support(self, support: Support) -> None:
        # print("Adding support")
        # this.supp_id = len(self.supports)
        supp_label = len(self.supports)
        self.supports[supp_label] = support
        # self.supports[this.supp_id] = this
        # self.support_nodes.append(this.coordinate)

        """ If the support is located between end nodes of a member,
            add a node to that member

        for member in self.members.values():                
            if member.point_intersection(this.coordinate):                    
                member.add_node_coord(this.coordinate)
        """

    def add(self,this) -> 'Raami':
        """
        Adding different things to the frame
        """
        # NODE
        if isinstance(this,FrameNode):
            self.add_node(this)
        # MEMBER
        elif isinstance(this, FrameMember):
            self.add_member(this)
        # LOADCASE
        elif isinstance(this,LoadCase):
            self.add_loadcase(this)
        # LOADCOMBINATION
        elif isinstance(this,LoadCombination):
            self.add_loadcombination(this)
        # POINTLOAD
        elif isinstance(this, PointLoad):
            self.add_pointload(this)
        # LINELOAD
        elif isinstance(this, LineLoad) or isinstance(this, PWLineLoad):
            self.add_lineload(this)
        # SUPPORT
        elif isinstance(this, Support):
            self.add_support(this)
        else:
            raise ValueError(f"Type : {type(this)} cannot be added to Raami -frame")

        return self

    def add_self_weight(self):
        """ Adds self weight to frame's members
        """
        if self.self_weight is False:
            self.self_weight = True
            for member in self.members.values():
                member.add_self_weight()
            # self.add_loads('self_weight')
            # self.generate_frame()
        else:
            self.remove_self_weight()
            self.add_self_weight()    
    
    def remove_self_weight(self):
        """ Removes self weight from frame's members
        """
        self.self_weight = False
        #del (self.line_loads["self_weight"])
        for member in self.members.values():
            member.remove_self_weight()
        
    
    def generate_fem(self):
        """ Generates FEM model 
            The idea is that the FEM model is generated once and subsequently
            the same model is modified.
        """
        if not self.fem_generated:
            self.fem_generated = True

        self.clear_fem()
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
    
    def clear_fem(self):
        """ Clears all data related to FEM and set up an empty FrameFEM object """        
        self.fem = fem.FrameFEM()
        self.fem_generated = False
        
        for mem in self.members.values():
            mem.clear_fem()
            
    
    def structural_analysis(self, load_id=None, support_method='ZERO', design=True):
        """ Calculates forces and displacements
            
            Parameters
            ----------
            :param load_id: Id of the loads to be added in the calculation model
            
            :type load_id: int / str
        """
        self.generate_fem()
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
            load_ids = []
            for lid in self.load_cases.keys():         
                load_ids.append(lid)
            
            #self.structural_analysis(load_id=load_ids, support_method=support_method)
            self.fem.linear_statics(lcase=load_ids, support_method=support_method)
            
            for lid in load_ids:
                self.calc_nodal_forces(lid)
                self.calc_nodal_displacements(lid)
                
                # Do member design
                if design:
                    self.design_members(lid) 
        else:
            #print('Calculate case:' + str(load_id))
            self.fem.linear_statics(support_method=support_method,lcase=load_id)
            # Using the FEM results, calculate nodal forces and displacements
            self.calc_nodal_forces(load_id)
            self.calc_nodal_displacements(load_id)
            #self.calc_nodal_displacements(load_id)
            # Assign forces to members
            #self.assign_forces(load_id)
            # Do member design
            if design:
                self.design_members(load_id)            
            # self.alpha_cr, _ = self.f.linear_buckling(k=4)
    
    def calculate(self, load_id='ALL', support_method='REM', design=True):
        """ This method is essentially the same as structural_analysis.
            It is for the purposes of structural optimization, see
            structopt.py -> OptimizationProblem.fea.
            
            OptimizationProblem was written for Frame2D object, which has
            'calculate' as a method for structural optimization.
            
            If OptimizationProblem.fea is changed such that it calls
            'structural_analysis' instead of 'calculate', then this
            method can be deleted.
        """
        
        self.structural_analysis(load_id,support_method,design)
    
    def design_members(self,load_id = LoadIDs['ULS']):
        """ Designs frame members (checks resistance)
        """
        
        for member in self.members.values():            
            member.design(load_id)
    
    def print_member_utilization(self,filename=None,details=False):
        """ Prints utilization ratios for all members """
        
        if filename is None:        
            print(f"** Member utilization ratios for {self.__name__} **\n")
        
            for member in self.members.values():
                member.print_utilization(details=details)
        else:
            with open(filename, 'w') as file:
                file.write(f"** Member utilization ratios for {self.__name__} **\n")
                for member in self.members.values():
                    member.print_utilization(file,details=details)    
                
    
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
            # TODO! IMPLEMENT calc_nodal_displacements for FrameMember
            member.calc_nodal_displacements(lcase)
            
    
    def plot(self, print_text=True, show=True,
             loads=True, color=False, axes=None, save=False, mem_dim=False,
             saveopts={'filename':'default.svg','format':'svg','orientation':'landscape', 'papertype':'a3'}):
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
            #plt.savefig('default.svg', format='svg')
            plt.savefig(saveopts['filename'],format=saveopts['format'],\
                        orientation=saveopts['orientation'],papertype=saveopts['papertype'])
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
            member.bmd(scale, load_id=load_id, show=False, plot_self=False)
        
        if save:
            plt.savefig('bending moment diagram.svg', format='svg')
        if show:
            plt.show()

    def plot_deflection(self, scale=1, prec=4, show=True, load_id=LoadIDs["ULS"], **kwargs):
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

        self.plot(print_text=False, show=False, **kwargs)
        self.calc_nodal_displacements(lcase=load_id)

        for member in self.members.values():
            X = []
            Y = []
            Z = []
            max_disp = 0
            coord = member.nodes[0].coords
            moved_node = member.nodes[0]
            """ Calculate the deflected location of each node of the member """
            for i, node in enumerate(member.fem_nodes):

                disp = node.u[load_id][:len(node.coord)]
                new_coord = node.coord + disp
                scaled_coord = node.coord.copy() + disp * scale

                #dist = np.sqrt((x0 - (x0+x1)) ** 2 + (y0 - (y0+y1)) ** 2)
                dist = np.linalg.norm(np.array(new_coord) - np.array(node.coord))
                if dist > max_disp:
                    max_disp = dist
                    coord = scaled_coord
                    moved_node = node
                """ Store deflected locations to X and Y """

                if len(coord) == 2:
                    x, y = scaled_coord
                    X.append(x)
                    Y.append(y)


                elif len(coord) == 3:
                    x, y, z = scaled_coord
                    X.append(x)
                    Y.append(y)
                    Z.append(z)

            if len(coord) == 2:
                """ Plot deflected locations """
                plt.plot([moved_node.coord[0], coord[0]], [moved_node.coord[1], coord[1]], color="grey", linestyle="--")
                plt.plot(X, Y, color='gray')
                plt.plot(*coord, 'ro')
                plt.text(*coord, "{0:5.{1}g} mm".format(max_disp, prec))

            elif len(coord) == 3:
                # TODO: Implement 3D plotting
                pass
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
                    
    def to_abaqus(self,target_dir='C:/Users/kmela/Data/',filename="Raami",partname='Frame',options=AbaqusOptions()):
        """ Writes the finite element model to ABAQUS 
        
            1. Create folder for the files
            2. Create main file using 'filename'
            3. Write data to main file
            4. For each member, write a file for cross-sectional data        
        """
        
        path = os.path.join(target_dir,filename + '_ABAQUS')
        
        print(path)
        
        try:
            os.mkdir(path)
        except OSError as error:
            print(f"Folder {path} already exists.")
        
        inp_file = path + '/' + filename + '.inp'
        
        secfiles = []
        matfiles = []
        
        # File for linear elastic material models
        elastic_mat_file = path + '/' + filename + '_Materials_Elastic.txt'
        
        # Create material data for elastic material model
        with open(elastic_mat_file,'w') as file:
            file.write(f'**\n** This file contains materials for {filename}\n**\n')
            
            for mem_name, mem in self.members.items():
                if len(options.included_members) > 0 and mem_name in options.included_members:
                    file.write(f'*Material, name={mem_name}_Mat\n')
                    file.write('*Density\n7850.,\n*Elastic\n2.10000e+11, 0.3\n**\n')
            
            for set_name, els in options.elsets.items():
                file.write(f'*Material, name={set_name}_Mat\n')
                file.write('*Density\n7850.,\n*Elastic\n2.10000e+11, 0.3\n**\n')
        
        # File for all section properties
        sections_file = path + '/' + filename + '_Sections.txt'
        
        # Write header for the file that contains cross-section data
        with open(sections_file,'w') as file:
            now = datetime.now()            
            date_str = now.strftime(" Date: %d/%m/%Y    %H:%M ")
            file.write(f'** Cross-sectional data for {filename}\n')
            file.write(f'** Created on {date_str}\n**\n')
            
        for mem_name, mem in self.members.items():
            # Create cross-section a file for each member
            if len(options.included_members) > 0 and mem_name in options.included_members:
                sec_name = mem.cross_section.__repr__().replace('.','_')
                sec_name = sec_name.replace(' ','_')
                sec_name += '_' + mem.material.__repr__()
                secfile = partname + '_' + mem_name + '_' + sec_name + '.txt'
                matfile = partname + '_' + mem_name + '_' + sec_name + '_Mat.txt'
                secfile_path = path + '/' + secfile
                
                #mem.cross_section.abaqus(secfile_path,matname=mem_name,setname=mem_name)
                #secfiles.append(secfile)
                mem.cross_section.abaqus(sections_file,matname=mem_name,setname=mem_name,open_type='a')
                matfiles.append(matfile)
        
        for set_name, els in options.elsets.items():
            # Create cross-section for each element set
            sec_name = els[0].section.__repr__().replace('.','_')
            sec_name = sec_name.replace(' ','_')
            sec_name += '_' + els[0].material.__repr__()
            secfile = partname + '_' + set_name + '_' + sec_name + '.txt'
            matfile = partname + '_' + set_name + '_' + sec_name + '_Mat.txt'
            secfile_path = path + '/' + secfile
            els[0].section.abaqus(sections_file,matname=set_name,setname=set_name,open_type='a')
            #els[0].section.abaqus(secfile_path,matname=set_name,setname=set_name)
            secfiles.append(secfile)
            matfiles.append(matfile)
            
        with open(sections_file,'a') as file:
            file.write('**')
        
        with open(inp_file,'w') as file:
            now = datetime.now()
            print(now)
            date_str = now.strftime(" Date: %d/%m/%Y    %H:%M ")
            file.write(f'** {filename}\n')
            file.write(f'** Created on {date_str}\n')
            file.write('**\n')
            file.write(f'*Part, name={partname}\n')
            # Write nodal coordinates
            file.write('*Node\n')
            for node in self.fem.nodes:
                file.write(f'{node.nid+1:>7g}, {node.x*1e-3:>14.10f},  0.0000000000, {node.y*1e-3:>14.10f}\n')
                    
            # Write elements.
            # Eccentricity elements are handled by multi-point constraints (MPC),
            # so they are excluded here.
            file.write('*Element, type=B31\n')
            for eid, el in enumerate(self.fem.elements):
                if not el in options.ecc_elements:
                    file.write(f'{eid+1:>5g}, {el.nodes[0].nid+1:>6g}, {el.nodes[1].nid+1:>6g}\n')
            
            # Element sets: the element set for each member
            # is for setting the cross-section for the elements
            s1_rel = []
            s2_rel = []
            for mem_name, mem in self.members.items():
                if len(options.included_members) > 0 and mem_name in options.included_members:
                    nel = len(mem.fem_elements)
                    
                    if mem.hinges[0]:
                        # This append does not make sense!
                        s1_rel.append
                    
                    if nel > 1:
                        #file.write(f'*Elset, elset={mem_name}, generate\n')
                        file.write(f'*Elset, elset={mem_name}\n')
                        for i, ele in enumerate(mem.fem_elements):
                            ndx = self.fem.elements.index(ele)+1
                            if i < nel-1:
                                file.write(f'{ndx:>5g}, ')
                            else:
                                file.write(f'{ndx:>5g}')
                            
                            # If member has releases add it to the corresponding
                            # release set. '2' is for rotational release at the first node
                            # and '5' is for rotational release at the second node
                            if 2 in ele.releases:
                                s1_rel.append(ndx)
                            elif 5 in ele.releases:
                                s2_rel.append(ndx)
                                                    
                        file.write('\n')
                    else:
                        file.write(f'*Elset, elset={mem_name}\n')
                        e1 = self.fem.elements.index(mem.fem_elements[0])+1
                        file.write(f'{e1:>5g}\n')
                            #e1 = self.fem.elements.index(mem.fem_elements[0])+1
                            #e2 = self.fem.elements.index(mem.fem_elements[-1])+1
                            #file.write(f'{e1:>5g}, {e2:>5g},   1\n')
            
            # Gap elements (if any)
            top_ndx = []            
            for el in options.top_gap_elements:
                top_ndx.append(self.fem.elements.index(el)+1)
            
            if len(top_ndx) > 1:
                write_elset(file, 'TopChordGaps', top_ndx)
                
            bottom_ndx = []
            for el in options.bottom_gap_elements:
                bottom_ndx.append(self.fem.elements.index(el)+1)
            
            if len(bottom_ndx) > 1:
                write_elset(file, 'BottomChordGaps', bottom_ndx)
                
            # If there are any other element sets, they are printed here
            for set_name, els in options.elsets.items():                
                el_ndx = [self.fem.elements.index(el)+1 for el in els]                
                write_elset(file, set_name, el_ndx)
                
            # Release sets
            file.write('**\n')
            
            write_elset(file,"S1Release",s1_rel)
            file.write('**\n')
            write_elset(file,"S2Release",s2_rel)
           
            # Include cross-sectional data from file 'sections_file'
            file.write(f'*Include, Input={sections_file}\n')
            #*Include, Input=S023C005M002_RedMat.txt            
            #for secfile in secfiles:
            #    file.write(f'*Include, Input={secfile}\n')
            
    
            # State element releases            
            file.write('**\n*RELEASE\n')
            file.write('S1Release, S1, M1\n')
            file.write('S2Release, S2, M1\n')

            #[elementti setin nimi alkupää], S1, M1
            #[elementti setin nimi loppupää], S2, M1
            
            file.write('*End Part\n')
            
            # Start Assembly part of the file. This includes
            # 1) Supports
            # 2) Multi-point constraints (hinges)
            # 3) Loads
            file.write('**\n**\n** ASSEMBLY \n**\n')
            file.write(f'*Assembly, name={partname}_Assembly \n**\n')
            
            # Name of the assembly instance
            insname = filename
            file.write(f'*Instance, name={insname}, part={partname}\n')
            file.write('*End Instance\n**\n')

            # Write element sets for loads
            for set_name, els in options.load_elsets.items():
                if len(els) > 0:
                    el_ndx = [self.fem.elements.index(el)+1 for el in els]                
                    write_elset(file, set_name, el_ndx, insname)
                    
            for set_name, nodes in options.load_nsets.items():
                if len(nodes) > 0:
                    node_ndx = [self.fem.nodes.index(node)+1 for node in nodes]
                    write_nset(file, set_name, node_ndx, insname)

            # Node sets for the following groups
            # 1) Symmetry conditions to prevent displacement in y direction
            #   (all nodes EXCEPT slave nodes of MPC)
            # 2) Supports: a) pinned support, b) slide support, where displacement in z direction
            #   is prevented.
            # 3) Monitored node: node, whose displacement will be traced.
            file.write('** Node set for out-of-plane support\n')
            file.write(f'*Nset, nset=BC_Y_symm_nodes, instance={insname}\n')
            
            # By default, all nodes are master nodes            
            y_symm_nodes = [n.nid+1 for n in self.fem.nodes]
            for mpc in options.mpc:
                # If MPCs have been defined in 'options', the first
                # element in each member of mpc list is the slave node index
                # (in abaqus numbering)
                y_symm_nodes.remove(mpc[0])
            
            i = 0
            while i < len(y_symm_nodes):
                file.write(', '.join(str(r) for r in y_symm_nodes[i:i+16]))
                file.write('\n')
                i += 16
            
            file.write('** Node set for monitored nodes.\n')
            # Find 'n_monitored' closed nodes with respect to the x coordinate in the FEM
            # nodes. These will be monitored
            monitored_nodes_set = 'Monitored_Nodes'
            monitored_nodes = list(np.argsort([abs(n.x-options.x_monitor) for n in self.fem.nodes])[:options.n_monitored])
            file.write(f'*Nset, nset={monitored_nodes_set}, instance={insname}\n')
            file.write(', '.join(str(mn) for mn in monitored_nodes))
            file.write('\n')
            
            # Supported nodes
            supp_fem_nodes = {'fixed':[], 'hinged':[], 'y_supp': [], 'x_supp': []}            
            for supp in self.supports.values():                
                if isinstance(supp,XYHingedSupport):
                    supp_fem_nodes['hinged'].append(supp.node.fem_node.nid)
                elif isinstance(supp,YHingedSupport):
                    supp_fem_nodes['y_supp'].append(supp.node.fem_node.nid)
                elif isinstance(supp,XHingedSupport):
                    supp_fem_nodes['x_supp'].append(supp.node.fem_node.nid)
                elif isinstance(supp,FixedSupport):
                    supp_fem_nodes['fixed'].append(supp.node.fem_node.nid)
            
            #print(supp_fem_nodes)
            
            file.write('** Node sets for supports\n')
            for supp_type, nids in supp_fem_nodes.items():
                nnodes = len(nids)
                if nnodes > 0:
                    file.write(f'*Nset, nset={supp_type}, instance={insname}\n')
                    for i, nid in enumerate(nids):                    
                        if i < nnodes-1:
                            file.write(f'{nid+1:>5g}, ')
                        else:
                            file.write(f'{nid+1:>5g}')
                    file.write('\n')
            
            # Eccentricity elements: write them as MPCs.
            for mpc in options.mpc:
                file.write('*MPC\n')
                file.write(f'BEAM, {insname}.{mpc[0]}, {insname}.{mpc[1]}\n')
            
            file.write('*End Assembly\n')
            
            # Write material files:
            file.write('**\n** MATERIALS\n**\n')
            file.write(f'*Include, Input={elastic_mat_file}\n')
            #for matfile in matfiles:
            #    file.write(f'*Include, Input={matfile}\n')
            
            # Write supports
            file.write('**\n** BOUNDARY CONDITIONS\n**\n')        
            for supp_type, nids in supp_fem_nodes.items():
                nnodes = len(nids)
                if nnodes > 0:
                    file.write('*Boundary\n')
                    if supp_type == 'hinged':
                        supp_code = 'PINNED'
                    elif supp_type == 'fixed':
                        supp_code = 'ENCASTRE'
                    elif supp_type == 'y_supp':
                        supp_code = '3' # pitäisikö olla ZSYMM?
                    elif supp_type == 'x_supp':
                        supp_code = 'XSYMM'
                    file.write(f'{supp_type}, {supp_code}\n')
            
            file.write('** --------------------------------\n')
            file.write('** STEP\n')
            file.write('*Step, name=Buckling, nlgeom=YES, extrapolation=NO, inc=200\n')
            file.write('*Static, riks\n')
            file.write('0.01, 1., 1e-16, 0.05, 1.,\n')

            # Write loads
            file.write('**\n** LOADS\n**\n')
            
            # Self-weight, imposed on all elements
            file.write('*Dload\n')
            file.write(', GRAV, 9.81, 0., 0., -1.\n')
            
            # Write output requests
            file.write('**\n** OUTPUT REQUESTS\n**\n')
            file.write('**FIELD OUTPUT:\n')
            file.write('*Output, field, variable=PRESELECT\n')
            file.write('*Output, field\n')
            file.write('*Node Output\n')
            file.write('U\n')
            file.write('*Element Output, directions=YES\n')
            file.write('E, PEEQ, S, MISES, SF\n')
            
            #file.write('*Output, field\n')
            #file.write(f'*Node Output, nset={monitored_nodes_set}\n')
            #file.write('U1, U2, U3, UR1, UR2, UR3\n')
            
            file.write('**HISTORY OUTPUT:\n')
            file.write('*Output, history, variable=PRESELECT\n')
            file.write('*Output, history\n')
            file.write(f'*Node Output, nset={monitored_nodes_set}\n')
            file.write('U1, U2, U3, UR1, UR2, UR3\n')

            file.write('*End Step\n')
            
if __name__ == '__main__':
    
    from metku.sections.steel.RHS import RHS, SHS
    import cProfile
    
    f = Raami()



    L = 5000.0
    H = 2000.0
    nodes = [[0.0,0.0],[0.0,H],[L,H],[L,0.0]]
    
    s = RHS(200,200,6)
    
    for node in nodes:
        f.add(FrameNode(node))
    
    f.add(SteelFrameMember(nodes=f.nodes[:2],section=s,mem_type='column',nel=20))
    f.add(SteelFrameMember(nodes=f.nodes[1:3],section=s,mem_type='beam',nel=20,hinges=[True,True]))
    f.add(SteelFrameMember(nodes=f.nodes[2:4],section=s,mem_type='column',nel=20))
    
    f.add(FixedSupport(f.nodes[0]))
    f.add(FixedSupport(f.nodes[-1]))
    
    #f.add(PointLoad(f.nodes[1], [0, -2000, 0],load_id=LoadIDs.ULS))
    
    f.add(LineLoad(f.members[1],[-20,-20],'y',coord_sys='local',load_id=0,ltype='snow',name='Lumi'))
    f.add(LineLoad(f.members[0],[1000,1000],'x',coord_sys='global',load_id=1,ltype='wind',name='Tuuli x+'))
    
    f.add(LoadCombination(comb_id=2,comb_type='ULS',load_cases=list(f.load_cases.values())))


    #f.plot(loads=True)
    f.generate_fem()
    #f.to_robot()
    f.structural_analysis(load_id=2,support_method="REM")
    #f.bmd(scale=10,load_id=2, loads=False)


    #f.members[0].bmd(load_id=2)
    # print(f.members[1].nodal_displacements[2])
    f.plot_deflection(load_id=2)
    #f.optimize_members("CURRENT")
    
            
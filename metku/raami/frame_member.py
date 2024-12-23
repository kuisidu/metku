# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Tue Dec 28 13:46:06 2021

Class for Frame Members, to be used with Raami

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.path as mpath

from metku.structures.steel.steel_member import SteelMember
from metku.raami.frame_node import FrameNode
import metku.framefem.framefem
from metku.framefem.elements.ebbeam import EBBeam, EBBeam3D
from metku.framefem.elements.rod import Rod
from metku.sections.steel.catalogue import make_section, ipe_profiles, h_profiles, rhs_profiles, shs_profiles, chs_profiles, hea_profiles, heb_profiles
from metku.sections.steel.CHS import CHS
from metku.sections.steel.RHS import RHS, SHS
from metku.sections.steel.ISection import IPE, HEA, HEB

#from loadIDs import LoadIDs

from .frame_loads import LoadIDs, LineLoad

def profile_list(cross_section,prof_type="CURRENT",Tmin=3.0,sec_class=3):
    """ Returns a list of profiles matching the type of 'cross_section.
        This is used in 'optimum_design' methods.    
       
        Tmin .. minimum wall thickness, for RHS/SHS profiles
        sec_class .. maximum section class (in compression)
    """
    
    prof_type = prof_type.upper()
    
    if prof_type == "CURRENT":
        """ Use the cross section family of the current profile """
        if isinstance(cross_section,IPE):
            prof_type = "IPE"
        elif isinstance(cross_section,HEA):
            prof_type = "HEA"
        elif isinstance(cross_section,HEB):
            prof_type = "HEB"
        elif isinstance(cross_section,RHS):
            prof_type = "RHS"
        elif isinstance(cross_section,SHS):
            prof_type = "SHS"
        elif isinstance(cross_section,CHS):
            prof_type = "CHS"
        else:
            prof_type = "IPE"
    
    if prof_type == "IPE":            
        profiles = ipe_profiles.keys()
    elif prof_type == "H":
        profiles = h_profiles.keys()
    elif prof_type == "HEA":
        profiles = hea_profiles.keys()
    elif prof_type == "HEB":
        profiles = heb_profiles.keys()
    elif prof_type == "RHS":
        profiles = []
        for key, value in rhs_profiles.items():
            pnew = make_section(key)
            pnew.Ned = -100e3
            if pnew.section_class() > sec_class:
                continue
            if value['t'] >= Tmin:
                profiles.append(key)
        #profiles = rhs_profiles.keys()
    elif prof_type == "CHS":
        profiles = chs_profiles.keys()
    elif prof_type == "SHS":
        profiles = shs_profiles.keys()
    else:
        raise ValueError(f'"{prof_type}" is not a valid profile type!')
    
    return profiles


class FrameMember:
    """ General frame member class
        
        Relevant member data includes:
            1) Cross-section and material
            2) Nodes
            3) Type: beam, column, brace, etc.
    
        Main operations:
            - Use member data to generate FEM model
            - From FEM results, read internal forces and displacements
            - Design member based on internal forces and displacements
        
    """

    def __init__(self,nodes,section,mem_type="beam",mem_id="",nel=2,\
                 hinges=[False,False],lcr=[1.0,1.0],sway=False):
        """
        Constructor

        Parameters
        ----------
        nodes : list of FrameNode objects
            End nodes of the member.
        section : varies
            Cross-sectional data.
        mem_type: string
            Type of member: 'beam', 'column', 'brace', etc.
        nel: int
            Number of elements for dividing the member
        hinges: list (boolean)
            States whether or not the ends of the member are hinged.
        Returns
        -------
        None.

        """
        
        self.active = True
        self.resistance_check = True
        self.allowable_utilization = 1.0
        
        self.fem_nodes = []
        self.fem_elements = []
        
        self.member = None
        
        # Buckling length factor
        self.__lcr = None
        self.lcr = lcr
                
        self.nodes = nodes
        
        for node in nodes:
            node.members.append(self)
        
        self.symmetry_pair = None # Potential symmetry member pair
        self.group = [] # List of members that are in the same group. These have
                          # have the same cross-section as the current member.
        
        self.cross_section = section
        self.mtype = mem_type
        self.mem_id = mem_id
        self.frame = None
        self.num_elements = nel
        self.local_node_coords = [] # Local coordinates of member nodes in FEM model
        self.global_node_coords = [] # Local coordinates of member nodes in FEM model
        
        self.loads = [] # List of loads acting on the member
        
        self.nodal_forces = {} # Nodal forces in different load cases
        self.nodal_displacements  = {} # Nodal displacements in different load cases
        self.NEd = {} # Design value of the axial force (maximum) for each load case
        self.VzEd = {} # Design value of the shear force (maximum) for each load case
        self.MyEd = {} # Design value of the bending moment (maximum) for each load case
        self.utilization = {} # Utilization ratio for different load cases
        
        # Create nodal coordinates
        self.calc_nodal_coordinates(nel)
        
        # Flags for stating whether or not the member ends are hinged.
        
        self.hinges = hinges
        
        # Flag stating whether or not the member can buckle in sway mode
        self.sway = sway
        
        self.costs = {}
        
        self.self_weight = False
        
        
        
    def __repr__(self):        
        return f"{type(self).__name__} [{self.mem_id}]"
    
    @property
    def cross_section(self):        
        return self.__cross_section
    
    @cross_section.setter
    def cross_section(self,val):
        
        self.__cross_section = val
        if not self.member is None:
            self.member.profile = val            
        
        self.update_element_sections()
        
        if not self.symmetry_pair is None:
            self.symmetry_pair.__cross_section = val
            if not self.symmetry_pair.member is None:
                self.symmetry_pair.member.profile = val
                self.symmetry_pair.update_element_sections()
            
    @property
    def material(self):
        """ Material of the member """
        return self.cross_section.material
    
    @material.setter
    def material(self,val):
        """ Set member material """
        self.cross_section.material = val
    
    @property
    def lcr(self):
        return self.__lcr
    
    @lcr.setter
    def lcr(self,val):
        self.__lcr = val
        
        if not self.member is None:
            self.member.lcr = val
    
    def coords(self):
        """ Return nodal coordinates """
        x1, x2 = np.asarray([self.nodes[0].coords,self.nodes[1].coords])
        return x1, x2
    
    def length(self):
        """ Calculate the length of the member """
        start_node, end_node = self.coords()
        return np.linalg.norm(end_node - start_node)
    
    def weight(self):
        """
        Function that calculates member's weight
        :return weight: float, weight of the member
        """

        return self.cross_section.A * self.length() * self.cross_section.density
    
    def cost(self):
        """ Total cost of the member 
            By default this is the material cost
        
        """
        
        if not self.member is None:
            c = self.member.cost()
        else:
            c = self.length()*self.cross_section.cost()
        
        return c
    
    @property
    def dir_vector(self):
        """
        Returns direction vector of the member.
        :return:
        """
        start, end = self.coords()
        unit = np.array(end - start) / self.length()
        return unit
    
    def local_to_global(self, loc):
        """
        Returns local coordinate in global coordinates
        :param loc: local coordinate
        :return: global coordinate
        """
        Lx = self.length() * loc
        return self.nodes[0].coords + Lx * self.dir_vector

    def global_to_local(self, glob) -> float:
        """
        Converts a global 3D coordinate to a local coordinate.
        :param glob: global 3D coordinate
        :return: local coordinate
        """

        displacement = glob - self.nodes[0].coords
        Lx = np.dot(displacement, self.dir_vector)
        loc = Lx / self.length()

        return loc

    @property
    def perpendicular(self):
        """
        Returns vector perpendicular to member.
        NOTE! This needs to be modified for 3D!        
        :return:
        """
        unit = self.dir_vector
        #if self.reverse:
        #    return np.array([-unit[1], unit[0]])

        return np.array([unit[1], -unit[0]])
    
    @property
    def angle(self):
        """
        Returns member's angle with respect to x-axis in radians
        """
        start_node, end_node = self.coords()
        x0, y0 = start_node
        x1, y1 = end_node
        if abs(x1 - x0) < 1e-6:
            # 90 degrees
            angle = 0.5*np.pi # np.radians(90)
        else:
            angle = np.arctan((y1 - y0) / (x1 - x0))
            if angle < 0:
                return np.pi + angle #np.radians(180) + angle
        return angle
    
    def add_self_weight(self):
        """ Adds self-weight to the member's loads
        """
        if not self.self_weight:
            self.self_weight = True
            load_id = LoadIDs['SELF-WEIGHT']
            
            # Default units: N/mm
            value = -self.cross_section.self_weight()
            self.frame.add(LineLoad(self, [value, value], 'y',load_id,ltype='dead',
                                    name='self-weight'))

    def remove_self_weight(self):
        """ Removes self-weight from loads
        """
        if self.self_weight:
            self.self_weight = False
        
    def recalculate_coordinates(self):
        """ Clear global node coordinates and recalculate those """
        self.global_node_coords.clear()
        self.calc_nodal_coordinates()

    def calc_nodal_coordinates(self, num_elements=0):
        """ Calculates nodal coordinates along the member 
            
        """
                        
        self.global_node_coords = []
        
        if num_elements > 0:
            self.num_elements = num_elements
        
        self.local_node_coords = np.linspace(0, 1, self.num_elements + 1)
        #for loc in np.arange(0, 1 + dloc, dloc):
        for loc in self.local_node_coords:
            coord = self.local_to_global(loc)
            self.global_node_coords.append(coord)



    def calc_nodal_displacements(self, load_id=LoadIDs['ULS']):
        """
        Collects the displacements at the nodes into
        nodal_displacements' attribute which is a dict.
        """
        local_coords = []
        displacements = []
        for node in self.fem_nodes:
            displacements.append(node.u[load_id])
            local_coords.append(self.global_to_local(node.coord))
        # Sort displacements starting from member start coordinate
        displacements = [x for _, x in sorted(zip(local_coords, displacements))]
        self.nodal_displacements[load_id] = displacements


    def calc_nodal_forces(self,load_id=LoadIDs['ULS']):
        """ Collects the internal forces at the nodes into
            'nodal_forces' attribute which is a dict.
        """
                
        def absmax(vals):
            min_val = min(vals)
            max_val = max(vals)
            abs_max = max(abs(min_val), abs(max_val))
            if abs_max == abs(min_val):
                return min_val
            else:
                return max_val

        """ Store also maximum bending moment, shear force and axial force """
        max_med = 0
        max_ved = 0
        max_ned = 0
        
        # set up the nodal_forces dict for this load case
        self.nodal_forces[load_id] = {'N':[],'Vz':[],'My':[]}
        
        for element in self.fem_elements:
            # The elements in fem_elements appear in the order from the first
            # end node to the last. Reading the forces in the first node of each element
            # generates a list of nodal forces that appear in the same order as
            # the local coordinates.
            self.nodal_forces[load_id]['N'].append(element.fint[load_id]['fx'][0])
            self.nodal_forces[load_id]['Vz'].append(element.fint[load_id]['fz'][0])
            self.nodal_forces[load_id]['My'].append(element.fint[load_id]['my'][0])

        # Read the forces in the last node
        self.nodal_forces[load_id]['N'].append(self.fem_elements[-1].fint[load_id]['fx'][1])
        self.nodal_forces[load_id]['Vz'].append(self.fem_elements[-1].fint[load_id]['fz'][1])
        self.nodal_forces[load_id]['My'].append(self.fem_elements[-1].fint[load_id]['my'][1])
        
        # Determine maximum (design) values of the internal forces
        self.NEd[load_id] = absmax(self.nodal_forces[load_id]['N'])
        self.VzEd[load_id] = absmax(self.nodal_forces[load_id]['Vz'])
        self.MyEd[load_id] = absmax(self.nodal_forces[load_id]['My'])
        
        """
        # i is node index
        i = 0        
        # node_ids is a list of node indices for the member
        node_ids = list(self.nodes.keys())
        node_forces = {}
        for element in self.elements.values():
            node_id = node_ids[i]
            #axial_force = absmax(element.axial_force)
            axial_force = absmax(element.fint[load_id]['fx'])
            if abs(axial_force) > abs(max_ned):
                max_ned = axial_force

            #shear_force = absmax(element.shear_force)
            shear_force = absmax(element.fint[load_id]['fz'])
            if abs(shear_force) > abs(max_ved):
                max_ved = shear_force

            #bending_moment = absmax(element.bending_moment)
            bending_moment = absmax(element.fint[load_id]['my'])
            if abs(bending_moment) > abs(max_med):
                max_med = bending_moment

            
            #self.nodal_forces[node_id] = [axial_force,
                                          shear_force,
                                          bending_moment]
            
            node_forces[node_id] = [axial_force,
                                    shear_force,
                                    bending_moment]
            i += 1

            if i == len(self.elements):

                node_id = node_ids[i]
                #axial_force = absmax(element.axial_force)
                axial_force = absmax(element.fint[load_id]['fx'])
                if abs(axial_force) > abs(max_ned):
                    max_ned = axial_force

                #shear_force = absmax(element.shear_force)
                shear_force = absmax(element.fint[load_id]['fz'])
                if abs(shear_force) > abs(max_ved):
                    max_ved = shear_force

                #bending_moment = absmax(element.bending_moment)
                bending_moment = absmax(element.fint[load_id]['my'])
                if abs(bending_moment) > abs(max_med):
                    max_med = bending_moment

                
                #self.nodal_forces[node_id] = [axial_force,
                                              shear_force,
                                              bending_moment]
                
                node_forces[node_id] = [axial_force,
                                    shear_force,
                                    bending_moment]

        self.nodal_forces[load_id] = node_forces
        self.med = max_med
        self.ved = max_ved
        self.ned = max_ned
        """

    def get_max_forces(self, load_ids=None):
        """
        Returns max and min forces for member.

        :return: dict of max and min forces
        example item:
        'min N': {'N':0,
                  'Vz': 100,
                  'My': 200,
                  'load_id': 300}

        """
        if load_ids is None:
            load_ids = self.frame.load_ids

        force_keys = list(next(iter(self.nodal_forces.values())).keys())
        functions = [min, max]
        combined = [f"{y.__name__} {x}" for x in force_keys for y in functions]
        empty_result = {k: 0 for k in force_keys}
        empty_result["load_id"] = None

        result_dict = {k: empty_result for k in combined}

        sliced_dict = {key: self.nodal_forces[key] for key in load_ids if key in self.nodal_forces}

        for load_id, forces in sliced_dict.items():
            for force in force_keys:
                for func in functions:
                    val = func(forces[force])
                    idx = forces[force].index(val)
                    key = f"{func.__name__} {force}"
                    if result_dict[key]["load_id"] is None:
                        result_dict[key] = {key: value[idx] for key, value in sliced_dict[load_id].items()}
                        result_dict[key]["load_id"] = load_id
                    elif abs(val) > abs(result_dict[key][force]):
                        result_dict[key] = {key: value[idx] for key, value in sliced_dict[load_id].items()}
                        result_dict[key]["load_id"] = load_id

        return result_dict


    def design(self,load_id = LoadIDs['ULS']):
        """ Designs the member (check resistance) """
        
        """ Here, the internal forces of the member at sections must
            somehow be inserted. Note that SteelMember does not yet
            provide several load cases.
        """

        rmax = self.member.design()
        
        if rmax > self.allowable_utilization:
            self.resistance_check = False
        else:
            self.resistance_check = True
            
        self.utilization[load_id] = rmax
    
    def print_utilization(self,file=None,details=False):
        """ Prints utilization ratios for different load cases """
        if file is None:
            print(f' Member {self.mem_id:3.0f}: {self.cross_section.__repr__():s}')
            if details:
                print(f'  Length = {self.length():6.0f} mm')
                print(f'  Buckling Length Factor, y-axis= {self.lcr[0]:4.2f} ')
                print(f'  Buckling Length Factor, z-axis= {self.lcr[1]:4.2f} ')
            for load_id, U in self.utilization.items():
                print(f'  Load case {load_id}: {U:4.3f}')    
        else:
            file.write(f' Member {self.mem_id:3.0f}: {self.cross_section.__repr__():s}\n')
            if details:
                file.write(f'  Length = {self.length():6.0f} mm\n')
                file.write(f'  Buckling Length Factor, y-axis= {self.lcr[0]:4.2f}\n')
                file.write(f'  Buckling Length Factor, z-axis= {self.lcr[1]:4.2f}\n')
            for load_id, U in self.utilization.items():
                file.write(f'  Load case {load_id}: {U:4.3f}\n')
                       
    def optimum_design(self,prof_type,verb):
        """ Searches for the minimum weight profile of 'prof_type' 
            This method must be implemented for each material separately
            because of the profile types.
        """
        
        return False
        
    def opposite_node(self,node):
        """ Returns the node of the member opposite to the 'node'
            in the input argument.
        """
        n = None
        
        for n in self.nodes:
            if n != node:
                break
        
        return n
    
    def generate_elements(self,fem):
        """ Creates internal nodes and elements between nodes 
            Assumes that the member end nodes have already been created
        """
    
        # First end node                    
        self.fem_nodes.append(self.nodes[0].fem_node)
        
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
        self.fem_nodes.append(self.nodes[1].fem_node)
        
        # Create members
        # The zip command allows to create two iterables, where the
        # items are actually consecutive items of the list fem_nodes.
        for n1, n2 in zip(self.fem_nodes,self.fem_nodes[1:]):            
            if self.mtype == 'beam' or self.mtype == 'column':
                #section = framefem.BeamSection(self.cross_section.A,self.cross_section.Iy)
                #material = self.make_framefem_material()
                if self.frame.dim == 2:        
                    newElement = EBBeam(n1,n2,self.cross_section,self.material)
                else:
                    newElement = EBBeam3D(n1,n2,self.cross_section,self.material)
            elif self.mtype == 'bar' or self.mtype == 'brace':
                """ Simple bar element that carries only axial forces """
                #material = self.make_framefem_material()
                newElement = Rod(n1,n2,self.cross_section,self.material)
            elif self.mtype == "top_chord" or self.mtype == "bottom_chord": 
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
    
    def clear_fem(self):
        """ Removes FEM nodes and elements """
        self.fem_nodes.clear()
        self.fem_elements.clear()
    
    def update_element_sections(self):
        """ Updates the cross sections of the finite elements of the member """
        if len(self.fem_elements) > 0:
            for element in self.fem_elements:
                element.section = self.cross_section
    
    def plot(self, print_text=True, c='k', axes=None, mem_dim=False, dy = 50):
        
        if self.active:
            if axes is None:
                fig, ax = plt.subplots(1)
            else:
                ax = axes
    
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
            
            # Plot text
            if self.mtype == 'beam':
                horzalign = 'center'                
                vertalign = 'bottom'
    
            elif self.mtype == 'column':
                horzalign = 'right'
                vertalign = 'center'
    
            else:
                horzalign = 'center'
                vertalign = 'center'
    
            x, y = self.local_to_global(0.5) - 0.5*self.perpendicular * dy
            rot = np.degrees(self.angle)
    
            if rot > 90:
                rot = -(180-rot)
    
            if print_text:                
                ax.text(x, y, str(self.mem_id) + ": " + str(self.cross_section),
                         rotation=rot, rotation_mode='anchor',
                         horizontalalignment=horzalign,
                         verticalalignment=vertalign)
            #else:
            #    ax.text(x, y, str(self.mem_id),
            #             rotation=rot, horizontalalignment=horzalign,
            #             verticalalignment=vertalign)
    
    def plot3d(self, print_text=True, c='k', axes=None, mem_dim=False, dy = 50):
        
        if self.active:
            if axes is None:
                fig, ax = plt.subplots(1,subplot_kw=dict(projection='3d'))
            else:
                ax = axes
    
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
                # zdir = 'y' indicates that mplot3D package considers the second coordinate
                # as the z coordinate. 
                ax.plot([X0[0], X1[0]], [X0[1], X1[1]], [X0[2], X1[2]], zdir='y', color=color)
            
            """
            # Plot text
            if self.mtype == 'beam':
                horzalign = 'center'
                vertalign = 'bottom'
    
            elif self.mtype == 'column':
                horzalign = 'right'
                vertalign = 'center'
    
            else:
                horzalign = 'center'
                vertalign = 'center'
    
            x, y = self.local_to_global(0.5) - self.perpendicular * dy
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
            """

    def bmd(self, scale=1, load_id = LoadIDs['ULS'], plot_self=True, show=True):
        """ Plots bending moment diagram """
        if plot_self:
            self.plot(  )
        # Scales Nmm to kNm
        unit_scaler = 1e-6

        # Member's position vector
        #v = self.dir_vector
        # Vector perpendicular to v
        u = -self.perpendicular
        X = []
        Y = []
        start_coord, end_coord = self.coords()
        x01, y01 = start_coord
        x02, y02 = end_coord
        X.append(x01)
        Y.append(y01)
        
        # Get bending moment values for the given load case
        moment_values = self.nodal_forces[load_id]['My']
        
        max_moment = max(moment_values)
        min_moment = min(moment_values)

        for elem in self.fem_elements:
            x0, y0 = elem.nodes[0].coord
            bending_moment = elem.fint[load_id]['my'][0]
            # Scale the bending moment
            val = bending_moment * unit_scaler * scale
            # Add point to the bending moment diagram
            x, y = np.array([x0, y0]) - val * u
            X.append(x)
            Y.append(y)
            horzalign = 'center'
            #vertalign = 'center'
            # Add the value of the maximum and minimum values as text in the figure
            if bending_moment == max_moment or bending_moment == min_moment:
                plt.text(x, y, f'{bending_moment * unit_scaler:.2f} kNm',
                         horizontalalignment=horzalign)
        
        # These points are the last in the list
        x0, y0 = elem.nodes[1].coord
        bending_moment = elem.fint[load_id]['my'][1]
        val = bending_moment * unit_scaler * scale
        x, y = np.array([x0, y0]) - val * u
        X.append(x)
        Y.append(y)
        X.append(x02)
        Y.append(y02)
        plt.text(x, y, f'{bending_moment * unit_scaler:.2f} kNm',
                 horizontalalignment=horzalign)
        plt.plot(X, Y, color='gray')

        if show:
            plt.show()

    def add_hinge(self, loc=0):
        """ Adds a hinge to the member
        input:
            loc .. location along the member (should be 0 or 1)
        """
        if loc == 0:
            self.hinges[0] = True
        else:
            self.hinges[1] = True
    
    def robot_releases(self):
        """ Returns a string indicating releases that Robot Structural Analysis
            can use in creating a str file.
        """
        if self.hinges[0] == True and self.hinges[1] == True:
            s = f'ELEments {self.mem_id+1} ORIgin RY END RY\n'
        elif self.hinges[0] == True and self.hinges[1] == False:
            s = f'ELEments {self.mem_id+1} ORIgin RY \n'
        elif self.hinges[0] == False and self.hinges[1] == True:
            s = f'ELEments {self.mem_id+1} END RY \n'
        else:
            s = None
        return s

class MultiSpanMember:
    """ Class for member with multiple spans """
    
    
    def __init__(self,nodes,section,mem_type="beam",mem_id="",nel=2,hinges=[False,False],lcr=[1.0,1.0],sway=False):
        """
        Constructor

        Parameters
        ----------
        nodes : list of FrameNode objects
            Provides the nodes of the member. Each pair of consecutive nodes constitutes
            one span of the member
        section : SteelSection
            Cross-section of the member.
        mem_type : string, optional
            member type. The default is "beam".
        mem_id : string, optional
            Member id or name. The default is "".
        nel : int, optional
            Number of finite elments per span. The default is 2.
        hinges : TYPE, optional
            DESCRIPTION. The default is [False,False].
        lcr: list, optional
            buckling length factors for the principal axes
            

        Returns
        -------
        None.

        """
        
        self.active = True
        self.resistance_check = True
        self.allowable_utilization = 1.0
        
        self.fem_nodes = []
        self.fem_elements = []
        
        self.symmetry_pair = None # Potential symmetry member pair
        
        self.members = None
        
        # Buckling length factor
        self.lcr = lcr
        
        
        self.nodes = nodes        
        self.cross_section = section
        self.mtype = mem_type
        self.mem_id = mem_id
        self.frame = None
        
        # num_elements is the number of elements per member segment.
        # If 'nel' is given as a scalar value, then a list is made
        # with each element having the value 'nel'.
        if isinstance(nel,list):
            self.num_elements = nel
        else:
            self.num_elements = [nel for i in range(len(nodes)-1)]
            
        self.local_node_coords = [] # Local coordinates of member nodes in FEM model
        self.global_node_coords = [] # Local coordinates of member nodes in FEM model
        
        self.loads = [] # List of loads acting on the member
        
        self.nodal_forces = {} # Nodal forces in different load cases
        self.NEd = {} # Design value of the axial force (maximum) for each load case
        self.VzEd = {} # Design value of the shear force (maximum) for each load case
        self.MyEd = {} # Design value of the bending moment (maximum) for each load case
        self.utilization = {} # Utilization ratio for different load cases
        
        # Create nodal coordinates
        #self.calc_nodal_coordinates(nel)
        
        # Flags for stating whether or not the member ends are hinged.
        
        self.hinges = hinges
        
        self.sway = sway
        
    def __repr__(self):        
        return f"{type(self).__name__} [{self.mem_id}]"
    
    @property
    def cross_section(self):        
        return self.__cross_section
    
    @cross_section.setter
    def cross_section(self,val):
        
        self.__cross_section = val
        if not self.members is None:
            for mem in self.members:
                mem.cross_section = val
        
            self.update_element_sections()
    
    @property
    def lcr(self):
        return self.__lcr
    
    @lcr.setter
    def lcr(self,val):
        self.__lcr = val
        
        if not self.members is None:
            for mem in self.members:
                mem.lcr = val
        
    @property
    def material(self):
        """ Material of the member """
        return self.cross_section.material
    
    def plot(self, print_text=True, c='k', axes=None, mem_dim=False, dy = 50):
        
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
                            
        for mem in self.members:
            mem.plot(print_text, c, ax, mem_dim, dy)
    
    def bmd(self, scale=1, load_id = LoadIDs['ULS']):
        """ Bending moment diagram """
        
        for mem in self.members:
            mem.bmd(scale,load_id)
        
    
    def coords(self):
        """ Returns list of nodal coordinates """
        x = [np.asarray(node.coords) for node in self.nodes]
        return x
    
    def end_coords(self):
        """ Returns the coordinates of the end nodes """
        x1, x2 = np.asarray([self.nodes[0].coords,self.nodes[-1].coords])
        return x1, x2

    def length(self):
        """ Calculate the length of the member """
        start_node, end_node = self.end_coords()
        return np.linalg.norm(end_node - start_node)
    
    def weight(self):
        """
        Function that calculates member's weight
        :return weight: float, weight of the member
        """        
        return self.cross_section.A * self.length() * self.cross_section.density
    
    @property
    def dir_vector(self):
        """ Direction vector of the member """
        return self.members[0].dir_vector

    @property
    def perpendicular(self):
        """ Unit vector perpendicular to the axis of the member """        
        return self.members[0].perpendicular
    
    @property
    def angle(self):
        """
        Returns member's angle with respect to x-axis in radians
        """
        return self.members[0].angle
    
    def local_to_global(self, loc):
        """
        Returns local coordinate in global coordinates
        :param loc: local coordinate
        :return: global coordinate
        """
        Lx = self.length() * loc
        return self.nodes[0].coords + Lx * self.dir_vector
    
    def global_to_local(self, x):
        """ Returns the local coordinate of the point x along the
            member axis.
            
            It is assumed that 'x' lies on the member axis, this not checked.
        """
        x0 = self.nodes[0].coords
        L = self.length()
        
        if isinstance(x,list):
            xloc = [np.linalg.norm(xi-x0)/L for xi in x]
        else:
            xloc = np.linalg.norm(x-x0)/L
        
        return xloc
            
    
    def calc_nodal_coordinates(self, num_elements=0):
        """ Calculates nodal coordinates along the member 
            
        """
                
        if num_elements > 0:
            self.num_elements = num_elements
        
        self.global_node_coords = []
        
        for mem in self.members:
            mem.calc_nodal_coordinates(num_elements)
            self.global_node_coords += mem.global_node_coords[:-1]
        
        self.global_node_coords.append(self.members[-1].global_node_coords[-1])
        
        self.local_node_coords = self.global_to_local(self.global_node_coords)
    
    def calc_nodal_displacements(self,load_id=LoadIDs['ULS']):
        """ Calculates the nodal displacements of each 'submember'.
        """
        for mem in self.members:
            mem.calc_nodal_displacements(load_id)
    
    def calc_nodal_forces(self,load_id=LoadIDs['ULS']):
        """ Collects the internal forces at the nodes into
            'nodal_forces' attribute which is a dict.
        """
        for mem in self.members:
            mem.calc_nodal_forces(load_id)
        
        
    def generate_elements(self,fem):
        """ Creates internal nodes and elements between nodes """
    
        for mem in self.members:
            mem.generate_elements(fem)
            self.fem_elements += mem.fem_elements
            
        if self.hinges[0]:
            
            """ There is a hinge in the first end """
            """ NOTE! This is wrong: the for loop above adds the fem elements and nodes
                for the individual members and not to the MultiSpan member
            """
            self.fem_elements[0].releases = [2]
        
        if self.hinges[1]:
            """ There is a hinge in the last end """
            self.fem_elements[-1].releases = [5]    
    
    def clear_fem(self):
        
        for mem in self.members:
            mem.clear_fem()
        
        self.fem_elements.clear()
    
    def update_element_sections(self):
        """ Updates the cross sections of the finite elements of the member """
        for mem in self.members:
            mem.update_element_sections()
    
    
    def robot_releases(self):
        """ Returns a string indicating releases that Robot Structural Analysis
            can use in creating a str file.
        """
        if self.hinges[0] == True and self.hinges[1] == True:
            s = f'ELEments {self.mem_id+1} ORIgin RY END RY\n'
        elif self.hinges[0] == True and self.hinges[1] == False:
            s = f'ELEments {self.mem_id+1} ORIgin RY \n'
        elif self.hinges[0] == False and self.hinges[1] == True:
            s = f'ELEments {self.mem_id+1} END RY \n'
        else:
            s = None
        return s
class SteelFrameMember(FrameMember):
    """ Class for steel frame members """
    
    def __init__(self,nodes,section,mem_type='beam',mem_id="",nel=2,hinges=[False,False],lcr=[1.0,1.0],sway=False):
        """

        Parameters
        ----------
        nodes : TYPE
            DESCRIPTION.
        section : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        
        super().__init__(nodes,section,mem_type,mem_id,nel,hinges,lcr,sway)
                
        self.member = SteelMember(self.cross_section, self.length(), Lcr=lcr, mtype=self.mtype)
    
        self.create_sections()
    
    def __repr__(self):
        
        s = super().__repr__()
        s += ': ' + self.cross_section.__repr__()
        return s
       
    def material_cost(self):
        """ Material cost of the member """
        
        return self.member.cost()
    
    def blasting_cost(self):
        
        return self.frame.workshop.cost_centres['blasting'].cost(self.length())
        
    def painting_cost(self):
    
        return self.frame.workshop.cost_centres['painting'].cost(self.length()*self.cross_section.Au)
    
    def make_framefem_material(self):
        """ Creates a Material class object that can be used to generate material for FrameFem elements """
        
        mat = framefem.Material(E=self.cross_section.E,rho=self.cross_section.density,nu=self.cross_section.material.nu)
        
        return mat
    
    def create_sections(self):
        """ Creates cross sections to be checked """
        
        for xloc in self.local_node_coords:
            self.member.add_section(loc=xloc)
    
    def is_section(self,profile):
        """ Compares the current cross section of the member with
            'profile'.
        """
        
        return (self.cross_section.A-profile.A)**2 + (self.cross_section.Iy-profile.Iy)**2 < 1e-6
        
    
    def design(self,load_id = LoadIDs['ULS']):
        """ Designs the member (check resistance) """
        
        # Assign forces to SteelMember sections
        #self.member.ned = self.nodal_forces[load_id]['N']
        #self.member.vzed = self.nodal_forces[load_id]['Vz']
        #self.member.myed = self.nodal_forces[load_id]['My']

        N = self.nodal_forces[load_id]['N']
        Vz = self.nodal_forces[load_id]['Vz']
        My = self.nodal_forces[load_id]['My']
        
        for loc, n, vz, my in zip(self.local_node_coords, N, Vz, My):
            self.member.sections[loc]['N'] = n
            self.member.sections[loc]['Vz'] = vz
            self.member.sections[loc]['My'] = my

        rmax = self.member.design(sway=self.sway)
        
        if rmax > self.allowable_utilization:
            self.resistance_check = False
        else:
            self.resistance_check = True
            
        self.utilization[load_id] = rmax
    
    def optimum_design(self, prof_type='CURRENT',verb=False,material="S355",
                       sec_class = 3, max_utility=1.0,
                       bmin=10, bmax=1e5, max_slenderness=5.0):
        """ Goes through all profiles in list and stops iterating when 
            cross-section can bear it's loads.
            If member is previously designed, iterating starts from 3 profiles
            before currrent profile.
        """        
        
        # Get list of profiles profiles of chosen type
        profiles = profile_list(self.cross_section,prof_type)
        
        initial_profile = self.cross_section

        for profile in profiles:
            # profile is a string containing the name of the profile
            # It is converted to a SteelSection object using the make_section
            # function of the 'catalogue.py' file.
            self.cross_section = make_section(profile)
            self.material = material
            self.cross_section.Ned = -1
            
            if self.cross_section.b < bmin or self.cross_section.b > bmax:
                continue
            
            if self.cross_section.section_class() > sec_class:
                continue
            
            if any(self.member.slenderness() > max_slenderness):
                continue
            
            for load_id in self.frame.load_ids:
                self.design(load_id)
                
                # Design also symmetry pair, if it exists
                if not self.symmetry_pair is None:
                    self.symmetry_pair.design(load_id)
                
                            
            if all(r<=max_utility for r in self.utilization.values()):# <= 1.0:                            
                # Check also symmetry pair, if it exists
                if not self.symmetry_pair is None:
                    if all(r<=max_utility for r in self.symmetry_pair.utilization.values()):# <= 1.0:                        
                        break
                else:
                    break
        
        """ If the profile changed during iteration, return 'True'. """
        if not self.is_section(initial_profile):
            if verb:
                print(f'Initial profile: {initial_profile}')
                print(f'New profile: {self.cross_section}')
            return True
        else:
            return False


class MultiSpanSteelMember(MultiSpanMember):
    """ Class for steel frame members """
    
    def __init__(self,nodes,section,mem_type='beam',mem_id="",nel=2,hinges=[False,False],lcr=[1.0,1.0],sway=False):
        """

        Parameters
        ----------
        nodes : TYPE
            DESCRIPTION.
        section : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        super().__init__(nodes,section,mem_type,mem_id,nel,hinges,lcr,sway)
        
        self.members = []
                        
        n = len(self.members)
        for nel, n1, n2 in zip(self.num_elements,nodes[:-1],nodes[1:]):
            new_mem_id = mem_id + str(n)
            n += 1
            self.members.append(SteelFrameMember([n1,n2],section,mem_type,new_mem_id,nel,lcr=lcr,sway=sway))
            self.members[-1].frame = self.frame
    
    def is_section(self,profile):
        """ Compares the current cross section of the member with
            'profile'.
        """
        
        return (self.cross_section.A-profile.A)**2 + (self.cross_section.Iy-profile.Iy)**2 < 1e-6
    
    def print_utilization(self,file=None, details=False):
        """ Prints utilization ratios for different load cases """
        if file is None:
            print(f' Member {self.mem_id:3.0f}: {self.cross_section.__repr__():s}')
            if details:
                print(f'  Length = {self.length():6.0f} mm')
                print(f'  Buckling Length Factor, y-axis= {self.lcr[0]:4.2f} ')
                print(f'  Buckling Length Factor, z-axis= {self.lcr[1]:4.2f} ')
            for load_id, U in self.utilization.items():
                print(f'  Load case {load_id}: {U:4.3f}')    
        else:
            file.write(f' Member {self.mem_id:3.0f}: {self.cross_section.__repr__():s}\n')
            if details:
                file.write(f'  Length = {self.length():6.0f} mm\n')
                file.write(f'  Buckling Length Factor, y-axis= {self.lcr[0]:4.2f}\n')
                file.write(f'  Buckling Length Factor, z-axis= {self.lcr[1]:4.2f}\n')
            for load_id, U in self.utilization.items():
                file.write(f'  Load case {load_id}: {U:4.3f}\n')
    
    def design(self,load_id = LoadIDs['ULS']):
        """ Designs the member (check resistance) """
        
        rmax = 0
        
        for mem in self.members:
            mem.design(load_id)
            rmax = max(rmax,mem.utilization[load_id])
                        
        if rmax > 1.0:
            self.resistance_check = False
        else:
            self.resistance_check = True
            
        self.utilization[load_id] = rmax
    
    def optimum_design(self,prof_type="CURRENT",verb=False,material="S355",
                       sec_class=3,max_utility=1.0):
        """ Goes through all profiles in list and stops iterating when 
            cross-section can bear it's loads.
            If member is previously designed, iterating starts from 3 profiles
            before currrent profile.
        """
        
        #Get list of profiles profiles of chosen type
        profiles = profile_list(self.cross_section,prof_type)        
        initial_profile = self.cross_section        

        for profile in profiles:
            # profile is a string containing the name of the profile
            # It is converted to a SteelSection object using the make_section
            # function of the 'catalogue.py' file.
            self.cross_section = make_section(profile,material)
            
            #print(self.cross_section,self.cross_section.material)
            
            sec = self.cross_section
            sec.Ned = -1
            
            if sec.section_class() > sec_class:
                continue
            
            for load_id in self.frame.load_ids:   
                self.design(load_id)
                            
            if all(r<=self.allowable_utilization for r in self.utilization.values()):# <= 1.0:
                break
            
        
        """ If the profile changed during iteration, return 'True'. """    
        if not self.is_section(initial_profile):
            if verb:
                print(f'Initial profile: {initial_profile}')
                print(f'New profile: {self.cross_section}')
            return True
        else:
            return False

class MemberGroup:
    """ Class for a group of members. Used primarily for setting the cross-section
        and optimization.
    """
    
    def __init__(self,members=[],name="Group"):
        """
        

        Parameters
        ----------
        members : List
            List of initial members.
        name : TYPE, optional
            DESCRIPTION. The default is "Group".

        Returns
        -------
        None.

        """

        self.members = members
        if len(members) == 0:
            self.cross_section = None
            self.frame = None
        else:
            self.cross_section = members[0].cross_section        
            self.frame = members[0].frame
        self.utilization = {}
        self.resistance_check = True
        
        

    @property
    def cross_section(self):        
        return self.__cross_section
    
    @cross_section.setter
    def cross_section(self,val):
        
        self.__cross_section = val
        if len(self.members) > 0:
            for mem in self.members:
                mem.cross_section = val
        
            mem.update_element_sections()
    
    @property
    def material(self):
        """ Material of the member """
        return self.__cross_section.material
    
    @material.setter
    def material(self,val):
        """ Set member material """
        self.__cross_section.material = val
        if len(self.members) > 0:
            for mem in self.members:
                mem.cross_section.material = val
            
    def add_member(self,new_member):
        
        self.members.append(new_member)
        
        if self.cross_section is None:
            self.cross_section = new_member.cross_section
        
        if self.frame is None:
            self.frame = new_member.frame
    
    def is_section(self,profile):
        """ Compares the current cross section of the member with
            'profile'.
        """
                        
        return (self.cross_section.A-profile.A)**2 + (self.cross_section.Iy-profile.Iy)**2 < 1e-6
    
    def design(self,load_id = LoadIDs['ULS']):
        """ Designs the member group (check resistance) """
        
        rmax = 0
        
        for mem in self.members:
            mem.design(load_id)
            rmax = max(rmax,mem.utilization[load_id])
                        
        if rmax > 1.0:
            self.resistance_check = False
        else:
            self.resistance_check = True
            
        self.utilization[load_id] = rmax
    
    def optimum_design(self,prof_type="CURRENT",verb=False,material="S355",
                       sec_class=3,max_utility=1.0,
                       bmin=10, bmax=1e5, max_slenderness=5.0):
        """ Goes through all profiles in list and stops iterating when 
            cross-section can bear it's loads.
            If member is previously designed, iterating starts from 3 profiles
            before currrent profile.
        """
        
        #print(self.cross_section)
        #print(prof_type)
        #Get list of profiles profiles of chosen type        
        profiles = profile_list(self.cross_section,prof_type)
        
        #print(profiles)

        initial_profile = self.cross_section

        for profile in profiles:
            # profile is a string containing the name of the profile
            # It is converted to a SteelSection object using the make_section
            # function of the 'catalogue.py' file.
            self.cross_section = make_section(profile)            
            self.material = material
            sec = self.cross_section
            sec.Ned = -1            
            
            #print(sec)
            
            if self.cross_section.b < bmin or self.cross_section.b > bmax:
                continue
            
            if sec.section_class() > sec_class:
                continue
            
            for member in self.members:
                # For multispan members, the member of a group contains
                # several members.
                if isinstance(member,MultiSpanMember):
                    for mem in member.members:
                        if any(mem.member.slenderness() > max_slenderness):
                            continue
                else:
                    if any(member.member.slenderness() > max_slenderness):
                        continue
            
            
            for load_id in self.frame.load_ids:   
                self.design(load_id)
                            
            if all(r<=max_utility for r in self.utilization.values()):# <= 1.0:
                break
            
        
        """ If the profile changed during iteration, return 'True'. """    
        if not self.is_section(initial_profile):
            if verb:
                print(f'Initial profile: {initial_profile}')
                print(f'New profile: {self.cross_section}')
            return True
        else:
            return False

if __name__ == "__main__":
    
    L = 8000
    
    X = [0, 0.6*L, L]
    Y = [0, 0 ,0]
    
    nodes = []
    for x, y in zip(X,Y):
        nodes.append(FrameNode([x,y]))
    
    section = HEA(260)
    
    m = MultiSpanSteelMember(nodes, section, mem_id="B")
    m.plot(dy=0)
    m.calc_nodal_coordinates()
    
        
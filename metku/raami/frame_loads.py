# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 14:52:58 2021

Loads for Raami objects

@author: kmela
"""

import numpy as np
from copy import copy

import framefem.framefem as fem

#from loadIDs import LoadIDs

LoadIDs = {'ULS': 0,
           'SLS_Characteristic': 1,
           'SLS_Frequent': 2,
           'SLS_Quasi_permanent': 3,
           'ACC': 4}

class LoadCase:
    """ Class for load cases. One load case represents a set of loads that
        are due to a single source. These loads can be of different type (PointLoad or LineLoad)
        and they can act on different parts of the structure, but they act
        simultaneously.
    """
    
    def __init__(self,load_id=0,loads=[]):
        """        

        Parameters
        ----------
        load_id : int, optional
            Unique load case identifier. The default is 0.
        loads : list, optional
            List of Load objects. The default is None.

        Returns
        -------
        None.

        """
        
        self.load_id = load_id
        self.loads = loads
    
    def __repr__(self):
        
        return f"Load case {self.load_id}. Total loads: {len(self.loads)}"
    
    @property
    def load_type(self):
        """ Returns the type of load """
        if len(self.loads) > 0:
            return self.loads[0].ltype
        else:
            return None
    
    def add(self,load):
        """ Adds a Load class object 'load' to the load case """
        load.load_id = self.load_id
        
        self.loads.append(load)

    def combine(self,comb_type='ULS',leading=True):
        """ Create factored loads for load combination of type comb_type """
        
        combined_loads = []
        
        for load in self.loads:
            comb_load = copy(load)
            comb_load.factored(comb_type,leading)
            combined_loads.append(comb_load)
        
        return combined_loads
            
        
class LoadCombination:
    """ Class for load combinations """
    
    
    def __init__(self,comb_id,comb_type='ULS',load_cases=[]):
        """
        
        Parameters
        ----------
        comb_id : float
            Unique identifier of the load combination.
        comb_type : string
            Type of combination. Can be 'ULS', 'SLS', or 'ACC'. This
            determines the combination factors
        load_cases : TYPE, optional
            DESCRIPTION. The default is [].

        Returns
        -------
        None.

        """
        
        self.comb_id = comb_id
        self.comb_type = comb_type
        self.load_cases = load_cases
        self.lead_variable_load = None
        
        # The first non-dead load is the leading variable load case.
        for lc in load_cases:
            if not lc.load_type == 'dead':
                self.lead_variable_load = lc
                break
        
        self.frame = None
    
    def combine_loads(self):
        """ Does the actual load combination:
            Each load in each load case is multiplied by its corresponding
            load combination factor that depends on the load type.
        """
        
        for lc in self.load_cases:
            # Do the load factorization for given load case
            # combo_loads is a list of factored loads.
            if lc == self.lead_variable_load:            
                combo_loads = lc.combine(self.comb_type,leading=True)
            else:
                combo_loads = lc.combine(self.comb_type,leading=False)
                
            for cl in combo_loads:
                cl.load_id = self.comb_id
                self.frame.add(cl)
        

class Load:
    """ General class for loads """

    def __init__(self, load_id, ltype='live', name='Load', f=1.0):
        """        

        Parameters
        ----------
        load_id : TYPE
            DESCRIPTION.
        ltype : string, optional
            Load type. Can be 'live', 'dead', 'wind', or 'snow'. The default is 'live'.
        name : string, optional
            Name of the load. The default is 'Load'.
        f : float, optional
            scaling factor for the load. The default is 1.0.

        Returns
        -------
        None.

        """
        
        
        self.load_id = load_id
        self.name = name
        self.ltype = ltype
        self.f = f
    
    def factored(self,comb_type='ULS',lead=True):
        
        factor = 1.0
        
        # variable loads
        if self.ltype == 'live' or self.ltype == 'wind' or self.ltype == 'snow':
            if comb_type == 'ULS':
                factor = 1.5
            
            if not lead:
                # load is not leading, so it must be multiplied by
                # corresponding factor. These factors apply also in
                # characteristic combation of the serviceability limit state
                if self.ltype == 'live':
                    factor *= 0.7
                elif self.ltype == 'snow':
                    factor *= 0.7
                elif self.ltype == 'wind':
                    factor *= 0.6
        # permanent loads
        elif self.ltype == 'dead':
            if comb_type == 'ULS':
                factor = 1.15
        
        self.f = factor
            

class PointLoad(Load):
    """ Class for point loads (forces and moments) """

    def __init__(self, node, v, load_id=LoadIDs['ULS'], ltype='live',
                 name='PointLoad',f=1.0):
        """        

        Parameters
        ----------
        node : FrameNode object
            Node of point of load application.
        v : numpy array or list
            Vector of the load.
        load_id : TYPE, optional
            Identifier for the load to be used in structural analysis. The default is LoadIDs.ULS.
        ltype : string, optional
            Load type. The default is 'live'.
        name : string, optional
            Name of the load. The default is 'PointLoad'.
        f : float, optional
            Scaling factor. The default is 1.0.

        Returns
        -------
        None.

        """
        
        super().__init__(load_id, ltype, name, f)
    
        
        self.node = node
        self.v = np.asarray(v)        

    @property
    def coordinate(self):
        return self.node.coords

    def add_load(self, fem_model):
        """ Adds load to a FEM model """        
        fem_model.add_load(fem.PointLoad(self.load_id, self.node.fem_node, self.v, self.f))


class LineLoad(Load):
    def __init__(self, member, values, direction, load_id=LoadIDs['ULS'], f=1.0,
                 ltype='live', name='LineLoad',coord_sys="global"):
        """ Class for line loads. Line load can be uniform or trapezoidal
        
            input:
                coordinates .. FrameMember type object that is the subject of the load
                values .. list type object. values[0] and values[1] are the magnitudes
                          of the load at the ends of the member
                direction .. 'x' or 'y'
                load_id .. integer type load identifier
                f .. scaling factor
                ltype .. type of load
                name .. name of load
                coord_sys .. "global" or 'local'. If 'local' the load is given in the
                             local coordinate system of the member.
        """
        super().__init__(load_id, ltype, name, f)
        
        self.member = member
        self.member.loads.append(self)
        
        self.member.has_load = True
        self.member.q0 = values[0]
        self.member.q1 = values[1]            
        
        
        self.values = np.asarray(values)
        self.direction = direction
        self.f = f
        self.element_ids = None
        self.coord_sys = coord_sys

    @property
    def coordinates(self):
        return self.member.coords()

    def calc_k(self):
        """ Evaluate the change in the value of a line load
            per element on a member
        """
        v0, v1 = self.values
        k = (v1 - v0) / self.member.num_elements
        return k

    def add_load(self, fem_model):
        """ Add the load to the FEM model """
        #k = self.calc_k()
        v0, v1 = self.values
        qvals = self.f*np.linspace(v0,v1,self.member.num_elements+1)        
        for i, elem_id in enumerate(self.member.fem_elements):
            #v1 = (i * k) + v0
            load = fem.LineLoad(self.load_id,
                                elem_id,
                                [0.0, 1.0],
                                [qvals[i],qvals[i+1]],                                
                                self.direction,
                                self.coord_sys)
            fem_model.add_load(load)
            #v0 = v1


# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 14:06:38 2021

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Union

class FrameNode:
    """ Class for frame nodes. Used mainly for storing coordinate data """
    
    def __init__(self, *coordinates: Union[list[float, float, float], float]):
        """
        

        Parameters
        ----------
        coordinates : list, optional
            nodal coordinates. The default is [0.0,0.0,0.0].

        Returns
        -------
        None.

        """
        if len(coordinates) == 1 and isinstance(coordinates[0], (list, tuple, np.ndarray)):
            coordinates = coordinates[0]
        
        self.node_id = ""
        self.coords = np.asarray(coordinates)
        self.fem_node = None
        
        self.members = [] # List of members to which the node is connected
        self.symmetry_node = None # Another frame node making a symmetry pair
    
    def __repr__(self):
        if len(self.coords) == 2:
            return f'FrameNode [{self.node_id}]: ({self.x:.3f},{self.y:.3f})'
        else:
            return f'FrameNode [{self.node_id}]: ({self.x:.3f},{self.y:.3f},{self.z:.3f})'
    
    @property
    def x(self):
        """ x coordinate """
        return self.coords[0]
    
    @x.setter
    def x(self,val):
        """ Set x coordinate value """
        self.coords[0] = val
        # Update also FE node, if it exists
        if not self.fem_node is None:
            self.fem_node.x = val
        
        # Update member lengths
        for mem in self.members:
            mem.member.length = mem.length()
        
    @property
    def y(self):
        """ y coordinate """
        return self.coords[1]
    
    @y.setter
    def y(self,val):
        """ Set y coordinate value """
        self.coords[1] = val
        # Update also FE node, if it exists
        if not self.fem_node is None:
            self.fem_node.y = val
        
        # Update connected member lengths
        for mem in self.members:
            mem.member.length = mem.length()
    
    @property
    def z(self):
        """ z coordinate """
        return self.coords[2]
    
    @z.setter
    def z(self,val):
        """ Set z coordinate value """
        self.coords[2] = val
        
        # Update also FE node, if it exists
        if not self.fem_node is None:
            self.fem_node.z = val
        
        # Update member lengths
        for mem in self.members:
            mem.member.length = mem.length()
    
    def plot(self, print_text=True, marker='ok', axes=None):
        """ Plots the node """
        
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        ax.plot(self.x,self.y,marker)
        
        if print_text:
            ax.text(self.x,self.y,self.node_id)


if __name__ == "__main__":
    n1 = FrameNode(0, 0, 0)
    n2 = FrameNode([0, 0, 0])
    n3 = FrameNode(1, 2)
    n4 = FrameNode([1, 2])
    n5 = FrameNode(np.asarray([1, 2]) + np.asarray([3, 5]))

    nodes = (n1, n2, n3, n4, n5)
    for n in nodes:
        print(n.x, n.y)
    print(n1, n2, n3, n4)
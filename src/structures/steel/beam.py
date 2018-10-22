# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 20:13:49 2018

@author: kmela

Class for steel beam

"""

import sys
import math
import numpy as np

import eurocode3
from cross_section import CrossSection

class SteelBeam:
    
    def __init__(self,profile,spans=[6.0e3],supports=["Hinge","X"]):
        """ Constructor
            
            profile -- member profile (of cross_section class)
            spans -- list of spans
            supports -- list of support conditions. Possible values are:
                    "Free" .. no support
                    "Hinge" .. x- and y- displacements fixed
                    "X" .. x-direction fixed
                    "Y" .. y-direction fixed
                    "Rigid" .. all degrees of freedom fixed
            
        """
        
        self.profile = profile
        self.spans = spans
        self.supports = supports
        
        self.forces = []
        self.moments = []
        self.loads = []
        
    
    def nspan(self):
        """ Number of spans """
        return len(self.spans)
    
    def length(self):
        """ Beam length """
        return sum(self.spans)
        
    def add_force(x,val):
        """ Adds a point force to the beam
            input: x .. x-coordinate
                val .. magnitude. The sign gives the direction
        """
        
        
    def mesh(self,nel=10):
        """ Generate FEM model 
            input: nel .. number of elements for each span
        """
        
        """ Generate nodes """
        ycoord = 0.0;
        x = 0.0;
        
        fen = FrameFEM()
        
        fem.add_node(x,ycoord)
        
        for span in self.spans:
            dx = span/nel
            for i in range(nel):
                x += dx
                fem.add_node(x,ycoord)

if __name__ == "__main__":
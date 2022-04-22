# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 08:16:41 2020

@author: kmela
"""
import numpy as np

class CrossSection():
    """ Base class for generic cross-sections
    
    Specific cross sections are subclasses of this class
    
    NOTE: The original purpose of this class is to act as a generic
            cross-section to be used in structural optimization
    
    """

    # Constructor
    def __init__(self,     
                 material=None,
                 A=0,
                 I=[0,0]):
        
        self.material = material
        self.A = A
        self.I = np.asarray(I)       

    def __repr__(self):
        return type(self).__name__
    
    @property
    def E(self):
        return self.material.E
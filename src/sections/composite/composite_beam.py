# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 19:53:07 2019

Composite beam sections

@author: kmela
"""

import math
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

try:
    from src.eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
    from src.sections.steel.steel_section import SteelSection
except:
    from eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
    from eurocodes.en1994 import en1994_1_1
    from sections.steel.steel_section import SteelSection
    

def CompositeIBeam:
    """ Class for composite beams with an I section as
        steel part and a slab (composite or concrete) on top
    """
    
    def __init__(self,steel_part,slab):
        """ Constructor
            Input:
                steel_part: ISection class profile or WI profile steel section
                slab: concrete or composite slab
        """
        
        self.steel = steel_part
        self.slab = slab
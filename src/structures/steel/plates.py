# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 20:01:30 2019

Classes for various steel plates

@author: kmela
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from materials.steel_data import Steel



class RectPlate:
    """ Rectangular plate """
    
    def __init__(self,width,depth,thickness,material="S355"):
        """ Constructor
            Input: Dimensions of the plate and material
        
        """
        
        self.b = width
        self.h = depth
        self.t = thickness
        self.material = Steel(material)
        
    def area(self):
        return self.h*self.b
    
    def end_area(self):
        return self.b*self.t
    
    def side_area(self):
        return self.h*self.t
    
    def volume(self):
        return self.h*self.b*self.t
    
    def weight(self):
        
        return self.density*self.volume()
    
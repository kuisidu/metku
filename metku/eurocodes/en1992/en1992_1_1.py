# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Sun Oct 6 2019

EN 1992-1-1 General rules for concrete structures

@author: kmela
"""

import math

from metku.eurocodes.en1992.constants import gammaC, gammaS, alpha_cc, concrete
from metku.materials.steel_data import Steel


try:
    from metku.eurocodes.en1992.constants import gammaC, gammaS, alpha_cc, concrete
except:
    from eurocodes.en1992.constants import gammaC, gammaS, alpha_cc, concrete
    from materials.steel_data import Steel



def fcd(fck,a_cc=alpha_cc,gC=gammaC):
    """ Design value of compression strength
    """
    return a_cc*fck/gC


class Rebar:
    """ Class for rebar """
    
    def __init__(self,diameter,material="B500B"):
        """ Constructor

            Parameters:
            -----------
                :param diameter: diameter of the rebar (mm)
                :param material: name of the material, to be found in

                Variables:
                ----------
                :ivar values: list of allowable discrete values
                :ivar name: name of the variable (string)
        """
        
        self.d = diameter
        self.mat = Steel(material)
        
        
    @property
    def r(self):
        """ Radius of rebar """
        return 0.5*self.d
    
    @property
    def E(self):
        """ Young's modulus """
        return self.mat.E
    
    @property
    def weight(self):
        """ Weight per unit length (kg/m) """
        return self.mat.rho*self.A*1e3
    
    @property
    def fy(self):
        """ Yield strength """
        return self.mat.fy
    
    @property
    def fyd(self):
        """ Design value of yield strength """
        return self.fy/gammaS
    
    @property
    def fu(self):
        """ Ultimate strength """
        return self.mat.fu
    
    @property
    def A(self):
        """ Cross-sectional area """
        return math.pi*self.d**2/4
    
    @property
    def I(self):
        """ Second moment of area with respect to centroid """
        return 0.25*math.pi*self.r**4
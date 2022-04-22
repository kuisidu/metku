# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Material data

Created on Wed Aug  7 09:23:08 2019

@author: kmela
"""

from math import sqrt

try:
    import metku.eurocodes.en1993.en1993_1_2 as en1993_1_2
except:
    import eurocodes.en1993.en1993_1_2 as en1993_1_2

steels = {"S235": {"fy": 235.0, "fu": 360.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S275": {"fy": 275.0, "fu": 430.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S355": {"fy": 355.0, "fu": 510.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S420": {"fy": 420.0, "fu": 500.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "B500B": {"fy": 500.0, "fu": 550.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S220GD": {"fy": 220.0, "fu": 300.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S250GD": {"fy": 250.0, "fu": 330.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S280GD": {"fy": 280.0, "fu": 360.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S320GD": {"fy": 320.0, "fu": 390.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S350GD": {"fy": 350.0, "fu": 420.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
        "S355MC": {"fy": 355.0, "fu": 430.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S500": {"fy": 500.0, "fu": 550.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S700": {"fy": 700.0, "fu": 750.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "DP1000": {"fy": 590.0, "fu": 980.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "HCT980X": {"fy": 590.0, "fu": 980.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       }

stainless = {'1.4301': {'cold_strip': {'fy': 230, 'fu': 540}, 
                        'hot_strip': {'fy': 210, 'fu': 520},
                        'hot_plate': {'fy': 210, 'fu': 520},
                        'sections': {'fy': 190, 'fu': 500},
                        'E': 200000.0, 'v':0.3, 'rho': 7850e-9}}

class Steel:
    
    def __init__(self,steel_grade="S355"):
        """ Constructor for Steel material class """

        self.name = steel_grade
        self.fy = steels[steel_grade]["fy"]
        self.fu = steels[steel_grade]["fu"]
        self.E = steels[steel_grade]["E"]
        self.nu = steels[steel_grade]["v"]
        self.rho = steels[steel_grade]["rho"]
        self.G = self.E/2/(1+self.nu)
    
    def __repr__(self):
        
        return f"{self.name}"
    
    @property
    def young(self):
        """ Young's modulus: alternative way of calling 'E' """
        return self.E
    
    def eps(self):
        return sqrt(235.0/self.fy)
    
    def Etemp(self,t):
        """ Modification of Young's modulus according to
            EN 1993-1-2
        """
        return en1993_1_2.EaT(t,self.E)
    
    def fy_temp(self,t):
        """ Modification of yield strength according to
            EN 1993-1-2
        """
        return en1993_1_2.fyT(t,self.fy)
    
    def fp_temp(self,t):
        """ Modification of proportional limit according to
            EN 1993-1-2
        """
        return en1993_1_2.fpT(t,self.fy)
    
    def fu_temp(self,t):
        """ Modification of ultimate strength according to
            EN 1993-1-2
        """
        return en1993_1_2.fuT(t,self.fy_temp(t))
        
    
    def stress_temp(self,strain,t,hardening=False):
        """ Calculate stress in steel for given strain at temperature
            't', possiblity with hardenin effect included.
        """
        return en1993_1_2.stress(strain,t,self.E,self.fy,hardening)
    
class StainlessSteel:
    
    def __init__(self,steel_grade="1.4301",product='sections'):
        """ Constructor """ 
        self.name = steel_grade
        self.product = product
        
        self.fy = stainless[steel_grade][product]["fy"]
        self.fu = stainless[steel_grade][product]["fu"]
        self.E = stainless[steel_grade]["E"]
        self.nu = stainless[steel_grade]["v"]
        self.rho = stainless[steel_grade]["rho"]
        self.G = self.E/2/(1+self.nu)
    
    def __repr__(self):
        
        return f"{self.name}"
    
    def eps(self):
        return sqrt(235.0/self.fy*self.E/200000.0)


s = Steel()
s.fy
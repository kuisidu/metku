# -*- coding: utf-8 -*-
"""
Material data

Created on Wed Aug  7 09:23:08 2019

@author: kmela
"""

from math import sqrt

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
       }

class Steel:
    
    def __init__(self,steel_grade="S355"):
        """ Constructor for Steel material class """

        self.name = steel_grade
        self.fy = steels[steel_grade]["fy"]
        self.fu = steels[steel_grade]["fu"]
        self.E = steels[steel_grade]["E"]
        self.nu = steels[steel_grade]["v"]
        self.rho = steels[steel_grade]["rho"]
    
    def __repr__(self):
        
        return f"{self.name}"
    
    def eps(self):
        return sqrt(235.0/self.fy)
        
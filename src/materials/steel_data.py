# -*- coding: utf-8 -*-
"""
Material data

Created on Wed Aug  7 09:23:08 2019

@author: kmela
"""

steels = {"S235": {"fy": 235.0, "fu": 360.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S275": {"fy": 275.0, "fu": 430.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S280GD": {"fy": 280.0, "fu": 360.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S355": {"fy": 355.0, "fu": 510.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
       "S420": {"fy": 420.0, "fu": 500.0, "E": 210000.0, "v": 0.3, "rho": 7850e-9},
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
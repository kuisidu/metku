# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Wed Mar 14 22:50:51 2018

Constants for EN 1992

@author: kmela
"""

# Constants: Young's modulus, Shear modulus, density
# Units: N, mm, i.e. MPa =N/mm^2

# Young's modulus Mpa
E = 210.0e3

# Shear modulus Mpa
G = 81.0e3

# Poisson's ratio MPa
NU = 0.3

# Density kg/mm^3
density = 7850.0e-9

# Partial safety factors
gammaC = 1.5    # Concrete
gammaS = 1.15   # Steel
gammaV = 1.25   # shear stud

# EN 1992-1-1 3.1.6
alpha_cc = 0.85 


""" Material data
fck .. lieriölujuus
fck_cube .. kuutiolujuus
fcm
Ecm
"""
concrete = {"C12/15": {"fck": 12.0, "fck_cube": 15.0, "fcm": 20, "Ecm": 27000.0},
            "C16/20": {"fck": 16.0, "fck_cube": 20.0, "fcm": 24, "Ecm": 29000.0},
            "C20/25": {"fck": 20.0, "fck_cube": 25.0, "fcm": 28, "Ecm": 30000.0},
            "C25/30": {"fck": 25.0, "fck_cube": 30.0, "fcm": 33, "Ecm": 31000.0},
            "C30/37": {"fck": 30.0, "fck_cube": 37.0, "fcm": 38, "Ecm": 33000.0},
            "C35/45": {"fck": 35.0, "fck_cube": 45.0, "fcm": 43, "Ecm": 34000.0},
            "C40/50": {"fck": 40.0, "fck_cube": 50.0, "fcm": 48, "Ecm": 35000.0},
            }

class Concrete:
    
    def __init__(self,material="C25/30"):
        """ Constructor for Concrete material class """

        self.name = material
        self.fck = concrete[material]["fck"]
        self.fck_cube = concrete[material]["fck_cube"]
        self.fcm = concrete[material]["fcm"]
        self.Ecm = concrete[material]["Ecm"]
        
    def fcd(self,a_cc=alpha_cc,gC=gammaC):
        """ Design value for compressive strength """
        return a_cc*self.fck/gC
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 2019

EN 1994-1-1 Design of composite steel and concrete structures: Part 1-1 
General rules and rules for buildings

@author: kmela
"""

import math

try:
    from src.eurocodes.en1992.constants import gammaC, gammaS, alpha_cc, concrete
except:
    from eurocodes.en1992.constants import gammaC, gammaS, alpha_cc, concrete
    

S3L_Studs = {19: {"D": 19.0, "A": 9.52, "H": 31.75},
            22: {"D": 22.0, "A": 9.52, "H": 34.925},
            25: {"D": 19.0, "A": 12.7, "H": 41.275},
            }


def ShearStud:
    """ Class for shear connectors """
    
    def __init__(self):
        """ Constructor """
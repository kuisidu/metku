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
    from eurocodes.en1994.constants import gammaV
    

S3L_Studs = {19: {"D": 19.0, "A": 9.52, "H": 31.75},
            22: {"D": 22.0, "A": 9.52, "H": 34.925},
            25: {"D": 19.0, "A": 12.7, "H": 41.275},
            }


class ShearStud:
    """ Class for shear connectors """
    
    def __init__(self,stud=S3L_Studs[19],length=105.0,fy=350,fu=450):
        """ Constructor """
        
        self.d = stud["D"]
        self.a = stud["A"]
        self.h = stud["H"]
        self.length = length
        self.fy = fy
        self.fu = fu
        
    def PRd(self,fck,Ecm,verb=False):
        """ Shear force resistance of a stud """
        
        PRd_steel = 0.8*self.fu*math.pi**2*self.d**2/4
        
        if self.length/self.d > 4:
            alpha = 1.0
        else:
            alpha = 0.2*(self.length/self.d+1)
        
        PRd_concrete = 0.29*alpha*self.d**2*math.sqrt(fck*Ecm)/gammaV
        
        if verb:
            print("PRd(steel) = {0:4.2f} [kN]".format(PRd_steel*1e-3))
            print("PRd(concrete) = {0:4.2f} [kN]".format(PRd_concrete*1e-3))
            print("d = {0:4.2f} [MPa]".format(self.d))
            print("fck = {0:4.2f} [MPa]".format(fck))
            print("Ecm = {0:4.2f} [MPa]".format(Ecm))
            print("alpha = {0:4.2f} [MPa]".format(alpha))
            
        
        
        return min(PRd_steel,PRd_concrete)
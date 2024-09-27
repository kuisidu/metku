# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 23:12:31 2018

Circular hollow sections

@author: kmela
"""

import math


from metku.eurocodes.en1993 import en1993_1_1
from metku.sections.steel.steel_section import SteelSection

import warnings
warnings.simplefilter("once", UserWarning)
class CHS(SteelSection):
    """ Circular hollow sections """

    def __init__(self,D,T,fy=355):
        """ Circular hollow section
            D -- diameter
            T -- wall thickness
            fy -- yield strength
        """
        self.D = D
        self.T = T
        
        A = self.area()
        Au = self.paint_area()
        I = self.second_moment()
        Ashear = self.shear_area()
        Wel = self.elastic_modulus()
        Wpl = self.plastic_modulus()
        self.It = self.torsion_modulus()
        self.Iw = 0.0

        super().__init__(fy,A,I,Au,Wpl,Wel,Ashear)
        

        self.imp_factor = [en1993_1_1.buckling_curve["c"],en1993_1_1.buckling_curve["c"]]

    def __repr__(self):
        return f"{type(self).__name__} {self.D:.0f}X{self.T:.1f}"

    def area(self):
        r = 0.5*self.D
        return math.pi*r**2 - math.pi*(r-self.T)**2

        
    def paint_area(self):
        return math.pi*self.D

    def shear_area(self):
        return 2*self.area()/math.pi


    def second_moment(self):
        r2 = 0.5*self.D
        r1 = r2-self.T
        return 0.25*math.pi*(r2**4-r1**4)

    def elastic_modulus(self):
        return self.second_moment()/(0.5*self.D)

    def plastic_modulus(self):
        D2 = self.D
        D1 = self.D-2*self.T
        return (D2**3-D1**3)/6

    
    def torsion_modulus(self):
        Iv2 = math.pi*self.D**4/32
        Iv1 = math.pi*(self.D-2*self.T)**4/32
        Iv = Iv2-Iv1
        return Iv

    def web_class_bend(self, verb=False):
        r = self.D/self.T
        cWeb = en1993_1_1.CHS_class(r,self.eps)
        return cWeb

    def web_class_comp(self, verb=False):
        r = self.D/self.T
        cWeb = en1993_1_1.CHS_class(r,self.eps)
        return cWeb
        
    def web_class_comp_bend(self, verb=False):
        r = self.D/self.T
        cWeb = en1993_1_1.CHS_class(r,self.eps)
        return cWeb

    def flange_class(self, verb=False):
        return 1
    
    def bending_resistance(self, C=0, verb=False):
        # Bending resistance, Nmm
        if C == 0:
            C = self.section_class()

        if C < 3:
            WRd = self.Wpl
        elif C == 3:
            WRd = self.Wel
        else:
            #raise NotImplementedError(f"Cross-section class 4 not implemented in: {self}")
            warnings.warn(f"!Cross-section class 4 not implemented for CHS profiles! "
                  f"\nUsing elastic bending resistance, check results accordingly!", UserWarning, stacklevel=2)

            WRd = self.Wel

        MRd = self.code.bending_resistance(WRd, self.fy)

        if verb:
            print("MRd = {0:4.2f} kNm".format(MRd * 1e-6))

        return MRd, C
    
    def MomentAxialForceInteract(self):
        """ This interaction formula is only for plastic design """
        MRd = self.BendingResistance(1);
        NRd = self.AxialForceResistance;
            
        """ EN 1993-1-1: 6.2.9.1(5)
            In class 1 or 2 use plastic design
        """
        n = abs(self.Ned)/NRd;

        if n > 1.0:
            MNRd = 0
        else:
            MNRd = MRd*min(1-n**1.7,1)

        return MNRd
    
    
def CHS_area(D,T):
    r = 0.5*D
    A = math.pi*r**2 - math.pi*(r-T)**2
    return A
        
def CHS_paint_area(D):
    Au = math.pi*D
    return Au

def CHS_shear_area(A):
    Ashear = 2*A/math.pi
    return Ashear

def CHS_second_moment(D,T):
    r2 = 0.5*D
    r1 = r2-T
    I = 0.25*math.pi*(r2**4-r1**4)
    return I

def CHS_elastic_modulus(D,I):
    Wel = I/(0.5*D)
    return Wel

def CHS_plastic_modulus(D,T):
    D2 = D
    D1 = D-2*T
    Wpl = (D2**3-D1**3)/6
    return Wpl

def CHS_torsion_modulus(D,T):
    Iv2 = math.pi*D**4/32
    Iv1 = math.pi*(D-2*T)**4/32
    Iv = Iv2-Iv1
    return Iv

class CustomCHS(CHS):
    """ Class for custom CHS profiles, mainly for optimization """
    
    def __init__(self,D,T,fy=355):
        
        super().__init__(D,T,fy)
    
    def __getattribute__(self, name):
        """ override the attribute access for those attributes
            that depend on the section dimensions and that are
            constant for rolled sections
        """
        if name == "A":
            if self.area() < 0:
                print(self.D, self.T)
                
            return self.area()
        elif name == "I":
            I = self.second_moment()
            #if I[1] < 0:
            #    print(self.h, self.b, self.tt, self.tb, self.tw)
            return self.second_moment()
        elif name == "Wel":
            return self.section_modulus()
        elif name == "Wpl":
            return self.plastic_section_modulus()
        elif name == "Ashear":
            return self.shear_area()
        elif name == "Au":
            return self.paint_area()
        elif name == "It":
            return self.torsion_modulus()
        else:
            return super().__getattribute__(name)
    
    
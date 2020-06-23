# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 23:12:31 2018

Rectangular hollow sections

@author: kmela
"""

import math
try:
    from src.eurocodes.en1993 import en1993_1_1
    from src.sections.steel.steel_section import SteelSection
except:
    from eurocodes.en1993 import en1993_1_1
    from sections.steel.steel_section import SteelSection

class CHS(SteelSection):
    """ Circular hollow sections """

    def __init__(self,D,T,fy=355):
        """ Circular hollow section
            D -- diameter
            T -- wall thickness
            fy -- yield strength
        """
        A = CHS_area(D,T)
        Au = CHS_paint_area(D)
        I = CHS_second_moment(D,T)
        Ashear = CHS_shear_area(A)
        Wel = CHS_elastic_modulus(D,I)
        Wpl = CHS_plastic_modulus(D,T)
        It = CHS_torsion_modulus(D,T)

        SteelSection.__init__(self,fy,A,I,Au,Wpl,Wel,Ashear)
        self.D = D
        self.T = T
        self.It = It
        self.Iw = 0.0

        self.imp_factor = [en1993_1_1.buckling_curve["c"],en1993_1_1.buckling_curve["c"]]

    def web_class_bend(self):
        r = self.D/self.T
        cWeb = en1993_1_1.CHS_class(r,self.eps)
        return cWeb

    def web_class_comp(self):
        r = self.D/self.T
        cWeb = en1993_1_1.CHS_class(r,self.eps)
        return cWeb
        
    def web_class_comp_bend(self):
        r = self.D/self.T
        cWeb = en1993_1_1.CHS_class(r,self.eps)
        return cWeb

    def flange_class(self):
        return 1
        
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

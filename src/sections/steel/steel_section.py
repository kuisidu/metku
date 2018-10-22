# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 22:42:48 2018

Class for general steel cross sections 
    NOTE: This could be a package with the following structure
        cross_section
            hollow_sections
                RHS
                SHS
                CHS
            I_sections


@author: kmela
"""

import math

from eurocodes.en1993 import en1993_1_1
from eurocodes.en1993 import constants

INFEASIBLE = 999

class SteelSection:
    """ General cross section class
    
    Specific cross sections are subclasses of this class
    
    """
    
    """
    methods (Abstract = true)
        cWeb = WebClassBend(obj)
        cWeb = WebClassComp(obj)
        cWeb = WebClassCompBend(obj)
        cFlange = FlangeClass(obj)
        MNRd = MomentAxialForceInteract(obj,UN,MRd)
    end
    """
    
    # Constructor
    def __init__(self,fy,A,I,Au,Wpl,Wel,Ashear,Ned=0.0,Med=0.0,Ved=0.0):
        self.fy = fy
        self.eps = en1993_1_1.epsilon(self.fy)
        self.E = constants.E
        self.A = A
        self.I = I
        self.Au = Au
        self.Wpl = Wpl
        self.Wel = Wel
        self.Ashear = Ashear
        self.Ned = Ned
        self.Med = Med
        self.Ved = Ved
        self.imp_factor = None
    
    def weight(self):
        """ Weight per unit length """
        w = self.A*constants.density
        return w
    
    def stress_ratio(self):
        """ computes the stress ratio Psi required
        in cross-section classification of parts in
        compression and bending
        See. EN 1993-1-1, Table 5.2
        NOTE: compression is positive, thus
        the '-' sign in front of Ned.
        """
        
        stmax = -self.Ned/self.A - self.Med/self.Wel[0]
        stmin = -self.Ned/self.A + self.Med/self.Wel[1]
        
        p = stmin/stmax
        
        return p
        
    def section_class(self):
        """ Determines cross-section class """
            
        C = 1
        type = "Tension"
        
        if abs(self.Ned) < 1e-4 and abs(self.Med) > 1e-4:
            # Pure bending
            type = "Bending"
            Cweb = self.web_class_bend()
            Cflange = self.flange_class()
            C = max(Cweb,Cflange)
        elif self.Ned < 0.0:
            # Compression
            if abs(self.Med) < 1e-4:
                # Pure compression
                type = "Compression"
                Cweb = self.web_class_comp()
                Cflange = self.flange_class()
                C = max(Cweb,Cflange)
            else:
                # Bending and compression
                type = "Bending and compression"
                
                # Flange is in compression
                Cflange = self.flange_class()
                # Classify web as a part in compression and bending
                Cweb = self.web_class_comp_bend(self.Ned)
                C = max(Cweb,Cflange)
        elif self.Ned > 0.0 and abs(self.Med) > 1e-4:
            # Axial tension and bending
            type = "Tension and bending"
            C = 1
        
        return C


    def axial_force_resistance(self):        
        if self.Ned >= 0:
            NRd = en1993_1_1.tension_resistance(self.A,self.fy)
        else:
            NRd = en1993_1_1.compression_resistance(self.A,self.fy)
            
        return NRd                    
            
    def bending_resistance(self,C=0,axis="y"):
        # Bending resistance, Nmm
        if C == 0:
            C = self.section_class()
            
        if axis == "y":
            n = 0
        else:
            n = 1
            
        if C < 3:                
            WRd = self.Wpl[n]
        elif C == 3:
            WRd = self.Wel[n]
        else:
            WRd = self.Wel[n]
                
        MRd = en1993_1_1.bending_resistance(WRd,self.fy)
            
        return MRd, C
    
    def elastic_bending_resistance(self,axis="y"):
        """ Elastic bending resistance """
        if axis == "y":
            n = 0
        else:
            n = 1
        return en1993_1_1.bending_resistance(self.Wel[n],self.fy)
    
    def plastic_bending_resistance(self,axis="y"):
        """ Plastic bending resistance """
        if axis == "y":
            n = 0
        else:
            n = 1
        return en1993_1_1.bending_resistance(self.Wpl[n],self.fy)
        
    def shear_force_resistance(self,C=0):
        if C == 0:
            C = self.section_class()
                
        if C < 3:
            VRd = en1993_1_1.shear_resistance(self.Ashear,self.fy)
        elif C == 3:
            VRd = en1993_1_1.shear_resistance(self.Ashear,self.fy)
        else:
            VRd = en1993_1_1.shear_resistance(self.Ashear,self.fy)
            
        return VRd
        
    def moment_axial_force_interact(obj,UN,MRd):        
        return None
        
        """
        def a = increase_M_and_N(self):
            "" Increases bending moment and axial force proportionally,
            until cross-section resistance is obtained. Used for classification
            
            Needs to be modified for Python!
            ""
            
            a0 = 1.0
            
            #%a = fzero(@MNInc,a0)
            aIter = []
            amax = 5
            while isempty(aIter)
                amax = 2*amax
                [f,aIter] = secant(@MNInc,[1.0,amax])
                   
            a = aIter(end)                        
        
        
        def res = MNInc(k):
            selfN.Med = k*self.Med
            selfN.Ned = k*self.Ned
                
            MNRd = selfN.MomentAxialForceInteract
            if MNRd > 0:
                res = selfN.Med/MNRd-1
            else
                res = 1
        """
            
    def section_resistance(self,class1or2=False, return_list=False):
        """ Calculates resistance of cross-section
            Checks the following:
                Axial force
                Shear force
                Bending moment
                Interaction of axial force and bending moment
            output: r .. maximum of different utilization ratios
        """
        
        """ Classify section """
        if class1or2:
            C = 1
        else:
            C = self.section_class()
            
        
        NRd = self.axial_force_resistance()
        UN = abs(self.Ned)/NRd
                    
        VRd = self.shear_force_resistance()
        UV = abs(self.Ved)/VRd
                
        MRd, C = self.bending_resistance(C)
        UM = abs(self.Med)/MRd
                         
        if C < 3:
            MNRd = self.moment_axial_force_interact(UN,MRd)
            if MNRd > 0.0:
                UMN = abs(self.Med)/MNRd
            else:
                UMN = INFEASIBLE
        else:
            UMN = UN+UM
        
        if return_list:
            r = [UN,UV,UM,UMN]
        else:
            r = max([UN,UV,UM,UMN])
        return r
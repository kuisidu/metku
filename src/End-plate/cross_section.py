""" Class for general steel cross sections """

import sys
import math

import eurocode3


class CrossSection:
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
    # Constructor
    def __init__(self,fy,A,I,Au,Wpl,Wel,Ashear,Ned=0,Med=0,Ved=0):
        self.fy = fy
        self.eps = eurocode3.epsilon(self.fy)
        self.E = eurocode3.young
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
        
        if self.Ned == 0 & self.Med != 0:
            # Pure bending
            type = "Bending"
            Cweb = self.web_class_bend
            Cflange = self.flange_class
            C = max(Cweb,Cflange)
        elif self.Ned < 0:
            # Compression
            if self.Med == 0:
                # Pure compression
                type = "Compression"
                Cweb = self.web_class_comp
                Cflange = self.flange_class
                C = max(Cweb,Cflange)
            else:
                # Bending and compression
                type = "Bending and compression"
                
                # Flange is in compression
                Cflange = self.flange_class
                # Classify web as a part in compression and bending
                Cweb = self.web_class_comp_bend
                C = max(Cweb,Cflange)
        elif self.Ned > 0 & self.Med != 0:
            # Axial tension and bending
            type = "Tension and bending"
            C = 1
        
        return C


        def axial_force_resistance(self):
            if self.Ned >= 0:
                NRd = self.tension_resistance
            else:
                NRd = self.compression_resistance

            return NRd

        def tension_resistance(self):
            return self.fy*self.A


        def compression_resistance(self):
            return self.fy*self.A

        def bending_resistance(self,C=0,axis=1):
            # Bending resistance, Nmm
            if C == 0:
                C = self.section_class

            if C < 3:
                WRd = self.Wpl(axis)
            elif C == 3:
                WRd = self.Wel(axis)
            else:
                WRd = self.Wel(axis)

            MRd = self.fy*WRd

            return MRd, C

        def shear_force_resistance(self,C=0):
            if C == 0:
                C = self.SectionClass

            if C < 3:
                VRd = self.fy*self.Ashear/math.sqrt(3)
            elif C == 3:
                VRd = self.fy*self.Ashear/math.sqrt(3)
            else:
                VRd = self.fy*self.Ashear/math.sqrt(3)

            return VRd
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

        def section_resistance(self,verb=0):
            """ In MATLAB: Check
            return r
            """

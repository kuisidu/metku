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
    # Original initial parameters (self,fy,A,I,Au,Wpl,Wel,Ashear)
    def __init__(self,fy,h,b,t_f,A,I,Au,Wpl,Wel,Iw,It,Ashear,Ned=0.0,Med=0.0,Ved=0.0):
        self.fy = fy
        self.eps = eurocode3.epsilon(self.fy)
        self.E = eurocode3.young
        self.h = h
        self.b = b
        self.t_f = t_f
        self.A = A
        self.I = I
        self.Au = Au
        self.Wpl = Wpl
        self.Wel = Wel
        self.Iw = Iw
        self.It = It
        self.Ashear = Ashear
        self.Ned = Ned
        self.Med = Med
        self.Ved = Ved
        self.imp_factor = [0,0]

    
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
    
    
    def define_imp_factor(self):
        """ Returns imperfection factor alpha for buckling
            Stores values in array [y-y, z-z]
        """         
        if self.h/self.b > 1.2:
            if self.t_f <= 40:
                if self.fy < 460:
                    self.imp_factor[0] = eurocode3.buckling_curve['a']
                    self.imp_factor[1] = eurocode3.buckling_curve['b']
                else:
                    self.imp_factor[0] = eurocode3.buckling_curve['a']
                    self.imp_factor[1] = eurocode3.buckling_curve['b']
            elif self.t_f > 40 and self.t_f <= 100:
                if self.fy < 460:
                    self.imp_factor[0] = eurocode3.buckling_curve['b']
                    self.imp_factor[1] = eurocode3.buckling_curve['c']
                else:
                    self.imp_factor[0] = eurocode3.buckling_curve['a']
                    self.imp_factor[1] = eurocode3.buckling_curve['a']
        
        else:
            if self.t_f <= 100:
                if self.fy < 460:
                    self.imp_factor[0] = eurocode3.buckling_curve['b']
                    self.imp_factor[1] = eurocode3.buckling_curve['c']
                else:
                    self.imp_factor[0] = eurocode3.buckling_curve['a']
                    self.imp_factor[1] = eurocode3.buckling_curve['a']
            else:
                if self.fy < 460:
                    self.imp_factor[0] = eurocode3.buckling_curve['d']
                    self.imp_factor[1] = eurocode3.buckling_curve['d']
                else:
                    self.imp_factor[0] = eurocode3.buckling_curve['c']
                    self.imp_factor[1] = eurocode3.buckling_curve['c']
      
    
    
    def define_imp_factor_LT(self):
        """
            Currently works only with hot-rolled sections
            TODO: Add coefficients for welded sections
        """        
        if self.h / self.b <= 2:
            self.ImperfectionLT = eurocode3.buckling_curve['a']
        else:
            self.ImperfectionLT = eurocode3.buckling_curve['b']
    
    def flange_class(self):
        """ TODO """
        return 3
    def web_class_bend(self):
        """ TODO """
        return 3
    def web_class_comp(self):
        """ TODO """
        return 3
    
    def web_class_comp_bend(self):
        """ TODO """
        return 3
    
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
                Cweb = self.web_class_comp_bend()
                C = max(Cweb,Cflange)
        elif self.Ned > 0.0 and abs(self.Med) > 1e-4:
            # Axial tension and bending
            type = "Tension and bending"
            C = 1
        
        return C


    def axial_force_resistance(self):        
        if self.Ned >= 0:
            NRd = self.tension_resistance()
        else:
            NRd = self.compression_resistance()
            
        return NRd
        
    def tension_resistance(self):
        return self.fy*self.A
            
        
    def compression_resistance(self):
            return self.fy*self.A
            
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
                
        MRd = self.fy*WRd
            
        return MRd, C
        
    def shear_force_resistance(self,C=0):
        if C == 0:
            C = self.section_class()
                
        if C < 3:
            VRd = self.fy*self.Ashear/math.sqrt(3)
        elif C == 3:
            VRd = self.fy*self.Ashear/math.sqrt(3)
        else:
            VRd = self.fy*self.Ashear/math.sqrt(3)
            
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
            
    def section_resistance(self,class1or2=False):
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
            MNRd, C = self.bending_resistance(C)
            UMN = abs(self.Med)/MNRd
        else:
            UMN = UN+UM
                
        r = max([UN,UV,UM,UMN])
        return r
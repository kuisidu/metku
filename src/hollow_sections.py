""" Hollow sections class """

import sys
import math

import eurocode3
from cross_section import CrossSection

class RHS(CrossSection):
    """ Rectangular hollow sections """

    def __init__(self,H,B,T,fy=355):
        """ Rectangular hollow sections
            
            H -- height
            B -- width
            T -- wall thickness
            fy -- yield strength
        """
        R = RHS_outer_radius(T)
        A = RHS_area(H,B,T,R)
        I = [0.0,0.0]
        I[0] = RHS_second_moment(H,B,T,R)
        I[1] = RHS_second_moment(B,H,T,R)
        Ashear = RHS_shear_area(A,H,B)
        Au = RHS_paint_area(H,B,R)
        Wel = [0.0,0.0]
        Wel[0] = RHS_elastic_modulus(H,I[0])
        Wel[1] = RHS_elastic_modulus(H,I[1])
        Wpl = [0.0,0.0]
        Wpl[0] = RHS_plastic_modulus(H,B,R,T)
        Wpl[1] = RHS_plastic_modulus(B,H,R,T)
        
        CrossSection.__init__(self,fy,A,I,Au,Wpl,Wel,Ashear)
        self.H = H
        self.B = B
        self.T = T
        self.R = R
        
        self.imp_factor = [eurocode3.buckling_curve["c"],eurocode3.buckling_curve["c"]]

    def flange_class(self):
        """ Determine class of compressed flange """
        cf = self.B-2*self.R
        rf = cf/self.T
        cFlange = eurocode3.internal_part_in_compression(rf,self.eps)
        
        return cFlange

    def web_class_comp(self):
        """ Determine class of compressed web """
        cw = self.H-2*self.R
        rw = cw/self.T
        cWeb = eurocode3.internal_part_in_compression(rw,self.eps)
        
        return cWeb
        
    def web_class_bend(self):
        """ Determine class of web in bending """
        cw = self.H-2*self.R
        rw = cw/self.T
        cWeb = eurocode3.internal_part_in_bending(rw,self.eps)

        return cWeb
	
    def moment_axial_force_interact(obj,UN,MRd):        
        """ Interaction rule for section resistance for combined
            axial force and bending
            
            input: UN .. NEd/NRd
                  MRd .. moment resistance
        """
        aw = min((self.A-2*self.B*self.T)/self.A,0.5)   
                        
        if UN > 1.0:
            MNRd = 0.0
        else:
            MNRd = MRd*min((1-UN)/(1-0.5*aw),1)
            
        return MNRd

class SHS(RHS):
    """ Square hollow sections
        Subclass of RHS
    """
    def __init__(self,H,T,fy=355):
        RHS.__init__(self,H,H,T,fy)


# Functions for computing RHS section properties
def RHS_area(H,B,T,R):
    """ cross-sectional area """
    A = 2*(H-2*R)*T+2*(B-2*R)*T + math.pi*(R**2-(R-T)**2)
    return A

def RHS_paint_area(H,B,R):
    """ circumference """
    Au = 2*(H+B)-8*R+2*math.pi*R
        
def RHS_shear_area(A,H,B):
    Ashear = A*H/(H+B)
    return Ashear

def RHS_second_moment(H,B,T,R):
    """ Second moment of area """
    I1 = 2*(1/12*T*(H-2*R)**3+1/12*(B-2*R)*T**3+(0.5*(H-T))**2*(B-2*R)*T)
    IR = 4*((math.pi/16-4/9/math.pi)*R**4+(0.5*H-R+4*R/3/math.pi)**2*math.pi*R**2/4)
    IRneg = 4*((math.pi/16-4/9/math.pi)*(R-T)**4+(0.5*H-R+4*(R-T)/3/math.pi)**2*math.pi*(R-T)**2/4)
    I1 = I1+IR-IRneg

    return I1
    
def RHS_elastic_modulus(H,I1):
    Wel = I1/(0.5*H)

    return Wel
        
def RHS_plastic_modulus(H,B,R,T):
    Wp = 2*((0.5*H-R)*T*(0.5*H-R)) + (B-2*R)*T*(H-T)
    WR = 2*math.pi/4*R**2*2*(0.5*H-R+4*R/3/math.pi)
    Wneg = 2*math.pi/4*(R-T)**2*2*(0.5*H-R+4*(R-T)/3/math.pi)
    Wp = Wp+WR-Wneg
    
    return Wp
        
def RHS_outer_radius(T):
    if T <= 6.0:
        R = 2*T
    elif T <= 10.0:
        R = 2.5*T
    else:
        R = 3*T

    return R

class CHS(CrossSection):
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

        CrossSection.__init__(self,fy,A,I,Au,Wpl,Wel,Ashear)
        self.D = D
        self.T = T
        self.It = It

        self.imp_factor = [eurocode3.buckling_curve["c"],eurocode3.buckling_curve["c"]]

    def web_class_bend(self):
        r = self.D/self.T
        cWeb = eurocode3.CHS_class(r,self.eps)
        return cWeb

    def web_class_comp(self):
        r = self.D/self.T
        cWeb = eurocode3.CHS_class(r,self.eps)
        return cWeb
        
    def web_class_comp_bend(self):
        r = self.D/self.T
        cWeb = eurocode3.CHS_class(r,self.eps)
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

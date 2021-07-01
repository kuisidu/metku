# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-1 General rules for steel structures

@author: kmela
"""

import math

try:
    from metku.eurocodes.en1993.constants import gammaM0, gammaM1, gammaM2
except:
    from eurocodes.en1993.constants import gammaM0, gammaM1, gammaM2

#                 a0     a    b      c     d 
buckling_curve = {"a0": 0.13, "a": 0.21, "b": 0.34, "c": 0.49, "d": 0.76}

def sway_imperfection(h,m):
    """ Global sway imperfection of 5.3.2 """
    phi0 = 1/200
    
    alpha_h = min(max(2/3,2/math.sqrt(h)),1.0)
    alpha_m = math.sqrt(0.5*(1+1/m))
    
    phi = phi0*alpha_h*alpha_m
    
    return phi

def epsilon(fy):
    """ Computes the epsilon required in various design rules 
        input: fy -- yield strength in MPa
    
    """

    e = math.sqrt(235.0/fy)
    return e

def internal_part_in_bending(r,e):
    """ Classification of internal parts of cross-sections in bending
    
    EN 1993-1-1, Table 5.2, part 1
    r = c/t
    e = epsilon
    """

    # Initialize class
    C = 4
      
    if r <= 72.0*e:
        C = 1
    elif r <= 83.0*e:
        C = 2
    elif r <= 124.0*e:
        C = 3
        
    return C

def internal_part_in_compression(r,e):
    """ Classification of internal parts of cross-sections in compression
    
    EN 1993-1-1, Table 5.2, part 1
    r = c/t
    e = epsilon
    
    """
    
    # If e > 1, it means that e is actually yield strength
    if e > 1.0:
       e = epsilon(e)

    C = 4

    if r <= 33.0*e:
        C = 1
    elif r <= 38.0*e:
        C = 2
    elif r <= 42.0*e:
        C = 3   
            
    return C 
               
def internal_part_comp_bend(r,e,a,p):
    """ Classification of internal parts of cross-sections in compression and bending

    EN 1993-1-1, Table 5.2, part 1
    r = c/t
    e = epsilon
    a = alpha: location of plastic neutral axis
    p = psi: location of elastic neutral axis
    
    """

    # If e > 1, it means that e is actually yield strength
    if e > 1.0:
        e = epsilon(e)

    C = 4

    # Check first for class 1 or 2
    if a == 0:
        C = 1
    elif a > 0.5:        
        if r <= 396.0*e/(13*a-1):
            C = 1
        elif r <= 456.0*e/(13*a-1):
            C = 2
    else:
        if r <= 36.0*e/a:
            C = 1
        elif r <= 41.5*e/a:
            C = 2

    if C > 2:
        # Check for class 3
        if p > -1:
            if r <= 42.0*e/(0.67+0.33*p):
                C = 3
        elif r <= 62.0*e*(1-p)*math.sqrt(-p):
            C = 3
            
    return C

def outstand_part_in_compression(r,e):
    """ Classification of outstand parts of cross-sections in compression

    EN 1993-1-1, Table 5.2, part 2
    r = c/t
    e = epsilon
    
    """
    
    # If e > 1, it means that e is actually yield strength
    if e > 1.0:
        e = epsilon(e)
    
    C = 4

    if r <= 9.0*e:
        C = 1
    elif r <= 10.0*e:
        C = 2
    elif r <= 14.0*e:
        C = 3
         
    return C

def CHS_class(r,e):
    """ Classification circular hollow sections
    
    EN 1993-1-1, Table 5.2, part 2
    r = d/t
    e = epsilon
    
    """

    C = 4

    if r <= 50*e**2:
        C = 1
    elif r <= 70*e**2:
        C = 2
    elif r <= 90*e**2:
        C = 3

    return C

def tension_resistance(A,fy):
    """ EN 1993-1-1, Eq. (6.6) 
    """
    NplRd = fy*A/gammaM0
    return NplRd

def net_section_resistance(Anet,fu):
    """ EN 1993-1-1, Eq (6.7)
    """
    NuRd = 0.9*Anet*fu/gammaM2
    return NuRd

def compression_resistance(A,fy):
    """ EN 1993-1-1, Eq (6.10) and (6.11)
        input:
            A .. gross cross section in class 1, 2, or 3
                 effective cross section in class 4
    """
    NcRd = A*fy/gammaM0
    return NcRd

def bending_resistance(WRd,fy):
    """ EN 1993-1-1, Eqs. (6.13) to (6.15)
        input:
            WRd .. section modulus
                   plastic modulus for class 1 or 2
                   elastic modulus for class 3
                   effective elastic modulus for class 4                   
    """
    McRd = WRd*fy/gammaM0
    return McRd

def shear_resistance(Ashear,fy):
    """ EN 1993-1-1, Eq. (6.18)
        input:
            Ashear .. shear area
    """
    VcRd = Ashear*fy/math.sqrt(3)/gammaM0
    return VcRd

def slenderness(A,fy,Ncr):
    """ Non-dimensional slenderness
        input:
            Ncr .. elastic buckling force
    """
    return math.sqrt(A*fy/Ncr)

def buckling_reduction_factor(slend,imp_factor=buckling_curve["c"]):
    """ Reduction factor for flexural buckling
        EN 1993-1-1, Eq. (6.49)
    """
    p = 0.5*(1+imp_factor*(slend-0.2)+slend**2)
    chi = min(1/(p+math.sqrt(p**2-slend**2)),1.0)            
    
    return chi

def buckling_strength(A,fy,chi):
    """ EN 1993-1-1, Eqs. (6.47) and (6.48)
    """
    NbRd = chi*A*fy/gammaM1
    return NbRd

#def beam_column_utility(NEd,NRd,chi,MyEd,MyRd,ky,MzEd=0.0,MzRd=1.0,chiLT=1.0):
def beam_column_utility(UN,UMy,ky,chiLT=1.0,UMz=0.0):
    """ EN 1993-1-1, Eqs. (6.61) and (6.62) 
        input:
            UN = NEd/NRd/chi
            UMy = MyEd/MyRd
            UMz = MzEd/MzRd
                
    """
    Ut = UN/gammaM1 + ky*UMy/chiLT/gammaM1 + UMz/gammaM1
    return Ut

def kyy(UN,slend,Cmy=1.0,section_class=2):
    """ EN 1993-1-1 Table B.1, "Method B"
    """
    if section_class < 3:
        ky = Cmy*(1+min(slend-0.2,0.8)*UN/gammaM1)
    else:
        ky = Cmy*(1+min(slend,0.6)*UN/gammaM1)
    
    return ky

def kzy(ky=0.0,UN=0.0,slend=1.0,CmLT=1.0,section_class=2,susceptible_to_torsion=False):
    """ EN 1993-1-1, Table B.1
    """
    if section_class < 3:
        if susceptible_to_torsion:
            if slend < 0.4:
                kzy = min(0.6+slend,1-0.1*slend/(CmLT-0.25)*UN/gammaM1)
            else:
                kzy = 1-0.1*min(slend,1)/(CmLT-0.25)*UN/gammaM1
        else:
            kzy = 0.6*ky
    else:
        if susceptible_to_torsion:
            kzy = 1-0.05*min(slend,1)/(CmLT-0.25)*UN/gammaM1
        else:
            kzy = 0.8*ky
    
    return kzy        

def kzz(UN,slend,Cmz,section_class=2,profile="hollow"):
    """ EN 1993-1-1 Table B.1
    """
    if section_class < 3:
        if profile == "hollow":
            kz =  Cmz*(1+min(slend-0.2,0.8)*UN/gammaM1)
        else:
            kz =  Cmz*(1+min(2*slend-0.6,1.4)*UN/gammaM1)
    else:
        kz = Cmz*(1+min(slend,0.6)*UN/gammaM1)
    
    return kz
    

def kyz(kz,section_class=2):
    """ EN 1993-1-1 Table B.1 and B.2
    """
    if section_class < 3:
        ky = 0.6*kz
    else:
        ky = kz
    return ky

def equivalent_moment_factor(M,load="uniform"):
    """ EN 1993-1-1, Table B.3
        input:
            M .. array
                 M[0] .. moment at end 1
                 M[1] .. moment at mid-span
                 M[2] .. moment at end 2
            load .. "uniform" or "point"
    """
        
                         
    
""" Eurocode 3: General expressions for steel structures

Functions:

epsilon
shear_eta
internal_part_in_bending
internal_part_in_compression
internal_part_comp_bend
outstand_part_in_compression

"""
import sys
import math

# Constants: Young's modulus, Shear modulus, density
young = 210e3
G = 81e3
density = 7850e-9

#                 a0     a    b      c     d 
buckling_curve = {"a0": 0.13, "a": 0.21, "b": 0.34, "c": 0.49, "d": 0.76}

def epsilon(fy):
    """ Computes the epsilon required in various design rules """

    e = math.sqrt(235.0/fy)
    return e
    
def shear_eta(fy):
    """ Computes the 'eta' according to EN 1993-1-5, clause 5.1(2) """

    eta = 1.0
    
    if fy <= 460.0:
        eta = 1.2
    
    return eta
   
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
    if a > 0.5:
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

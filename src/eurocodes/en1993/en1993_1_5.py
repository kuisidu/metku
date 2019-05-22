# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-5 Rules of plated structures

@author: kmela
"""

import math

from eurocodes.en1993.constants import E, gammaM1


def shear_eta(fy):
    """ Computes the 'eta' according to EN 1993-1-5, clause 5.1(2) """

    eta = 1.0
    
    if fy <= 460.0:
        eta = 1.2
    
    return eta

def reduction_factor_internal(lp,psi):
    """ Reduction factor for internal compression elements
        4.4(2), and AC:2009
        input: lp .. slenderness (lambda_p)
              psi .. stress ratio
    """
    
    if lp <= 0.5+math.sqrt(0.085-0.055*psi):
        rho = 1.0
    else:
        rho = min((lp-0.055*(3+psi))/lp**2,1.0)
    
    return rho

def reduction_factor_outstand(lp):
    """ Reduction factor for outstand compression elements
        4.4(2), and AC:2009
        
        input: lp .. slenderness (lambda_p)
    """
    
    if lp <= 0.748:
        rho = 1.0
    else:
        rho = min((lp-0.188)/lp**2,1.0)
    
    return rho

def lambda_p(b,t,epsilon,ksigma):
    """ Slenderness according to 4.4(2)
        input:
            b .. width of the element
            t .. thickness of the element
            epsilon = sqrt(235/fy)
            ksigma .. buckling factor
    """
    return b/t/28.4/epsilon/math.sqrt(ksigma)
    
    
def buckling_factor_internal(psi=1.0):
    """ Table 4.1 """
    
    if psi == 1.0:
        ksigma = 4.0
    elif psi > 0:
        ksigma = 8.2/(1.05+psi)
    elif psi == 0:
        ksigma = 7.81
    elif psi > -1:
        ksigma = 7.81-6.29*psi+9.78*psi**2
    elif psi == -1.0:
        ksigma = 23.9
    elif psi > -3:
        ksigma = 5.98*(1-psi**2)
        
    return ksigma

def buckling_factor_outstand(psi=1.0,sigma1="tip"):
    """ Table 4.2 
        input:
            sigma1 .. location of the maximum compressive stress
                      "tip" = unsupported end
                      "support" = supported end
    """
    
    if sigma1 == "tip":
        if psi == 1.0:
            ksigma = 0.43
        elif psi == 0.0:
            ksigma = 0.57
        elif psi == -1:
            ksigma = 0.85
        else:
            ksigma = 0.57-0.21*psi+0.07*psi**2
    else:
        if psi == 1.0:
            ksigma = 0.43
        if psi < 1.0 and psi > 0.0:
            ksigma = 0.578/(psi+0.34)
        elif psi == 0.0:
            ksigma = 1.7
        elif psi < 0 and psi > -1:
            ksigma = 1.7-5*psi+17.1*psi**2
        else:
            ksigma = 23.8
    
    return ksigma

def transverse_force_resistance(fyw,hw,tw,fyf,bf,tf,ss,a=0,type="a",c=0):
    """ Clause 6 
        input:
            fyw .. yield strength of the web [MPa]
            hw .. height of the web (mm)
            tw .. web thickness (mm)
            fyf .. yield strength of the flange
            bf .. width of flange
            tf .. thicknes of flange
            a .. distance between transverse stiffeners (a=0 for no stiffeners)
            ss .. length of stiff bearing
            type .. type of load application
                    "a" .. through one flange at the middle of the member
                    "b" .. through two flanges at the middle of the member
                    "c" .. through one side of the flange at the end of the member
            c .. distance of the edge of the load application surface from the end
                of the member (for type "c")
    """
    if a == 0:
        a = 1000*hw
    
    ss = min(ss,hw)
    
    if type == "a":
        kF = 6 + 2*(hw/a)**2
    elif type == "b":
        kF = 3.5+2*(hw/a)**2
    elif type == "c":
        kF = min(2+6*(ss+c)/hw,6.0)
    
    # Next, the reduction factor for effective length (Clause 6.4) and
    # Effective loaded length are determined. This may require iteration
    m1 = fyf*bf/fyw/tw
    m2 = 0
    
    Fcr = 0.9*kF*E*tw**3/hw
            
    if type == "a" or type == "b":
        ly = min(ss + 2*tf*(1+math.sqrt(m1+m2)),a)
    elif type == "c":
        le = min(0.5*kF*E*tw**2/fyw/hw,ss+c)
        ly = min(le + tf*math.sqrt(0.5*m1+(le/tf)**2+m2),le+tf*math.sqrt(m1+m2))
    
    lambdaF = math.sqrt(ly*tw*fyw/Fcr)
    
    if lambdaF > 0.5:
        # in this case, the initial assumption of m2=0 was wrong, and
        # ly needs to be corrected:
        m2 = 0.02*(hw/tf)**2
    
    if type == "a" or type == "b":
        ly = min(ss + 2*tf*(1+math.sqrt(m1+m2)),a)
    elif type == "c":
        le = min(0.5*kF*E*tw**2/fyw/hw,ss+c)
        ly = min(le + tf*math.sqrt(0.5*m1+(le/tf)**2+m2),le+tf*math.sqrt(m1+m2))
    
    lambdaF = math.sqrt(ly*tw*fyw/Fcr)
        
    # Reduction factor (Eq. (6.3))
    chiF = min(0.5/lambdaF,1.0)
    
    Leff = chiF*lambdaF
    
    FRd = fyw*Leff*tw/gammaM1
    
    return FRd
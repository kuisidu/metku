# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-5 Rules of plated structures

@author: kmela
"""

import math



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
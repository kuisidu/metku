# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-3 Supplementary rules for cold-formed members and sheeting

@author: kmela
"""

import math
import numpy as np

from eurocodes.en1993.constants import gammaM0, gammaM1, gammaM2

""" Chapter 3: Materials """
def fya(fyb,fu,t,Ag,n,k=7):
    """ Average yield strength, Eq. (3.1) 
        input:
            fyb .. basic yield strength [MPa]
            fu .. ultimate strength [MPa]
            k .. numerical coefficient depending on the type of forming
            n .. number of 90 degree bends
            t .. design core thickness of steel before colf-forming [mm]
    """

    return min(fyb+(fu-fyb)*k*n*t**2/Ag,0.5*(fu+fyb))

""" Chapter 5: Structural analysis """

def gr(r,t,phi):
    """ Calculate the length gr of Fig. 5.1 """
    rm = r + 0.5*t
    v = math.radians(phi)
    return rm*(math.tan(v) - math-sin(v))

def notional_width(b,rm,phi=np.array([90,90])):
    """ Notional width of a part
        input:
            b .. length of straight part of the part
            rm .. list or numpy array of midline radii of the corners
                 can be also a single number, if the part has only
                 one corner
            phi .. list or numpy array of corner angles.
                   can also be a single number
    """
    rm = np.array(rm)
    phi = np.array(phi)
    #print(sum(rm*np.sin(np.radians(0.5*phi))))
    return b + sum(rm*np.sin(np.radians(0.5*phi)))

def rounding_factor(r,phi,bp):
    """ Factor 'delta' of Eq. (5.1d) 
        input:
            r .. inner radius of different parts (list)
            phi .. angles (list)
            bp .. notional widths (list)
    """
    r = np.array(r)
    phi = np.array(phi)
    bp = np.array(bp)
    
    return 0.43*sum(r*phi/90)/sum(bp)
    

""" Chapter 8: Design of Joints """

""" Table 8.2: Design resistances for self-tapping screws """
def screw_validity(e1,e2,p1,p2,d,t,t1):
    """ Check the conditions on the range of validity """
    if e1 >= 3*d and e2 >= 1.5 and p1 >= 3*d and p2 >= 3*d and d <= 8.0 and d >= 3.0:
        return True
    else:
        return False

def screw_bearing_resistance(fu,d,t,t1,e1):
    """ Bearing resistance of a screw loaded in shear
        input:
            fu .. ultimate strength of the plate
            d .. diameter of the screw
            t .. thickness of the thinner plate
            t1 .. thickness of the thicker plate
            e1 .. the end distance from the centre of the fastener
                    to the adjacent end of the connected part, in 
                    the direction of load transfer
    """
    if t == t1:
        alfa = min(3.6*math.sqrt(t/d),2.1)
    elif t1 >= 2.5*t:
        alfa = 2.1
    else:
        # Linear interpolation not implemented!
        alfa = 2.1
    
    FbRdMax = fu*e1*t/1.2/gammaM2
    
    FbRd = min(alfa*fu*d*t/gammaM2,FbRdMax)
    return FbRd

def screw_net_secton_resistance(Anet,fu):
    """ Net-section resistance
        input:
            Anet .. net cross-sectional area of the connected part
            fu ..  ultimate tensile strength of the supporting member 
                    into which a screw is fixed
    """
    FnRd = Anet*fu/gammaM2
    return FnRd
    
def screw_shear_resistance_fin(diameter,material="hardened"):
    """ Screw shear resistance according to the Finnish
        National Annex
        
        input:
            diameter .. diameter of the screw [mm]
            material .. either "hardened" or "stainless"
        
        output:
            shear strength of the screw [N]
    """
    
    """ Values from the table """
    FvRd_hardened = {4.8:5200, 5.5:7200, 6.3:9800, 8.0:16300}
    FvRd_stainless = {4.8:4600, 5.5:6500, 6.3:8500, 8.0:14300}
    
    if material == "hardened":        
        FvRd = FvRd_hardened[diameter]
    elif material == "stainless":
        FvRd = FvRd_stainless[diameter]
        
    return FvRd
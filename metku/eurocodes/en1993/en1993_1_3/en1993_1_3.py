# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-3 Supplementary rules for cold-formed members and sheeting

@author: kmela
"""

import math
import numpy as np

from eurocodes.en1993.constants import gammaM0, gammaM1, gammaM2

def linear_interpolation(p1,p2,x):
    """ Performs linear interpolation.
        p1 = (x1,y1), p2 = (x2,y2)
        
        The function finds the value y corresponding to x1 <= x <= x2:
        Such that y(x) = k*x + b
    """
    
    y1 = p1[1]
    x1 = p1[0]
    y2 = p2[1]
    x2 = p2[0]
    
    k = (y2-y1)/(x2-x1)
    b = y1-k*x1
    
    return k*x+b

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
    return rm*(math.tan(v) - math.sin(v))

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

def distortional_buckling_reduction_factor(lambda_d):
    """ EN 1993-1-3, Eq. (5.12)
        Reduction factor for distortional buckling
        
        input:
            lambda_d .. slenderness = sqrt(fyb/sigma_cr,s)
    """
    
    if lambda_d <= 0.65:
        chi_d = 1.0
    elif lambda_d <= 1.38:
        chi_d = 1.47-0.723*lambda_d
    else:
        chi_d = 0.66/lambda_d
        
    return chi_d

def buckling_factor_edge_stiffener(bratio):
    """ EN 1993-1-3, Eq. (5.13b) and (5.13c)
        Buckling factor for single edge stiffener
        
        input:
            bratio = b_p,c/b_p, where
            b_p,c = notional width of the edge stiffener
            b_p = notional width of the flange etc.
    """
    
    if bratio <= 0.35:
        ksigma = 0.5
    elif bratio <= 0.6:
        ksigma = 0.5 + 0.83*((bratio-0.35)**2)**(1/3)
    else:
        print("Error: b ratio must be at most 0.6.")
        
    return ksigma

def distortional_buckling_stress(Is,As,E,K):
    """ EN 1993-1-3, Eq. (5.15)
        Critical stress for distortional buckling
    """    
    return 2*math.sqrt(K*E*Is)/As

""" Chapter 8: Design of Joints """

""" Table 8.2: Design resistances for self-tapping screws """
def screw_validity(e1,e2,p1,p2,d,t,t1):
    """ Check the conditions on the range of validity """
    if e1 >= 3*d and e2 >= 1.5 and p1 >= 3*d and p2 >= 3*d and d <= 8.0 and d >= 3.0:
        return True
    else:
        return False

def screw_bearing_resistance(fu,d,t,t1,verb=False):
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
        alfa = min(3.2*math.sqrt(t/d),2.1)
    elif t1 >= 2.5*t and t < 1.0:
        alfa = min(3.2*math.sqrt(t/d),2.1)
    elif t1 >= 2.5*t:
        alfa = 2.1
    else:
        # Linear interpolation not implemented!
        if verb:
            print("Linear interpolation for alpha")
        p0 = (t,min(3.2*math.sqrt(t/d),2.1))
        if t < 1.0:
            a1 = min(3.2*math.sqrt(t/d),2.1)
        else:
            a1 = 2.1
        p1 = (2.5*t,a1)
        
        print(p0,p1)
        alfa = linear_interpolation(p0,p1,t1)
        #alfa = 2.1
    # This is for rivets!
    #FbRdMax = fu*e1*t/1.2/gammaM2
        
    FbRd = alfa*fu*d*t/gammaM2 
    if verb:
        print("Self-tapping screw:")
        print("  t = {0:4.2f} mm".format(t))
        print("  t1 = {0:4.2f} mm".format(t1))
        print(" alpha = {0:4.2f}".format(alfa))
        print(" FbRd = {0:4.3f} kN".format(FbRd*1e-3))
    
    return FbRd, alfa

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

def sheeting_shear_stiffness(t,s,hw,broof):
    """ Eq. (10.1b) 
        input:
            t .. design thickness of sheeting
            s .. distance between purlins
            hw .. height of sheeting (between centerlines)
            broof .. overall length of the roof diaphragm in the direction of the
                    span of the sheeting in mm.
    """
    return 1000*math.sqrt(t**3)*(50 + 10*broof**(1/3))*s/hw

def CD_A(ba,tnom,position,A,bT,bR,load='gravity'):
    """ EN 1993-1-3, Eq. (10.17) """
    if ba <= 125:
        kba = (ba/100)**2
    elif ba < 200:
        kba = 1.25*(ba/100)
    else:
        print("Warning: width of purlin flange must not exceed 200mm!")
        kba = 1.25*(ba/100)
        
    if tnom > 0.75:
        if position == 'pos':
            p = 1.1
        else:
            p = 1.5
    else:
        p = 1.5
    
    kt = (tnom/0.75)**p
        
    if bR <= 185:
        kbR = 1.0
    else:
        kbR = 185/bR
    
    if load == 'gravity':
        if position == 'pos':
            if tnom == 0.75:
                kA = 1.0+(A-1.0)*0.08
            elif tnom >= 1.0:
                kA = 1.0+(A-1.0)*0.095
        else:
            if tnom == 0.75:
                kA = 1.0+(A-1.0)*0.16
            elif tnom >= 1.0:
                kA = 1.0+(A-1.0)*0.095
    else:
        kA = 1.0
    
    """ Table 10.3
        Assume that pitch of fasteners is e = bR and sheet is fastened
        through troughs.
    """
    if load == 'gravity':    
        if position == 'pos':    
            bT_max = 40
            C100 = 5.2
        else:
            bT_max = 120
            C10 = 3.1
    else:
        C100 = 2.6
        bT_max = 40
        
    if bT > bT_max:
        kbT = math.sqrt(bT_max/bT)
    else:
        kbT = 1.0
    
    CDA = C100*kba*kt*kbR*kA*kbT
    
    return CDA

class GrooveStiffener():
    """ Class for a plane element with intermediate stiffener
        EN 1993-1-3 (2006): 5.5.3.3
    """
    
    def __init__(self,nodes,t):
        """ Contructor
            :param nodes: list of [y,z] coordinates of the points
            :param t: calculation thickness
        """
        self.nodes = nodes
        self.t = t
    
    def bp1(self):
        """ Determine the first notional width """

    def bp2(self):
        """ Determine the second notional width """

if __name__ == "__main__":
    
    from materials.steel_data import Steel    
    
    plate = Steel("S350GD")
    
    screw_bearing_resistance(fu=plate.fu,d=4.2,t=0.5,t1=3.0,verb=True)
    
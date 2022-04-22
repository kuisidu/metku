# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-3 Supplementary rules for cold-formed members and sheeting

@author: kmela
"""

import math
import numpy as np

from eurocodes.en1993.constants import gammaM0, gammaM1, gammaM2, gammaMfi
from eurocodes.en1993 import en1993_1_2 

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

""" Chapter 6: Design of members """
def axial_force_resistance(self,Ned,A,**kwargs):
    try:
        if Ned >= 0:
            NRd = kwargs['fya']*A/gammaM0        
        else:
            Aeff = kwargs['Aeff']        
            if Aeff < A:
                NRd = kwargs['fyb']*Aeff/gammaM0
            else:
                lambda_ratio = kwargs['lambda_ratio']
                fyb = kwargs['fyb']
                fya = kwargs['fya']
                NcRd = A*(fyb+4*(fya-fyb)*(1-lambda_ratio))
                    
                NRd = min(NcRd,fya*A)/gammaM0
    
    except:
        print("axial_force_resistance: invalid keyword!")
        NRd = 0.0
          
    return NRd

def RwRd(fyb,t,hw,ss,c,r,phi=90,e=0.0,flange_stiffeners=False,web_stiffeners=False,verb=False):
    """ Resistance to transverse forces, cross sections with a single
        web. EN 1993-1-3, 6.1.7.2
    

    Parameters
    ----------
    fyb : float
        Yield strength of the base material.
    t : float
        Design value of wall thickness.

    Returns
    -------
    RwRd : Resistance value, in N.

    """
    
    def k1(k):
        return 1.33-0.33*k
    
    def k2(r,t):
        return min(max(1.15-0.15*r/t,0.5),1.0)
    
    def k3(phi):
        return 0.7 + 0.3*(phi/90**2)

    def k4(k):
        return 1.22-0.22*k
    
    def k5(r,t):
        return min(1.06-0.06*r/t,1.0)
    
    def k5_star(k):
        return max(1.49-0.53*k,0.6)
    
    def k6(t):
        return 0.88-0.12*t/1.9
    
    def k7(hw,t):
        hw_t = hw/t
        
        if hw_t < 150:
            return 1+hw_t/750
        else:
            return 1.2
    
    def k8(hw,t,k):
        hw_t = hw/t
        
        if hw_t < 66.5:
            return 1/k
        else:
            return (1.1-hw_t/665)*k
    
    def k9(t):
        return 0.82 + 0.15*t/0.19
    
    def k10(hw,t,k):
        return (0.98-hw/t/865)*k
    
    def k11(t):
        return 0.64+0.31*t/1.9
    
    # Base value of the resistance
    Rw0 = t**2*fyb/gammaM1
    
    k = fyb/228
    
    if abs(e) < 1e-6:
        # In this case, the transverse load is acting at the end
        # of the member.
        if c <= 1.5*hw:
            if web_stiffeners:
                C0 = k7(hw,t)*(8.8+1.1*np.sqrt(ss/t))
            else:
                # Load is close to the end of the member
                C0 = k1(k)*k2(r,t)*k3(phi)
                ss_t = ss/t
                if flange_stiffeners:
                    C = C0*(9.04-hw/t/60)*(1+0.01*ss_t)
                else:
                    if ss_t <= 60:
                        C = C0*(5.92-hw/t/132)*(1+0.01*ss_t)
                    else:
                        C = C0*(5.92-hw/t/132)*(0.71+0.015*ss_t)
        else:
            if web_stiffeners:
                C0 = k5_star(k)*k6(t)*(13.2+2.87*np.sqrt(ss/t))                
            else:
                C0 = k3(phi)*k4(k)*k5(r,t)
                ss_t = ss/t
                if ss_t <= 60:
                    C = C0*(14.7-hw/t/49.5)*(1+0.007*ss_t)
                else:
                    C = C0*(14.7-hw/t/49.5)*(0.75+0.011*ss_t)
    else:
        # Two opposite forces.
        # Assume for now that the distance between forces is
        # e < 1.5*hw.
        # If this is not the case, then we go for case 1 and examine
        # the two forces separately.
        if c <= 1.5*hw:
            if web_stiffeners:
                C0 = k10(hw,t,k)*k11(t)*(8.8+1.1*np.sqrt(ss/t))
            else:
                # Load is close to the end of the member
                C0 = k1(k)*k2(r,t)*k3(phi)
                C = C0*(6.66-hw/t/64)*(1+0.01*ss/t)
        else:
            if web_stiffeners:
                C0 = k8(hw,t,k)*k9(t)*(13.2+2.87*np.sqrt(ss/t))
            else:
                C0 = k3(phi)*k4(k)*k5(r,t)
                C = C0*(21.0-hw/t/16.3)*(1+0.0013*ss/t)
            
    RwRd = C*Rw0
    
    return RwRd
    
def RwRd_two_webs(fyb,r,t,phi,hw,ss,verb=False,code='old',**kwargs):
    """    

    Parameters
    ----------
    fyb : float
        Yield strength of the base material [MPa].
    r : float
        internal corner radius.
    t : float
        thickness of the steel core
    phi : float
        angle of the web relative to flanges (degrees), 45 <= phi <= 90
    hw : float
        Height of the web.
    ss : float
        Support width.
    category : integer, optional
        Loading category. value is either 1 or 2. The default is 1.
    force_type : string
        Type of loading. Either 'load' for local load or 'support' for support
    **kwargs : 
        either empty, or 'e = ' or 'c = '.

    Returns
    -------
    Resistance [N]

    """
    # Young's modulus
    E = 210000
    
    e = -1
    c = -1
    
    for key, value in kwargs.items():        
        if key == 'e':
            e = value            
        elif key == 'c':
            c = value        
        elif key == 'VEd':
            # list of shear forces of transverse shear forces on each side
            # of the local load or support reaction, with
            # abs(VEd[0]) > abs(VEd[1])
            VEd = value
        else:
            raise ValueError('Incorrect input value, give either c or e.')
            
    category = 1
        
    # The last condition corresponds to internal support reaction, which means
    # that neither 'c' or 'e' is given as input parameter    
    if e > 1.5*hw or c > 1.5*hw or (e < 0 and c < 0):
        
        category = 2

        
    if category == 1:
        # This a value is for trapezoidal sheeting
        # for liner trays and hats, use a = 0.057
        if code == 'old':
            a = 0.075
        else:
            if e > 0:
                # Sheeting with opposite loading
                a = 0.075
            else:
                # Sheeting near clear end
                a = min(0.18*np.sqrt((c+ss)/1.5/hw),0.15)
        # Effective bearing length
        la = 10.0
    else:
        # This a value is for trapezoidal sheeting
        # for liner trays and hats, use a = 0.115
        a = 0.15
        
        if code == 'old':
            VEd = [abs(V) for V in VEd]
            
            beta_v = (VEd[0]-VEd[1])/(VEd[0]+VEd[1])
            
            if beta_v <= 0.2:
                la = ss
            elif beta_v >= 0.3:
                la = 10.0
            else:
                # linear intepolation
                la = ss - (beta_v-0.2)/0.1*(ss-10)
        else:
            la = ss
        
        la = min(la,200)
    
    RwRd = a*t**2*np.sqrt(fyb*E)*(1-0.1*np.sqrt(r/t))*(0.5+np.sqrt(0.02*la/t))*(2.4+(phi/90)**2)/gammaM1   
    
    if verb:
        print("Resistance for transverse loads:")
        print(f"Category: {category}")
        print(f"Web height: hw = {hw:{4}.{2}f} mm")
        print(f"Support width: ss = {ss:{4}.{2}f} mm")
        print(f"alpha = {a:{4}.{2}f}")
        
        if category == 2 and code == "old":
            print(f"beta_v = {beta_v:{4}.{2}f}")
        
        print(f"Effective bearing length: la = {la:{4}.{2}f} mm")
        
        print(f" Rw_Rd = {RwRd*1e-3:{4}.{2}f} kN")
        
        
    
    return RwRd

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

def screw_bearing_resistance_elevated_temp(fu,d,t,fub,e1,e2,T,verb=False):
    """ Bearing resistance at elevated temperatures """
    
    fuT = en1993_1_2.fyT(T,fu)
    kb = en1993_1_2.kbT(T)
    fubT = kb*fub

    ad = e1/3/d
    abT = min([ad,fubT/fuT,1.0])

    k1 = min(2*e2/d-1,2.5)
    
    FbRd = k1*abT*fubT*d*t/gammaMfi
    
    if verb:
        print("Self-tapping screw:")
        print(f"fuT = {fuT:.2f} MPa")
        print(f"kbT = {kb:.3f}")
        print(f"fubT = {fubT:.2f} MPa")
        print(f"d = {d} mm")
        print(f"t = {t} mm")
        print(f"alpha_(b,T) = {abT:.2f}")
        print(f"k1 = {k1:.2f}")
        print(" FbRd = {0:4.3f} kN".format(FbRd*1e-3))
    
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
    
    #plate = Steel("S350GD")
    
    #R = RwRd_two_webs(fyb=350, r=2.4, t=1.16, phi=65, hw=130, ss=150,verb=True,code='new',
    #                  c=1.6*130,VEd=[300,100])
    #print(R)
    
    d = 5.5
    t = 1.0-0.04
    FbRd = screw_bearing_resistance_elevated_temp(fu=420,d=d,t=t,fub=800,e1=109*d,e2=2.1*d,T=639,verb=True)    
    
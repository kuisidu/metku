# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Thu Oct 17 19:53:07 2019

Composite beam sections

@author: kmela
"""

import math
import numpy as np

from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.patches as patches


from metku.eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
from metku.eurocodes.en1994 import en1994_1_1
    

class ConcreteFilledTube:
    """ Class for composite columns that are concrete filled tubes
    """
    
    def __init__(self,length,steel_part,concrete,rebar=None,rebar_pos=None):
        """ Constructor
            Input:
                steel_part: RHS or SHS class profile of the steel section
                concrete: Concrete class material
                rebar: 
                
            parameters
                length .. length of the beam
                steel .. steel profile
                concrete .. Concrete class material
                rebar .. Rebar class. Rebar is placed symmetrically by default,
                        using a group of four rebars
                rebar_pos .. coordinates of the top right rebar, assumed to lay
                            in the positive orthant
                nbars .. number of rebars
                ned .. axial force (assumed compression)
                nged .. axial force due to dead loads
                med .. med[0] .. bending moment at the top of the column
                       med[1] .. bending moment at the bottom of the column
                
        """
        
        self.length = length
        self.steel = steel_part
        self.concrete = concrete
        self.rebar = rebar
        self.rebar_pos = [0,0]
        self.nbars = 4
        self.ned = 0.0
        self.nged = 0.0
        self.med = [0.0,0.0]
        
        self.lcr = [1.0,1.0]

        
        if rebar_pos is None:
            """ Set default position """
            self.u = 30
            self.v = 30
        else:    
            self.rebar_pos = rebar_pos
    
    
    @property
    def weight_per_length(self):
        """ Weight per unit length (kN/m) """
        Wsteel = self.steel.weight()*1e3*9.81e-3
        Wconcrete = self.Ac*1e-6*self.concrete.density 
        Wrebar = self.rebar.weight*self.nbars*9.81e-3
        
        #print(Wsteel,Wconcrete,Wrebar)
        return Wsteel + Wconcrete + Wrebar
    
    @property
    def NEd(self):
        """ Design value of axial force """
        return abs(self.ned)
    
    @property
    def M1(self):
        """ Minimum of bending end moments (absolute value) """
        if abs(self.med[0]) <= abs(self.med[1]):
            return self.med[0]
        else:
            return self.med[1]
        
    
    @property
    def M2(self):
        """ Maximum of bending end moments (absolute value) """
        if abs(self.med[0]) <= abs(self.med[1]):
            return self.med[1]
        else:
            return self.med[0]
    
    @property
    def fyd(self):
        """ Design value of yield strength of steel"""
        return self.steel.fy
    
    @property
    def fcd(self):
        """ Design value of yield strength of concrete """
        return self.concrete.fcd(a_cc=1.0)
    
    @property
    def fsd(self):
        """ Design value of yield strength of reinforcement """
        return self.rebar.fyd
    
    @property
    def lcry(self):
        """ Buckling length with respect to major axis """
        return self.length*self.lcr[0]
    
    @property
    def lcrz(self):
        """ Buckling length with respect to minor axis """
        return self.length*self.lcr[1]
    
    @property
    def u(self):
        """ Distance in z-direction (vertical) of the rebar edge from
            inner surface of the steel profile
        """
        return 0.5*self.steel.H - self.steel.T -  (self.rebar_pos[1] + self.rebar.r)
        
    @u.setter
    def u(self,u):
        self.rebar_pos[1] = 0.5*self.steel.H - self.steel.T - u - self.rebar.r
    
    @property
    def v(self):
        """ Distance in y-direction (horizontal) of the rebar edge from
            inner surface of the steel profile
        """
        return 0.5*self.steel.B - self.steel.T - self.rebar_pos[0] + self.rebar.r
        
    @v.setter
    def v(self,v):
        self.rebar_pos[0] = 0.5*self.steel.B - self.steel.T - v - self.rebar.r
    
    @property
    def z_rebar(self):
        """ Vertical position of the rebar (center) from the center of the profile """
        return self.rebar_pos[1]
    
    @property
    def y_rebar(self):
        """ Horizontal position of the rebar (center) from the center of the profile """
        return self.rebar_pos[0]
    
    @property
    def As(self):
        """ Total area of rebar
        """
        return 4*self.rebar.A
    
    @property
    def Aa(self):
        """ Area of steel part """
        return self.steel.A
    
    @property
    def Ia(self):
        """ Second moment of area of steel part with respect to
            the major axis at centroid of steel part
        """
        return self.steel.I[0]
    
    @property
    def Wpa(self):
        
        return self.steel.Wpl[0]
    
    @property
    def Ea(self):
        """ Young's modulus of steel part """
        return self.steel.E
    
    @property
    def EAa(self):
        """ Axial stiffness of steel part """
        return self.Ea*self.Aa
    
    @property
    def EIa(self):
        """ Bending stiffness of steel part 
            with respect to the neutral axis of steel part
        """
        return self.Ea*self.Ia
    
    @property
    def Is(self):
        """ Second moment of area of reinforcement """
        Is0 = self.rebar.I
        Is1 = self.nbars*(Is0 + self.rebar.A*self.z_rebar**2)
        
        return Is1
    
    @property
    def Wps(self):
        """ Bending resistance of reinforcement """
        return self.nbars*self.rebar.A*self.z_rebar
    
    @property
    def EIs(self):
        """ Bending stiffness of reinforcement """
        return self.rebar.E*self.Is
    
    @property
    def Eceff(self):
        """ Effective Young's modulus for concrete 
            EN 1994-1-1, Eq. (6.41)
        """
        phi = self.concrete.creep_coefficient(h0=1e11,t=1e11,t0=28)
        
        print(phi)
        return self.Ecm*1/(1+self.nged/self.ned*phi)
        
    
    @property
    def Ecm(self):
        """ Young's modulus of concrete """
        return self.concrete.Ecm
    
    @property
    def EIc(self):
        """ Bending stiffness of concrete part 
            with respect to the neutral axis of concrete part
        """
        return self.Ecm*self.Ic
    
    @property
    def Ac(self):
        """ Area of concrete part """
        Ac = (self.steel.H-2*self.steel.R)*(self.steel.B-2*self.steel.T) + \
                2*(self.steel.B-2*self.steel.R)*(self.steel.R-self.steel.T) + \
                math.pi*(self.steel.R-self.steel.T)**2 -\
                self.As
        return Ac
        #print(Ac)
        #return (self.steel.H-2*self.steel.T)*(self.steel.B-2*self.steel.T)-self.As
    
    @property
    def Ic(self):
        """ Second moment of area of concrete part with respect to
            major axis at the centroid
        """
        b = self.steel.B
        h = self.steel.H
        t = self.steel.T
        # Inner radius
        r = self.steel.R-t
        
        Icy = (b-2*t)*(h-2*t)**3/12 - r**4/3 - r**2*(h-2*t-r)**2 + \
                +r**4*(9*math.pi**2-64)/36/math.pi + \
                +math.pi*r**2*(0.5*h-t+4*r/3/math.pi)**2 - self.Is
        
        return Icy
    
    
    @property
    def Wpc(self):
        """ Section modulus of concrete """
        b = self.steel.B
        h = self.steel.H
        t = self.steel.T
        # Inner radius
        r = self.steel.R-t
        
        return 0.25*(b-2*t)*(h-2*t)**2 - 2*r**3/3 - r**2*(4-math.pi)*(0.5*h-t-r)-self.Wps

    
    @property
    def EAc(self):
        """ Axial stiffness of concrete """
        return self.Ecm*self.Ac
        

    def EAcom(self):
        """ axial stiffness composite section """
        return self.EAa + self.EAc

    def EIcom(self):
        """ sending stiffness of composite section """        
        return self.Ea*(self.Ia+self.ya_com()**2*self.Aa) + self.Ecm*(self.Ic+self.yc_com()**2*self.Ac)

    
    def EIeff(self,axis="y"):
        """ Effective flexural stiffness for determining Ncr 
            EN 1994-1-1, Eq. (6.40)
        """
        Ke = 0.6
        
        if axis == "y":
            Ic = self.Ic
    
        return self.EIa + self.EIs + Ke*self.Eceff*Ic
    
    
    def EIeffII(self,axis="y"):
        """ Effective flexural stiffness
            EN 1994-1-1, Eq. (6.42)
        """
        Ko = 0.9
        Ke = 0.5
        
        if axis == "y":
            Ic = self.Ic
    
        return Ko*(self.EIa + self.EIs + Ke*self.Eceff*Ic)
    
    def e0(self):
        """ Member imperfection 
            EN 1994-1-1, table 6.5
        """
        
        if self.reinforcement_ratio() <= 0.03:
            r = 300
        else:
            r = 200
            
        return self.length/r
    
    def Ra(self):
        """ Full resultant of steel part """
        
        return self.steel.fy*self.Aa
    
    def Rc(self):
        """ Full resultant of concrete part """
        
        return 0.85*self.slab.material.fcd()*self.Ac
    
    def Rf(self):
        """ Resultant of top flange of the steel beam """
        
        Af = self.steel.tf * self.steel.b        
        
        return self.steel.fy*Af
    
    def Rw(self):
        """ Resultant of web of the steel beam """
        
        Aw = self.steel.hw * self.steel.tw        
        
        return self.steel.fy*Aw
    
    def steel_contribution_ratio(self,verb=False):
        
        delta_a = self.NaRd()/self.NplRd()
        
        if verb:
            print("Steel contribution ratio, delta_a = {0:4.2f}".format(delta_a))
        
        return delta_a
    
    
    def reinforcement_ratio(self,verb=False):
        
        rho_s = self.As/self.Ac
        
        if verb:
            print("Reinforcement ratio, rho_s = {0:4.2f}".format(rho_s))
        
        return rho_s
    
    def NaRd(self,verb=False):
        """ Axial resistance of the steel part"""
        NaRd = self.steel.fy*self.Aa
        
        if verb:
            print("Na,Rd = {0:4.2f} kN".format(NaRd*1e-3))
        return NaRd
    
    def NplRd(self,verb=False):
        """ Axial resistance """
        NplRd = self.steel.fy*self.Aa + self.rebar.fyd*self.As + self.fcd*self.Ac
        
        if verb:
            print("Npl,Rd = {0:4.2f} kN".format(NplRd*1e-3))
        return NplRd
    
    def NplRk(self,verb=False):
        """ Axial resistance with characteristic strength values"""
        NplRk = self.steel.fy*self.Aa + self.rebar.fy*self.As + self.concrete.fck*self.Ac
        
        if verb:
            print("Npl,Rk = {0:4.2f} kN".format(NplRk*1e-3))
            
        return NplRk
    
    def Ncr(self,axis="y"):
        """ Elastic critical force """
        
        if axis == "y":
            Lcr = self.lcry
        else:
            Lcr = self.lcrz
            
        return math.pi**2*self.EIeff(axis)/Lcr**2
    
    def Ncreff(self,axis="y"):
        """ Elastic critical force """
        
        if axis == "y":
            Lcr = self.lcry
        else:
            Lcr = self.lcrz
            
        return math.pi**2*self.EIeffII(axis)/Lcr**2
    
    
    def slenderness(self,axis="y"):
        
        Ncr = self.Ncr(axis)
        
        return math.sqrt(self.NplRk()/Ncr)
    
    def NbRd(self,axis="y",verb=False):
        """ Buckling strength """
        l = self.slenderness(axis)
        if self.reinforcement_ratio() <= 0.03:
            imp_factor = en1993_1_1.buckling_curve["a"]
        else:
            imp_factor = en1993_1_1.buckling_curve["b"]
            
        chi = en1993_1_1.buckling_reduction_factor(l,imp_factor)
    
        NbRd = chi*self.NplRd()
    
        if verb:
            print("**Column buckling: **")
            print("Slenderness: {0:4.3f}".format(l))
            print("Imperfection factor: {0:4.2f}".format(imp_factor))
            print("Buckling reduction factor: {0:4.3f}".format(chi))
            print("Buckling strength: {0:4.2f} kN".format(NbRd*1e-3))
            print("Utilisation: {0:4.2f} ".format(self.NEd/NbRd))
        
        return NbRd 
    
    
    def NpmRd(self):
        
        return self.fcd*self.Ac
    
    def Mpl_max_Rd(self):
        """ Maximum moment resistance """        
        
        return self.fyd*self.Wpa + self.fsd*self.Wps + 0.5*self.fcd*self.Wpc
    
    def MplRd(self,verb=False):
        """ Plastic moment resistance """
        b = self.steel.B
        t = self.steel.T
        fsd = self.fsd
        fcd = self.fcd
        fyd = self.fyd
        
        """ Assume first that the neutral axis is between rebars """
        Asn = 0
        
        hn = (self.NpmRd() - Asn*(2*fsd-fcd))/(2*b*fcd + 4*t*(2*fyd-fcd))
        
        
        """ If neutral axis is below top rebar row, there is no negative
            part for the plastic resistance
        """
        if hn < self.z_rebar - self.rebar.r:
            Wps_n = 0
        else:
            Asn = self.As
            hn = (self.NpmRd() - Asn*(2*fsd-fcd))/(2*b*fcd + 4*t*(2*fyd-fcd))
            Wps_n = self.Wps
            
        Wpc_n = (b-2*t)*hn**2 - Wps_n
        Wpa_n = b*hn**2 - Wpc_n - Wps_n

        MnRd = fyd*Wpa_n + fsd*Wps_n + 0.5*fcd*Wpc_n

        MplRd = self.Mpl_max_Rd() - MnRd 

        if verb:
            print("Plastic moment resistance:")
            print("hn = {0:4.2f} mm (z_rebar = {1:4.2f} mm)".format(hn,self.z_rebar))
            print("Wps_n = {0:4.2f} mm3".format(Wps_n))
            print("Wpc_n = {0:4.2f} mm3".format(Wpc_n))
            print("Wpa_n = {0:4.2f} mm3".format(Wpa_n))
            print("MplRd = {0:4.2f} kNm".format(MplRd*1e-6))

        return MplRd 
    
    def alpha_m(self):
        
        if self.steel.fy <= 355:
            alpha_m = 0.9
        else:
            alpha_m = 0.8
            
        return alpha_m
    
    
    
    def VRd(self):
        """ Shear resistance
            This equals shear resistance of the steel part
        """
        
        return self.steel.shear_force_resistance()
    
    def design(self,verb=False):
        """ Carry out design according to EN 1994-1-1, simplified method """
        
        """ Check validity """
        side_ratio = self.steel.H/self.steel.B
        
        A_ratio = self.reinforcement_ratio()
        
        if side_ratio < 0.2 or side_ratio > 5.0:
            """ steel part dimensions out of scope """
            print("Warning: the ratio h/b = {0:4.2f} is out of scope.".format(side_ratio))
        
        if A_ratio > 0.06:
            """ too much rebar """
            print("Warning: too much rebar - {0:4.2f} > 0.06".format(A_ratio))
        
        self.steel_contribution_ratio(verb)
        
        slend = self.slenderness()
        
        if verb:
            print("Slenderness = {0:4.3f}".format(slend))
            
        if slend > 2.0:
            """ too slender member """
            print("Warning: slenderness too large - {0:4.2f} > 2.0".format(A_ratio))
        
        
        if self.M2 == 0:
            """ Only axial force """
            self.NbRd(verb=True)
        
        """ Design for combined bending moment and axial force """
        e0 = self.e0()
        Med0 = self.NEd*e0
        
        if verb:
            print("Initial imperfection e0 = {0:4.2f} mm".format(e0))
            print("Bending moment due to Initial imperfection MEd0 = {0:4.2f} kNm".format(Med0*1e-6))
            
        """ Check for need to include second-order effects """
        Ncr_eff = self.Ncreff()
        
        if verb:
            print("Check for second-order effects:")
            print("Ncr,eff = {0:4.2f} kN".format(Ncr_eff*1e-3))
            print("10*NEd = {0:4.2f} kN".format(10*self.NEd*1e-3))
        
        if Ncr_eff < 10*self.NEd:
            if verb:
                print("Include second-order effects.")
                
            if self.M2 > 0:
                bm = max(0.66+0.44*self.M1/self.M2,0.44)
            else:
                bm = 0
                
            MEdI = bm*self.M2 + Med0
            km = max(1/(1-self.NEd/Ncr_eff),1.0)
            
            MEd = max(km*MEdI,self.M2)
            
            if verb:
                print("bm = {0:4.3f}".format(bm))
                print("MEd,I = {0:4.2f} kNm".format(MEdI*1e-6))
                print("km = {0:4.3f}".format(km))
                print("MEd = {0:4.2f} kNm".format(MEd*1e-6))
                            
        else:
            MEd = Med0 + self.M2
            if verb:
                print("No need to include second-order effects.")
                
        
        """ Do M+N -interaction """
        MplRd = self.MplRd(True)
        NplRd = self.NplRd()
        chiC = self.NpmRd()/NplRd
        chiD = self.NEd/NplRd
        muD = min((chiD-1)/(chiC-1),1.0)
        
        Mpl_N_Rd = muD*MplRd
        
        UMN = MEd/Mpl_N_Rd/self.alpha_m()
        
        if verb:
            print("chi_C = Npm_Rd/Npl_Rd = {0:4.3f}".format(chiC))
            print("chi_d = NEd/Npl_Rd = {0:4.3f}".format(chiD))
            print("mu_d = (chi_d-1)/(chi_C-1) = {0:4.3f}".format(muD))
            print("MplRd = {0:4.3f} kNm".format(MplRd*1e-6))
            print("Mpl_N_Rd = mu_d*Mpl_Rd = {0:4.3f} kNm".format(Mpl_N_Rd*1e-6))
            print("MEd/Mpl_N_Rd = {0:4.3f}".format(MEd/Mpl_N_Rd))
            print("Utilization ratio: {0:4.3f}".format(UMN))
            
        
        
        
    
    def draw(self):
        """ Draw the profile             
        """
        
        fig, ax = plt.subplots(1)

        self.steel.draw(axes=ax)        
                    
        # Draw rebar
        r = self.rebar.r        
        
        xy = np.array((self.rebar_pos))
        
        dxy = np.array([[1,1],[-1,1],[-1,-1],[1,-1]])
        
        for xyloc in dxy:
            ax.add_patch(patches.Circle(tuple(xyloc*xy),r))
        
        
        ax.set_aspect('equal')

if __name__ == '__main__':
    
    from sections.steel.RHS import SHS
    from eurocodes.en1992 import en1992_1_1, constants
    """
    steel = SHS(300,8)
    con = constants.Concrete("C40/50")
    rebar = en1992_1_1.Rebar(20,"B500B")
    ned = 2.4e6
    nged = 1.7e6
    med = [0,0]
    
    """
    p = ConcreteFilledTube(6500,steel,con,rebar)
    p.u = 30
    p.ned = 2.118e6
    p.nged = 767.696e3
    p.design(True)
    #p.NplRd(True)
    #p.NaRd(True)
    #p.NplRk(True)
    p.steel_contribution_ratio(True)
    #p.draw()
    """
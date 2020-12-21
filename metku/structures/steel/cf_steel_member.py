# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 21:22:03 2020

@author: kmela
"""

import math
import numpy as np
from metku.structures.steel.steel_member import SteelMember

class CFSteelMember(SteelMember):
    """ Class for cold-formed steel members """
    
    def __init__(self, profile, length, Lcr=[1.0, 1.0], mtype="beam",
                 LT_buckling=False):
        
        super().__init__(profile,length,Lcr,mtype,LT_buckling)
    
    def slenderness(self, verb=False):
        """ Non-dimensional slenderness according to EN 1993-1-1 """
        NRd = self.profile.Aeff * self.fy()
        Ncr = self.ncrit(verb)
        # if NRd <= 1e-6:
        # print(self.profile.A)
        # print(self.profile.h)
        # print(NRd, Ncr)

        slend = np.sqrt(NRd / Ncr)

        if verb:
            print("lambda,y = {0:4.2f}".format(slend[0]))
            print("lambda,z = {0:4.2f}".format(slend[1]))

        return slend
    
    def buckling_strength(self, verb=False):
        """ Member buckling strength according to EN 1993-1-3 """

        if verb:
            print("** Buckling strength** ")

        slend = self.slenderness(verb)
        NbRd = []
        NRd = self.profile.Aeff * self.fy()
        # self.profile.imp_factor()
        # print(slend)
        # p = 0.5*(1+np.array(self.profile.imp_factor)*(slend-0.2)+slend**2)
        # r = 1/(p+np.sqrt(p**2-slend**2))


        for i in range(len(slend)):
            p = 0.5 * (1 + self.profile.imp_factor[i] * (slend[i] - 0.2) +
                       slend[i] ** 2)
            r = min(1 / (p + math.sqrt(p ** 2 - slend[i] ** 2)), 1.0)
            if verb:
                print("Psi,{0:1} = {1:4.2f}".format(i, p))
                print("chi,{0:1} = {1:4.2f}".format(i, r))
                print("NbRd,{0:1} = {1:4.2f}".format(i, NRd * r * 1e-3))

            NbRd.append(NRd * r)

        # NbRd = self.profile.A*self.fy()*r;
        return NbRd    
    
    def Sreq(self):
        """ Required shear stiffness for sheeting to provide full
            lateral support to the purlin
        """
        
        Sw = self.profile.E*self.profile.Iw*math.pi**2/self.length**2
        St = self.profile.G*self.profile.It
        Sz = self.profile.E*self.profile.Iz*math.pi**2/self.length**2*0.25*self.profile.h**2
        
        return (Sw + St + Sz)*70/self.profile.h**2
    
    
    
if __name__ == '__main__':
    
    from eurocodes.en1993.en1993_1_3 import cf_profs
    from eurocodes.en1993.en1993_1_3.en1993_1_3 import sheeting_shear_stiffness
    
    c = cf_profs.CSection(t_nom=2.0,h=120,a=60,b=60,ca=20,cb=20,r=2.0)
    z = cf_profs.ruukki_z(300,2.5)
    m = CFSteelMember(z,length=5000)
    
    print(m.kh())
    
    #print(c.kh())
    #print(m.Sreq())
    #print(sheeting_shear_stiffness(t=0.96,s=1500,hw=(70-1),broof=24000))
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 21:22:03 2020

@author: kmela
"""

import math
from steel_member import SteelMember

class CFSteelMember(SteelMember):
    """ Class for cold-formed steel members """
    
    def __init__(self, profile, length, Lcr=[1.0, 1.0], mtype="beam",
                 LT_buckling=False):
        
        super().__init__(profile,length,Lcr,mtype,LT_buckling)
        
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
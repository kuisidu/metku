# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:26:02 2019

@author: kmela
"""

import numpy as np

from eurocodes.en1993 import en1993_1_1
from structures.steel.steel_member import SteelMember
from sections.steel.ISection import IPE

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def kiepahdus():
    p = IPE(300)
    m = SteelMember(profile=p, length=6000, Lcr=[1.0, 1.0], mtype="beam", LT_buckling=True)
    m.add_section(myed=30.0e6)
    za = 0.5*m.profile.h
    Mcr = m.mcrit(C=[1.132, 0.459], za=za, k=[1, 1])
    
    MbRd = m.LT_buckling_strength(Mcr, axis='y',method='I',verb=True)
    
    print("Iz = {0:4.2f} mm4".format(p.I[1]))
    print("It = {0:4.2f} mm4".format(p.It))
    print("Iw = {0:4.2f} mm6".format(p.Iw))
    print("Mcr = {0:4.4f} kNm".format(Mcr*1e-6))
    
    return m, p
    


if __name__ == '__main__':
    
    m, p = kiepahdus()
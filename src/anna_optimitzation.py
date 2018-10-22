# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 10:33:13 2018

@author: mercader
"""

import sys
# path needs to be added to use tables_and_tuples and I_sections
sys.path.append('S:\91202_METKU\Kristo\Python\src\End-plate')


from anna_semi_rigid_frame import SemiRigidFrame

storey_height = 3.65
bay_length = 7.3
storeys = 2
bays = 1
num_elements = 2

SRF = SemiRigidFrame(storeys=storeys,
                     bays=bays,
                     storey_height=storey_height,
                     bay_length=bay_length,
                     num_elements=num_elements)

class optimitzation:
    def __init__(self,SRF,r,h):
        self.r=r
        self.h=h
        
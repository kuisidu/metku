# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 13:02:03 2017

@author: kmela
"""
import numpy as np

from hollow_sections import SHS
from roof_truss import RoofTruss

p = SHS(150,6.0)

t = RoofTruss(24000,2200,1/20,"KT",divX=2)

"""
t = truss.Truss()
t.add_node(0,0)
t.add_node(0,1)
t.add_node(1,0)

t.add_member(0,1,p)
t.add_member(0,2,p)
t.add_member(1,2,p)
"""

t.draw()

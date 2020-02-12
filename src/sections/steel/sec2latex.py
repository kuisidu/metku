# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 19:46:57 2019

@author: kmela
"""

from sections.steel.catalogue import ipe_profiles
from sections.steel.ISection import IPE

for prof_name, dims in ipe_profiles.items():
    p = IPE(int(dims["h"]))
    
    
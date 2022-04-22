# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 19:46:57 2019

@author: kmela
"""

from sections.steel.catalogue import ipe_profiles
from sections.steel.ISection import IPE

for prof_name, dims in ipe_profiles.items():
    p = IPE(int(dims["h"]))
    
    
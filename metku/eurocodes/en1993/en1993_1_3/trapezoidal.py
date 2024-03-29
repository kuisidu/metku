# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 19:34:15 2021

Cold-formed trapezoidal sheeting profiles

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from metku.materials.steel_data import Steel
from metku.eurocodes.en1993.en1993_1_3.open_prof import OpenProf
from metku.eurocodes.en1993.en1993_1_3 import en1993_1_3
from metku.eurocodes.en1993 import en1993_1_5
from metku.eurocodes.en1993.constants import gammaM0, gammaM1


class TrapezoidalSheet:
    """ Class for trapezoidal sheeting """
    
    def __init__(self,t,h,bt,bb,angle,r=3,material=Steel("S350GD"),ntroughs=3):
        """
        Constructor

        Parameters
        ----------
        t : float
            Nominal thickness.
        h : float
            height of profile.
        bt : float
            width of top flange.
        bb : float
            width of bottom flange.
        angle : float
            angle between bottom flange and chord (degrees).
        r : float
            rounding (interal)
        material : Steel, optional
            material. The default is Steel("350GD").
        ntroughs : float
            number of troughs in one sheet

        Returns
        -------
        None.

        """
        
        self.tnom = t
        self.tcoat = 0.04
        self.t = self.tnom-self.tcoat
        self.h = h
        self.bt = bt
        self.bb = bb
        self.angle = angle
        self.angle_rad = np.radians(angle)
        self.material = material
        self.ntrough = ntroughs
        
        
    @property
    def E(self):
        return self.material.E
    
    @property
    def nu(self):
        return self.material.nu
    
    @property
    def fy(self):
        return self.material.fy
    
    @property
    def fu(self):
        return self.material.fu
    
    @property
    def hw(self):
        """ Profile height between center lines """
        return self.h-self.tnom
    
    @property
    def sw(self):
        """ Length of one web """
        return self.hw*np.sin(self.angle_rad)
    
    @property
    def pitch(self):
        """ Pitch of corrugations """
        return self.bb+self.bt+2*self.h*np.cos(self.angle_rad)
    
    @property
    def width(self):
        """ Sheet width """
        return self.pitch*self.ntrough
    
    @property
    def A(self):
        """ Cross-sectional area """
        return self.ntrough*self.t*(self.bt+self.bb+2*self.sw)
    
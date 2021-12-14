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
from materials.steel_data import Steel
from eurocodes.en1993.en1993_1_3.open_prof import OpenProf
from eurocodes.en1993.en1993_1_3 import en1993_1_3
from eurocodes.en1993 import en1993_1_5
from eurocodes.en1993.constants import gammaM0, gammaM1


class TrapezoidalSheet:
    """ Class for trapezoidal sheeting """
    
    def __init__(self,t,h,bt,bb,angle,r=3,material=Steel("350GD"),ntroughs=3):
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
    def hw(self):
        """ Profile height between center lines """
        return self.h-self.tnom
    
    @property
    def sw(self):
        """ Length of one web """
        return self.hw*np.sin(self.angle_rad)
    
    @property
    def A(self):
        """ Cross-sectional area """
        return self.ntrough*self.t*(self.bt+self.bb+2*self.sw)
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:26:02 2019

@author: kmela
"""

import numpy as np

from eurocodes.en1993 import en1993_1_1

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def beamcolumn():
    """ Beam-column interaction """

    xLT = 1.0
    UN = 0.8
    UM = 0.8
    Cmy = 0.9
    
    l = np.linspace(0.2,2)
    
    kyy = []
    
    for i in range(l):
        kyy.append(en1993_1_1.kyy(UN,l[i],Cmy))
    
    


if __name__ == '__main__':
    
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-5 Rules of plated structures

@author: kmela
"""

import math



def shear_eta(fy):
    """ Computes the 'eta' according to EN 1993-1-5, clause 5.1(2) """

    eta = 1.0
    
    if fy <= 460.0:
        eta = 1.2
    
    return eta
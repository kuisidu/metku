# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 15:03:13 2022

Various functions for Raami package

@author: kmela
"""

import numpy as np


def divide_line_segment(X1,X2,nseg):
    """ Divides line segment between points X1 and X2 into 'nseg'
        segments.
        
        Returns the points, including X1 and X2
    """
    
    # Determine direction vector from X1 to X2
    #v = (X2-X1)/np.linalg.norm(X2-X1)
    v = X2-X1
    
    T = np.linspace(0,1,nseg+1)
    
    X = [X1 + t*v for t in T]
    
    return X

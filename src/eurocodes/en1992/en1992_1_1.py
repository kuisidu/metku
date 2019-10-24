# -*- coding: utf-8 -*-
"""
Created on Sun Oct 6 2019

EN 1992-1-1 General rules for concrete structures

@author: kmela
"""

import math

try:
    from src.eurocodes.en1992.constants import gammaC, gammaS, alpha_cc, concrete
except:
    from eurocodes.en1992.constants import gammaC, gammaS, alpha_cc, concrete




def fcd(fck,a_cc=alpha_cc,gC=gammaC):
    """ Design value of compression strength
    """
    return a_cc*fck/gC

# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:59:36 2022

Constants to be used in optimization

@author: kmela
"""

CACHE_BOUND = 2**8
INT_TOL = 1e-4
""" Tolerance for discreteness violation of a discrete variable """
DISC_TOL = 1e-3
""" Tolerance for comparing variable vectors """
X_TOL = 1e-10

""" Artificial variable upper and lower bounds """
XUB = 1e8
XLB = -1e8

""" Tolerances for solvers """
GRAD_TOL = 1e-8  # Tolerance for the norm of gradient
ABS_F_TOL = 1e-8 # Tolerance for objective function difference, absolute
REL_F_TOL = 1e-4 # Tolerance for objective function difference, relative
G_TOL = 1e-8     # Tolerance for constraint feasibility
X_TOL = 1e-8     # Tolerance for change in variable values
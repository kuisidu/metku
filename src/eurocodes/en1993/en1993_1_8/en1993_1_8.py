# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-8 Rules for connections

@author: kmela
"""

# Bolt materials [MPa], SFS EN 1993-1-8, table 3.1
mat_bolt = {4.6: {"f_yb": 240.0, "f_ub": 400.0},
            4.8: {"f_yb": 320.0, "f_ub": 400.0},
            5.6: {"f_yb": 300.0, "f_ub": 500.0},
            5.8: {"f_yb": 400.0, "f_ub": 500.0},
            6.8: {"f_yb": 480.0, "f_ub": 600.0},
            8.8: {"f_yb": 640.0, "f_ub": 800.0},
            10.9: {"f_yb": 900.0, "f_ub": 1000.0}}

# Standard bolt sizes
bolt_sizes = [12, 16, 20, 24, 30, 36]

# Standard bolt lengths
bolt_lengths = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 90.0, 100.0, 110.0,
                120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0]

# Bolt sizes, stress area [mm^2]
bolt_size = {12: {"A_s": 84.3},
             16: {"A_s": 157.0},
             20: {"A_s": 245.0},
             24: {"A_s": 353.0},
             30: {"A_s": 561.0},
             36: {"A_s": 817.0},
             }

import math

from eurocodes.en1993.constants import gammaM5

""" CHAPTER 3: Bolts """
def bolt_shear_resistance(fub,A,bolt_class=8.8,threads_in_plane=False):
    """ EN 1993-1-8, Table 3.4
        Input:
            fub .. ultimate strength of bolt [MPa]
            A .. cross-sectional area [mm^2]
            threads_at_plane .. True, if shear plane passes through threaded
                                portion of the bolt
                                False, otherwise
        Output:
            FvRd .. bolt shear resistance [N]
    """
    
    if threads_in_plane:
        av = 0.6
    else:
        if bolt_class = {10.9,6.8,5.8,4.8}:
            av = 0.5
        elif bolt_class = {4.6,5.6,8.8}:
            av = 0.6
    
    FvRd = av*fub*A/gammaM2
    return FvRd

def bolt_bearing_resistance(fub,fu,d,t,e,p,d0,pos_perp="edge",pos_load="edge"):
    """ EN 1993-1-8, Table 3.4
        Input:
            fub .. ultimate strength of bolt [MPa]
            fu .. ultimate strength of plate [MPa]
            d .. diameter of bolt hole [mm]
            t .. thickness of plate [mm]
            pos_perp .. position of bolt in direction perpendicular to
                        load transfer ("edge" or "innter")
            pos_load .. position of bolt in direction of load transfer
                        ("edge" or "inner")
        Output:
            FbRd .. bearing resistance [N]
    """
    e1 = e[0]
    e2 = e[1]
    p1 = p[0]
    p2 = p[1]
    
    if pos_load == "edge":
        ad = e1/3/d0
    else:
        ad = p1/3/d0 - 0.25    

    if pos_perp == "edge":
        k1 = min(2.8*e2/d0-1.7,2.5)
    else:
        k1 = min(1.4*p2/d0-1.7,2.5)
        
    
    ab = min(ad,fub/fu,1.0)
    
    FbRd = k1*ab*fu*d*t
    
    return FbRd

def bolt_tension_resistance(fub,As,k2=0.9):
    """ EN 1993-1-8, Table 3.4
        Input:
            fub .. ultimate strength of bolt [MPa]
            As .. stress area of bolt [mm^2]
            k2 = 0.63 for countersunk bolts
               = 0.9 otherwise (default)
        Output:
            FtRd .. bolt tension resistance [N]
    """
    FtRd = k2*fub*As/gammaM2
    return FtRd

def bolt_punching_shear_resistance(dm,tp,fu):
    """ EN 1993-1-8, Table 3.4
        Input:
            dm .. mean of the across points and across flats 
                dimensions of the bolt head or the nut, whichever is smaller
            tp .. thickness of the plate under the bolt or nut
        Output:
            BpRd .. punching shear resistance [N]
    """
    BpRd = 0.6*math.pi*dm*tp*fu/gammaM2
    return BpRd

def shear_and_tension_resistance(FvEd,FvRd,FtEd,FtRd):
    """ EN 1993-1-8, Table 3.4
        Input:
            FvEd .. shear force in bolt
            FvRd .. shear resistance of bolt
            FtEd .. tension force in bolt
            FtRd .. tension resistance of bolt
        Output:
            U .. utilization ratio
    """
    U = FvEd/FvRd + FtEd/FtRd/1.4
    return U
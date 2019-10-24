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

from eurocodes.en1993.constants import gammaM0, gammaM2, gammaM5

""" CHAPTER 3: Bolts """
class Bolt:
    """ Class for bolts """
    
    def __init__(self,size,bolt_class=8.8,length=50.0):
        """ Constructor
             input: size .. bolt size
                    bolt_class .. bolt class
                    length .. length of bolt
        """
        self.d = size # diameter of unthreaded shank
        
        # bolt hole: EN 1090-2 (2018): Table 11
        if size < 16:
            self.d0 = size+1
        elif size < 27:
            self.d0 = size+2
        else:
            self.d0 = size+3
        
        self.A = math.pi*size**2/4
        self.As = bolt_size[size]["A_s"]
        
        self.length = length
        self.fyb = mat_bolt[bolt_class]["f_yb"]
        self.fub = mat_bolt[bolt_class]["f_ub"]
        self.bolt_class = bolt_class
        
    def shear_resistance(self,threads_in_plane=False,verb=False):
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
            As = self.A
        else:
            As = self.As
            if self.bolt_class in {10.9,6.8,5.8,4.8}:
                av = 0.5
            elif self.bolt_class in {4.6,5.6,8.8}:
                av = 0.6
    
        
        FvRd = av*self.fub*As/gammaM2
        
        if verb:
            print("Shear resistance:")
            if threads_in_plane:
                print("Bolt threads in shear plane")
            else:
                print("Bolt threads not in shear plane")
                
            print("Area: As = {0:4.2f}".format(As))            
            print("alpha_v = {0:4.2f}".format(av))            
            print("fub = {0:4.2f}".format(self.fub))            
            print("FvRd = {0:4.2f} kN".format(FvRd*1e-3))    
        
        
        return FvRd
    
    def tension_resistance(self,k2=0.9):
        """ EN 1993-1-8, Table 3.4
            Input:
                k2 = 0.63 for countersunk bolts
                   = 0.9 otherwise (default)
            Output:
                FtRd .. bolt tension resistance [N]
        """
        FtRd = k2*self.fub*self.As/gammaM2
        return FtRd
    
 
    def bearing_resistance(self,fu,t,e,p,pos_perp="edge",pos_load="edge",\
                           verb=False):
        """ EN 1993-1-8, Table 3.4
            Input:
                t .. thickness of plate [mm]
                pos_perp .. position of bolt in direction perpendicular to
                            load transfer ("edge" or "inner")
                pos_load .. position of bolt in direction of load transfer
                            ("edge" or "inner")
                e .. array for edge distances e[0] = e1, e[1] = e2
                p .. array for bolt row distances p[0] = p1, p[1] = p2
                fu .. ultimate strength of plate [MPa]
            Output:
                FbRd .. bearing resistance [N]
        """
        e1 = e[0]
        e2 = e[1]
        p1 = p[0]
        p2 = p[1]
        
        if pos_load == "edge":
            ad = e1/3/self.d0
        else:
            ad = p1/3/self.d0 - 0.25    
    
        if pos_perp == "edge":
            k1 = min(2.8*e2/self.d0-1.7,2.5)
        else:
            k1 = min(1.4*p2/self.d0-1.7,2.5)
            
        
        ab = min(ad,self.fub/fu,1.0)
        
        FbRd = k1*ab*fu*self.d*t/gammaM2
        
        if verb:
            print("Bearing resistance:")
            print("Edge distances: e1 = {0:4.2f}, e2 = {1:4.2f}".format(e1,e2))
            print("Bolt row distances: p1 = {0:4.2f}, p2 = {1:4.2f}".format(p1,p2))
            print("alpha_d = {0:4.2f}".format(ad))
            print("alpha_b = {0:4.2f}".format(ab))
            print("k1 = {0:4.2f}".format(k1))
            print("d = {0:4.2f}".format(self.d))
            print("t = {0:4.2f}".format(t))
            print("t = {0:4.2f}".format(t))
            print("FbRd = {0:4.2f} kN".format(FbRd*1e-3))
        
        return FbRd

class BoltRow:
    """ Class for bolt rows """
    
    def __init__(self,bolt,p,z):
        """ Constructor
            input:
                bolt .. Bolt class object defining the bolt
                p .. array of distances between adjacent bolt centroids
                z .. position of bolt row
        """
        
        self.bolt = bolt
        self.p = p
        self.z = z
        self.bolts = len(p)+1
        
        
class BoltRowGroup:
    """ Class for group of bolt rows """
    
    def __init__(self,bolt_rows):
        """ Constructor
            input:
                bolt_rows .. array of BoltRow objects
        """
        
        self.rows = bolt_rows
        self.nrows = len(bolt_rows)
        self.p = []
        
        """ Evaluate distance between adjacent groups
            it is assumed that the rows as inserted in ordered list
        """
        for i in range(len(bolt_rows)-1):
            self.p.append(bolt_rows[i+1].z-bolt_rows[i].z)
    

def block_tearing(fy,fu,Ant,Anv,concentric_load=True,verb=False):
    """ Check block tearing resistance 
        fy .. yield strength of the plate [MPa]
        fu .. ultimate strength of the plate [MPa]
        Ant .. net area subject to tension [mm2]
        Anv .. net area subject to shear [mm2]
    """
    
    if concentric_load:
        Veff = fu*Ant/gammaM2 + fy/math.sqrt(3)*Anv/gammaM0
    else:
        Veff = 0.5*fu*Ant/gammaM2 + fy/math.sqrt(3)*Anv/gammaM0
    
    if verb:
        print("Block tearing:")
        print("fy = {0:4.2f} MPa, fu = {1:4.2f} MPa".format(fy,fu))
        print("Ant = {0:4.2f} mm^2, Anv = {1:4.2f} mm^2".format(Ant,Anv))
        print("Veff,Rd = {0:4.2f} kN".format(Veff*1e-3))
    
    return Veff
    

    
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
        if bolt_class in {10.9,6.8,5.8,4.8}:
            av = 0.5
        elif bolt_class in {4.6,5.6,8.8}:
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
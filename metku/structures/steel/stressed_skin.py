# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 19:48:26 2021

Stressed skin design

@author: kmela
"""

import numpy as np

from metku.eurocodes.en1993.en1993_1_3.trapezoidal import TrapezoidalSheet
from metku.eurocodes.en1993.en1993_1_3 import en1993_1_3

from metku.sections.steel.RHS import RHS, SHS
from metku.structures.steel.steel_member import SteelMember

K1 = {0: [[0.013, 0.03, 0.041, 0.041, 0.046, 0.05, 0.066, 0.103, 0.193],
          [0.042, 0.096, 0.131, 0.142, 0.142, 0.153, 0.199, 0.311, 0.602],
          [0.086, 0.194, 0.264, 0.285, 0.283, 0.302, 0.388, 0.601, 1.188],
          [0.144, 0.323, 0.438, 0.473, 0.468, 0.494, 0.629, 0.972, 1.935],
          [0.216, 0.438, 0.654, 0.703, 0.695, 0.729, 0.922, 1.42, 2.837],
          [0.302, 0.674, 0.911, 0.98, 0.965, 1.008, 1.266, 1.938, 3.892],
          [0.402, 0.895, 1.208, 1.3, 1.277, 1.329, 1.661, 2.536, 5.098],
          [0.516, 1.146, 1.546, 1.662, 1.631, 1.692, 2.107, 3.208, 6.453]],
      5: [[0.014, 0.031, 0.041, 0.044, 0.044, 0.049, 0.066, 0.107, 0.205],
          [0.05, 0.099, 0.128, 0.134, 0.132, 0.146, 0.198, 0.336, 0.652],
          [0.107, 0.202, 0.253, 0.26, 0.254, 0.28, 0.386, 0.681, 1.548],
          [0.188, 0.338, 0.413, 0.417, 0.404, 0.448, 0.629, 1.158, 2.639],
          [0.295, 0.507, 0.604, 0.601, 0.578, 0.648, 0.934, 1.783], 
          [0.429, 0.706, 0.823, 0.806, 0.772, 0.877, 1.306, 2.586], 
          [0.591, 0.935, 1.066, 1.028, 0.983, 1.135, 1.756, 3.605], 
          [0.78, 1.191, 1.328, 1.264, 1.208, 1.423, 2.299, 4.838]],
      10: [[0.016, 0.031, 0.04, 0.042, 0.042, 0.048, 0.065, 0.111, 0.221],
           [0.056, 0.101, 0.123, 0.125, 0.123, 0.139, 0.2, 0.366, 0.873],
           [0.125, 0.204, 0.238, 0.233, 0.226, 0.264, 0.402, 0.786], 
           [0.222, 0.338, 0.375, 0.356, 0.345, 0.418, 0.689, 1.445], 
           [0.349, 0.494, 0.526, 0.486, 0.473, 0.605, 1.082, 2.428], 
           [0.502, 0.668, 0.682, 0.615, 0.608, 0.837, 1.607], 	
           [0.677, 0.851, 0.834, 0.736, 0.752, 1.128, 2.308], 	
           [0.869, 1.035, 0.975, 0.844, 0.907, 1.494, 3.2]],
      15: [[0.017, 0.031, 0.04, 0.041, 0.041, 0.047, 0.066, 0.115],
           [0.062, 0.102, 0.118, 0.115, 0.113, 0.134, 0.209, 0.403],
           [0.139, 0.202, 0.218, 0.204, 0.2, 0.254, 0.44, 0.945],
           [0.244, 0.321, 0.325, 0.293, 0.294, 0.414, 0.796], 
           [0.37, 0.448, 0.426, 0.371, 0.396, 0.636, 1.329], 
           [0.506, 0.568, 0.508, 0.434, 0.513, 0.941], 	
           [0.646, 0.668, 0.561, 0.483, 0.664, 1.349], 	
           [0.766, 0.735, 0.578, 0.527, 0.861]],
      20: [[0.018, 0.032, 0.039, 0.039, 0.039, 0.046, 0.066, 0.111],
           [0.068, 0.101, 0.111, 0.106, 0.104, 0.131, 0.221, 0.452],
           [0.148, 0.193, 0.194, 0.174, 0.177, 0.255, 0.492], 
           [0.249, 0.289, 0.267, 0.23, 0.259, 0.444, 0.931], 
           [0.356, 0.372, 0.315, 0.27, 0.364, 0.725], 	
           [0.448, 0.42, 0.326, 0.303, 0.512],
           [0.509, 0.423, 0.301, 0.346],
           [0.521, 0.372, 0.259, 0.413]],      
      25: [[0.019, 0.032, 0.038, 0.038, 0.038, 0.045, 0.068, 0.126],
           [0.072, 0.099, 0.103, 0.095, 0.095, 0.129, 0.236, 0.513],
           [0.151, 0.178, 0.166, 0.144, 0.16, 0.268, 0.557], 
           [0.238, 0.244, 0.204, 0.176, 0.247, 0.494],
           [0.306, 0.272, 0.203, 0.204, 0.376],
           [0.333, 0.248, 0.172, 0.241],
           [0.3, 0.174, 0.142],
           [0.204, 0.081]],
      30: [[0.02, 0.032, 0.037, 0.036, 0.036, 0.044, 0.07, 0.133],
           [0.075, 0.095, 0.094, 0.084, 0.087, 0.132, 0.256], 
           [0.148, 0.157, 0.135, 0.116, 0.152, 0.291], 	
           [0.208, 0.186, 0.139, 0.139, 0.253], 
           [0.226, 0.161, 0.112, 0.176], 
           [0.18, 0.089, 0.093],
           [0.077]],
      35: [[0.021, 0.032, 0.036, 0.034, 0.034, 0.043, 0.072, 0.142],
           [0.076, 0.089, 0.083, 0.072, 0.082, 0.137, 0.281], 
           [0.137, 0.13, 0.102, 0.093, 0.151], 
           [0.162, 0.119, 0.082, 0.12],
           [0.123, 0.059],
           [0.032]],
      40: [[0.023, 0.032, 0.034, 0.032, 0.032, 0.043, 0.075, 0.155],
           [0.075, 0.081, 0.07, 0.06, 0.077, 0.146], 	
           [0.116, 0.096, 0.068, 0.078], 	
           [0.1, 0.053, 0.048],
           [0.024]],
      45: [[0.024, 0.031, 0.032, 0.029, 0.03, 0.043, 0.079], 
           [0.071, 0.069, 0.056, 0.05, 0.073], 
           [0.086, 0.057, 0.041],  
           [0.032]]
      }

class StressedSkin:
    """ General class for stressed skins """
    
    
    def __init__(self,sheet,a,b,n,seam_fastener={'d':4.8, 'spacing':500, 'loc':'trough'}):
        """
        

        Parameters
        ----------
        sheet : TYPE
            DESCRIPTION.
        seam_fastener : TYPE
            DESCRIPTION.
        a : TYPE
            DESCRIPTION.
        b : TYPE
            DESCRIPTION.
        n : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
    
        self.sheet = sheet
        self.seam_fastener = seam_fastener
        
        self.a = a
        self.b = b
        self.n = n
                
        self.p = sheet.pitch
        
        self.fastener_calc = 'en1993'
        
        # These values are from Table 5.1 of ECCS (1995)
        # slip in the fastening between sheeting and rafter
        self.sp = 0.35
        # slip in the seam fasteners
        self.ss = 0.25
        # slip in edge members/shear connectors
        self.ssc = 0.35
    
    def Fs(self):
        """ Seam fastener strength """
        
        fu = self.sheet.fu
        d = self.seam_fastener['d']
        t = self.sheet.t
        
        if self.fastener_calc == 'en1993':        
            # Calculate resistance as bearing resistance according to
            # EN 1993-1-3
            t1 = self.sheet.t
            Fs = en1993_1_3.screw_bearing_resistance(fu,d,t,t1)[0]
        else:
            # Calculation according to ECCS guide
            Fs = min(2.9*np.sqrt(t/d)*fu*d*t,3.8e3)
            
        
        return Fs

    def beta3(self):
        """ Value of parameter beta3 which takes into account
            the location of seam fasteners
        """
        
        if self.seam_fastener['loc'] == 'trough':
            beta3 = 1.0
        else:
            beta3 = (self.nf-1)/self.nf

        return beta3

    def beta1(self):
        """ Factor to allow for the number of sheet/purlin
            fasteners per sheet width.
            
            ECCS (1995), Annex C.1
        """
        
        nf = self.nf()
        
        if (nf % 2) == 0:
            # Even number of fasteners
            if self.seam_fastener['loc'] == 'trough':
                # Decking
                beta1 = sum([((2*(i+1)-1)/(nf-1))**3 for i in range(0,int(nf/2))])
            else:
                # Sheeting
                beta1 = sum([((2*(i+1)-1)/nf)**3 for i in range(0,int(nf/2))])
        else:
            # Odd number of fasteners
            if self.seam_fastener['loc'] == 'trough':
                # Decking
                beta1 = sum([(2*(i+1)/(nf-1))**3 for i in range(0,int((nf-1)/2))])
            else:
                # Sheeting
                beta1 = sum([(2*(i+1)/nf)**3 for i in range(0,int((nf-1)/2))])
        
        return beta1

    def nf(self):
        """ Number of sheet/purlin fasteners per sheet width
            including those at the overlaps.
            
            This is equal to the number of corrugations
        """
        
        return self.sheet.ntrough
    
    def K(self,fasteners_in_every_trough=True,verb=False):
        """ Table 5.9 """
        
        d = self.sheet.pitch
        h = self.sheet.h
        l = self.sheet.bt
        
        l2d = l/d
        h2d = h/d
        
        lnd = int(l2d/0.1)-1
        hnd = int(h2d/0.1)-1
        
        angle = 90-self.sheet.angle
        
        int_angle = int(5*np.floor(angle/5))
        
        if fasteners_in_every_trough:
            K = K1[int_angle][hnd][lnd]
        else:
            raise ValueError("Fasteners in every other trough is not included yet.")
        
        if verb:
            print(f"l/d = {l2d:.2f}")
            print(f"h/d = {h2d:.2f}")
            print(f"angle = {angle:.2f}")
            print(f"rounded angle = {int_angle:.2f}")
            print(f"K = {K:.4f}")
        
        return K

class StressedSkinRafters(StressedSkin):
    """ Class for stressed skin design using rafters.
        The corrugations are spanning in the direction of the building.
    """
    
    def __init__(self,sheet,a,b,n,seam_fastener,edge_fastener,rafter,edge_member):
        """

        Parameters
        ----------
        sheet : TYPE
            DESCRIPTION.
        seam_fastener : TYPE
            DESCRIPTION.
        a : TYPE
            DESCRIPTION.
        b : TYPE
            DESCRIPTION.
        n : TYPE
            DESCRIPTION.
        rafter : TYPE
            DESCRIPTION.
        edge_member : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        super().__init__(sheet,a,b,n,seam_fastener)
        
        self.rafter = rafter
        self.edge_member = edge_member
        self.edge_fastener = edge_fastener
        
    @property
    def ns(self):
        """ Number of seam fasteners per seam in one panel,
            excluding the seam-to-rafter fasteners.
        """
        return int(self.b/self.seam_fastener['spacing'])-1
    
    @property
    def nsc(self):
        """ Number of edge member fasteners in one panel,
            excluding the seam-to-rafter fasteners.
        """
        return int(self.b/self.edge_fastener['spacing'])
    
    @property
    def nsh(self):
        """ Number of sheet widths per shear panel """
        return round(self.a/self.sheet.width)
    
    @property
    def rafter_t(self):
        """ Return rafter wall thickness """
        
        if hasattr(self.rafter,'tf'):
            t = self.rafter.tf
        elif hasattr(self.rafter,'T'):
            t = self.rafter.T
        elif hasattr(self.rafter,'tt'):
            t = self.rafter.tt
        
        return t
    
    @property
    def edge_member_t(self):
        """ Return edge member wall thickness """
        
        if hasattr(self.edge_member,'tf'):
            t = self.edge_member.tf
        elif hasattr(self.edge_member,'T'):
            t = self.edge_member.T
        elif hasattr(self.edge_member,'tt'):
            t = self.edge_member.tt
        
        return t
    
    def Fp(self):
        """ Shear resistance of fastener between rafter and sheet """
        
        fu = self.sheet.fu
        d = self.seam_fastener['d']
        t = self.sheet.t
        
        if self.fastener_calc == 'en1993':        
            # Calculate resistance as bearing resistance according to
            # EN 1993-1-3            
            t1 = self.rafter_t
                        
            Fp = en1993_1_3.screw_bearing_resistance(fu,d,t,t1)[0]
        else:
            # Calculation according to ECCS guide
            if d == 5.5:
                Fpmax = 6.5e3
            elif d == 6.3:
                Fpmax = 8.0e3
            else:
                raise ValueError("ECCS Guide only allows 5.5 mm or 6.3 mm sheet-to-rafter fasteners.")
                
            Fp = min(1.9*fu*d*t,Fpmax)
        
        return Fp
    
    def Fsc(self):
        """ Shear resistance of fastener between edge member and sheet """
        
        fu = self.sheet.fu
        d = self.edge_fastener['d']
        t = self.sheet.t
        
        if self.fastener_calc == 'en1993':        
            # Calculate resistance as bearing resistance according to
            # EN 1993-1-3            
            t1 = self.edge_member_t
                        
            F = en1993_1_3.screw_bearing_resistance(fu,d,t,t1)[0]
        else:
            # Calculation according to ECCS guide
            if d == 5.5:
                Fmax = 6.5e3
            elif d == 6.3:
                Fmax = 8.0e3
            else:
                raise ValueError("ECCS Guide only allows 5.5 mm or 6.3 mm sheet-to-edge member fasteners.")
                
            F = min(1.9*fu*d*t,Fmax)
        
        return F
        
    def seam_strength(self,verb=False):
        """ Resistance of the seam """
        
        Fs = self.Fs()
        Fp = self.Fp()
        
        ns = self.ns
        beta1 = self.beta1()
        beta3 = self.beta3()
        
        Vs = self.a/self.b*(ns*Fs + beta1/beta3*Fp)
        
        if verb:
            print("Seam strength:")
            print(f"a = {self.a:.2f} mm")
            print(f"b = {self.b:.2f} mm")
            print(f"Fs = {Fs*1e-3:.2f} kN")
            print(f"Fp = {Fp*1e-3:.2f} kN")
            print(f"beta_1 = {beta1:.0f} ")
            print(f"beta_3 = {beta3:.0f} ")
            print(f"ns = {ns:.0f} ")
            print(f"Vs = {Vs*1e-3:.2f} kN")
                    
        return Vs

    def edge_member_fastener_strength(self,verb=False):
        """ ECCS (1995), 5.8.1.2 """
        
        Fsc = self.Fsc()
        nsc = self.nsc
        
        V = self.a/self.b*Fsc*nsc
        
        if verb:
            print("Edge member fastener strength:")
            print(f"Fsc = {Fsc*1e-3:.2f} kN")
            print(f"nsc = {nsc:.0f} ")
            print(f"V = {V*1e-3:.2f} kN")
        
    def sheet_rafter_fastener_strength(self,verb=False):
        """ ECCS (1995), 5.8.3.1 """
        
        V = 0.6*self.a*self.Fp()/self.p

        if verb:
            print("Sheet-to-rafter fastener strength:")
            print(f"a = {self.a:.2f} mm")
            print(f"p = {self.p:.2f} mm ")
            print(f"V = {V*1e-3:.2f} kN")
    
        return V

    def end_collapse(self,verb=False):
        """ ECCS (1995), 5.8.3.2 """
        
        t = self.sheet.t
        fy = self.sheet.fy
        a = self.a
        d = self.sheet.pitch
        
        every_corrugation_fastened = True
        
        if every_corrugation_fastened:
            K = 0.9
        else:
            K = 0.3
        
        V = K*t**1.5*a*fy/np.sqrt(d)
        
        if verb:
            print("End collapse resistance:")
            print(f"a = {a:.2f} mm")
            print(f"t = {t:.2f} mm ")
            print(f"fy = {fy:.2f} MPa ")
            print(f"d = {d:.2f} mm ")
            print(f"V = {V*1e-3:.2f} kN")

    def profile_distortion(self,verb=True):
        """ Shear flexibility for profile distortion,
            Table 5.9 of ECCS (1995)
        """
        
        a = self.a
        b = self.b
        t = self.sheet.t
        E = self.sheet.E*1e-3
        d = self.sheet.pitch
        
        K = self.K()
        
        # Factor to take into account sheet continuity
        alpha5 = 1.0
        
        c1_1 = a*d**2.5*alpha5*K/E/t**2.5/b**2
        
        if verb:
            print("Shear flexibility: profile distortion")
            print(f"a = {a:.2f} mm")
            print(f"b = {b:.2f} mm")
            print(f"t = {t:.2f} mm ")            
            print(f"d = {d:.2f} mm ")
            print(f"K = {K:.4f} ")
            print(f"alpha_5 = {alpha5:.1f} ")
            print(f"c1_1 = {c1_1:.4f} mm/kN")
        
        return c1_1

    def shear_strain(self,verb=True):
        """ Shear flexibility for shear strain,
            Table 5.9 of ECCS (1995)
        """
        
        a = self.a
        b = self.b
        t = self.sheet.t
        h = self.sheet.h
        E = self.sheet.E*1e-3
        nu = self.sheet.nu
        d = self.sheet.pitch
                        
        c1_2 = 2*a*(1+nu)*(1+2*h/d)/E/t/b
        
        if verb:
            print("Shear flexibility: shear strain")            
            print(f"c1_2 = {c1_2:.4f} mm/kN")

        return c1_2

    def fastener_deformation(self,verb=True):
        """ Fastener deformation components """
               
        sp = self.sp
        ss = self.ss
        ssc = self.ssc
        nsh = self.nsh
        ns = self.ns
        nsc = self.nsc
        beta1 = self.beta1()
        
        # Sheet-to-rafter
        c2_1 = 2*self.a*sp*self.p/self.b**2
        
        # Seam fasteners
        c2_2 = ss*sp*(nsh-1)/(ns*sp+beta1*ss)
        
        # Edge members
        c2_3 = 2*ssc/nsc
        
        if verb:
            print("Shear flexibility: sheet-to-rafter")            
            print(f"c2_1 = {c2_1:.4f} mm/kN")
            print("Shear flexibility: seam fasteners")            
            print(f"c2_2 = {c2_2:.4f} mm/kN")
            print("Shear flexibility: edge members")  
            print(f"c2_3 = {c2_3:.4f} mm/kN")
        
        return [c2_1,c2_2,c2_3]
        
    def edge_member_strain(self,verb=True):
        """ Axial strain in edge members """
        a = self.a
        b = self.b
        E = self.edge_member.E*1e-3
        A = self.edge_member.A
        n = self.n
                
        c3 = n**2*b**3/E/A/a**2
        
        if verb:
            print("Shear flexibility: axial strain of edge member")            
            print(f"c3 = {c3:.4f} mm/kN")
        
        return c3
    
    def shear_flexibility(self,verb):
        """ Shear flexibility of one panel """
        
        c1_1 = self.profile_distortion(verb)
        c1_2 = self.shear_strain(verb)
        c2 = self.fastener_deformation(verb)
        c3 = self.edge_member_strain(verb)
        
        c_prime = self.b**2/self.a**2*(c1_1+c1_2+sum(c2))
        
        c = c_prime + c3
        
        if verb:
            print("Shear flexibility: ")
            print(f"c_prime = {c_prime:.4f} mm/kN")
            print(f"c = {c:.4f} mm/kN")
        
        return c
        
        
        
if __name__ == "__main__":
    
    sheet = TrapezoidalSheet(1.0,130,111,75,61)
    rafter = RHS(160,160,6)
    edge_member = RHS(100,100,5)
    seam_fastener = {'d':4.8, 'spacing':500, 'loc':'trough'}
    edge_fastener = {'d':5.5, 'spacing':250}
    a = 24000
    b = 6000
    n = 12
    
    s = StressedSkinRafters(sheet, a, b, n, seam_fastener, edge_fastener, rafter, edge_member)    
    
    #s.seam_strength(True)
    #s.edge_member_fastener_strength(True)
    #s.sheet_rafter_fastener_strength(True)
    #s.end_collapse(True)
    
    #s.profile_distortion(True)
    #s.shear_strain(True)
    #s.fastener_deformation(True)
    
    s.shear_flexibility(True)
    
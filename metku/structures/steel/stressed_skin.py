# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 19:48:26 2021

Stressed skin design

@author: kmela
"""

import numpy as np

from eurocodes.en1993.en1993_1_3.trapezoidal import TrapezoidalSheet
from eurocodes.en1993.en1993_1_3 import en1993_1_3

from sections.steel.RHS import RHS, SHS
from structures.steel.steel_member import SteelMember



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
                
        self.fastener_calc = 'en1993'
    
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
    

class StressedSkinRafters(StressedSkin):
    """ Class for stressed skin design using rafters.
        The corrugations are spanning in the direction of the building.
    """
    
    def __init__(self,sheet,a,b,n,seam_fastener,rafter,edge_member):
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
        
    @property
    def ns(self):
        """ Number of seam fasteners per seam in one panel,
            excluding the seam-to-rafter fasteners.
        """
        return int(self.b/self.seam_fastener['spacing'])-1
    
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


if __name__ == "__main__":
    
    sheet = TrapezoidalSheet(1.0,130,70,50,60)
    rafter = RHS(160,160,6)
    edge_member = RHS(100,100,5)
    seam_fastener = {'d':4.8, 'spacing':500, 'loc':'trough'}
    a = 24000
    b = 6000
    n = 12
    
    s = StressedSkinRafters(sheet, a, b, n, seam_fastener, rafter, edge_member)
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 19:48:26 2021

Stressed skin design

@author: kmela
"""





class StressedSkin:
    """ General class for stressed skins """
    
    
    def __init__(self,sheet,seam_fastener,a,b,n):
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
        


class StressedSkinRafters(StressedSkin):
    """ Class for stressed skin design using rafters.
        The corrugations are spanning in the direction of the building.
    """
    
    def __init__(self,sheet,seam_fastener,a,b,n,rafter,edge_member):
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
        
        super().__init__(sheet,seam_fastener,a,b,n)
        
        self.rafter = rafter
        self.edge_member = edge_member
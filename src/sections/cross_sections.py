# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 08:16:41 2020

@author: kmela
"""

class CrossSection(metaclass=ABCMeta):
    """ Base class for generic cross-sections
    
    Specific cross sections are subclasses of this class
    
    NOTE: The original purpose of this class is to act as a generic
            cross-section to be used in structural optimization
    
    """

    # Constructor
    def __init__(self,     
                 material,
                 A=0,
                 I=[0,0]):
        
        if isinstance(material,str):
            self.material = Steel(material)
        elif isinstance(material,(float,int)):
            self.material = Steel("S" + str(int(material)))

        self.E = self.material.E
        self.A = A
        self.I = np.asarray(I)       

    def __repr__(self):
        return type(self).__name__
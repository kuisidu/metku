# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 19:53:07 2019

Composite beam sections

@author: kmela
"""

import math
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

try:
    from src.eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
    from src.sections.steel.steel_section import SteelSection
except:
    from eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
    from eurocodes.en1994 import en1994_1_1
    from sections.steel.steel_section import SteelSection

class ConcreteSlab:
    """ Simple concrete slab with rectangular cross-section
        
    
    """
    
    def __init__(self,thickness,material):
        """ Constructor
            input:
                thickness [mm] thickness of slab
                material .. Concrete class material
        """
        
        self.hc = thickness
        self.material = material
    
class CompositeSlab:
    """ Simple composite slab with rectangular cross-section
        
    
    """
    
    def __init__(self,h_concrete,concrete,h_steel):
        """ Constructor
            input:
                thickness [mm] thickness of slab
                material .. Concrete class material
        """
        
        self.hc = h_concrete
        self.material = concrete
        self.ha = h_steel
        
    @property
    def h(self):
        """ Total height of the slab """
        return self.hc + self.ha

class CompositeIBeam:
    """ Class for composite beams with an I section as
        steel part and a slab (composite or concrete) on top
    """
    
    def __init__(self,length,steel_part,slab,studs=None):
        """ Constructor
            Input:
                steel_part: ISection class profile or WI profile steel section
                slab: concrete or composite slab
        """
        
        self.length = length
        self.steel = steel_part
        self.slab = slab
        self.studs = studs
        self.b0 = 0        
    
    @property 
    def beff(self):
        """ Effective width of concrete flange """
        return self.b0 + 2.0*self.length/8.0
        
    @property
    def yc(self):
        """ Distance of centroid of concrete slab from the top
            of the composite section
        """
        return 0.5*self.slab.hc
    
    @property
    def ya(self):
        """ Distance of centroid of steel part from the top
            of the composite section
        """
        return self.slab.h + 0.5*self.steel.h
    
    @property
    def Aa(self):
        """ Area of steel part """
        return self.steel.A
    
    @property
    def Ia(self):
        """ Second moment of area of steel part with respect to
            the major axis at centroid of steel part
        """
        return self.steel.I[0]
    
    @property
    def Ea(self):
        """ Young's modulus of steel part """
        return self.steel.E
    
    @property
    def EAa(self):
        """ Axial stiffness of steel part """
        return self.Ea*self.Aa
    
    @property
    def EIa(self):
        """ Bending stiffness of steel part 
            with respect to the neutral axis of steel part
        """
        return self.Ea*self.Ia
    
    @property
    def Ecm(self):
        """ Young's modulus of concrete """
        return self.slab.material.Ecm
    
    @property
    def EIc(self):
        """ Bending stiffness of concrete part 
            with respect to the neutral axis of concrete part
        """
        return self.Ecm*self.Ic
    
    @property
    def Ac(self):
        """ Area of concrete part """
        return self.slab.hc*self.beff
    
    @property
    def Ic(self):
        """ Second moment of area of concrete part with respect to
            the major axis at centroid of concrete part
        """
        return 1/12*self.beff*self.slab.hc**3
    
    @property
    def EAc(self):
        """ Axial stiffness of concrete """
        return self.Ecm*self.Ac
        
    
    def centroid(self):
        """ Centroid of the composite section 
            Assume full interaction
        """
        return (self.EAa*self.ya + self.EAc*self.yc)/self.EAcom()
        
    
    def EAcom(self):
        """ axial stiffness composite section """
        return self.EAa + self.EAc


    def ya_com(self):
        """ distance of centroid of steel part from the neutral
            axis of the composite section
        """
        
        return abs(self.ya-self.centroid())
    
    def yc_com(self):
        """ distance of centroid of concrete part from the neutral
            axis of the composite section
        """
        
        return abs(self.yc-self.centroid())

    def EIcom(self):
        """ sending stiffness of composite section """        
        return self.Ea*(self.Ia+self.ya_com()**2*self.Aa) + self.Ecm*(self.Ic+self.yc_com()**2*self.Ac)

    def distance_centroids(self):
        """ distance between centroids of steel and concrete parts """
        return self.ya - self.yc
    
    def composite_coefficient(self):
        """ Coefficint alpha that indicates the degree of composite action 
            TRY/BY, pp. 24        
        """
        a = self.distance_centroids()
        return a**2*self.EAa*self.EAc/(self.EAa+self.EAc)/(self.EIa + self.EIc)
        

if __name__ == '__main__':
    
    from sections.steel.ISection import IPE
    from eurocodes.en1992 import en1992_1_1, constants
    
    steel = IPE(300)
    #slab = ConcreteSlab(140,constants.Concrete("C25/30"))
    slab = CompositeSlab(100,constants.Concrete("C25/30"),50)
    
    p = CompositeIBeam(12000,IPE(450),slab)
    
    print("Effective width of concrete flange {0:4.2f} mm".format(p.beff))    
    print("Centroid of concrete {0:4.2f} mm".format(p.yc))
    print("Centroid of steel part {0:4.2f} mm".format(p.ya))
    print("Centroid {0:4.2f} mm".format(p.centroid()))
    print("Bending stiffness of steel part: {0:4.2f} MNm".format(p.EIa*1e-12))
    print("Bending stiffness of concrete part: ",p.EIc*1e-12)
    print("Distance between centroids [mm]: ", p.distance_centroids())
    print("Liittovaikutuskerroin: ", p.composite_coefficient())
    print("Bending stiffness of composite section: ",p.EIcom()*1e-12)
    
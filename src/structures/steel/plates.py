# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 20:01:30 2019

Classes for various steel plates

@author: kmela
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

from cost.workshop import Workshop
from cost.cost_data import BASIC_STEEL_PRICE, steel_grade_add_on, thickness_add_on

from materials.steel_data import Steel

""" Default thicknesses from SSAB """
THICKNESSES = [4, 5, 6, 8, 10, 12, 14, 15, 16, 18, 20, 20, 22, 25, 30]

class RectPlate:
    """ Rectangular plate """
    
    def __init__(self,width,depth,thickness,material="S355"):
        """ Constructor
            Input: Dimensions of the plate and material
        
        """
        
        self.b = width
        self.h = depth
        self.t = thickness
        self.material = Steel(material)
        
        self.weld_size = 0.0
        
        # List of hole diameters
        self.holes = []
        
        self._cost = {'material':0.0,
                     'cutting':0.0,
                     'painting':0.0,
                     'drilling':0.0,
                     'welding':0.0,
                     'holes':0.0}
    
    @property
    def fy(self):
        return self.material.fy
    
    @property
    def fu(self):
        return self.material.fu
    
    def area(self):
        return self.h*self.b
    
    def end_area(self):
        return self.b*self.t
    
    def side_area(self):
        return self.h*self.t
    
    def volume(self):
        return self.h*self.b*self.t
    
    def circumference(self):
        return 2*(self.h+self.b)
    
    def paint_area(self):
        """ Area of the plate to be painted. 
            Assume that both sides are painted.
        """
        
        return 2*self.area() + self.t*self.circumference()
    
    def weight(self):
        
        return self.material.rho*self.volume()
    
    def cost(self,workshop,material_cost=BASIC_STEEL_PRICE,verb=False):
        """ Evaluate the cost of the plate
            units: e
        """
        
        self._cost['material'] = self.weight()*(material_cost + 
                                                steel_grade_add_on(self.material.name) + 
                                                thickness_add_on(self.t))*1e-3
        
        self._cost['cutting'] = workshop.cost_centres['cutting'].cost(self.t,self.circumference())
        self._cost['painting'] = workshop.cost_centres['painting'].cost(self.paint_area())
        self._cost['blasting'] = workshop.cost_centres['blasting'].cost(max(self.b,self.h))
        if self.weld_size > 0:
            self._cost['welding'] = workshop.cost_centres['assembly_welding'].cost(self.weld_size,self.circumference())
        
        
        if len(self.holes) > 0:
            total_holes = sum([d*math.pi for d in self.holes]) 
            self._cost['holes'] = workshop.cost_centres['cutting'].cost(self.t,total_holes)
        
        total_cost = 0.0
        for c in self._cost.values():
            total_cost += c
        
        self._total_cost = total_cost
        
        if verb:
            self.cost_distribution()
        
        return total_cost
    
    def cost_distribution(self,pie=False):
        """ Print cost distribution """
        
        if pie:
            pass
        else:
            print(""" Cost distribution of the plate """)
            p_tot = 0.0
            for key, value in self._cost.items():
                p =  value/self._total_cost*100
                p_tot += p
                print("   {0:10s}: {1:5.2f} € [{2:5.2f} %]".format(key,value,p))
            
            print("-----------------------")
            print("   {0:10s}: {1:5.2f} € [{2:5.2f} %]".format("Total",self._total_cost,p_tot))
                

if __name__ == '__main__':
    
    ws = Workshop()
    plate = RectPlate(230,380,25,'S355')
    
    plate.holes = [26,26,26,26,26,26]
    
    plate.cost(ws,verb=True)
    #print("Plate cost: {0:4.2f} €".format(plate.cost(ws)))
    #plate.cost_distribution()
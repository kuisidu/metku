# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Wed Oct  2 20:01:30 2019

Classes for various steel plates

@author: kmela
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

from metku.cost.workshop import Workshop
from metku.cost.cost_data import BASIC_STEEL_PRICE, steel_grade_add_on, thickness_add_on

from metku.materials.steel_data import Steel

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
    
    def draw(self,x0=(0,0),axes=None):
        """ Draw the plate """
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes

        plate = patches.Rectangle(x0, self.b, self.h,edgecolor='black',facecolor='white')
        
        ax.add_patch(plate)
        ax.axis('equal')
        
        plt.show()
        
        return ax

class RectPlateWithHoles(RectPlate):
    """ Class for rectangular steel plates with a regular pattern of circular holes """
    
    def __init__(self,width,depth,thickness,d0,x0,px=0,py=0,n1=3,n2=1,material="S355"):
        """

        Parameters
        ----------
        width : double
            width of the plate.
        depth : double
            height of the plate.
        thickness : double
            plate thickness.
        d0 : double
            bolt hole diameter
        x0 : numpy array or list
            location of the first hole centroid.
            x0[0] .. horizontal coordinate
            x0[1] .. vertical coordinate
        px : double
            distance between hole centroids, horizontal direction
        py:
            distance between hole centroids, vertical direction 
        n1 : integer, optional
            Number of bolts in vertical direction. The default is 3.
        n2 : integer, optional
            number of vertical bolt rows. The default is 1.
        material : string, optional
            Steel grade of the plate. The default is "S355".

        Returns
        -------
        None.

        """
        super().__init__(width,depth,thickness,material)
        
        self.x0 = x0
        self.d0 = d0
        self.n1 = n1
        self.n2 = n2
        self.px = px
        self.py = py
    
    

    @property
    def n(self):
        """ Number of bolts """
        return self.n1*self.n2
        
    
    def eLeft(self):
        """ Edge distance to the left edge """
        return self.x0[0]
    
    def eRight(self):
        """ Edge distance to the right edge """
        return self.b - self.x0[0] - (self.n2-1)*self.px
    def eTop(self):
        """ Edge distance to the top edge """
        return self.h - self.x0[1] - (self.n1-1)*self.py
    
    
    
    def draw(self,x0=(0,0),axes=None):
        
        ax = super().draw(x0,axes)
        
        # Draw holes
        for j in range(self.n1):
            print(j)
            for i in range(self.n2):
                x = x0[0] + self.x0[0] + i*self.px
                y = x0[1] + self.x0[1] + j*self.py
                hole = patches.Circle((x,y), radius= 0.5*self.d0, edgecolor='black',facecolor='white')
                ax.add_patch(hole)
                

if __name__ == '__main__':
    
    ws = Workshop()
    plate = RectPlate(230,380,25,'S355')
    
    plate.holes = [26,26,26,26,26,26]
    
    #plate.draw()
    
    plate = RectPlateWithHoles(170,230,10,22, [60,45],55,70,3,2)
    plate.draw()
    
    #plate.cost(ws,verb=True)
    #print("Plate cost: {0:4.2f} €".format(plate.cost(ws)))
    #plate.cost_distribution()
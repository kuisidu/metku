# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:23:00 2018

@author: huuskoj
"""
import copy
from decimal import Decimal
from truss2d import Truss2D, TopChord, BottomChord, TrussWeb

def test1():
    truss = Truss2D()
    
    coord1 = [[0,1], [5,2]]
    coord2 = [[0,0], [5,0]]
    
    top_chord = TopChord(coord1)
    bottom_chord = BottomChord(coord2)
    
    members = 3
    
    
    dx = coord1[1][0]/(members-1)
    y = lambda x: x*(coord1[1][1]- coord1[0][1]) / (coord1[1][0]- coord1[0][0])+coord1[0][1]
    c1 = copy.copy(coord1[0])
    c2 = copy.copy(coord2[0])
    for i in range(members):       
        if i%2 == 0:
            if i!= 0:               
                c1[0] += float(Decimal(2*dx))
                c1[1] = y(c1[0])
        else:
            c2[0] += float(Decimal(2*dx))
        if c1[0] > coord1[1][0] or c1[1] > coord1[1][1]:
            c1 = copy.copy(coord1[1])
        if c2[0] > coord2[1][0] or c2[1] > coord2[1][1]:
            c2 = copy.copy(coord2[1])
            

        member = TopChord([c1, c2])
        truss.add(member)
    
    
    
    truss.add(top_chord)
    truss.add(bottom_chord)
    truss.plot()
    truss.generate()
    truss.f.draw()
    
if __name__ == '__main__':
    test1()
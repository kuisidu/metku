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
    web1 = TrussWeb([top_chord, 0.9], [bottom_chord, 0.5])

    truss.add(top_chord)
    truss.add(bottom_chord)
    truss.add(web1)
    truss.plot()
    truss.generate()
    truss.f.draw()
    
if __name__ == '__main__':
    test1()
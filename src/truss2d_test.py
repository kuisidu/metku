# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:23:00 2018

@author: huuskoj
"""
import copy
from decimal import Decimal
from truss2d import Truss2D, TopChord, BottomChord, TrussWeb, TrussJoint
from src.frame2d.frame2d import LineLoad, XYHingedSupport, Frame2D, SteelColumn, FixedSupport

def test1():
    truss = Truss2D(num_elements=1)
    
    coord1 = [[0,10], [8,4]]
    coord2 = [[-2,-2], [8,0]]
    
    top_chord = TopChord(coord1)
    bottom_chord = BottomChord(coord2)
    truss.add(top_chord)
    truss.add(bottom_chord)

    members = 5

    c1 = 0
    c2 = 0

    for i in range(members+2):
        if i%2 == 0:
            c1 = i/members
            if c1 > 1:
                c1 = 1
        elif i!=0:
            c2 = i/members
            if c2 > 1:
                c2 = 1
        truss.add(TrussWeb(c1, c2))
    truss.add(LineLoad(top_chord, [-50, -50], 'y'))
    truss.add(XYHingedSupport(bottom_chord.local(0)))
    truss.add(XYHingedSupport(bottom_chord.local(1)))
    truss.generate()
    truss.calculate()
    truss.plot()

def test2():
    frame = Frame2D()

    col1 = SteelColumn([[0,0], [0,2]])
    col2 = SteelColumn([[4,0], [4,2]])

    frame.add(col1)
    frame.add(col2)
    frame.add(FixedSupport([0,0]))
    frame.add(FixedSupport([4,0]))

    frame.generate()

    truss = Truss2D(fem=frame.f)
    coord1 = [[-2, 6], [4, 4]]
    coord2 = [[0, 2], [4, 2]]

    top_chord = TopChord(coord1)
    bottom_chord = BottomChord(coord2)
    truss.add(top_chord)
    truss.add(bottom_chord)

    members = 7

    c1 = 0
    c2 = 0

    for i in range(members + 2):
        if i % 2 == 0:
            c1 = i / members
            if c1 > 1:
                c1 = 1
        elif i != 0:
            c2 = i / members
            if c2 > 1:
                c2 = 1
        truss.add(TrussWeb(c1, c2))
    truss.add(LineLoad(top_chord, [-50, -50], 'y'))
    truss.generate()
    frame.add(truss)
    frame.plot()


if __name__ == '__main__':
    test2()
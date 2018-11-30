# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:23:00 2018

@author: huuskoj
"""
from truss2d import Truss2D, TopChord, BottomChord, TrussWeb, TrussJoint
from frame2d.frame2d import LineLoad,PointLoad, XYHingedSupport, Frame2D, SteelColumn, FixedSupport, YHingedSupport
import matplotlib.pyplot as plt


def test1():
    truss = Truss2D(num_elements=1)
    
    coord1 = [[1,2], [2,2]]
    coord2 = [[0,0], [2,0]]
    
    coord3 = [[2,2], [4,1]]
    coord4 = [[2,0], [4,0]]
    
    top_chord = TopChord(coord3)
    bottom_chord = BottomChord(coord4)
    
    
    top_chord2 = TopChord(coord3)
    bottom_chord2 = BottomChord(coord4)
    
    truss.add(top_chord)
    truss.add(bottom_chord)
    truss.plot()
    members = 2
    c1 = 0
    c2 = 0
    flip = False
    
    
    for i in range(members+2):
        if i%2 == 0:
            c1 = i /members
            if c1 > 1:
                c1 = 1
        elif i!=0:
            c2 = i /members
            if c2 > 1:
                c2 = 1
        if flip:
            truss.add(TrussWeb(c1, c2))
        else:
            truss.add(TrussWeb(c2, c1))
    truss.plot()
    

    #truss.add(TrussWeb(0,0))
    #truss.add(LineLoad(top_chord, [-100, -100], 'y'))
    #truss.add(PointLoad([1,0.5], [0,-50,0]))
    truss.add(YHingedSupport(bottom_chord.local(0)))
    truss.add(XYHingedSupport(bottom_chord.local(1)))
    #truss.plot()
    truss.generate()
    truss.f.draw()
    #truss.calculate()
    #truss.plot_normal_force()
    #top_chord.coordinates = [[0,2],[2,1]]

    truss.f.draw()
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

    truss = Truss2D(fem=frame.f, num_elements=1)
    coord1 = [[0, 2], [4, 2]]
    coord2 = [[0.2, 1.5], [3.8, 1.5]]

    top_chord = TopChord(coord1)
    bottom_chord = BottomChord(coord2)
    truss.add(top_chord)
    truss.add(bottom_chord)

    members = 6

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
    truss.add(LineLoad(top_chord, [-150, -150], 'y'))   
    truss.generate()
    frame.add(truss)
    frame.calculate()
    frame.plot()
    frame.f.draw()
    #frame.plot_normal_force(show=False)
    #truss.plot_normal_force()


if __name__ == '__main__':

    test1()
"""
truss = Truss2D(num_elements=1)

coord1 = [[0,4], [2,1]]
coord2 = [[0,0], [2,-0]]

top_chord = TopChord(coord1)
bottom_chord = BottomChord(coord2)
truss.add(top_chord)
truss.add(bottom_chord)
"""


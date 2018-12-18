# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:23:00 2018

@author: huuskoj
"""
from truss2d import Truss2D, TopChord, BottomChord, TrussWeb, TrussJoint
from frame2d.frame2d import LineLoad,PointLoad, XYHingedSupport, Frame2D, SteelColumn, FixedSupport, YHingedSupport
from eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint
from sections.steel.RHS import RHS
import matplotlib.pyplot as plt
import numpy as np

def test1():
    truss = Truss2D(num_elements=1)
    
    coord1 = [[0,2], [2,4]]
    coord2 = [[0,0], [4,0]]
    
    coord3 = [[2,4], [4,2]]
    coord4 = [[2,0], [4,0]]
    
    top_chord = TopChord(coord1)
    bottom_chord = BottomChord(coord2)
    
    
    top_chord2 = TopChord(coord3)
    bottom_chord2 = BottomChord(coord4)
    
    truss.add(top_chord)
    truss.add(bottom_chord)
    truss.add(top_chord2)
    truss.add(TrussWeb([0, 0], [1,3], 'global'))
    truss.add(TrussWeb([1, 0], [1,3], 'global'))
    truss.add(TrussWeb([3, 0], [4,2], 'global'))
    
    members = 4
    c1 = 0
    c2 = 0
    flip = False
    
    """
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
    """
    #truss.add(TrussWeb(0,0))
    #truss.add(LineLoad(top_chord, [-100, -100], 'y'))
    #truss.add(PointLoad([1,0.5], [0,-50,0]))
    truss.add(YHingedSupport(bottom_chord.local(0)))
    truss.add(XYHingedSupport(bottom_chord.local(1)))
    #truss.plot()
    truss.generate()
    truss.calculate()
    truss.plot_normal_force()
    #top_chord.coordinates = [[0,2],[2,1]]
    truss.f.draw()
    truss.plot()
    
def test2():
    """
    Joint test
    example 3.4 from SSAB domex tube rakenneputket p.178
    """
    N0Ed = -1364 # kN
    N1Ed = -600 # kN
    N2Ed = 600 # kN

    profile0 = RHS(200,200,8, 420)
    profile1 = RHS(150,150,6, 420)
    profile2 = RHS(150,150,6, 420)
    
    profile0.Ned = N0Ed
    profile1.Ned = N1Ed
    profile2.Ned = N2Ed
     
    theta1 = 45
    theta2 = 45
    gap = 28 # mm
    
    joint = RHSKGapJoint(profile0, [profile1, profile2], [theta1, theta2], gap)
    
    # Correct answers
    beta = 0.75
    gamma = 12.5
    n = 0.5482
    kn = 1.0
    e = 20.1
    AV0 = 3584
    
    print("Checking results: ")
    print("beta: ", joint.beta() , " = ", beta)
    print("gamma: ", joint.gamma(), " = ", gamma)
    print("n: ", joint.eval_N(), " = ", n)
    print("kn: ",joint.eval_KN(), " = ", kn)
    print("e: ", joint.eccentricity(), " = ", e)
    print("Av0: ", joint.a_eta() , " = ", AV0)
    print("Chord face failure: ", joint.chord_face_failure(), " = 807.4 kN")
    print("Punching shear: ", joint.punching_shear(), " = 1566 kN")
    #print("Chord shear: ", joint.chord_shear())
    print("Brace failure: ", joint.brace_failure())
    
    
    
def test3():
    """
    Portal truss-frame
    """
    X = np.array([0,0,0.1714,0.3429, 0.5143, 0.6857, 0.8571, 1.0286, 1.2,
                  1.3714, 1.5429, 1.7143, 1.8857, 2.0571, 2.2286, 2.4, 2.4, 0, 2.4])
    Y = np.array([0, -0.18, 0.0086, -0.18, 0.0257, -0.18, 0.0429, -0.18, 0.06, -0.18,
                  0.0428, -0.18, 0.0257, -0.18, 0.0086, -0.18, 0, -0.98, -0.98])
    
    X = np.array([0,01.714,03.429, 05.143, 06.857, 08.571, 10.286, 12.,
                  13.714, 15.429, 17.143, 18.857, 20.571, 22.286, 24])
    Y = np.array([-01.8, 00.086, -01.8, 00.257, -01.8, 00.429, -01.8, 00.6, -01.8,
                  00.428, -01.8, 00.257, -01.8, 00.086, -01.8])

    frame = Frame2D(num_elements=10)
    col1 = SteelColumn([[0,0],[0, -9.8]])
    col2 = SteelColumn([[24,0],[24, -9.8]])
    
    frame.add(col1)
    frame.add(col2)      
    
    truss = Truss2D(num_elements=10, fem=frame.f)
    
    bottom_chord = BottomChord([[0, -1.8],[24,-1.8]], material="S420", profile="RHS 140x140x6")
    top_chord1 = TopChord([[0,0], [12,0.6]], material="S420", profile="RHS 120x120x6")
    top_chord2 = TopChord([[12,0.6], [24,0]], material="S420", profile="RHS 120x120x6")
    truss.add(bottom_chord)
    truss.add(top_chord1)
    truss.add(top_chord2)
    for i in range(len(X)-1):
        truss.add(TrussWeb([X[i], Y[i]], [X[i+1], Y[i+1]], 'global', profile="Rhs 90x90x4"))
        
    truss.plot(print_text=False)
    #truss.generate()
    truss.f.draw()
    
    frame.add(truss)
    # SUPPORTS
    frame.add(FixedSupport([0, -9.8]))
    frame.add(FixedSupport([24, -9.8]))
    # LOADS
    frame.add(LineLoad(bottom_chord, [-0.69, -0.69], 'y'))
    frame.add(LineLoad(top_chord1, [-21.45, -21.45], 'y'))
    frame.add(LineLoad(top_chord2, [-21.45, -21.45], 'y'))
    frame.add(LineLoad(col1, [3.51, 3.51], 'x'))
    frame.add(LineLoad(col2, [0.17, 0.17], 'x'))
    
    
    frame.generate()
    frame.calculate()
    frame.plot_deflection()
    frame.plot_normal_force()
    
    
    
def test4():
    truss = Truss2D(num_elements=3)
    coord1 = [[0,2], [4,2]]
    coord2 = [[0,0], [4,0]]  
    
    top_chord = TopChord(coord1, profile="RHS 200X200X8")
    bottom_chord = BottomChord(coord2, profile="RHS 200X200X8")  
    
    truss.add(top_chord)
    truss.add(bottom_chord)
    
    members = 4
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
            truss.add(TrussWeb(c1, c2, profile="RHS 150x150x6"))
        else:
            truss.add(TrussWeb(c2, c1, profile="RHS 150x150x6"))
    
    #truss.add(TrussWeb(0,0))
    truss.add(LineLoad(top_chord, [-250, -250], 'y'))
    truss.add(LineLoad(bottom_chord, [-250, -250], 'y'))
    #truss.add(PointLoad([1,0.5], [0,-50,0]))
    truss.add(YHingedSupport(bottom_chord.local(0)))
    truss.add(XYHingedSupport(bottom_chord.local(1)))
    truss.generate()
    truss.calculate()
    truss.plot()
    truss.plot_normal_force()
    
    for joint in truss.joints.values():
        if joint.joint_type in ["K", "N"]:
            print(joint.jid)
            print(joint.calculate())
    
    
    
    
if __name__ == '__main__':

    test3()
"""
truss = Truss2D(num_elements=1)

coord1 = [[0,4], [2,1]]
coord2 = [[0,0], [2,-0]]

top_chord = TopChord(coord1)
bottom_chord = BottomChord(coord2)
truss.add(top_chord)
truss.add(bottom_chord)
"""


# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:23:00 2018

@author: huuskoj
"""
from truss2d import Truss2D, TopChord, BottomChord, TrussWeb, TrussJoint, TrussColumn
from frame2d.frame2d import SteelBeam, LineLoad,PointLoad, XYHingedSupport, Frame2D, SteelColumn, FixedSupport, YHingedSupport
#from frame2d.frameoptimizer import FrameOptimizer
from eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint
from sections.steel.RHS import RHS
import matplotlib.pyplot as plt
import numpy as np
import math
import robot3D

from framefem import FrameFEM

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
    #print("Chord shear: ", joint.chord_shear(), " = 1106 kN")
    print("Brace failure: ", joint.brace_failure(), " = 1148 kN")
    
    
    
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
    
    truss = Truss2D(num_elements=10)
    
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
    truss.generate()
    truss.plot()
    """
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
    """
    
    
def test4():
    truss = Truss2D(num_elements=3)
    coord1 = [[0,2], [2,4]]
    coord2 = [[0,0], [2,0]]  
    
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
            truss.add(TrussWeb(c1, c2, profile="RHS 100x100x5"))
        else:
            truss.add(TrussWeb(c2, c1, profile="RHS 100x100x5"))
    
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
    
    truss.joints[4].coordinate = [1.6,0]
    truss.joints[4].plot_joint(0.05)
    truss.joints[4].cnode.Fy
    print(truss.joints[4].coordinate)
    truss.plot()


def portal_frame():
    X = np.array([0, 1.714, 3.429, 5.143, 6.857, 8.571, 10.286, 12.,
                  13.714, 15.429, 17.143, 18.857, 20.571, 22.286, 24])
    Y = np.array([-01.8, 0.086, -1.8, 0.257, -1.8, 0.429, -1.8, 0.6,
                  -1.8, 0.429, -1.8, 0.257, -1.8, 0.086, -1.8])
    
    # INITIALIZE EMPTY FRAME
    frame = Frame2D(num_elements=4)
    # COLUMNS
    col1 = SteelColumn([[0,0],[0, -9.8]], profile="RHS 220x220x7.1")
    col2 = SteelColumn([[24,0],[24, -9.8]], profile="RHS 220x220x7.1")
    frame.add(col1)
    frame.add(col2)
    # SUPPORTS
    frame.add(FixedSupport([0, -9.8]))
    frame.add(FixedSupport([24, -9.8]))

    truss = Truss2D(num_elements=40, fem_model=frame.f)
    bottom_chord = BottomChord([[0, -1.8],[24,-1.8]], material="S420", profile="RHS 140x140x6")
    top_chord1 = TopChord([[0,0], [12,0.6]], material="S420", profile="rhs 120x120x6")
    top_chord2 = TopChord([[12,0.6], [24,0]], material="S420", profile="rhs 120x120x6")
    truss.add(bottom_chord)

    truss.add(top_chord1)
    truss.add(top_chord2)
    for i in range(len(X)-1):
        truss.add(TrussWeb([X[i], Y[i]], [X[i+1], Y[i+1]], 'global', profile="Rhs 90x90x4"))

    # ADD TRUSS TO FRAME
    frame.add(truss)
    # LOADS
    frame.add(LineLoad(bottom_chord, [-0.69, -0.69], 'y'))
    frame.add(LineLoad(top_chord1, [-21.45, -21.45], 'y'))
    frame.add(LineLoad(top_chord2, [-21.45, -21.45], 'y'))
    frame.add(LineLoad(col1, [3.51, 3.51], 'x'))
    frame.add(LineLoad(col2, [0.17, 0.17], 'x'))
   
    # GENERATE
    frame.generate()
    #frame.plot(print_text=False)
    frame.f.draw()
    # Calculate
    frame.calculate()
    frame.plot_normal_force()
    frame.bmd(30)
    
    frame.to_robot('portal_frame')
    
def portal_frame2():
    
    print("PORTAL FRAME 2")
    # INITIALIZE EMPTY FRAME
    frame = Frame2D(num_elements=2, fem_model=FrameFEM())
    # COLUMNS
    col1 = SteelColumn([[0,0],[0, 5]], profile="RHS 200x200x7.1")
    col2 = SteelColumn([[10,0],[10, 5]], profile="RHS 200x200x7.1")
    frame.add(col1)
    frame.add(col2)
    # SUPPORTS
    frame.add(FixedSupport([0, 0]))
    frame.add(FixedSupport([10, 0]))
    
    truss = Truss2D(num_elements=4, fem_model=frame.f)
    
    coord1 = [[0,5], [10,5]]
    coord2 = [[0,2], [10,2]]   

    top_chord = TopChord(coord1)
    bottom_chord = BottomChord(coord2)
    
    truss.add(top_chord)
    truss.add(bottom_chord)
    
    
    members = 1
    c1 = 0.0
    c2 = 0.
    flip = False

    
    for i in range(1,members+1):
        if i%2 == 0:
            c1 = round(i/2 /members, 4)
            if c1 > 1:
                c1 = 1
        elif i!=0:
            c2 = round(i/2 /members, 4)
            if c2 > 1:
                c2 = 1
        if flip:
            truss.add(TrussWeb(c1, c2))
        else:
            truss.add(TrussWeb(c2, c1, Sj1=np.inf, Sj2=np.inf, num_elements=2))
            truss.add(TrussWeb(1-c2, 1-c1, Sj1=np.inf, Sj2=np.inf, num_elements=2))
    """

    lista = [[0.25, 0],
             [0.25, 0.5],
             [0.75, 0.5],
             [0.75, 1]]
        
    
    lista = [[0.2, 0.1],
            [0.8, 0.9],
            [0.2, 0.3],
            [0.5, 0.3],
            [0.5, 0.7],
            [0.8, 0.7]]
    
    
    random.shuffle(lista)
    print("LISTA: ", lista)
    
    for val in lista:
        a, b = val
        truss.add(TrussWeb(a, b))
    """
    
    frame.add(truss)
    
    #frame.add(PointLoad([5, 2], [0, -1000, 0]))
    frame.add(LineLoad(top_chord, [-10, -10], 'y'))
    #frame.plot()
    
    frame.generate()
    frame.calculate()
    frame.f.draw()
    #frame.plot_normal_force()
    
    #print(K)
    
    frame.bmd(50)
    frame.to_robot('portal_test')

    """
    print("ECC")
    ecc_elems = []
    for member in frame.members.values():
        for elem in member.ecc_elements.values():
            print([node.x for node in elem.nodes])
    
    #print("JOINTS")
    coords = []
    for joint in truss.joints.values():
        #print(joint.jid, "  ", joint.coordinate, "  ", joint.loc)
        coords.extend(joint.nodal_coordinates)
    coords.sort()
    #print(coords)
    elements = []
    print(" \n ELEMENTS \n")
    for i, elem in enumerate(frame.f.elements):
        print(i, elem.bending_moment)
        elements.extend(elem.axial_force)

    return frame.f.nodal_coords, coords
    """
    for mem in frame.members.values():
        print(mem.profile, "  ", mem.mem_id)
    for i, elem in enumerate(frame.f.elements):
        print(i, elem.section.Iy)
def frame_test():
    
    frame = Frame2D()
    c1 = [[0,0], [0,2]]
    c2 = [[5,0], [5,2]]
    c3 = [[0,2], [0.1, 2]]
    c4 = [[0.1, 2], [4.9, 2]]
    c5 = [[4.9,2], [5,2]]
    
    col1 = SteelColumn(c1)
    col2 = SteelColumn(c2)
    
    ecc_beam1 = SteelBeam(c3, profile="he 1000 A", num_elements=3)
    #ecc_beam1.Sj2 = 0
    ecc_beam2 = SteelBeam(c5, profile="he 1000 A", num_elements=3)
    #ecc_beam2.Sj1 = 0
    
    beam = SteelBeam(c4)
    beam.Sj1 = 0
    beam.Sj2 = 0
    
    
    frame.add(col1)
    frame.add(col2)
    frame.add(beam)
    frame.add(ecc_beam1)
    frame.add(ecc_beam2)
    
    frame.add(LineLoad(beam, [-10, -10], 'y'))
    frame.generate()
    frame.plot()
    
    frame.f.draw()
    frame.calculate()
    frame.bmd(10)
    
    
def truss_test():

    frame = Frame2D()
    c1 = [[0,0], [0,5]]
    c2 = [[10,0], [10,5]]

    
    col1 = SteelColumn(c1)
    col2 = SteelColumn(c2)
    frame.add(col1)
    frame.add(col2)
    
    
    truss = Truss2D(num_elements=4, fem_model=frame.f)
    
    coord1 = [[0,5], [10,5]]
    coord2 = [[0,2], [10,2]]   

    top_chord = TopChord(coord1)
    bottom_chord = BottomChord(coord2)
 
    
    truss.add(top_chord)
    truss.add(bottom_chord)
    #truss.add(top_chord2)
    
    members = 2
    c1 = 0.0
    c2 = 0.0
    flip = True
    
    web1 = TrussWeb(0.5, 0.5, Sj1=np.inf, Sj2=np.inf)
    web2 = TrussWeb(.5, 0.5, Sj1=np.inf, Sj2=np.inf)
    
    
    
    
    #truss.add(web2)
    #truss.add(web1)
    
    
    for i in range(1,members+1):
        if i%2 == 0:
            c1 = round(i/2 /members, 4)
            if c1 > 1:
                c1 = 1
        elif i!=0:
            c2 = round(i/2 /members, 4)
            if c2 > 1:
                c2 = 1
        if flip:
            truss.add(TrussWeb(c1, c2, Sj1=np.inf, Sj2=np.inf))
            truss.add(TrussWeb(1-c1, 1-c2, Sj1=np.inf, Sj2=np.inf))
        else:
            truss.add(TrussWeb(c2, c1, Sj1=np.inf, Sj2=np.inf))
            truss.add(TrussWeb(1-c2, 1-c1, Sj1=np.inf, Sj2=np.inf))
            
    
    """
    frame.add(FixedSupport([0,0]))
    frame.add(FixedSupport([10,0]))
    
    frame.add(LineLoad(top_chord, [-10, -10], 'y'))
    
    
    frame.add(truss)
    frame.generate()
    frame.calculate()
    frame.bmd(10)
    frame.to_robot('truss_test')
    """
    truss.add(FixedSupport([0,2]))
    truss.add(FixedSupport([10,2]))
    
    truss.add(LineLoad(top_chord, [-10, -10], 'y'))
    truss.generate()
 
    truss.calculate()
    truss.f.draw()
    truss.bmd(10)
    truss.plot_normal_force()
    truss.to_robot('truss_test')
    
    
    
   
    

def frame_truss():
    
    frame = Frame2D(num_elements=1)
    
    top_coord = [[0,5], [10,5]]
    bot_coord = [[0.1, 2], [9.9, 2]]
    
    col1_coord = [[0,0], [0, 4.95]]
    col2_coord = [[10, 0], [10, 4.95]]
    
    # Ecc to frame
    ecc1_coord = [[0,4.95], [0, 5]]
    ecc2_coord = [[0, 2], [0.1, 2]]
    ecc3_coord = [[10, 4.95], [10, 5]]
    ecc4_coord = [[9.9, 2], [10, 2]]
    
    web1_coord = [[0.05, 4.95], [4.95, 2.05]]
    web2_coord = [[5.05, 2.05], [9.95, 4.95]]
    # Web's ecc
    ecc5_coord = [[0.05, 5], [0.05, 4.95]]
    ecc6_coord = [[4.95, 2], [4.95, 2.05]]
    ecc7_coord = [[5.05, 2], [5.05, 2.05]]
    ecc8_coord = [[9.95, 4.95], [9.95, 5]]
    

    
    col1 = SteelColumn(col1_coord, profile="RHS 200x200x7.1")
    col2 = SteelColumn(col2_coord, profile="RHS 200x200x7.1")
    
    top_chord = SteelBeam(top_coord, profile="RHS 100x100x5")
    bottom_chord = SteelBeam(bot_coord, profile="RHS 100x100x5")
    
    web1 = SteelBeam(web1_coord, profile="RHS 50x50x2")
    web2 = SteelBeam(web2_coord, profile="RHS 50x50x2")
    
    ecc1 = SteelColumn(ecc1_coord, profile="he 1000 a")
    ecc2 = SteelColumn(ecc2_coord, profile="he 1000 a")
    ecc3 = SteelColumn(ecc3_coord, profile="he 1000 a")
    ecc4 = SteelColumn(ecc4_coord, profile="he 1000 a")
    ecc5 = SteelColumn(ecc5_coord, profile="he 1000 a")
    ecc6 = SteelColumn(ecc6_coord, profile="he 1000 a")
    ecc7 = SteelColumn(ecc7_coord, profile="he 1000 a")
    ecc8 = SteelColumn(ecc8_coord, profile="he 1000 a")
    
    
    #web3 = SteelBeam([[1, 2], [1.5, 1]], profile="RHS 50x50x2")
    #web4 = SteelBeam([[1.5, 1], [2, 2]], profile="RHS 50x50x2")
    
    frame.add(col1)
    frame.add(col2)
    frame.add(top_chord)
    frame.add(bottom_chord)
    frame.add(web1)
    frame.add(web2)
    #frame.add(web3)
    #frame.add(web4)
    frame.add(ecc1)
    frame.add(ecc2)
    frame.add(ecc3)
    frame.add(ecc4)
    frame.add(ecc5)
    frame.add(ecc6)
    frame.add(ecc7)
    frame.add(ecc8)

    frame.add(LineLoad(top_chord, [-10, -10], 'y'))
    
    frame.add(FixedSupport([0,0]))
    frame.add(FixedSupport([10, 0]))
    
    frame.generate()
    frame.calculate()
    frame.plot()
    for mem in frame.members.values():
        print(mem.mem_id, "  ",mem.elements.keys())
    frame.f.draw()
    frame.bmd(50)
    #frame.plot_normal_force()

    #frame.to_robot("testi")



def simple_truss():
    
    H0 = 8.5
    H1 = 2
    H2 = 6
    H3 = 2
    L1 = 24
    L2 = L1
    n = 30
    
    frame = Frame2D(num_elements=20, fem_model=FrameFEM())
    # COLUMNS
    c1 = [0,0]
    col1 = SteelColumn([c1,[0, H0+H1]], profile="RHS 200x200x7.1")
    col2 = SteelColumn([[L1+L2,0],[L1+L2, H0+H3]], profile="RHS 200x200x7.1")
    frame.add(col1)
    frame.add(col2)
    # SUPPORTS
    frame.add(FixedSupport(c1))
    frame.add(FixedSupport([L1+L2, 0]))
    
    truss = Truss2D(simple=[H0, H1, H2,H3, L1,L2, n],num_elements=10, fem_model=frame.f)
    #bchord = truss.bottom_chords[0]
    #bchord.coordinates = [[0.5, H0], [L1 + L2 - 0.5, H0]]
    frame.add(truss)
    frame.add(LineLoad(truss.top_chords[0], [-100, -100], 'y'))
    frame.add(LineLoad(truss.top_chords[1], [-100, -100], 'y'))
    frame.add(PointLoad([0, H0+H1], [100, 0, 0]))
    #for web in truss.webs.values():
    #    web.Sj1 = 1e10
    #    web.Sj2 = 1e10
    #frame.add(PointLoad([L/2,H0+H], [0,-100,0]))
    frame.generate()
    #frame.f.draw()
    frame.calculate()
    frame.plot_deflection(10)
    
    #frame.bmd(30)
    #frame.plot_normal_force()
    #frame.plot_buckling(50, 10)
    frame.to_robot("simple_truss", num_frames=10, s=10)
    """
    # num elements affects results
    # the order in which web's are added affects the result
    truss = Truss2D(simple=[H0, H, L, n], num_elements=4)
    truss.add(LineLoad(truss.top_chords[0], [-100, -100], 'y'))
    truss.add(FixedSupport([0, H0]))
    truss.add(FixedSupport([L, H0]))
    truss.generate()
    truss.calculate()
    truss.f.draw()
    truss.bmd(1)
    truss.plot_normal_force()
    truss.to_robot("simple_truss")
    """
    
    
def optimization_test():
    
    H0 = 2
    H1 = 0.5
    H2 = 0.8
    H3 = H1
    L1 = 4
    L2 = 4
    n = 6
    
    frame = Frame2D()
    frame.add(SteelColumn([[0,0], [0,H0 + H1]]))
    frame.add(SteelColumn([[L1 + L2,0], [L1+L2,H0 + H1]]))
    frame.add(FixedSupport([0,0]))
    frame.add(FixedSupport([L1 + L2, 0]))

        
    
    truss = Truss2D(simple=[H0, H1, H2,H3, L1,L2, n], num_elements=2)
    
    truss.add(XYHingedSupport([0,H0+H1]))
    truss.add(XYHingedSupport([L1+L2,H0+H1]))
    
    truss.add(LineLoad(truss.top_chords[0], [-20, -20], 'y'))
    truss.add(LineLoad(truss.top_chords[1], [-20, -20], 'y'))
    
    frame.add(truss)
    frame.generate()
    
    lb = [25] * len(truss.members)
    ub = [300] * len(truss.members)

    optimizer = FrameOptimizer(frame,
                             opt_algorithm='basinhopping',
                             lb=lb,
                             ub=ub,
                             maxiter=100)
                
    frame.plot()               
    #optimizer.optimize()
    frame.plot()


def truss_column_test():

    coord = [0,0]
    L1 = 0.3
    L2 = 1
    H1 = 5.5
    H2 = 5
    n = 2
    
    
    col = TrussColumn(coord, H1, H2, L1, L2, n)
    col.plot()
    

    col.generate()
    col.f.draw()
    print(col.top_chords[0].nodal_coordinates)
    col.to_robot("truss_column")
    
def strength_test():
    
    col = SteelColumn([[0,0], [0,3.65]], material='S420',profile="RHS 200x200x8")
    print(col.cross_section.MRd)
    print(col.steel_member.NbRd)
    
 
def robot_test_3D():
    
    H0 = 4
    H1 = 1
    H2 = 1.5
    H3 = H1
    L1 = 10
    L2 = L1
    n = 20
    n2 = 6
    s = 6
    num_frames = 4
    col1 = SteelColumn([[0,0], [0, H0 + H1]], profile="He 300 a")
    col2 = SteelColumn([[L1, 0], [L1, H0 + H2]], profile="He 300 a")
    col3 = SteelColumn([[L1 + L2, 0], [L1 + L2, H0 + H3]], profile="He 300 a")
    beam1 = SteelBeam([[0, H0+H1], [L1, H0+H2]])
    beam2 = SteelBeam([[L1, H0+H2], [L1 + L2, H0+H3]])
    
    #col4 = SteelColumn()
    col5 = SteelColumn([[0,0], [0, H0 + H1]], profile="He 300 a")
    col6 = SteelColumn([[L1 + L2, 0], [L1 + L2, H0 + H3]], profile="He 300 a")
    
    end_frame1 = Frame2D()
    end_frame1.add(col1)
    end_frame1.add(col2)
    end_frame1.add(FixedSupport([L1, 0]))
    end_frame1.add(col3)
    end_frame1.add(beam1)
    end_frame1.add(beam2)
    end_frame1.add(SteelColumn([[L1/2, 0], [L1/2, beam1.shape(L1/2)]]))
    end_frame1.add(SteelColumn([[L1 + L2/2, 0], [L1 + L2/2, beam2.shape(L1 + L2/2)]]))
    end_frame1.add(FixedSupport([L1/2, 0]))
    end_frame1.add(FixedSupport([L1 + L2/2, 0]))
    
    
    end_frame2 = Frame2D()
    long_frame1 = Frame2D(simple=[1, num_frames+1, H0 + H1, s], supports='fixed')
    long_frame1.add(SteelBeam([[0, 0], [s, H0+H1]]))
    long_frame1.add(SteelBeam([[s*(num_frames+1), 0], [s*num_frames, H0+H1]]))
    long_frame2 = Frame2D(simple=[1, num_frames+1, H0 + H3, s], supports='fixed')
    long_frame2.add(SteelBeam([[50, (H0+H1)/2], [60, (H0+H1)/2]]))
    long_frame2.add(SteelBeam([[60, (H0+H1)/2], [70, (H0+H1)/2]]))
    main_frame = Frame2D()
    main_frame.add(col5)
    main_frame.add(col6)
    truss = Truss2D(simple=[H0, H1, H2, H3, L1, L2, n], fem_model=main_frame.f)
    main_frame.add(truss)
    stiff_truss = Truss2D(simple=[0, s, s, s, L1, L2, n2])
    stiff_truss.plot()
    
    end_frame1.plot()
    long_frame1.plot()
    long_frame2.plot()
    main_frame.plot()
    
    robot3D.to_robot_3D("testi",
                        end_frame1,
                        long_frame1,
                        main_frame,
                        num_frames=num_frames-1,
                        s=s,
                        L=L1 + L2,
                        stiff_truss=stiff_truss)
    
if __name__ == '__main__':
    #truss_test()
    #frame_truss()
    portal_frame()
    #frame_test()
    #simple_truss()
    #optimization_test()
    #truss_column_test()
    #strength_test()
    #robot_test_3D()

"""
truss = Truss2D(num_elements=1)

coord1 = [[0,4], [2,1]]
coord2 = [[0,0], [2,-0]]

top_chord = TopChord(coord1)
bottom_chord = BottomChord(coord2)
truss.add(top_chord)
truss.add(bottom_chord)
"""


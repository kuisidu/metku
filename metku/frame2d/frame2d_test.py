# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
from metku.frame2d.frame2d import Frame2D, SteelBeam, SteelColumn, PointLoad,\
     FixedSupport, XHingedSupport, LineLoad, YHingedSupport, XYHingedSupport,\
     Hinge

import random
from decimal import Decimal
import numpy as np
import time
import matplotlib.pyplot as plt
import copy

def test1():
    """
    Test for creating frame
    """
    # Coordinates of members, loads and supports
    coord1 = [[0,1], [1,1]]
    coord2 = [[1,0], [1,1]]
    coord3 = [[0.0,0], [0.0, 2]]
    supp_coord1 = [0.0,0]
    supp_coord2 = [1,0]
    #hinge_coord = [0.5, 1]

    # Loads and members
    col1 = SteelColumn(coord2)
    col2 = SteelColumn(coord3)
    beam1 = SteelBeam(coord1, num_elements=10)
    
    load1 = PointLoad(coord3[1], [50, 0,0])
    load2 = LineLoad(beam1, [-50, -50], 'y')
    # Create empty frame 'envelope'
    frame = Frame2D(num_elements=5)

    # Add members
    frame.add(col1)
    frame.add(col2)
    frame.add(beam1)

    # Add loads
    frame.add(load1)
    frame.add(load2)

    # Add supports
    frame.add(FixedSupport(supp_coord1, supp_id=1))
    frame.add(FixedSupport(supp_coord2, supp_id=2))
    
    # Add hinge
    #frame.add(Hinge(hinge_coord))

    # Generate frame
    frame.generate()
    frame.plot()    
    # Calculate result
    frame.calculate()
    frame.f.draw()
    frame.bmd(10)
    #frame.plot_deflection()
"""
def test2():
    frame = Frame2D(num_elements=2)
    frame.load_robot_csv()
    for elem in frame.f.elements:
        if elem.length() == 0:
            print(frame.f.elements.index(elem))
    frame.f.draw()
    frame.plot()
    frame.calculate()
    frame.bmd(5)
"""

def test3():
    # SIMPLE = [storeys, bays, storey height, bay length]
    frame = Frame2D(simple=[2,1,2,2], supports='fixed', num_elements=4)
    for member in frame.members.values():
        if member.mtype == 'beam':
            frame.add(LineLoad(member, [-50, -50], 'y'))
       # if member.mtype =='column' and member.mem_id  == 1:
        #    frame.add(PointLoad(member.coordinates[1], [100,0,0]))
    frame.generate()
    #frame.plot()
    frame.hinge_joints()
    frame.calculate()
    for member in frame.members.values():
        member.plot_results()
    frame.plot()
    frame.bmd(20)

def test4():
   
    # Coordinates of members, loads and supports
    coord1 = [[0, 1], [1, 1]]
    coord2 = [[1, 0], [1, 1]]
    coord3 = [[0.0, 0], [0.0, 1.5]]
    supp_coord1 = [0.0, 0]
    supp_coord2 = [1, 0]
    # Loads and members
    beam1 = SteelBeam(coord1)
    col1 = SteelColumn(coord3)
    col2 = SteelColumn(coord2)
    load1 = PointLoad(coord3[1], [10, 0, 0])
    load2 = LineLoad(beam1, [-10, -10], 'y')
    # Create empty frame 'envelope'
    frame = Frame2D(num_elements=3)
    # Add members
    frame.add(col1)
    frame.add(col2)
    frame.add(beam1)
    # Add loads
    frame.add(load1)
    frame.add(load2)
    # Add supports
    frame.add(FixedSupport(supp_coord1))
    frame.add(FixedSupport(supp_coord2))
    # Generate frame
    frame.generate()
    frame.calculate()
    #for mem in frame.members.values():
    #    print(mem.nodal_coordinates)
    frame.f.draw()
    frame.bmd()
    
def test5():
    """
    Function for testing how the number of elements affects calculation times
    and calculation accuracy
    """
    bays = 1
    storeys = 1
    num_elem = 0
    beams = bays * storeys
    cols = (bays +1) * storeys
    num_mem = beams + cols
    prev =[0] *num_mem
    cum = 0
    Y_acc = []
    Y_time = []
    X = []
    step = 1
    for i in range(10):
        time.sleep(0.1)
        num_elem += step
        # SIMPLE = [storeys, bays, storey height, bay length]
        frame = Frame2D(simple=[storeys,bays,2,2],
                        supports='fixed',
                        num_elements=num_elem)
        for member in frame.members.values():
            if member.mtype == 'beam':
                frame.add(LineLoad(member, [-50, -50], 'y'))
            elif member.mem_id < storeys:
                frame.add(PointLoad(member.coordinates[1], [50, 0, 0]))
        
        start = time.time()
        frame.generate()
        frame.calculate()
        frame.bmd(10)
        end = time.time()
        for member in frame.members.values():
            change = abs((abs(member.med) - abs(prev[member.mem_id])) / member.med)
            prev[member.mem_id] = member.med
            #print(member.med)
        
        if num_elem > step:
            cum += change
            X.append(num_elem)
            Y_acc.append(change*100)
            Y_time.append((end-start)*1000)
        print(f'elems: {num_elem},change: {change*100:.4f} %, cumulative: {cum*100:.4f} %, time: {(end - start)*1000:.4f} ms')
    print()
    print("MEMBERS: ", num_mem)
    plt.plot(X, Y_acc)
    plt.xlabel('num_elements / member')
    plt.ylabel('difference in values [%]')
    plt.show()
    
    plt.plot(X, Y_time)
    plt.xlabel('num_elements / member')
    plt.ylabel('calculation time [ms]')
    plt.show()
    
    
def test6():
    frame = Frame2D(num_elements=6)
    coord1 = [[0,0],[0,1]]
    coord2 = [[0,1], [1,1]]
    coord3 = [[1,0],[1,1]]

    supp_coord1 = [0,0]
    supp_coord2 = [1,0]

    col1 = SteelColumn(coord1, num_elements=2)
    col2 = SteelColumn(coord3, num_elements=3)
    beam = SteelBeam(coord2)

    load2 = LineLoad(beam, [-50, -50], 'y')

    frame.add(col1)
    frame.add(col2)
    frame.add(beam)
    frame.add(load2)
    frame.add(FixedSupport(supp_coord1))
    frame.add(FixedSupport(supp_coord2))
    frame.generate()
    frame.calculate()
    #frame.plot()
    """
    frame = Frame2D(simple=[1,1,1,1], supports='fixed', num_elements=3)
    for member in frame.members.values():
            if member.mtype == 'beam':
                frame.add(LineLoad(member, [-50, -50], 'y'))
            elif member.mem_id < 1:
                frame.add(PointLoad(member.coordinates[1], [50, 0, 0]))
    frame.generate()
    frame.calculate()
    """

    frame.plot_deflection(25)
    frame.bmd(25)
    frame.plot_buckling(2)
    frame.to_robot("result_test")

def test7():
    frame = Frame2D(simple=[1,1,1,1], supports='fixed', num_elements=4)
    frame.add(LineLoad(frame.members[2], [-50, -50], 'y'))
    print(frame.members[2].loads)
    print(frame.line_loads)
    frame.generate()
    frame.plot()
    frame.delete_member(2)
    frame.plot()
    frame.add(SteelBeam([[0,1],[1,2]]))
    frame.add(SteelColumn([[1,1], [1,2]]))
    frame.add(PointLoad([1,2], [-50, 0, 0]))

def test8():
    frame = Frame2D(simple=[3, 3, 10, 14], supports='xyhinged', num_elements=4)
    for member in frame.members.values():
        if member.mtype == 'beam':
            member.profile = 'he 550 a'
            frame.add(LineLoad(member, [-33, -33], 'y'))
        elif member.mtype == 'column':
            member.profile = 'he 300 b'
            if member.mem_id == 0:
                frame.add(LineLoad(member, [4, 4], 'x'))
            elif member.mem_id == 1:
                frame.add(LineLoad(member, [2, 2], 'x'))

    #frame.add(LineLoad(frame.members[2], [-10, -50], 'y'))
    #frame.hinge_joints()
    frame.generate()
    frame.calculate()
    #for elem in frame.members[0].elements.values():
    #    print(elem.bending_moment[0])
    #print(elem.bending_moment[1])
    frame.f.draw()
    frame.plot_deflection(25)
    frame.bmd(25)
    frame.plot_buckling(20)
    frame.to_robot("result_test")
    

def test9():
    """
    Test for buckling shape plotting
    """
    frame = Frame2D(simple=[1,10,4,4], supports='fixed', num_elements=2)
    for member in frame.members.values():
        if member.mtype == "beam":
            member.profile = 'ipe 80'
        else:
            member.profile = 'HE 300 A'
    for beam in frame.beams:
        frame.add(LineLoad(beam, [-10, -10], 'y'))
    #frame.add(PointLoad([4,2], [70, 0, 0]))
    #frame.add(PointLoad([2,4], [0, -80, 0]))
    #frame.members[2].Sj2 = 1e-6
    frame.generate()
    frame.calculate()

    frame.bmd(10)
    
def test10():
    frame = Frame2D(num_elements=2)
    
    coord1 = [[0,1], [5,1]]
    coord2 = [[0,0], [5,0]]
    
    top_chord = SteelBeam(coord1)
    bottom_chord = SteelBeam(coord2)
    
    members = 7
    

    dx = coord1[1][0]/(members-1)
    y = lambda x: 1
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
            

        member = SteelBeam([c1, c2])
        frame.add(member)
        
    
    
    frame.add(top_chord)
    frame.add(bottom_chord)
    frame.add(LineLoad(top_chord, [-50, -50], 'y'))
    frame.add(FixedSupport(coord2[0]))
    frame.add(FixedSupport(coord2[1]))
    frame.generate()
    #frame.hinge_joints()
    frame.plot()
    frame.calculate()
    frame.bmd(5)
    frame.f.draw()

    frame.to_robot("frame-truss-test")

def test_all():
    test1()
    #test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()
    test10()

if __name__ == '__main__':

    test9()
    # from time import time
    #
    # frame = Frame2D(simple=[1 ,1, 5000, 5000], supports='fixed', num_elements=10)
    # frame.add(LineLoad(frame.members[2], [-10, -10], 'y'))
    # for i in range(10):
    #     start = time()
    #     frame.generate()
    #     frame.calculate()
    #     end = time()
    #     print("Time elapsed: ", end - start)

    """
    for member in frame.members.values():
            if member.mtype == 'beam':
                frame.add(LineLoad(member, [-50, -50], 'y'))
            elif member.mem_id < 1:
                frame.add(PointLoad(member.coordinates[1], [50, 0, 0]))
    frame.generate()
    frame.calculate()
    """





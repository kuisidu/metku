import sys
sys.path.append(".\src")

from frame2d import Frame2D, SteelBeam, SteelColumn, PointLoad,\
     FixedSupport, XHingedSupport, LineLoad, YHingedSupport, Hinge

import random
import time
import matplotlib.pyplot as plt

def test1():

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
    beam1 = SteelBeam(coord1)
    
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

    #col1.profile = 'HE 200 A'
    frame.plot()
    
    # Calculate result
    frame.calculate()
    #print(frame.nodal_forces[2])

    #frame.f.draw()
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

    col1 = SteelColumn(coord1)
    col2 = SteelColumn(coord3)
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

    print(frame.members[2].nodal_forces)
    frame.bmd(25)

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
    frame = Frame2D(simple=[1, 2, 10, 14], supports='xyhinged', num_elements=10)
    for member in frame.members.values():
        if member.mtype == 'beam':
            member.profile = 'he 550 a'
            frame.add(LineLoad(member, [-33, -33], 'y'))
        elif member.mtype == 'column':
            member.profile = 'he 300 b'
            if member.mem_id == 0:
                frame.add(LineLoad(member, [4, 4], 'x'))
            elif member.mem_id == 2:
                frame.add(LineLoad(member, [2, 2], 'x'))

    #frame.add(LineLoad(frame.members[2], [-10, -50], 'y'))
    frame.generate()
    frame.calculate()
    frame.plot_buckling(50, 2)

def test9():
    frame = Frame2D(simple=[1,1,4,4], supports='fixed', num_elements=2)
    for member in frame.members.values():
        member.profile = 'ipe 80'
    frame.add(PointLoad([4,2], [70, 0, 0]))
    frame.add(PointLoad([2,5.5], [0, -80, 0]))
    #frame.members[2].Sj2 = 1e-6
    frame.generate()
    frame.members[0].coordinates = [[0, 0], [0,7]]
    frame.calculate()
    #frame.f.draw()

    frame.plot_buckling(20,2)
if __name__ == '__main__':
    test8()
    """
    frame = Frame2D(simple=[1,1,2,2], supports='fixed', num_elements=10)
    for member in frame.members.values():
            if member.mtype == 'beam':
                frame.add(LineLoad(member, [-50, -50], 'y'))
            elif member.mem_id < 1:
                frame.add(PointLoad(member.coordinates[1], [50, 0, 0]))
    frame.generate()
    frame.calculate()
    """




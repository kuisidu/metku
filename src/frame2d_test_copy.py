from frame2d_copy import Frame2D, SteelBeam, SteelColumn, PointLoad,\
     FixedSupport, XYHingedSupport, LineLoad, YHingedSupport, Hinge

from anna_optimization_frame2d import optimization

def test1():

    # Coordinates of members, loads and supports
    coord1 = [[0,1], [1,1]]
    coord2 = [[1,0], [1,1]]
    coord3 = [[0.0,0], [0.0, 2]]
    supp_coord1 = [0.0,0]
    supp_coord2 = [1,0]
    hinge_coord = [0.5, 1]

    # Loads and members
    load1 = PointLoad(coord3[1], [50, 0,0])
    col1 = SteelColumn(coord1)
    col2 = SteelColumn(coord3)
    beam1 = SteelBeam(coord2)

    load2 = LineLoad(col1, [-50, -50], 'y')
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
    frame.add(Hinge(hinge_coord))

    # Generate frame
    frame.generate()

    #col1.profile = 'HE 200 A'
    frame.plot()

    # Calculate result
    frame.calculate()
    #print(frame.nodal_forces[2])
    for member in frame.members.values():
        for element in member.elements.values():
            print(element.bending_moment)
    #frame.f.draw()
    frame.bmd(20)
    frame.design_members()
    frame.plot()
    
    frame.members[0].bmd_test()
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
        for element in member.elements.values():
            print(element.bending_moment)
   
    #print()
    frame.f.draw()
    frame.bmd(50)
    
    frame.plot_deflection()
    print ("weight =", frame.weight,"kg")
    print ("area =", frame.weight,"kg")
    
    #frame.plot()
     
    
    
    optimization.optimize_frame(frame,optimizer='slsqp',maxiter=50, swarmsize=150)

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
    
    frame = Frame2D(simple=[5,5,5,5], supports='fixed')
    for member in frame.members.values():
        if member.mtype == 'beam':
            frame.add(LineLoad(member, [-50, -50], 'y'))
    
    frame.generate()
    frame.calculate()
    frame.plot(loads=False)
    frame.design_members()
    frame.plot(loads=False)
    
def test6():
    frame = Frame2D(num_elements=2)
    beam = SteelBeam([[0,0], [10,0]])
    frame.add(beam)
    frame.add(YHingedSupport([0,0]))
    frame.add(XYHingedSupport([10,0]))
    frame.add(LineLoad(beam, [-10,-10], 'y'))
    frame.generate()
    frame.f.draw()
    frame.calculate()
    frame.bmd(20)
    frame.members[0].bmd_test()
    
if __name__ == '__main__':
    
    test6()




def split(elem, parts):
    n1, n2 = elem.nodes
    x1, y1 = n1.x
    x2, y2 = n2.x

    dx = (x2 - x1) / parts
    dy = (y2 - y1) / parts
    for i in range(parts-1):
        x = x1 + (i+1) * dx
        y = y1 + (i+1) * dy
        fem.add_node(x, y)
        n2 = fem.nodes[-1]
        sect = fem.sections[0]
        mat = fem.materials[0]
        split_elem = EBBeam(n1, n2, sect, mat)
        fem.add_element(split_elem)
        n1 = n2
        if ll1.elem == elem:
            ll = LineLoad(2, split_elem, [0, 1], [-2, -2], 1, 1)
            fem.add_load(ll)
        


# Import dependencies
from fem.frame_fem import FrameFEM, FEMNode, Material, BeamSection,\
     LoadCase, LineLoad, PointLoad, Support
from fem.elements import EBBeam

# Initialize fem
fem = FrameFEM()

# Nodal coordinates
c1 = [0, 0]
c2 = [0, 2]
c3 = [0, 3]
c4 = [5, 3]
c5 = [8, 3]
c6 = [8, 1]

# Add nodes
fem.add_node(*c1)
fem.add_node(*c2)
fem.add_node(*c3)
fem.add_node(*c4)
fem.add_node(*c5)
fem.add_node(*c6)

# Add section
A = 1e3
Iy = 1e6
sect = BeamSection(A, Iy)
fem.add_section(sect)

# Add material
E = 210e3
nu = 0.3
rho = 7850e-9
fem.add_material(E, nu, rho)

# Create and add elements
for i in range(fem.nnodes()-1):
    n1 = fem.nodes[i]
    n2 = fem.nodes[i + 1]
    # using same section and material for every element
    sect = fem.sections[0]
    mat = fem.materials[0]
    elem = EBBeam(n1, n2, sect, mat)
    fem.add_element(elem)

# Draw the frame to see if it's correct
#fem.draw()

# Add supports
fem.add_support(1, 0, [0, 1], 0)
fem.add_support(1, -1, [1], 0)

# Add loads
pl1 = PointLoad(2, fem.nodes[1], [8, 0, 0], 1)
pl2 = PointLoad(2, fem.nodes[4], [0, -5, 0], 1)
ll1 = LineLoad(2, fem.elements[2], [0, 1], [-2, -2], 1, 1)

fem.add_load(pl1)
fem.add_load(pl2)
fem.add_load(ll1)

# Add load case
#fem.add_loadcase(1, 2)

### Calculate results
##fem.nodal_dofs()
##fem.linear_statics()

# Print results
print(f'Reaction at {c1}: Fx: {-fem.elements[0].shear_force[0]:.3f} kN,\
Fy: {-fem.elements[0].axial_force[0]:.3f} kN') 
print(f'Reaction at {c6}: Fx: {-fem.elements[-1].shear_force[0]:.3f} kN,\
Fy: {-fem.elements[-1].axial_force[0]:.3f} kN')


# Shear force diagram
import numpy as np
import matplotlib.pyplot as plt
import math

##for elem in fem.elements:
##    scale = 10
##    s1, s2 = np.asarray(elem.shear_force) / 10
##    n1, n2 = elem.nodes
##    x1, y1 = n1.x
##    x2, y2 = n2.x
##    if y2 - y1 == 0:
##        plt.plot([x1, x1, x2, x2], [y1, y1 - s1, y2 - s2, y2], c='k')
##    else:
##        plt.plot([x1, x1 - s1, x2 - s2, x2], [y1, y1, y2, y2], c='k')
##
##fem.draw()

# Bending moment diagram
for i in range(fem.nels()):
    elem = fem.elements[i]
    split(elem, 10)
    
fem.add_loadcase(1, 2)
fem.nodal_dofs()
fem.linear_statics()
for elem in fem.elements:
    scale = 10
    s1, s2 = np.asarray(elem.shear_force) / scale
    m1, m2 = np.asarray(elem.bending_moment) / scale
    n1, n2 = elem.nodes
    x1, y1 = n1.x
    x2, y2 = n2.x
    if y2 - y1 == 0:
        plt.plot([x1, x2], [y1 - m1, y2 - m2], c='k')
    else:
        plt.plot([x1 - m1,  x2 - m2], [y1, y2], c='k')
    
fem.draw()


        


    

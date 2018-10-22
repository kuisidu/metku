# -*- coding: utf-8 -*-
"""
Created on Tue May  9 08:51:34 2017

Plane frame from Arbabi, Example 12.5

@author: kmela
"""

import frame_fem as fem

f = fem.FrameFEM()

""" Units: kN,mm """

E = 29000 # ksi
A = 10 # in2
I = 250 # in4
L = 15.0*12 # in
L2 = 27.0*12 # in
h = 10.0*12 # in

f.add_node(0.0,0.0)
f.add_node(0.0,h)
f.add_node(L,h)
f.add_node(L2,0)
    
# Add material    
f.add_material(E,0.3,7850e-9)
    
# Add section properties
s = fem.BeamSection(A,I)
f.add_section(s)
    
# Add first element
el = fem.EBBeam(f.nodes[0],f.nodes[1],f.sections[0],f.materials[0])
f.add_element(el)
    
# Add second element    
el2 = fem.EBBeam(f.nodes[1],f.nodes[2],f.sections[0],f.materials[0])    
f.add_element(el2)

# Add third element    
el3 = fem.EBBeam(f.nodes[2],f.nodes[3],f.sections[0],f.materials[0])    
f.add_element(el3)
    
# Add supports
f.add_support(1,0,[0,1,2],0)
f.add_support(1,3,[0,1,2],0)
    
# Assign nodal dofs
f.nodal_dofs()
    
"""
for n in f.nodes:
    print(n.dofs)
"""

# Add point load
load1 = fem.PointLoad(2,f.nodes[1],[15.0,0,0],1.0)
f.add_load(load1)

# Add point load
load2 = fem.PointLoad(2,f.nodes[2],[0.0,41.0,0],1.0)
f.add_load(load2)

# Add line load
q = -2.2/12.0
load2 = fem.LineLoad(2,f.elements[1],[0.0,1.0],[q,q],"y")
f.add_load(load2)
    
# Add load case    
f.add_loadcase(supp_id=1,load_id=2)

f.draw()

K = f.global_stiffness_matrix()
p = f.global_load_vector(2)

print(K)
print(p)

u = f.linear_statics()
print(u)

for n in f.nodes:
    print(n.u)

print(" ** Calculate element forces ** ")

for elem in f.elements:
    print("** Element **")
    elem.internal_forces()
    
    print("** Axial forces **")
    print(elem.axial_force)
    print("** Bending moments **")
    print(elem.bending_moment)
    print("** Shear forces **")
    print(elem.shear_force)
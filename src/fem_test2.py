# -*- coding: utf-8 -*-
"""
Created on Tue May  9 08:51:34 2017

Plane frame from MCGuire et al, Example 5.7

@author: kmela
"""

import frame_fem as fem

f = fem.FrameFEM()

""" Units: kN,mm """

E = 200.0 # GPa
A = 6000.0 # mm2
I = 200.0e6 # mm4
L = 7416.0 # mm
L2 = 8000.0 # mm
h = 3000.0 # mm

f.add_node(0.0,0.0)
f.add_node(L,h)
f.add_node(L+L2,h)
    
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
    
# Add supports
f.add_support(1,0,[0,1,2],0)
f.add_support(1,2,[0,1,2],0)
    
# Assign nodal dofs
f.nodal_dofs()
    
"""
for n in f.nodes:
    print(n.dofs)
"""

# Add point load
load1 = fem.PointLoad(2,f.nodes[1],[18.75,-46.35,0],1.0)
f.add_load(load1)

# Add line load
q = -4e-3
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

for elem in f.elements:
    elem.internal_forces()
    print("Element")
    print(elem.axial_force)
    print(elem.bending_moment)
    print(elem.shear_force)
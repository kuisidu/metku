# -*- coding: utf-8 -*-
"""
Created on Tue May  9 08:51:34 2017

Simple portal frame

@author: kmela
"""

import frame_fem as fem
from sections.steel.ISection import IPE, HEA
from eurocodes.en1993 import constants

f = fem.FrameFEM()

""" Units: N,mm """

pcol = HEA(100,catalogue=True)
pbeam = IPE(100,catalogue=True)

# Load: N/mm
q = -50 # kN/m

E = 210000 # MPa
#A = 10 # in2
#I = 250 # in4

h = 1000 # mm
L = 1000 # mm

elements_per_member = 3

# Add nodes for column 1
dy = h/elements_per_member
y = 0.0
x = 0.0
dx = L/elements_per_member

# Add material    
f.add_material(E,0.3,7850e-9)
    
# Add section properties
sbeam = fem.BeamSection(pbeam.A,pbeam.I[0])
f.add_section(sbeam)
scol = fem.BeamSection(pcol.A,pcol.I[0])
f.add_section(scol)

nodes = 1

f.add_node(x,y)
for i in range(elements_per_member):
    y += dy
    f.add_node(x,y)
    el = fem.EBBeam(f.nodes[nodes-1],f.nodes[nodes],f.sections[0],f.materials[0])
    f.add_element(el)
    nodes += 1
    #y += dy
    

# add beam
for i in range(elements_per_member):
    x += dx
    f.add_node(x,y)    
    el = fem.EBBeam(f.nodes[nodes-1],f.nodes[nodes],f.sections[1],f.materials[0])
    f.add_element(el)
    # Add line load    
    load = fem.LineLoad(2,el,[0.0,1.0],[q,q],"y")
    f.add_load(load)
    nodes += 1
    
    
# add the second column     
for i in range(elements_per_member):
    y -= dy
    f.add_node(x,y)
    el = fem.EBBeam(f.nodes[nodes-1],f.nodes[nodes],f.sections[0],f.materials[0])
    f.add_element(el)
    nodes += 1


f.draw()

print(nodes)
    
# Add supports
f.add_support(1,0,[0,1,2],0)
f.add_support(1,nodes-1,[0,1,2],0)
    
# Assign nodal dofs
f.nodal_dofs()
    
# Add load case    
f.add_loadcase(supp_id=1,load_id=2)

u = f.linear_statics()
#print(u)

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

# -*- coding: utf-8 -*-
"""
Created on Tue May  9 08:51:34 2017

Simple portal frame

@author: kmela
"""

import frame_fem
from fem.frame.elements.ebbeam import EBBeam
from sections.steel.ISection import IPE, HEA
from eurocodes.en1993 import constants

f = frame_fem.FrameFEM()

""" Units: N,mm """

pcol = HEA(100,catalogue=True)
pbeam = IPE(100,catalogue=True)

# Load: N/mm
q = -50 # kN/m

E = 210000 # MPa
#A = 10 # in2
#I = 250 # in4

h = 1000 # mm
L = 2000 # mm

elements_per_member = 4

# Add nodes for column 1
dy = h/elements_per_member
y = 0.0
x = 0.0
dx = L/elements_per_member

# Add material    
f.add_material(E,0.3,7850e-9)
    
# Add section properties
sbeam = frame_fem.BeamSection(pbeam.A,pbeam.I[0])
f.add_section(sbeam)
scol = frame_fem.BeamSection(pcol.A,pcol.I[0])
f.add_section(scol)

nodes = 1

# Add column 1
f.add_node(x,y)
for i in range(elements_per_member):
    y += dy
    f.add_node(x,y)
    el = EBBeam(f.nodes[nodes-1],f.nodes[nodes],f.sections[1],f.materials[0])
    f.add_element(el)
    nodes += 1
    #y += dy
    

# add beam
for i in range(elements_per_member):
    x += dx
    f.add_node(x,y)    
    el = EBBeam(f.nodes[nodes-1],f.nodes[nodes],f.sections[0],f.materials[0])
    f.add_element(el)
    # Add line load    
    load = frame_fem.LineLoad(2,el,[0.0,1.0],[q,q],"y")
    f.add_load(load)
    nodes += 1
    
    
# add the second column     
for i in range(elements_per_member):
    y -= dy
    f.add_node(x,y)
    el = EBBeam(f.nodes[nodes-1],f.nodes[nodes],f.sections[1],f.materials[0])
    f.add_element(el)
    nodes += 1


f.draw()

#print(nodes)
    
# Add supports
f.add_support(1,0,[0,1,2],0)
f.add_support(1,nodes-1,[0,1,2],0)
    
# Assign nodal dofs
f.nodal_dofs()
    
# Add load case    
f.add_loadcase(supp_id=1,load_id=2)

u = f.linear_statics()
#print(u)

print("Results of analysis:")
print("Nodal displacements:")
    
f.print_displacements()

print("Member forces:")
mem_cnt = 0
for el in f.elements:
    print('*** Element {:d} ***'.format(mem_cnt))
    mem_cnt += 1
    #print(el.floc)        
    print('Axial force (node 1) = {0:5.3f} kN'.format(el.axial_force[0]*1e-3))
    print('Axial force (node 2) = {0:5.3f} kN'.format(el.axial_force[1]*1e-3))
    print('Bending moment (node 1) = {0:5.3f} kNm'.format(el.bending_moment[0]*1e-6))
    print('Bending moment (node 2) = {0:5.3f} kNm'.format(el.bending_moment[1]*1e-6))
    print('Shear force (node 1) = {0:5.3f} kN'.format(el.shear_force[0]*1e-3))
    print('Shear force (node 2) = {0:5.3f} kN'.format(el.shear_force[1]*1e-3))

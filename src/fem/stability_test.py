# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 12:49:34 2018

Testing linear stability analysis for plane frames

@author: kmela
"""

import numpy as np
import frame_fem
from sections.steel import RHS
from fem.frame.elements.ebbeam import EBBeam
from timeit import default_timer as timer
    
p = RHS.SHS(100,4)
    
f = frame_fem.FrameFEM()

# Add material
f.add_material(210e3,0.3,7850e-9)
    
# Add section properties    
s = frame_fem.BeamSection(p.A,p.I[0])
f.add_section(s)
    
# Create nodes   
L = 5000.0
h = 5000.0

nels = 20
nnodes = nels+1
dx = L/nels

# Create nodes
for i in range(nels+1):
    f.add_node(0,i*dx)    
    # Create column elements
    if i > 0:
        el = EBBeam(f.nodes[i-1],f.nodes[i],s,f.materials[0])  
        f.add_element(el)
        
load1 = frame_fem.PointLoad(2,f.nodes[nnodes-1],[0.0,-25000.0,0.0],1.0)
f.add_load(load1)

f.draw()

# Add supports
f.add_support(sid=1,nid=0,dof=[0,1,2],val=0.0)

# Assign nodal dofs
f.nodal_dofs()

# Add load case
f.add_loadcase(supp_id=1,load_id=2)

start = timer()
u = f.linear_statics()
w,v = f.linear_buckling()
end = timer()
print("Elapsed time on buckling analysis: " + str(end - start) + " [s]\n")
print(w)
#print(v)
f.print_displacements()

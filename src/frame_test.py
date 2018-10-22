# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:31:04 2018

Testing semi-rigid frame with two storeys

@author: kmela
"""

import numpy as np
import frame_fem
from hollow_sections import SHS
#from fem.frame.elements.ebbeam import EBBeam
#from fem.frame.elements.eb_semi_rigid_beam import EBSemiRigidBeam
    
p = SHS(100,4)
    
f = frame_fem.FrameFEM()

# Add material
f.add_material(210e3,0.3,7850e-9)
    
# Add section properties    
s = frame_fem.BeamSection(p.A,p.I[0])
f.add_section(s)
    
# Create nodes   
L = 5000.0
h = 5000.0

storeys = 2
nbeam_els = 4
dx = L/nbeam_els

xends = [0,L]
    
# Create nodes
for j in range(2):
    for i in range(storeys+1):
        f.add_node(xends[j],h*i)
        # Create column elements
        if i > 0:
            el = frame_fem.EBBeam(f.nodes[j*(storeys+1)+i-1],f.nodes[j*(storeys+1)+i],s,f.materials[0])  
            f.add_element(el)

# Create beam elements
sr_frame = True
line_load = True
X = 0.0
N = f.nnodes()
DN = N-1
Sj1 = 1e9
Sj2 = np.inf
for i in range(nbeam_els-1):
    for j in range(storeys):
        f.add_node((i+1)*dx,(j+1)*h)
        if i == 0:
            if sr_frame:
                el = frame_fem.EBSemiRigidBeam(f.nodes[N-DN],f.nodes[N],s,f.materials[0],rot_stiff=[Sj1,Sj2])            
            else:
                el = frame_fem.EBSemiRigidBeam(f.nodes[N-DN],f.nodes[N],s,f.materials[0])            
        elif i == nbeam_els-2:
            el = frame_fem.EBBeam(f.nodes[N-storeys],f.nodes[N],s,f.materials[0])
            f.add_element(el)
            if line_load:
                load = frame_fem.LineLoad(2,el,[0,1],-5*np.array([1,1]),1)
                f.add_load(load)
            if sr_frame:
                el = frame_fem.EBSemiRigidBeam(f.nodes[N],f.nodes[storeys+2+j],s,f.materials[0],rot_stiff=[Sj2,Sj1])            
            else:
                el = frame_fem.EBBeam(f.nodes[N],f.nodes[storeys+2+j],s,f.materials[0])
        else:
            el = frame_fem.EBBeam(f.nodes[N-storeys],f.nodes[N],s,f.materials[0])
            
        f.add_element(el)
        if line_load:
            load = frame_fem.LineLoad(2,el,[0,1],-5*np.array([1,1]),1)
            f.add_load(load)
        elif i == int(nbeam_els/2):            
            load1 = frame_fem.PointLoad(2,f.nodes[N],[0.0,-25000.0,0.0],1.0)
            f.add_load(load1)

        N = N+1
        print(N)
        
f.draw()

# Add supports
f.add_support(1,0,[0,1,2],0.0)
f.add_support(1,storeys+1,[0,1,2],0.0)

# Assign nodal dofs
f.nodal_dofs()

# Add load case
f.add_loadcase(supp_id=1,load_id=2)
    
p = f.global_load_vector(2)
"""
    
    
    # Add load    
    
    load1 = LineLoad(2,f.elements[1],[0,1],-5*np.array([1,1]),1)
    load2 = LineLoad(2,f.elements[2],[0,1],-5*np.array([1,1]),1)
    f.add_load(load1)
    f.add_load(load2)
    
    nload = f.nodes[-(2+int(nbeam_els/2))]
    print(nload.x)
    
    #f.add_load(load2)
    
    f.draw()
    # Add point load
    
    load1 = PointLoad(2,f.nodes[0],[0.0,0.0,100.0],1.0)
    f.add_load(load1)
    load2 = PointLoad(2,f.nodes[1],[0.0,0.0,-100.0],1.0)
    f.add_load(load2)
    
    
    
    #print(p)
     Two-storey frame
"""
    
K = f.global_stiffness_matrix()
u = f.linear_statics()
    
print("Results of analysis:")
print("Nodal displacements:")
    
#f.print_displacements()

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
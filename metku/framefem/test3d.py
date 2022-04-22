# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Sat Aug 17 14:01:25 2019

@author: kmela
"""

from framefem import FrameFEM, Section, PointLoad
from framefem.elements.rod import Rod

def truss3dTest():
    truss3d = FrameFEM()
    # Units: kN,mm
    
    truss3d.add_material(200,0.3,7850e-9)
    
    truss3d.add_node(0.0,0.0,-4e3)
    truss3d.add_node(-3e3,0.0,0.0)
    truss3d.add_node(0.0,0.0,4e3)
    truss3d.add_node(0.0,5e3,0.0)
    
    A1 = Section(A=0.001e6)
    A2 = Section(A=0.002e6)
    
    truss3d.add_section(A1)
    truss3d.add_section(A2)
    
    truss3d.add_load(PointLoad(2,truss3d.nodes[3],[12,0,0,0,0,0],1))
    
    truss3d.add_element(Rod(truss3d.nodes[0],truss3d.nodes[3],A1,truss3d.materials[0]))
    truss3d.add_element(Rod(truss3d.nodes[1],truss3d.nodes[3],A2,truss3d.materials[0]))
    truss3d.add_element(Rod(truss3d.nodes[2],truss3d.nodes[3],A1,truss3d.materials[0]))
    
    truss3d.add_support(1,0,[0,1,2,3,4,5],val=0.0)
    truss3d.add_support(1,1,[0,1,2,3,4,5],val=0.0)
    truss3d.add_support(1,2,[0,1,2,3,4,5],val=0.0)
    
    #print(truss3d.nodes[1].dofs)
    
    truss3d.add_loadcase(supp_id=1,load_id=2)
    
    truss3d.nodal_dofs()
    K = truss3d.global_stiffness_matrix()
    p = truss3d.global_load_vector(2)
    u, K0 = truss3d.linear_statics()

    return truss3d

    #truss3d.draw()
    
def frame3dTest():
    f = FrameFEM()
    
    f.add_material(200,0.3,7850e-9)
    
    n0 = [0.0,0.0,0.0]
    
    
    truss3d.add_node(0.0,0.0,-4e3)
    truss3d.add_node(-3e3,0.0,0.0)
    truss3d.add_node(0.0,0.0,4e3)
    truss3d.add_node(0.0,5e3,0.0)

if __name__ == "__main__":
    
    t = truss3dTest()
    


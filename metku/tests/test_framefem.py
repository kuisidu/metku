# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
import unittest

#from .framefem import FrameFEM, BeamSection, LineLoad, PointLoad
#from .elements import EBBeam

import numpy as np
import metku.framefem as ff
from metku.framefem.elements.ebbeam import EBBeam, EBBeam3D
from metku.framefem.elements.springs import LinearSpring
from metku.sections.steel.ISection import IPE, HEA

class TestFrameFEM(unittest.TestCase):
    def test_adding_material(self):
        fem = FrameFEM()
        fem.add_material(100, 200, 300)
        mat = fem.materials[0]
        self.assertEqual(100, mat.young)
        self.assertEqual(200, mat.nu)
        self.assertEqual(300, mat.density)

    def test_adding_node(self):
        fem = FrameFEM()
        fem.add_node(0, 0)
        node = fem.nodes[0]
        self.assertEqual([0, 0], [node.x, node.y])

    def test_linear_statics(self):
        fem = FrameFEM()
        # Nodes
        fem.add_node(0, 0)
        fem.add_node(0, 1)
        fem.add_node(1, 1)
        fem.add_node(1, 0)
        # Material
        fem.add_material(210e3, 0.3, 7850e-9)
        # Supports
        fem.add_support(1, 0, [0, 1, 2])
        fem.add_support(1, 3, [0, 1, 2])
        # Sections
        sect = BeamSection(1e3, 2e6)
        fem.add_section(sect)
        # Elements
        for nid in range(fem.nnodes() - 1):
            n1 = fem.nodes[nid]
            n2 = fem.nodes[nid + 1]
            mat = fem.materials[0]
            sect = fem.sections[0]
            ele = EBBeam(n1, n2, sect, mat)
            fem.add_element(ele)
        # Loads
        pointload = PointLoad(1, fem.nodes[1], [10, 0, 0], f=1)
        lineload = LineLoad(1, fem.elements[1], [0, 1], [-10, -10], 1)
        fem.add_load(pointload)
        fem.add_load(lineload)
        # Loadcase
        fem.add_loadcase(1, 1)
        # Linear statics
        fem.nodal_dofs()
        fem.linear_statics()
        

def simple_portal_fem(L,H,ne):
    """ Yksinkertainen portaalikehä """
    pf = ff.FrameFEM()
    
    col_sect = IPE(220)
    beam_sect = IPE(220)
    
    #col = ff.BeamSection(col_sect.A,col_sect.Iy)
    #beam = ff.BeamSection(beam_sect.A,beam_sect.Iy)
    
    pf.add_section(col_sect)
    pf.add_section(beam_sect)
    
    pf.add_material(col_sect.E,col_sect.material.nu,col_sect.density)
    # Generoidaan solmut
    
    Y = np.linspace(0,H,ne+1)
    X = np.linspace(0,L,ne+1)
    
    for y in Y[:-1]:
        pf.add_node(X[0],y)
    
    for x in X[:-1]:
        pf.add_node(x,Y[-1])
    
    for y in reversed(Y):
        pf.add_node(X[-1],y)
    
    for (k,node) in enumerate(pf.nodes[:-1]):
        if (node.x == 0 and pf.nodes[k+1].x == 0) or (node.x == L and pf.nodes[k+1].x == L):
                el_sec = col_sect
        else:
            el_sec = beam_sect
            
            
        pf.add_element(EBBeam(node,pf.nodes[k+1],el_sec,pf.materials[0]))
        
    pf.add_support(sid=0,nid=0,dof=[0,1,2],val=0.0)
    pf.add_support(sid=0,nid=pf.nnodes()-1,dof=[0,1,2],val=0.0)
    
    
    #F = -200e3
    F = -5e3
    
    pl1 = ff.PointLoad(sid=1,node=pf.nodes[ne],v=[F,0,0],f=1.0)
    pl2 = ff.PointLoad(sid=1,node=pf.nodes[2*ne],v=[F,0,0],f=1.0)
    pf.add_load(pl1)
    pf.add_load(pl2)
    
    for i in range(ne,2*ne):
        #print(i)
        pf.add_load(ff.LineLoad(sid=1,eid=pf.elements[i],xloc=[0,1],qval=[-20,-20],direction='y'))
    #pf.add_load(ll1)    
    
    #pl3 = ff.PointLoad(sid=1,node=pf.nodes[ne],v=[-F,0,0],f=1.0)
    #pf.add_load(pl3)
    
    pf.add_loadcase(supp_id=0,load_id=1)
    
    pf.nodal_dofs()
    
    KG = []
    
    pf.linear_statics(lcase=1,support_method="REM")
    #(w,v,KG) = pf.linear_buckling(k=10)
    #print(w)
    #nd = np.argmin(abs(pf.load_factors))
    #nd_sort = np.argsort(abs(pf.load_factors))
    #pf.draw(buckling_mode=nd_sort[0],scale=40000.0)
    pf.draw(deformed=1,scale=50)
    
    return pf, KG

def simple_column(L,ne=6):
    """ Mastopilari """
    pf = ff.FrameFEM()
    
    col_sect = IPE(220)
    
    pf.add_section(col_sect)
    
    pf.add_material(col_sect.E,col_sect.material.nu,col_sect.density)
    # Generoidaan solmut
    
    Y = np.linspace(0,L,ne+1)
    
    for y in Y:
        pf.add_node(0,y)
    
    for (k,node) in enumerate(pf.nodes[:-1]):                    
        pf.add_element(EBBeam(node,pf.nodes[k+1],col_sect,pf.materials[0]))
        
    pf.add_support(sid=0,nid=0,dof=[0,1,2],val=0.0)
    
    F = -200e3
    
    pl1 = ff.PointLoad(sid=1,node=pf.nodes[ne],v=[0,F,0],f=1.0)
    pf.add_load(pl1)
    
    pf.add_loadcase(supp_id=0,load_id=1)
    
    pf.nodal_dofs()
    pf.linear_statics(support_method="REM")
    """
    (w,v, KG) = pf.linear_buckling(k=3)
    #print(w)
    nd_sort = np.argsort(abs(pf.load_factors))
    nd = np.argmin(abs(pf.load_factors))
    print(nd)
    pf.draw(buckling_mode=nd_sort[0],scale=1.0)
    """
    return pf #, KG
    
def test_hinge():
    
    a = 5000
    b = 2000
    
    f = ff.FrameFEM()
    
    sect = IPE(200)
    
    f.add_section(sect)
    
    f.add_material(sect.E,sect.material.nu,sect.density)
    
    f.add_node(0,0)
    f.add_node(a,0)
    f.add_node(a+b,0)
    
    f.add_element(EBBeam(f.nodes[0],f.nodes[1],sect,f.materials[0]))
    f.add_element(EBBeam(f.nodes[1],f.nodes[2],sect,f.materials[0]))

    f.add_support(sid=0,nid=0,dof=[0,1,2],val=0.0)
    f.add_support(sid=0,nid=2,dof=[0,1,2],val=0.0)
    
    F = -200e3
    
    pl1 = ff.PointLoad(sid=1,node=f.nodes[1],v=[0,F,0],f=1.0)
    f.add_load(pl1)
    
    ll1 = ff.LineLoad(sid=1,eid=f.elements[0],xloc=[0,1],qval=[-10,-10],direction='y')
    #f.add_load(ll1)    
    f.add_loadcase(supp_id=0,load_id=1)
    
    # Vapautetaan momentti
    f.add_release(0,[5])
    #f.add_release(1,[2])
    
    f.nodal_dofs()
    f.linear_statics(lcase=1,support_method="REM")

    f.draw()
    
    return f

def test_springs():
    
    f = ff.FrameFEM()
    
    a = 1.0
    
    f.add_node(0,0)
    f.add_node(a,0)
    f.add_node(2*a,0)
    f.add_node(3*a,0)
    
    f.add_element(LinearSpring(f.nodes[0],f.nodes[1],k=1000))
    f.add_element(LinearSpring(f.nodes[1],f.nodes[2],k=2000))
    f.add_element(LinearSpring(f.nodes[2],f.nodes[3],k=3000))

    f.add_support(sid=0,nid=0,dof=[0,1,2],val=0.0)
    f.add_support(sid=0,nid=3,dof=[0,1,2],val=0.0)
    
    F = 5000
    
    pl1 = ff.PointLoad(sid=1,node=f.nodes[2],v=[F,0,0],f=1.0)
    f.add_load(pl1)
    
    f.add_loadcase(supp_id=0,load_id=1)
    
    
    f.nodal_dofs()
    f.linear_statics(lcase=1,support_method="REM")

    return f

def test_beam_and_spring():
    
    a = 3
    
    f = ff.FrameFEM()
    
    sect = ff.BeamSection(A=0.01,Iy=2e-4)
    
    f.add_section(sect)
    
    f.add_material(E=210e6,nu=.3,density=7850)
    
    f.add_node(0,0)
    f.add_node(a,0)
    f.add_node(2*a,0)
    f.add_node(2*a,-1e-3)
    
    f.add_element(EBBeam(f.nodes[0],f.nodes[1],sect,f.materials[0]))
    f.add_element(EBBeam(f.nodes[1],f.nodes[2],sect,f.materials[0]))
    f.add_element(LinearSpring(f.nodes[2],f.nodes[3],k=200))

    f.add_support(sid=0,nid=0,dof=[0,1,2],val=0.0)
    f.add_support(sid=0,nid=1,dof=[1],val=0.0)
    f.add_support(sid=0,nid=3,dof=[0,1],val=0.0)
    
    F = -50
    
    pl1 = ff.PointLoad(sid=1,node=f.nodes[2],v=[0,F,0],f=1.0)
    f.add_load(pl1)
    
    #ll1 = ff.LineLoad(sid=1,eid=f.elements[0],xloc=[0,1],qval=[-10,-10],direction='y')
    #f.add_load(ll1)    
    f.add_loadcase(supp_id=0,load_id=1)
    
    # Vapautetaan momentti
    #f.add_release(0,[5])
    #f.add_release(1,[2])
    
    f.nodal_dofs()
    f.linear_statics(lcase=1,support_method="REM")

    #f.draw()
    
    return f

def test3d_frame():
    
    L1 = 3000
    L2 = 6000
    h = 4500
    X = []
    X.append([0,0,0])
    X.append([0,L2,0])
    X.append([L1,L2,0])
    X.append([L1,0,0])
    X.append([0,0,h])
    X.append([0,L2,h])
    X.append([L1,L2,h])
    X.append([L1,0,h])
        
    col_sect = HEA(240)
    beam_sect = IPE(220)
    
    #col = ff.BeamSection(col_sect.A,col_sect.Iy)
    #beam = ff.BeamSection(beam_sect.A,beam_sect.Iy)
    
    f = ff.FrameFEM()
    
    for x in X:
        f.add_node(x[0],x[1],x[2])

    f.add_section(col_sect)
    f.add_section(beam_sect)
    
    f.add_material(col_sect.E,col_sect.material.nu,col_sect.density)

    for i in range(4):
        n1 = f.nodes[i]
        n2 = f.nodes[i+4]
        f.add_element(EBBeam3D(n1, n2, col_sect,f.materials[0]))

    for i in range(4,7):
        n1 = f.nodes[i]
        n2 = f.nodes[i+1]
        f.add_element(EBBeam3D(n1, n2, beam_sect,f.materials[0]))
        
        if i+1 == 7:
            n1 = f.nodes[i+1]
            n2 = f.nodes[4]
            f.add_element(EBBeam3D(n1, n2, beam_sect,f.materials[0]))
        
    for i in range(4):        
        f.add_support(sid=0,nid=i,dof=[0,1,2,3,4,5])
    
    F = -50e3
    
    pl1 = ff.PointLoad(sid=1,node=f.nodes[7],v=[0,0,F,0,0,0],f=1.0)
    f.add_load(pl1)
    
    #ll1 = ff.LineLoad(sid=1,eid=f.elements[0],xloc=[0,1],qval=[-10,-10],direction='y')
    #f.add_load(ll1)    
    f.add_loadcase(supp_id=0,load_id=1)
    
    # Vapautetaan momentti
    #f.add_release(0,[5])
    #f.add_release(1,[2])
    
    f.nodal_dofs()
    f.linear_statics(lcase=1,support_method="REM")
    
    f.draw()
    

    return f

def test3d_frame2():
    
    L1 = 3000
    L2 = 6000
    h = 4500
    X = []

    X.append([L2,L1,0])
    X.append([L2,L1,h])
    X.append([0,L1,h])    
    #X.append([L2,0,h])
    
    col_sect = HEA(240)
    beam_sect = IPE(220)
    
    #col = ff.BeamSection(col_sect.A,col_sect.Iy)
    #beam = ff.BeamSection(beam_sect.A,beam_sect.Iy)
    
    f = ff.FrameFEM()
    
    for x in X:
        f.add_node(x[0],x[1],x[2])

    f.add_section(col_sect)
    f.add_section(beam_sect)
    
    f.add_material(col_sect.E,col_sect.material.nu,col_sect.density)

    f.add_element(EBBeam3D(f.nodes[0], f.nodes[1], col_sect,f.materials[0]))
    f.add_element(EBBeam3D(f.nodes[2], f.nodes[1], beam_sect,f.materials[0]))
    #f.add_element(EBBeam3D(f.nodes[2], f.nodes[1], beam_sect,f.materials[0]))

          
    f.add_support(sid=0,nid=0,dof=[0,1,2,3,4,5])
    f.add_support(sid=0,nid=2,dof=[0,1,2,3,4,5])
    #f.add_support(sid=0,nid=3,dof=[0,1,2,3,4,5])
    
    F = -100e3
    
    pl1 = ff.PointLoad(sid=1,node=f.nodes[1],v=[0,0,F,0,0,0],f=1.0)
    f.add_load(pl1)
    
    #ll1 = ff.LineLoad(sid=1,eid=f.elements[0],xloc=[0,1],qval=[-10,-10],direction='y')
    #f.add_load(ll1)    
    f.add_loadcase(supp_id=0,load_id=1)
    
    # Vapautetaan momentti
    #f.add_release(0,[5])
    #f.add_release(1,[2])
    
    f.draw()
    
    f.nodal_dofs()
    f.linear_statics(lcase=1,support_method="REM")
    
    
    

    return f

def test3d_beam():
    
    L1 = 3000
    L2 = 6000
    h = 4500
    X = []
    X.append([0,0,0])
    X.append([0,0,h])
        
    col_sect = HEA(240)
    beam_sect = IPE(220)
    
    #col = ff.BeamSection(col_sect.A,col_sect.Iy)
    #beam = ff.BeamSection(beam_sect.A,beam_sect.Iy)
    
    f = ff.FrameFEM()
    

    for x in X:
        f.add_node(x[0],x[1],x[2])

    f.add_section(col_sect)
    #f.add_section(beam_sect)
    
    f.add_material(col_sect.E,col_sect.material.nu,col_sect.density)

    f.add_element(EBBeam3D(f.nodes[0], f.nodes[1], col_sect,f.materials[0]))

    f.add_support(sid=0,nid=0,dof=[0,1,2,3,4,5])
    
    F = 100e3
    
    pl1 = ff.PointLoad(sid=1,node=f.nodes[1],v=[0,0,-F,0,0,0],f=1.0)
    f.add_load(pl1)
    
    #ll1 = ff.LineLoad(sid=1,eid=f.elements[0],xloc=[0,1],qval=[-10,-10],direction='z',coords='local')
    #f.add_load(ll1)    
    f.add_loadcase(supp_id=0,load_id=1)
    
    # Vapautetaan momentti
    #f.add_release(0,[5])
    #f.add_release(1,[2])
    
    f.draw()
    
    f.nodal_dofs()
    f.linear_statics(lcase=1,support_method="REM")
    
    
    

    return f

if __name__ == '__main__':
    
    #pf, KG = simple_portal_fem(L=4000,H=4000,ne=4)
    
    #pf = simple_column(L=8000,ne=6)
    
    #f = test_hinge()

    #f = test_springs()
    
    #f = test_beam_and_spring()
    
    f = test3d_frame()
    #f = test3d_frame2()
    #f = test3d_beam()
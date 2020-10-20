import unittest

#from .framefem import FrameFEM, BeamSection, LineLoad, PointLoad
#from .elements import EBBeam

import numpy as np
import metku.framefem as ff
from metku.framefem.elements.ebbeam import EBBeam 
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
    
    
    F = -200e3
    
    pl1 = ff.PointLoad(sid=1,node=pf.nodes[ne],v=[0,F,0],f=1.0)
    pl2 = ff.PointLoad(sid=1,node=pf.nodes[2*ne],v=[0,F,0],f=1.0)
    pf.add_load(pl1)
    pf.add_load(pl2)
    
    #pl3 = ff.PointLoad(sid=1,node=pf.nodes[ne],v=[-F,0,0],f=1.0)
    #pf.add_load(pl3)
    
    pf.add_loadcase(supp_id=0,load_id=1)
    
    pf.nodal_dofs()
    pf.linear_statics(support_method="REM")
    (w,v,KG) = pf.linear_buckling(k=10)
    print(w)
    nd = np.argmin(abs(pf.load_factors))
    nd_sort = np.argsort(abs(pf.load_factors))
    pf.draw(buckling_mode=nd_sort[0],scale=40000.0)
    
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
    
    

if __name__ == '__main__':
    
    pf, KG = simple_portal_fem(L=4000,H=4000,ne=20)
    
    #pf = simple_column(L=8000,ne=6)
    


# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 16:21:20 2022

Plane truss example using raami package

@author: Kristo Mela
"""

from raami.raami_plane_truss import Ktruss_example

if __name__ == "__main__":
    
    t =  Ktruss_example(h2=2000,h1=1500,dx1=1000,dx2=1000,first=False,edges=True)
    t.symmetry()
    t.generate_supports()
    t.generate_joints()
    t.generate_uniform_load(q=-25)
    t.generate_fem(model='ecc_elements')
    #t.generate_fem(model='en1993')
    #t.generate_fem(model='no_eccentricity')
    #t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    x, xc = t.joints[4].brace_chord_face_x()
    #t.bmd(scale=10,load_id=t.load_ids[0])

t.plot(geometry=True,loads=False)

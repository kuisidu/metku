# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Sun Apr 10 16:21:20 2022

Plane truss example using raami package

@author: Kristo Mela
"""

from raami.raami_plane_truss import Ktruss_example, Ntruss_example
from raami.raami_plane_truss import SlopedTruss
from raami.exports import AbaqusOptions
from raami.raami_truss_opt import minimize_eccentricity
from optimization.solvers.slsqp import SLSQP
from optimization.solvers.trust_region import TrustRegionConstr

def LauriKTruss(span,h2,h1,dx,nel_chord=4,nel_brace=4,ndiv=4):
    # Create K truss for Lauri
    t = SlopedTruss(L1=0.5*span,L2=0.5*span,h2=h2,h1=h1,dx1=dx,dx2=dx)
    t.braces_as_beams = True
    t.generate_topology('K',ndiv,nel_chord=nel_chord,nel_brace=nel_brace)
    t.symmetry()
    t.generate_supports()
    t.generate_joints()
    t.generate_uniform_load(q=-25)
    t.generate_fem(model='no_eccentricity')
        
    print("Structural analysis..")
    t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    print("Done.")
    
    top = {'material': 'S700', 'class': 2, 'utility': 1.0}
    bottom = {'material': 'S700', 'class': 2, 'utility': 1.0}
    braces = {'material': 'S700', 'class': 2, 'utility_tens': 1.0, 'utility_comp':1.0}
    
    t.optimize_members(verb=True,top=top,bottom=bottom,braces=braces)
     
    """
    P, x0 = minimize_eccentricity(t,min_gap=20)   
    
        
    solver = SLSQP()
    #solver = TrustRegionConstr()
    min_ecc, xmin = solver.solve(P,x0=x0,verb=True)
    #t.plot(geometry=True,loads=False)
    
    t.clear_fem()
    t.generate_fem(model="ecc_elements")
    #t.fem.draw()
    
    #t.generate_fem(model='ecc_elements')
    opts = AbaqusOptions(x_monitor = 0.5*t.span, n_monitored = 2)
    t.to_abaqus(filename='K-ristikko',partname="K-ristikko",options=opts)
    """
    
    return t
    

if __name__ == "__main__":
    
    t = LauriKTruss(span=24000,h2=2400,h1=1800,dx=1000,nel_chord=4,nel_brace=4,ndiv=4)
    #t =  Ktruss_example(h2=2000,h1=1500,dx1=1000,dx2=1000,first=False,edges=True)
    #t =  Ntruss_example(h2=2000,h1=1500,dx1=1500,dx2=1500,first=False,edges=True)
    
    #t.generate_fem(model='en1993')
    #t.generate_fem(model='no_eccentricity')
    #t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    #x, xc = t.joints[4].brace_chord_face_x()
    #t.bmd(scale=10,load_id=t.load_ids[0])

#t.plot(geometry=True,loads=False)

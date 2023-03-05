# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 16:48:20 2022

Optimization tools for Trusses using Raami framework

@author: kmela
"""

import numpy as np


from metku.eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint
from metku.raami.raami_plane_truss import TubularKGapJoint
from metku.optimization.structopt import OptimizationProblem
from metku.optimization.variables import Variable
from metku.optimization.constraints import NonLinearConstraint, LinearConstraint
from metku.optimization.objective import ObjectiveFunction
from metku.optimization.solvers.slsqp import SLSQP
from metku.optimization.solvers.trust_region import TrustRegionConstr


def minimize_eccentricity(truss,min_gap=20):
    """ Minimizes eccentrictity of the joints in a tubular truss """
    
    def gap_constraint(joint,min_gap=20):
        """ Generates constraints for the gap in the joint 
            The gap must satisfy:
                0.5*b0*(1-beta) <= gap <= 1.5*b0*(1-beta)
                gap >= t1 + t2
                gap >= min_gap
        """
        beta = joint.joint.beta()
        b0 = joint.joint.b0
        gub = 1.5*b0*(1-beta)
        
        glb = max([0.5*b0*(1-beta),min_gap,sum(joint.t)])
        
        def gap_ub(x):
            """ Returns gap upper bound constraint """
            gap = joint.calc_gap()
            
            return gap - gub
            #return gap/gub-1
        
        def gap_lb(x):
            """ Returns gap lower bound constraint """
            gap = joint.calc_gap()                
            
            return -gap + glb
            #return 1-gap/glb
        
        return gap_ub, gap_lb
    
    def total_eccentricity(opt_vars):
        # Objective function is the sum of eccentricities at K joints.
        C = 100
        yvars = []
        n = len(opt_vars)
        for i, var in enumerate(opt_vars):
            if var.target['property'] == 'yloc':
                yvars.append(i)
        
        def obj_value(x):
            #return 0
            return 1/C*sum(np.array(x)[yvars]**2)
        
        def obj_grad(x):
            """ Gradient of the objective function """
            df = np.zeros(n)
            df[yvars] = 2/C*np.array(x)[yvars]
            
            return df
        
        return obj_value, obj_grad
        
    
    xlb = -200
    xub = 200
    # Create empty optimization problem
    P = OptimizationProblem(structure=truss)
    
    # Create joint variables and constraints
    for key, joint in truss.joints.items():
        # Assume symmetry of the truss, so only need to consider
        # half of the joints.
        if joint.node.x <= 0.5*truss.span:            
            if joint == truss.ridge_joint:
                new_var = Variable(f'y_{key}',lb=-0.55*joint.h0,ub=0.25*joint.h0,value=0.1*joint.h0,
                                   target={'property':'yloc','objects':[joint]})
                P.add(new_var)
                
                # Add constraints for the gap
                gap_ub, gap_lb = gap_constraint(joint,min_gap)
                gap_ub_con = NonLinearConstraint(gap_ub,con_type="<",name=f'Gap upper bound, joint {key}',problem=P)
                gap_lb_con = NonLinearConstraint(gap_lb,con_type="<",name=f'Gap lower bound, joint {key}',problem=P)
                P.add(gap_ub_con)
                P.add(gap_lb_con)
            elif isinstance(joint,TubularKGapJoint):
                
                new_y_var = Variable(f'y_{key}',lb=-0.55*joint.h0,ub=0.25*joint.h0,value=0.1*joint.h0,
                                         target={'property':'yloc','objects':[joint]})
                new_x_var = Variable(f'x_{key}',lb=xlb,ub=xub,value=0,
                                         target={'property':'xloc','objects':[joint]})
                
                P.add(new_y_var)                
                P.add(new_x_var)
                
                # Add constraints for the gap
                gap_ub, gap_lb = gap_constraint(joint,min_gap)
                gap_ub_con = NonLinearConstraint(gap_ub,con_type="<",name=f'Gap upper bound, joint {key}',problem=P)
                gap_lb_con = NonLinearConstraint(gap_lb,con_type="<",name=f'Gap lower bound, joint {key}',problem=P)
                P.add(gap_ub_con)
                P.add(gap_lb_con)
    
    x0 = np.array([var.value for var in P.vars])
    
    # Generate constraints:
    # For vertical members, add constraints that keep the members vertical
    n = len(P.vars)
    for key, brace in truss.braces.items():        
        if brace.vertical:
            
            a = np.zeros(n)
            # We need to identify the joints of the member
            # and the corresponding variables in order to write constraints.
            
            for i, var in enumerate(P.vars):
                if var.target['objects'][0].node == brace.nodes[0]:
                    print("Variable connected to node 1")
                    if var.target['property'] == 'xloc':
                        print("Variable related to local x coordinate")
                        a[i] = var.target['objects'][0].chords[0].dir_vector[0]
                    elif var.target['property'] == 'yloc':
                        print("Variable related to local y coordinate")
                        a[i] = var.target['objects'][0].chords[0].perpendicular[0]
                elif var.target['objects'][0].node == brace.nodes[1]:
                    print("Variable connected to node 2")
                    if var.target['property'] == 'xloc':
                        print("Variable related to local x coordinate")
                        a[i] = -var.target['objects'][0].chords[0].dir_vector[0]
                    elif var.target['property'] == 'yloc':
                        print("Variable related to local y coordinate")
                        a[i] = -var.target['objects'][0].chords[0].perpendicular[0]
            
            if any(a):            
                P.add(LinearConstraint(a,b=0,con_type="=",name=f"Vertical brace {key}"))
            
    #
    
    # Objective function
    obj_fun, obj_grad = total_eccentricity(P.vars)
    P.add(ObjectiveFunction("Eccentricity", obj_fun))
    P.grad = obj_grad
    
    return P, x0

if __name__ == "__main__":
    
    from metku.raami.raami_plane_truss import Ktruss_example, Ntruss_example
    
    #t = Ktruss_example(h2=1800,h1=1500,dx1=1000, dx2=1000,first=True, edges=False)
    t = Ntruss_example(h2=1700,dx1=2000, dx2=2000, first=False,edges=True)
    t.symmetry()
    t.generate_supports()
    t.generate_joints()
    t.determine_symmetry_joints()
    t.generate_uniform_load(q=-25)
    t.generate_fem(model='en1993')
    #t.generate_fem(model='no_eccentricity')
    
    t.optimize_members(verb=True)
    t.plot(geometry=True,loads=False)
    
    P, x0 = minimize_eccentricity(t)
    #x0[0] = 20
    
    solver = SLSQP()
    #solver = TrustRegionConstr()
    min_ecc, xmin = solver.solve(P,x0=x0,verb=True)
    
    t.plot(geometry=True,loads=False)
    
    #P(x0)
    

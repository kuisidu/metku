# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Tue Oct 6 12:17 2020

Augmented Lagrangian method

For unconstrained optimization

@author: kmela
"""

import numpy as np

from scipy.optimize import minimize
import  matplotlib.pyplot as plt


import metku.optimization.structopt as sopt
from metku.optimization.solvers.optsolver import OptSolver



class ALM(OptSolver):
    """ Augmented Lagrangian method """

    def __init__(self, problem):
        """ Input:
            problem .. OptimizationProblem type object
        """
        super().__init__()
        self.result = None
        self.problem = problem
        
    def solve(self, x0, maxiter=200, Rineq=1.0,Req=1.0):
        """ Augmented Lagrangian method iteration """
        
        """ Display commands:
                'latex' .. print iteration results as a LaTeX table
                'multpliers' .. print Lagrange multipliers over iterations
                'iteration' .. print iteration data
        """
        #disp = 'latex'
        #disp = 'multipliers'
        disp = 'iteration'
        
        """ Flag for displaying results of each minimization of the
            augmented Lagrangian
        """
        disp_min = True
        
        k = 0
        x = [x0]
        
        self.xvals.append(x0)
        self.fvals.append(self.problem.obj(x0))
        
        """ Initialize augmented Lagrangian """
        aug_fun = self.problem.augmented_lagrangian(Rineq,Req)
        
        f_stop = 0
        
        """ Cineq: multiplier for penalty term for inequality constraints
            Ceq: multiplier for penalty term for equality constraints
        """
        Cineq = 1.5
        Ceq = 1.5
        
        print("*** Augmented Lagrangian method ***")
        
        while k < maxiter:
            
            if disp == "iteration":
                print("Iteration {0:4g}".format(k))
                print("Rineq = {0:5g}".format(Rineq))
                print("Lagrange multipliers (estimates):")
            
                for (i,con) in enumerate(self.problem.cons):
                    if con.type == "<":
                        s = 'g'.join(str(i))
                    else:
                        s = 'h'.join(str(i))
                    print("{0:s}: {1:6.4f}".format(s,con.mult))
                    
                for (i,var) in enumerate(self.problem.vars):
                    print("x{0:g} lower: {1:6.4f}".format(i,var.lb_mult))
                    print("x{0:g} upper: {1:6.4f}".format(i,var.ub_mult))
            elif disp == "multipliers":
                smult = str(k+1)
                for con in self.problem.cons:
                    smult += ' & ' + '{0:6.4f}'.format(con.mult)
                
                for var in self.problem.vars:
                    smult += ' & ' + '{0:6.4f}'.format(var.lb_mult)
                    smult += ' & ' + '{0:6.4f}'.format(var.ub_mult)
            
                smult += ' \\\\'
                print(smult)
            
            """ Minimize augmented Lagrangian function 
                By default, gradient is evaluated by differences
                and conjugate gradient method is used
            """
            res = minimize(aug_fun,x[k],method='CG',jac='2-point',
                           options={'disp':disp_min,'finite_diff_rel_step':1e-8})
            
            
            x.append(res.x)
            self.xvals.append(res.x)
            self.fvals.append(self.problem.obj(res.x))
            
            if disp == 'latex':
                xs = '& '
            else:
                xs = ','
                
            xs = xs.join('{0:6.4f} '.format(xi) for xi in x[k+1])
            
            if disp == 'iteration':
                print("Objective: {0:6.3f}, x = [{1:s}]".format(res.fun,xs))
                print("f(x) = {0:6.4f}".format(self.fvals[k+1]))
            elif disp == 'latex':
                """ Print stuff to a LaTeX table """
                siter = str(k+1) + ' & ' + str(Rineq) + ' & '
                
                siter += xs
                siter += ' & {0:6.3f} '.format(res.fun)
                siter += '& {0:6.3f}'.format(self.fvals[k+1])
                

            if disp is not None:            
                for (i,con) in enumerate(self.problem.cons):
                    if disp == 'iteration':
                        s = 'g' + str(i+1)
                        print('{0:s}(x) = {1:6.4f}'.format(s,con(x[k+1])))
                    elif disp == 'latex':
                        siter += ' & {0:6.4f}'.format(con(x[k+1]))
            
            if disp == 'latex':
                siter += ' \\\\'
                print(siter)
            
            """ Check stopping criterion """
            if self.stopping_criterion('x',x=[x[k+1],x[k]]):
                break
                                       
            #if np.linalg.norm(x[k+1]-x[k]) <= X_TOL*np.linalg.norm(x[k]):
            #    break            
            
            if self.stopping_criterion('f',f=self.fvals[k:k+2]):
                f_stop += 1
            #if abs(self.fvals[k+1]-self.fvals[k]) <= ABS_F_TOL + REL_F_TOL*abs(self.fvals[k]):
            #    f_stop += 1
                
            if f_stop == 2:
                break
            
            """ Update Lagrange multipliers """
            self.problem.update_lag_mult(x[k+1],Rineq,Req)            
            
            """ Update penalty factors """
            Rineq *= Cineq
            Req *= Ceq
            
            aug_fun = self.problem.augmented_lagrangian(Rineq,Req)
            
            k += 1

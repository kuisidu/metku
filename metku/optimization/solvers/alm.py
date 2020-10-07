# -*- coding: utf-8 -*-
"""
Created on Tue Oct 6 12:17 2020

Augmented Lagrangian method

For unconstrained optimization

@author: kmela
"""

import numpy as np

from scipy.optimize import minimize
import  matplotlib.pyplot as plt

try:
    import metku.optimization.structopt as sopt
    from metku.optimization.solvers.optsolver import OptSolver
except:
    import optimization.structopt as sopt
    from optimization.solvers.optsolver import OptSolver


GRAD_TOL = 1e-8
ABS_F_TOL = 1e-8
REL_F_TOL = 1e-4
X_TOL = 1e-8

def line_search_obj(fun,xk,dk):
    """ Objectve function for line search """
    
    def line_fun(a):
        #print(xk + a*dk)
        return fun(xk + a*dk)
    
    return line_fun

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
        
        #disp = 'latex'
        #disp = 'multipliers'
        disp = 'iteration'
        disp_min = True
        
        k = 0
        x = [x0]
        
        self.xvals.append(x0)
        self.fvals.append(self.problem.obj(x0))
        
        aug_fun = self.problem.augmented_lagrangian(Rineq,Req)
        
        f_stop = 0
        Cineq = 2
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
            
            """ Minimize augmented Lagrangian function """
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
            if np.linalg.norm(x[k+1]-x[k]) <= X_TOL*np.linalg.norm(x[k]):
                break            
            
            #print(x)
            if abs(self.fvals[k+1]-self.fvals[k]) <= ABS_F_TOL + REL_F_TOL*abs(self.fvals[k]):
                f_stop += 1
                
            if f_stop == 2:
                break
            
            """ Update Lagrange multipliers """
            self.problem.update_lag_mult(x[k+1],Rineq,Req)            
            
            """ Update penalty factors """
            Rineq *= Cineq
            Req *= Ceq
            
            aug_fun = self.problem.augmented_lagrangian(Rineq,Req)
            
            k += 1

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 17:12:02 2020

Steepest descent method

For unconstrained optimization

@author: kmela
"""

import numpy as np

from scipy.optimize import minimize_scalar
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

def line_search_obj(fun,xk,dk):
    """ Objectve function for line search """
    
    def line_fun(a):
        #print(xk + a*dk)
        return fun(xk + a*dk)
    
    return line_fun

class SteepestDescent(OptSolver):
    """ Steepest descent method """

    def __init__(self, fun, grad):
        """ Input:
            fun .. function to be minimized.
            grad .. gradient of the function. If None is given, gradient is
                    evaluated numerically
                    
            TODO:
                implement numerical gradient evaluation!
        """
        super().__init__()
        self.result = None
        self.obj = fun
        self.grad = grad
        
    def solve(self, x0, maxiter=200):
        """ Perform steepest descent iteration """
        
        k = 0
        x = [x0]
        d = []
        
        self.xvals.append(x0)
        self.fvals.append(self.obj(x0))
        
        f_stop = 0
        
        while k < maxiter:
            
            df = self.grad(x[k])
            
            ndf = np.linalg.norm(df,2)
            
            """ Check gradient stopping criterion """
            if ndf < GRAD_TOL:
                break
            else:
                d.append(-df/ndf)
            
            """ Perform line search """
            
            """
            ls_fun = line_search_obj(self.obj,x[k],d[k])
            a = np.linspace(0,5)
            y = np.zeros(a.shape)
            for i in range(a.shape[0]):
                y[i] = ls_fun(a[i])
            
            plt.plot(a,y)
            """
            
            ls_res = minimize_scalar(line_search_obj(self.obj,x[k],d[k]),
                                     method='brent')
            
            #print(d[k])
            #print("Line search result:", ls_res.x)
            #print(ls_res.success)
            if ls_res.success:
                x.append(x[k] + ls_res.x*d[k])
                self.xvals.append(x[-1])
                self.fvals.append(self.obj(x[-1]))
                
            
            #print(x)
            if abs(self.fvals[-1]-self.fvals[-2]) <= ABS_F_TOL + REL_F_TOL*abs(self.fvals[-2]):
                f_stop += 1
                
            if f_stop == 2:
                break
            
            k += 1

def steep_ex_fun(x):
    
    return (x[0]-4)**2+ (4*x[1]-2)**2

def steep_ex_grad(x):
    
    return np.array([2*(x[0]-4),8*(4*x[1]-2)])


if __name__ == "__main__":
    
    steep_sol = SteepestDescent(steep_ex_fun,steep_ex_grad)
    
    steep_sol.solve(np.array([10.0,10.0]),maxiter=15)

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:39:18 2022

Method of Moving Asymptotes

@author: kmela
"""

from copy import copy
import numpy as np

from metku.optimization.opt_functions import mma_fun, NumGrad
from metku.optimization.solvers.optsolver import OptSolver
from metku.optimization.solvers.slsqp import SLSQP
from metku.optimization.variables import Variable
from metku.optimization.structopt import OptimizationProblem
import metku.optimization.constraints as cons

class MMAProblem(OptimizationProblem):
    """ Class for MMA subproblems """
    
    def __init__(self,P,x,U,L):
        """
        Generates MMA approximation of the problem at point 'x' using the asymptotes 'U' and 'L'


        Parameters
        ----------
        x : numpy array
            Approximation point.
        U : numpy array
            Upper asymptote.
        L : numpy array
            Lower asymptote.

        Returns
        -------
        Optimization problem with nonlinear function replaced by their MMA approximations.

        """
        
        super().__init__(name="MMA")
                
        self.xk = x
        self.U = [U]
        self.L = [L]        
                
        for i, var in enumerate(P.vars):
            # Create variables for the mma problem
            # Rules for variable bounds are from Svanberg (2007)
            self.add(Variable(name=var.name,
                              lb=max(var.lb,L[i]+0.1*(x[i]-L[i]),x[i]-0.5*(var.ub-var.lb)),
                              ub=min(var.ub,U[i]-0.1*(U[i]-x[i]),x[i]+0.5*(var.ub-var.lb)),
                              value=x[i]))
        
        # Approximate objective
        if P.grad is None:
            dfx = NumGrad(P.obj,x)
        else:
            dfx = P.grad(x)
            
        self.obj = mma_fun(x,P.obj(x),dfx,L,U)
        
        # Approximate constraints
        for con in P.cons:
            if isinstance(con,cons.LinearConstraint):
                # Linear constraints are simply added
                self.add(con)
            else:
                # For nonlinear constraints, MMA approximation is made                
                mma_con = copy(con)
                gx = con.con(x)
                dgx = con.df(x)
                mma_con.con = mma_fun(x,gx,dgx,L,U)
                
                self.add(mma_con)            

class MMA(OptSolver):
    
    def __init__(self):
        super().__init__()
        
        self.L = []
        self.U = []
        self.mma_probs = []
    
    def gamma(self,xk,xk1,xk2):
        """ Parameter for modifying the asymptotes 
            xk .. current iterate
            xk1 .. x(k-1)
            xk2 .. x(k-2)
        """
        xtol = 1e-6
        xprod = (xk-xk1)*(xk1-xk2)
        return np.array([0.7 if x < -xtol else 1.2 if x > xtol else 1 for x in xprod])
        
    
    def take_action(self):
        """ Make the MMA sub-problem and solve it. To create the sub-problem
            means the following:
                
                1. Make MMA approximations for each function at the current
                   iterate
                2. Define asymptotes L and U for each variable. Here,
                   variable bounds, asymptotes from previous iterations
                   and previous iterates are needed.
        """
        #print("MMA ACTION")
        #print(self.xvals)
        iteration = len(self.xvals)-1
        xk = self.xvals[-1]
        xlb = self.problem.xlb()
        xub = self.problem.xub()
        Dx = xub-xlb
                                
        #print(iteration,xk,Dx)
        
        # Update asymptotes (Svanberg (2007))
        if iteration <= 1:            
            Lnew = xk-0.5*Dx
            Unew = xk+0.5*Dx         
            #print("L = ",Lnew)
            #print("U = ",Unew)
        else:
            xk1 = self.xvals[-2]
            g = self.gamma(xk,xk1,self.xvals[-3])
            Lnew = xk-g*(xk1-self.L[-1])
            Unew = xk+g*(self.U[-1]-xk1)
            
        Lnew = np.minimum(np.maximum(xk-10*Dx,Lnew),xk-0.1*Dx)
        Unew = np.minimum(np.maximum(xk+0.01*Dx,Unew),xk+10*Dx)
        
        #print("Lnew = ",Lnew)
        #print("Unew = ",Unew)
        
        self.L.append(Lnew)
        self.U.append(Unew)        
        
        # Create MMA sub-problem
        Pmma = MMAProblem(self.problem,xk,self.U[-1],self.L[-1])
        self.mma_probs.append(Pmma)
        
        # Solve the MMA sub-problem with SQP.
        # Here, a separate solver for the MMA problem could be implemented
        sqp_solver = SLSQP()
        fmma, xmma = sqp_solver.solve(Pmma,x0=xk,verb=False)
        
        #Pmma(xmma)        
        #print(Pmma.cons)
        #print(Pmma.xlb(),Pmma.xub())
        
        return xmma
    
    def step(self,action):
        """ Update design varible values """
        self.X = action
        
        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'
        
 
        
    
    
    
    
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:00:48 2019

Sequential programming with quadratic objective function.

At each iteration, the problem

min grad f(x_k)'*(x-xk) + 0.5*||x-x_k||**2
such that grad g_j(x_k)'*(x-xk) + g_j(x_k) <= 0

is solved.

@author: kmela
"""


import scipy.optimize as sciop
import numpy as np
import time

import gurobipy as grb

try:
    from metku.optimization.solvers.optsolver import OptSolver
    from metku.optimization.structopt import *
    from metku.optimization.benchmarks import *
except:
    from optimization.solvers.optsolver import OptSolver
    from optimization.structopt import *
    from optimization.benchmarks import *
    

#from ortools.linear_solver import pywraplp

class SLQP(OptSolver):

    def __init__(self):
        super().__init__()

    def take_action(self):
        """
        Defines action to take
        :return:
        """

        """ linearize problem at X
            A .. Jacobian matrix of the constraints
            B .. RHS of the linearized constraints
            df .. gradient of the objective function
            fx .. f(x_k)
        """
        #print(self.X)
        A, B, df, fx = self.problem.linearize(*self.X)

        # bounds = [x*self.move_limits for x in self.X]
        #
        # res = linprog(df, A, B, bounds=bounds)

        
        qp_fun = lambda x: 0.5*(x-self.X).dot(x-self.X) + df.dot(x-self.X)
        qp_jac = lambda x: (x-self.X) + df
        
        qp_con = {"type":'ineq',"fun":lambda x: A.dot(x) - B,"jac":lambda x: A}

    
        x_bnd = []
        for var in self.problem.vars:
            x_bnd.append((var.lb,var.ub))
            
        print(x_bnd)

        print(qp_fun(self.X))
        print(qp_jac(self.X))
        print(qp_con["fun"](self.X))
        print(qp_con["jac"](self.X))

        qp_res = sciop.minimize(qp_fun,self.X,method='SLSQP',constraints=qp_con,bounds=x_bnd)        
        print(qp_res)

        if qp_res.success == True:
            xnew = qp_res.x
        else:
            print('Something happened with qp problem.')
            xnew = qp_res.x
        
        #print(self.X)
        #return xnew - self.X
        """    
        print(A,B)
        qp_fun = lambda d: 0.5*d.dot(d) + df.dot(d)
        qp_jac = lambda d: d + df
        
        qp_con = {"type":'ineq',"fun":lambda d: A.dot(d) - B}#,"jac":lambda x: A.T}

        d0 = -df
        d_bnd = []
        for i, var in enumerate(self.problem.vars):
            d_bnd.append((var.lb-self.X[i],var.ub-self.X[i]))
            
        #print(x_bnd)

        #print(qp_fun(self.X))
        #print(qp_jac(self.X))
        #print(qp_con["fun"](self.X))
        #print(qp_con["jac"](self.X))

        qp_res = sciop.minimize(qp_fun,d0,method='SLSQP',jac=qp_jac,constraints=qp_con,bounds=d_bnd)    
        print(qp_res)

        if qp_res.success == True:
            dnew = qp_res.x
        else:
            print('Something happened with qp problem.')
            dnew = qp_res.x
            
        return dnew
        """

        
        qp = grb.Model("qp")
        
        # Number of constraints
        n = A.shape[0]
        # Number of variables
        m = A.shape[1]

        x = {}
        #d = {}
        # Variables
        qpvars = []
        qpobj = []
        for i, var in enumerate(self.problem.vars):
            # create variable
            x[i] = qp.addVar(var.lb,var.ub)
            #d[i] = qp.addVar(lb=-200,ub=200)
            # add terms to the objective function
            qpobj.append(0.5*(x[i]-self.X[i])*(x[i]-self.X[i])+df[i]*x[i])
            #qpobj.append(0.5*d[i]*d[i]+df[i]*d[i])
        
        """ Prepare quadratic part of the objective function:
            square terms
        """
        #qpobj = [0.5*va*va for va in qpvars]
        
        
        qp.setObjective(grb.quicksum(qpobj))
        
        # Constraints
        for i in range(n):
            qp.addConstr(grb.quicksum([A[i,j]*x[j] for j in range(m)]) <= B[i])
            #qp.addConstr(grb.quicksum([A[i,j]*d[j] for j in range(m)]) <= B[i])

        
        qp.setParam("LogToConsole",0)
        #qp.update()
        qp.optimize()
        

        X = []
        for v in qp.getVars():
            print('%s %g' % (v.varName, v.x))
            X.append(v.x)
        
        #time.sleep(10)

        X = np.asarray(X)
        
        #self.qp = qp
        
        d = X-self.X
        print(d)
        #return X
        return X-self.X
        
       
        
    def step(self, action):
        
        self.X += action
        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                                self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'

if __name__ == '__main__':

    problem = FifteenBarTruss(prob_type='continuous')
    x0 = [var.ub for var in problem.vars]
    solver = SLQP()
    solver.solve(problem, x0=x0, maxiter=100, log=True, verb=True)
    problem(solver.X)
    problem.structure.plot()
    
    
    
    
    
    # #
    # # import matplotlib.pyplot as plt
    # #
    # # plt.plot(solver.fvals)
    # # plt.show()

    # problem = OptimizationProblem("Quadratic Problem")
    # problem.prob_type = 'continuous'
    # problem.obj = lambda x: x[0] ** 2 + x[1] ** 2
    # var1 = Variable("X1", 0, 5)
    # var2 = Variable("X2", 0, 5)
    # problem.add_variables([var1, var2])
    # con1 = NonLinearConstraint(lambda x: x[0] ** 2 / 20 - x[1] + 1)
    # con2 = NonLinearConstraint(lambda x: x[1] ** 2 / 20 - x[0] + 1)
    # problem.add_constraints([con1, con2])
    # solver = SLP()
    # x0 = [0.1, 0.1]
    # solver.solve(problem, x0=x0, maxiter=30)


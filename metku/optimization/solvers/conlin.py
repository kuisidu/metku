# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 00:06:01 2020

CONLIN method

Convex-linear approximation of an optimization problem

@author: kmela
"""

import numpy as np
import time

try:
    from metku.optimization.solvers.optsolver import OptSolver
    from metku.optimization.structopt import *
    from metku.optimization.benchmarks import *
except:
    from optimization.solvers.optsolver import OptSolver
    from optimization.structopt import *
    from optimization.benchmarks import *
    

from ortools.linear_solver import pywraplp

class CONLIN(OptSolver):
    """ Solver class for Convex-Linear approximation Method (CONLIN) """
    
    def __init__(self, ):
        """ Constructor 
            :param move_limits: (tuple) relative lower and upper move limits
            :param gamma: parameter for reducing the relative move limits at each iteration
            :param C: penalty factor for safe-guarding againts potential infeasibility of LP problems
            :param beta: upper bound for the feasibility variable 'beta'. If beta = 0, then the feasibility variable
                        is ignored
            :param move_limit_type: type of move limits. 
                        Possible values: 
                            'range' .. move limits are taken with respect to the entire range of variable values. 
                            'relative' .. move limits are taken with respect to the current variable values.
        """
        super().__init__()
        self.move_limits = np.asarray(move_limits)
        self.move_limits_hist = [np.asarray(move_limits)]
        self.move_limit_type = move_limit_type
        self.gamma = gamma
        self.C = C
        self.beta = beta
        
        self.include_beta = True
        
        if beta == 0:
            self.include_beta = False

    def take_action(self):
        """
        Defines action to take
        Here, this means constructing the CONLIN approximation and solving it at x
        :return:
        """
        conlin_prob = self.problem.conlin(*self.X)
        # bounds = [x*self.move_limits for x in self.X]
        #
        # res = linprog(df, A, B, bounds=bounds)

        solver = pywraplp.Solver('SLP',
                                 pywraplp.Solver.CLP_LINEAR_PROGRAMMING)

        # Number of constraints
        n = A.shape[0]
        # Number of variables
        m = A.shape[1]

        x = {}
        # Variables
        for i, var in enumerate(self.problem.vars):
            # Delta is the range of variable values
            # First, lb and ub are the allowable changes in the variable value
            if self.move_limit_type == "range":
                delta = var.ub - var.lb        
                #lb, ub = self.move_limits * delta
            elif self.move_limit_type == 'relative':
                delta = var.value
                #lb = self.move_limits[0]*var.value
                #ub = self.move_limits[1]*var.value
            
            lb, ub = self.move_limits * delta
            """ Set variable bounds, taking into account true upper and lower bounds """
            lb = max(var.lb, var.value - lb)
            ub = min(var.ub, var.value + ub)
            x[i] = solver.NumVar(lb,
                                 ub,
                                 var.name)                     
        # Beta        
        if self.include_beta:
            beta = solver.NumVar(0, self.beta, 'beta')
        # Constraints
        for i in range(n):
            if self.include_beta:
                solver.Add(solver.Sum([A[i, j] * x[j] for j in range(m)]) <= B[i] + beta)
            else:
                solver.Add(solver.Sum([A[i, j] * x[j] for j in range(m)]) <= B[i])

        # Objective
        if self.include_beta:
            solver.Minimize(solver.Sum([df[j] * x[j] for j in range(m)]) + self.C * beta)
        else:
            solver.Minimize(solver.Sum([df[j] * x[j] for j in range(m)]))


        sol = solver.Solve()


        # If solution if infeasible
        if sol == 2:
            print("Solution found was infeasible!")
            # Return something != 0, otherwise iteration will
            # stop because two consecutive iterations produce
            # too similar results
            return np.zeros_like(self.X)

        if self.include_beta and beta.solution_value() < 1:
            self.beta = beta.solution_value()

        X = []
        for i in range(len(x)):
            X.append(x[i].solution_value())
        X = np.asarray(X)

        return X - self.X

# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela


import numpy as np

from ortools.linear_solver import pywraplp

from metku.optimization.solvers.optsolver import OptSolver

class SLP(OptSolver):
    """ Solver class for Sequential Linear Programming Method (SLP) """
    
    def __init__(self, move_limits=(0.05, 0.05), gamma=1e-2, C=2e5, beta=10, move_limit_type='range'):
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
        :return:
        """
        A, B, df, fx = self.problem.linearize(*self.X)
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

    def step(self, action):

        self.C += self.C * self.gamma
        """ Reduce move limits by 'gamma'x100 percent"""        
        self.move_limits -= self.move_limits * self.gamma        
        self.move_limits_hist.append(list(self.move_limits))
        #print(self.X,action)
        self.X += action
        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                                self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'

if __name__ == '__main__':
    from metku.optimization.benchmarks import *
    problem = ThreeBarTruss("continuous")
    #problem = FifteenBarTruss(prob_type='continuous')
    x0 = [var.ub for var in problem.vars]
    solver = SLP(move_limits=[0.25, 0.25], gamma=1e-3)
    fopt, xopt = solver.solve(problem, x0=x0, maxiter=100, log=False, verb=True, plot=False)
    problem(xopt)
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


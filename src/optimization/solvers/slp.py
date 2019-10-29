
from scipy.optimize import linprog
import numpy as np
import time

try:
    from src.optimization.solvers.optsolver import OptSolver
    from src.optimization.structopt import *
    from src.optimization.benchmarks import *
except:
    from optimization.solvers.optsolver import OptSolver
    from optimization.structopt import *
    from optimization.benchmarks import *
    

from ortools.linear_solver import pywraplp

class SLP(OptSolver):

    def __init__(self, move_limits=(0.05, 0.05), gamma=1e-2):
        super().__init__()
        self.move_limits = np.asarray(move_limits)
        self.gamma = gamma
        self.alpha = 1.0

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
            #lb = max(var.value - 0.5*self.alpha*(var.ub-var.lb),var.lb)
            #ub = min(var.value + 0.5*self.alpha*(var.ub-var.lb),var.ub)
            delta = var.ub - var.lb

            lb, ub = self.move_limits * delta
            lb = max(var.lb, var.value - lb)
            ub = min(var.ub, var.value + ub)
            x[i] = solver.NumVar(lb,
                                 ub,
                                 var.name)

        # Constraints
        for i in range(n):
            solver.Add(solver.Sum([A[i, j] * x[j] for j in range(m)]) <= B[i])

        # Objective
        solver.Minimize(solver.Sum([df[j] * x[j] for j in range(m)]))


        sol = solver.Solve()
        # If solution if infeasible
        if sol == 2:
            print("Solution found was infeasible!")
            self.move_limits += np.array([-self.gamma, self.gamma])
            # Return something != 0, otherwise iteration will
            # stop because two consecutive iterations produce
            # too similar results
            return np.ones_like(self.X)
        else:
            X = []
            for i in range(len(x)):
                if self.problem.prob_type == 'discrete':
                    X.append(int(x[i].solution_value()))
                else:
                    X.append(x[i].solution_value())
            X = np.asarray(X)

            return X - self.X
        # return res.x - self.X

    def step(self, action):

        self.move_limits -= self.move_limits * self.gamma
        self.alpha = self.alpha/(1+self.alpha)
        self.X += action
        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                                self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'

if __name__ == '__main__':

    problem = FifteenBarTruss(prob_type='continuous')
    x0 = [var.ub for var in problem.vars]
    solver = SLP(move_limits=[0.1, 0.1], gamma=1e-3)
    fopt, xopt = solver.solve(problem, x0=x0, maxiter=100, log=True, verb=True, plot=True)
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


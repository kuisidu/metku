
from src.optimization.solvers.optsolver import OptSolver
from src.optimization.structopt import *

from scipy.optimize import linprog
import numpy as np
import time

from src.optimization.benchmarks import *



class SLP(OptSolver):

    def __init__(self, move_limits=[0.75, 1.25], gamma=1e-2):
        super().__init__()
        self.move_limits = np.asarray(move_limits)
        self.gamma = gamma

    def take_action(self):
        """
        Defines action to take
        :return:
        """

        A, B, df, fx = self.problem.linearize(self.X)

        #bounds = [(x -self.step_length, x + self.step_length) for x in self.X]

        #res = linprog(df, A, B, bounds=bounds)


        from ortools.linear_solver import pywraplp

        solver = pywraplp.Solver('SLP',
                                 pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)

        # Number of constraints
        n = A.shape[0]
        # Number of variables
        m = A.shape[1]

        x = {}
        # Variables
        for i, var in enumerate(self.problem.vars):
            lb, ub = self.move_limits * var.value
            x[i] = solver.NumVar(lb,
                                 ub,
                                 var.name)

        # Constraints
        for i in range(n):
            solver.Add(solver.Sum([A[i, j] * x[j] for j in range(m)]) <= B[i])

        # Objective
        solver.Minimize(solver.Sum([df[j] * x[j] for j in range(m)]))

        sol = solver.Solve()
        X = []
        for i in range(len(x)):
            if self.problem.prob_type == 'discrete':
                X.append(int(x[i].solution_value()))
            else:
                X.append(x[i].solution_value())
        X = np.asarray(X)

        return X - self.X
        #return res.x - self.X




    def step(self, action):


        self.move_limits += (1 - self.move_limits) * self.gamma
        self.X += action
        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                           self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'

    def solve(self, problem, x0=None, maxiter=100, maxtime=60, log=False):
        """
        Solves given problem

        :param problem:
        :param maxiter:
        :param maxtime:
        :return:
        """
        self.problem = problem


        if x0:
            self.X = x0
        else:
            self.X = self.random_feasible_point()
        problem.substitute_variables(self.X)
        done = False
        for i in range(maxiter):

            if time.process_time() >= maxtime or done:
                break


            action = self.take_action()
            state, reward, done, info = self.step(action)
            self.X = state
            problem.substitute_variables(state)
            self.calc_constraints(self.X)
            print(state)

            if log:
                self.fvals.append(self.fval)


if __name__ == '__main__':

    problem = TenBarTruss(prob_type='continuous')
    solver = SLP(move_limits=[0.9, 1.1], gamma=2e-2)
    solver.solve(problem, maxiter=150, maxtime=300)
    problem(solver.X)
    print("Move limits: ", solver.move_limits)
    #problem.structure.plot_normal_force()

    # problem = OptimizationProblem("Quadratic Problem")
    # problem.obj = lambda x: x[0] ** 2 + x[1] ** 2
    # var1 = Variable("X1", 0, 5)
    # var2 = Variable("X2", 0, 5)
    # problem.add_variables([var1, var2])
    # con1 = NonLinearConstraint(lambda x: x[0] ** 2 / 20 - x[1] + 1)
    # con2 = NonLinearConstraint(lambda x: x[1] ** 2 / 20 - x[0] + 1)
    # problem.add_constraints([con1, con2])
    # solver = SLP(step_length=1)
    # solver.solve(problem, maxiter=30)


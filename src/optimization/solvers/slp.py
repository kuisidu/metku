
from src.optimization.solvers.optsolver import OptSolver
from src.optimization.structopt import *

from scipy.optimize import linprog
import numpy as np
import time

from src.optimization.benchmarks import *



class SLP(OptSolver):

    def __init__(self, step_length=100, step_factor=1e-3):
        super().__init__()
        self.step_length = step_length
        self.step_factor = step_factor

    def take_action(self):
        """
        Defines action to take
        :return:
        """

        A, B, df, fx = self.problem.linearize(self.X)

        bounds = [(x -self.step_length, x + self.step_length) for x in self.X]

        res = linprog(df, A, B, bounds=bounds)


        return res.x - self.X


    def random_feasible_point(self):
        """
        Creates a random feasible starting point for optimization

        :return: random feasible starting point
        """
        print("Creating random feasible starting point!")

        X = [0] * self.problem.nvars()
        for i, var in enumerate(self.problem.vars):
            X[i] = var.lb + (var.ub - var.lb) * np.random.rand()
        while np.any(self.calc_constraints(X) > 0):

            for i, var in enumerate(self.problem.vars):
                X[i] = var.lb + (var.ub - var.lb) * np.random.rand()

            self.problem.substitute_variables(X)
            self.problem.fea()

        print("Starting point created!")

        return np.asarray(X)

    def step(self, action):

        self.step_length -= self.step_length * self.step_factor
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
    solver = SLP(step_length=100)
    solver.solve(problem, maxiter=30, maxtime=300)
    problem(solver.X)
    problem.structure.plot_normal_force()

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


from optimization.solvers.optsolver import OptSolver
from optimization.structopt import LinearConstraint
from scipy.optimize import linprog
import numpy as np
import time

from optimization.benchmarks import *


def numeric_gradient(fun, h, x):
    """ Gradient by finite difference """
    fx = fun(x)
    n = len(x)
    fh = fun(x + h)
    df = (fh - fx) / h

    return df


def linearize(fun, x, grad=None):
    """ Linearization of a function

        input:
            fun .. function that takes x as input
            x .. point of linearization (array)
            grad .. returns the gradient of fun at x

        output:
            a .. grad(x)
            b .. fun(x)
            c .. grad(x)*x
    """

    b = fun(x)

    if grad == None:
        h = np.ones_like(x) * 1e-2
        a = numeric_gradient(fun, h, x)
    else:
        a = grad(x)

    c = a.dot(x)

    return a, b, c

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

        A = np.zeros((len(self.problem.cons), len(self.X)))
        B = np.zeros(len(self.problem.cons))
        # This needs to be defined more generally!
        C = np.zeros(len(self.X))
        for i, var in enumerate(self.problem.vars):
            mem = self.problem.structure.members[i]
            C[i] = mem.A * mem.length


        for i, con in enumerate(self.problem.cons):
            # TODO: Check if constraint is already linear
            print(self.X)
            a, b, c = linearize(con, self.X)
            print(self.X)
            A[i] = a
            B[i] = c - b


        bounds = [(x -self.step_length, x + self.step_length) for x in self.X]

        res = linprog(C, A, B, bounds=bounds, method='revised simplex')


        return res.x - self.X


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
            self.X = np.zeros(len(problem.vars))
            self.X += problem.vars[0].ub
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
            #print(state)

            if log:
                self.fvals.append(self.fval)


if __name__ == '__main__':

    problem = TenBarTruss(prob_type='continuous')
    solver = SLP()
    solver.solve(problem, maxiter=200)
    problem(solver.X)


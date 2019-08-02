from scipy.optimize import differential_evolution
from .optsolver import OptSolver
import numpy as np


class GA(OptSolver):

    def __init__(self, popsize=15, strategy='best1bin'):

        self.popsize = popsize
        self.strategy = strategy

    def obj(self, X, penalties=[]):
        X = [int(x) for x in X]
        if np.any(X != self.problem.X):
            self.problem.substitute_variables(X)
        obj_val = self.problem.obj(X)
        penalty = self.penalty(X)
        penalties.append(penalty)
        print(len(penalties))

        return obj_val + penalty


    def penalty(self, X):
        constr_vals = self.calc_constraints(X).clip(0, np.inf)
        penalty = 1e5 * np.sum(constr_vals ** 2)
        return penalty

    def solve(self, problem, maxiter=100):

        self.problem = problem
        bounds = self._create_bounds()




        result = differential_evolution(self.obj,
                                        bounds,
                                        popsize=self.popsize,
                                        strategy=self.strategy,
                                        maxiter=maxiter,
                                        updating='deferred')

        X = [int(x) for x in result.x]

        print(X, problem.obj(X))
        problem(X)
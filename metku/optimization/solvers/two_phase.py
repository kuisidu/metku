# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
import numpy as np


from metku.optimization.solvers import OptSolver
from metku.optimization.structopt import DiscreteVariable, Variable

"""
phase 1:
- continuous / discrete variables

phase 2:
- discrete / index variables

TwoPhase(solver1, solver2, limits=(1, 3), var_type1='', var_type2="", maxiter1=100, maxiter2=100)

fopt, xopt = solver.solve(problem, x0=x0, maxiter=100, verb=True, log=True)

"""

class TwoPhase:


    def __init__(self, solver1, solver2, limits=(1, 3),
                 var_type1="continuous", var_type2="discrete",
                 maxiter1=None, maxiter2=None):
        super().__init__()

        if not isinstance(solver1, OptSolver)\
                and isinstance(solver2, OptSolver):
            raise ValueError("Solvers must be 'OptSolver' class solvers!")

        self.solver1 = solver1
        self.solver2 = solver2
        self.limits = np.asarray(limits)
        self.var_type1 = var_type1
        self.var_type2 = var_type2
        self.maxiter1 = maxiter1
        self.maxiter2 = maxiter2

        self.problem = None

    def discretize(self):
        """
        Discretizes relaxed problem
        :return:
        """
        # Continuous solution
        best_x = self.solver1.X
        # Initialize problem as discrete problem
        self.problem.__init__(prob_type=self.var_type2)
        for i, var in enumerate(self.problem.vars):
            # if isinstance(var, DiscreteVariable):
            #     pass
            if type(var) is Variable:
                pass
            else:
                values = np.asarray(var.values)
                idx = np.argmin((values - best_x[i])**2)
                var.substitute(var.values[idx])
                ll, ul = self.limits
                lb = max(0, idx - ll)
                ub = min(len(values), idx + ul)
                var.__init__(var.name, var.values[lb:ub], var.target)
            
    def solve(self, problem, **kwargs):
        """
        Solves given problem
        :param problem:
        :param kwargs:
        :return:
        """
        self.problem = problem
        if problem.prob_type != self.var_type1:
            problem.__init__(prob_type=self.var_type1)
        print(f"Solving first phase using {type(self.solver1).__name__}")
        if self.maxiter1 is not None:
            kwargs['maxiter'] = self.maxiter1
        self.solver1.solve(problem, **kwargs)
        self.discretize()
        print(f"Discretizing first phase solution using {type(self.solver2).__name__}")
        x0 = [var.lb for var in self.problem.vars]
        kwargs['x0'] = x0
        if self.maxiter2 is not None:
            kwargs['maxiter'] = self.maxiter2
        xopt, fopt = self.solver2.solve(problem, **kwargs)
            
        return xopt, fopt


if __name__ == '__main__':
    from metku.optimization.solvers import *
    from metku.optimization.benchmarks import *
    from metku.optimization.problems import WIColumn
    import deap

    #solver1 = SLP(move_limits=[0.9, 1.5])
    solver1 = MISLP(move_limits=[0.85, 1.5])
    solver2 = GA(pop_size=100, mutate=deap.tools.mutShuffleIndexes,
                 mutation_kwargs={'indpb':0.15})

    #solver = TwoPhase(solver1, solver2, limits=[3, 3], var_type2='index', maxiter1=40, maxiter2=50)
    problem = WIColumn(prob_type='discrete')
    x0 = [var.ub for var in problem.vars]
    problem(x0)
    #fopt, xopt = solver1.solve(problem, x0=x0, maxiter=50, verb=True)
    #problem(xopt)
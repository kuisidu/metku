# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
import numpy as np

from metku.optimization.structopt import OptimizationProblem, Variable, \
    NonLinearConstraint, ObjectiveFunction


class ConstrainedRosenbrockProblem(OptimizationProblem):
    """
    https://en.wikipedia.org/wiki/Test_functions_for_optimization
    Rosenbrock function constrained with a cubic and a line
    """

    def __init__(self):
        super().__init__(name="Constrained Rosenbrock")

        self.prob_type = 'continuous'
        self.create_variables()
        self.create_constraints()
        self.create_objective()

    def create_variables(self):
        self.vars = []

        x1 = Variable(name="x1",
                      lb=-1.5,
                      ub=1.5)
        x2 = Variable(name="x2",
                      lb=-0.5,
                      ub=2.5)

        self.vars = [x1, x2]

    def create_constraints(self):
        self.cons = []

        def confun1(x):
            x1, x2 = x

            return (x1 - 1) ** 3 - x2 + 1

        def confun2(x):
            x1, x2 = x

            return x1 + x2 - 2

        con1 = NonLinearConstraint(con_fun=confun1,
                                   name="Constraint 1",
                                   con_type="<",
                                   problem=self)

        con2 = NonLinearConstraint(con_fun=confun2,
                                   name="Constraint 2",
                                   con_type="<",
                                   problem=self)

        self.cons = [con1, con2]

    def create_objective(self):
        def obj_fun(x):
            x1, x2 = x
            return (1 - x1) ** 2 + 100 * (x2 - x1 ** 2) ** 2

        obj = ObjectiveFunction(name="Objective ",
                                obj_fun=obj_fun)
        obj.problem = self

        self.obj = obj


class ConstrainedTownsendProblem(OptimizationProblem):
    """
    https: // en.wikipedia.org / wiki / Test_functions_for_optimization
    Townsend function (modified)
    """
    def __init__(self):
        super().__init__(name="Constrained Townsend")

        self.prob_type = 'continuous'
        self.create_variables()
        self.create_constraints()
        self.create_objective()

    def create_variables(self):
        self.vars = []

        x1 = Variable(name="x1",
                      lb=-2.25,
                      ub=2.5)
        x2 = Variable(name="x2",
                      lb=-2.5,
                      ub=1.75)

        self.vars = [x1, x2]

    def create_constraints(self):
        self.cons = []

        def confun(x):
            x1, x2 = x
            t = np.arctan2(x1, x2)

            b1 = (2 * np.cos(t) - 0.5 * np.cos(2 * t) - 0.25 * np.cos(
                3 * t) - 0.125 * np.cos(4 * t)) ** 2
            b2 = (2 * np.sin(t)) ** 2

            g = x1 ** 2 + x2 ** 2 - b1 - b2
            return g

        constr = NonLinearConstraint(con_fun=confun,
                                     name="Constraint",
                                     con_type='<',
                                     problem=self)

        self.cons = [constr]

    def create_objective(self):
        def obj_fun(x):
            x1, x2 = x

            val = -(np.cos((x1 - 0.1) * x2)) ** 2 - x1 * np.sin(3 * x1 + x2)
            return val

        obj = ObjectiveFunction(name="Objective ",
                                obj_fun=obj_fun)
        obj.problem = self

        self.obj = obj

if __name__ == '__main__':
    from metku.optimization.solvers import *

    problem = ConstrainedRosenbrockProblem()
    solver = SLP()
    x0 = [1,1]
    problem.substitute_variables(x0)
    problem.vars[0].lock()
    x0 = [var.ub for var in problem.vars]

    fopt, xopt = solver.solve(problem, x0=x0, maxiter=50)
    problem(xopt)

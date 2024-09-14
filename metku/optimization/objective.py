# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:21:39 2022

Objective function class for optimization problems

@author: kmela
"""

import numpy as np

from metku.optimization.optimzation_enums import ObjectiveTypeEnum


class ObjectiveFunction:

    def __init__(self,
                 name: str,
                 obj_fun: callable,
                 obj_type: ObjectiveTypeEnum = ObjectiveTypeEnum.MIN,
                 problem: 'OptimizationProblem' = None,
                 fea_required: bool = False):

        self.name = name
        if not callable(obj_fun):
            raise ValueError("obj_fun must be a function!")
        self.obj_fun = obj_fun
        self.obj_type = obj_type[:3].upper()
        self.problem = problem
        self.fea_required = fea_required

    def __call__(self, x):
        """ The method returns f(x) where x is an array with length
            of the number of free variables. The array will be appended
            for any fixed variables with corresponding fixed values.
            This way, if some variables are fixed, the original objective
            function can still be called.

            Read values from 'x' and reverse them.
            Reversing is done so that 'pop' method can be employed.
        """
        if self.fea_required and not self.problem.fea_done:
            self.problem.fea()

        X = [val for val in x]
        X.reverse()

        """ Fixed values values is an array
            with the values
            fixed_vals[i] is None if the variable vars[i] is free and
            fixed_vals[i] is val if the variable vars[i] has a fixed value val
        """
        fixed_vals = self.problem.fixed_vals.copy()

        for i, val in enumerate(fixed_vals):
            if val is None:
                """ If variable vars[i] is free, take its value from X
                    Because X was reversed, the pop method provides the correct
                    value here.
                """
                fixed_vals[i] = X.pop()
            if not X:
                """ If variable vars[i] is fixed, then use the corresponding
                    value from 'fixed_vals'
                """
                break

        self.problem.substitute_variables(fixed_vals)

        if self.obj_type == ObjectiveTypeEnum.MIN:
            return self.obj_fun(fixed_vals)
        else:
            return -self.obj_fun(fixed_vals)

    def neg_call(self, x):
        """
        Calculates negative value for objective function

        :return: negative value of objective function
        """
        return -self(x)


class LinearObjective(ObjectiveFunction):
    """ Class for linear objective functions """

    def __init__(self, name, c, obj_type="MIN", problem=None):
        """ Constructor:
            input:
                c .. constant array for performing the evaluation
                     c @ x
        """
        obj_fun = self.evaluate

        self.c = np.array(c)

        super().__init__(name, obj_fun, obj_type, problem)

    def evaluate(self, x):
        """ Evaluate function value """
        return self.c.dot(np.array(x))

# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:14:17 2022

@author: kmela
"""

import numpy as np
from collections.abc import Iterable

from metku.optimization.constants import X_TOL, POT_ACT_TOL
from metku.optimization.variables import Variable, IndexVariable
import metku.optimization.opt_functions as opt_funs


class Constraint:
    """ General class for constraints """

    def __init__(self,
                 name:str="",
                 con_type="<",
                 problem:'OptimizationProblem'=None,
                 fea_required:bool=False,
                 load_ids:list[int]=None,
                 vars:list[Variable]=None,
                 grad=None):
        """ Constructor
            Input:
                name: Name of the constraint (string)
                con_type: '<', '>' or '='
                problem: OptimizationProblem class object, to which the
                         constraint belongs
                fea_required: (Boolean) to indicate whether or not carry
                              out structural analysis, when the constraint
                              is evaluated
                vars: list of Variable type objects
                grad: callable function that returns the gradient of constraint function
        """
        self.name = name
        self.type = con_type
        self.problem = problem
        self.fea_required = fea_required
        if not isinstance(vars, Iterable):
            vars = [vars]
        self.vars = vars


        if not isinstance(load_ids, Iterable):
            load_ids = [load_ids]

        self.load_ids = load_ids

        # Lagrange multiplier
        self.mult = 1.0
        
        self.grad = grad
        self.potential = True # Flag for stating whether the constraint is potentially active or not
        self.pot_tol = POT_ACT_TOL

    def __call__(self, x):
        """ Evaluate the constraint using the syntax con(x)
            
        Runs finite element analysis on problem's structure if needed
        """
        if self.problem is not None:
            """ Collect current variable values stored in each
                optimization variable to X
            """
            
            X = np.array([var.value for var in self.problem.vars])            
            x = np.asarray(x)
            
            """ If the new values deviate from the old values
                in terms of L2 norm by more than X_TOL, substitute
                new values.

                NOTE: X_TOL should be small
            """
            if np.linalg.norm(X - x) > X_TOL:
                #print('variablet asetettiin')
                self.problem.substitute_variables(x)

        """ If structural analysis is part of the consraint evaluation
            AND
            if FEM analysis has not been carried out for X, perform FEM
        """
        if self.fea_required:
            for load_id in self.load_ids:
                if load_id not in self.problem.fea_done:
                    self.problem.fea(load_id)
                elif not self.problem.fea_done[load_id]:
                    self.problem.fea(load_id)

    def __repr__(self):
        return f"{type(self).__name__}: {self.name}"

    def neg_call(self, x):
        """ Evaluate -g(x) """
        return -self(x)


class LinearConstraint(Constraint):
    """ Class for linear constraints """

    def __init__(self,
                 a=None,
                 b=None,
                 con_type="<",
                 name="",
                 problem=None,
                 **kwargs):
        """ Constructor 
            Constraint is of the form: a'*x con_type b
            
            **kwargs: can be 'fea_required' and 'vars'.
            Probably fea_required is always False for linear constraints
        """

        self.a = np.array(a)
        self.b = b
        Constraint.__init__(self, name, con_type, problem, **kwargs)

    def __call__(self, x):
        """ Evaluate constraint at x 
        """
        
        """ This call substitutes variable values (if needed) and
            performs structural analysis (if needed)
        """
        super().__call__(x)

        X = []
        """ By iterating over 'all_vars', also the fixed variables
            are included in the function evaluation.
        """
        for var in self.problem.all_vars:
            if isinstance(var, IndexVariable):
                X.append(var.idx)
            else:
                X.append(var.value)

        if isinstance(self.b, Variable):
            b = self.b.value
        else:
            b = self.b

        gval = self.a.dot(np.array(X)) - b
        
        if self.type == '<':
            if gval < -self.pot_tol:
                self.potential = False
            else:
                self.potential = True
        elif self.type == '>':
            if gval > self.pot_tol:
                self.potential = False
            else:
                self.potential = True

        return gval
        


class NonLinearConstraint(Constraint):
    """ Class for linear constraints """

    def __init__(self, con_fun, con_type="<", name="", problem=None, **kwargs):
        """ Constructor
            con_fun: function that returns the value of constraint function
                     g(x) <= 0 (or >= 0 or = 0)
        """

        self.con = con_fun
        self.grad = None        
        Constraint.__init__(self, name, con_type, problem, **kwargs)

    def __call__(self, x):
        """ Evaluate constraint at x """
        
        """ This call substitutes variable values (if needed) and
            performs structural analysis (if needed)
        """        
        super().__call__(x)

        
        """ Make vector 'fixed_vals', corresponding to the vector of
            design variables that includes those in 'x' and the fixed
            values
        """
        X = [val for val in x]
        X.reverse()
        
        """ fixed_vals is a list with length equalling the
            number of original variable, including fixed and free.
            In place of free variables, 'fixed_vals' has None, and to this
            position, the value from 'X' should be inserted.
        """
        # IS COPY ENOUGH, OR SHOULD IT BE deepcopy?
        fixed_vals = self.problem.fixed_vals.copy() 
        for i, val in enumerate(fixed_vals):
            if val is None:
                fixed_vals[i] = X.pop()
            if not X:
                break

        gval = self.con(fixed_vals)    

        if self.type == '<':
            if gval < -self.pot_tol:
                self.potential = False
            else:
                self.potential = True
        elif self.type == '>':
            if gval > self.pot_tol:
                self.potential = False
            else:
                self.potential = True

        return gval
    
    def df(self,x):
        """ Evaluate gradient at x """
        
        if self.grad is None:
            # Evaluate gradient numerically
            df = opt_funs.NumGrad(self.con,x)
        else:
            df = self.grad(x)
        
        return df
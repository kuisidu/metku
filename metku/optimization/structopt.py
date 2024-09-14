# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Tue May 21 11:00:48 2019

Constructing and solving structural optimization problems

Idea: all the data is stored in OptimizationProblem object.
As different optimization solvers are applied, the methods of
OptimizationProblem class will be used to write the problem data in
a form suitable for each solver.

Classes:
    Variable
    IntegerVariable(Variable)
    BinaryVariable(IntegerVariable)
    IndexVariable(IntegerVariable)
    DiscreteVariable(Variable)


@author: kmela
"""

import copy
import time
from itertools import product

import numpy as np
import scipy.optimize as sciop

from metku.optimization.constants import XLB, XUB
from metku.optimization.constraints import Constraint, LinearConstraint, NonLinearConstraint
from metku.optimization.objective import ObjectiveFunction
from metku.optimization.opt_functions import NumGrad, conlin_fun, mma_fun
from metku.optimization.variables import BinaryVariable, IndexVariable
from metku.optimization.variables import Variable, DiscreteVariable, IntegerVariable


# from metku.eurocodes.en1991.loadIDs import LoadIDs

class OptimizationProblem:
    """ Class for defining and operating with optimization problems """

    def __init__(self, name="", variables=None, constraints=None, objective=None, \
                 gradient=None, hess=None, structure=None, profiles=None):
        """ Constructor

            Parameters:
            -----------
                :param name: name of the problem (string)
                :param variables: list of variables (Variable)
                :param constraints: list of constraints (Constraint)
                :param objective: objective function
                :param gradient: gradient of the objective function
                :param hess: Hessian of the objective function
                :param structure: structure to be optimized (optional).
                    If structural analysis is required in optimization, 'structure'
                    must have the method 'calculate' that performs structural analysis
                :param profiles: list of available profiles (optional)

                Variables:
                ----------
                :ivar name: name of the problem
                :ivar vars: variables (list)
                :ivar cons: constraints (list)
                :ivar obj: objective function
                :ivar grad: gradient of the objective function
                :ivar hess: Hessian of the objective function
                :ivar structure: structure to be optimized
                :ivar profiles: list of available profiles (optional)
                :ivar fea_done: indicates if FEM analysis has been done (bool)
                :ivar X: stored design variable values
                :ivar x0: initial point
                :ivar num_iters: number of iterations (int)
                :ivar num_fem_analyses: number of FEM analyses (int)
                :ivar fvals: number of function evaluations (int)
                :ivar states:
                :ivar gvals:
        """
        self.con_tol = 1e-4
        self.name = name
        if variables:
            self.__vars = []
            self.vars = variables
        else:
            self.__vars = []
        if constraints:
            self.cons = constraints
        else:
            self.cons = []
        self.obj = objective
        # Gradient and Hessian of the objective function
        # should also include the possibility that some of the variables
        # can be fixed!
        self.grad = gradient
        self.hess = hess
        self.structure = structure
        self.profiles = profiles
        self.fea_done = False
        self.X = None
        self.x0 = None
        self.num_iters = 0
        self.num_fem_analyses = 0
        self.fvals = []
        self.states = []
        self.gvals = []
    
    def xlb(self):
        """ Returns a numpy array of variable lower bounds """
        return np.array([x.lb for x in self.vars])

    def xub(self):
        """ Returns a numpy array of variable upper bounds """
        return np.array([x.ub for x in self.vars])

    def add(self, this):
        """
        Adds given object to the problem
        :param this: object to be added
        
        :type: Variable, Constraint, ObjectiveFunction
        """
        # VARIABLES
        if isinstance(this, Variable):
            if this not in self.__vars:
                self.__vars.append(this)

        # CONSTRAINTS
        elif isinstance(this, Constraint):
            if this not in self.cons:
            #if this not in self.cons and \
            #        (isinstance(this, NonLinearConstraint) and this.con not in [con.con for con in self.cons]):
                self.cons.append(this)
                this.problem = self

        # OBJECTIVES
        elif isinstance(this, ObjectiveFunction):
            # TODO: Multiple objectives
            self.obj = this
            this.problem = self
        else:
            raise ValueError(f"{this} must be either Variable, Constraint or ObjectiveFunction")

    @property
    def locked_vars(self):
        """
        Returns locked variables

        :return: list of locked variables (Variable objects, not values!)
        """
        return [var for var in self.__vars if var.locked]

    @property
    def all_vars(self):
        """
        Returns all variables

        :return: list of all variables
        """
        return self.__vars

    @property
    def vars(self):
        """
        Returns all unlocked variables

        :return: list of unlocked variables
        """
        return [var for var in self.__vars if not var.locked]

    @vars.setter
    def vars(self, vals):
        """
            Sets variables, to 'vals' (list of Variable objects)            
        """
        self.clear_vars()
        for val in vals:
            self.add(val)

    @property
    def fixed_vals(self):
        vals = []
        for var in self.__vars:
            if var.locked:
                vals.append(var.value)
            else:
                vals.append(None)
        return vals

    #@property
    def feasible(self,x=None):
        """
        Returns problem's feasibility at x
        :return: True/False
        """
        
        res = True

        if x is None:
            x = self.X

           
        for con in self.cons:
            g = con(x)             
            if con.type == '<':
                if g > self.con_tol:
                    #print('<= type inequality constraint infeasible')
                    #print(con,g)
                    res = False
                    break
            elif con.type == '>':
                if g < -self.con_tol:
                    #print('>= type inequality constraint infeasible')
                    #print(con,g)
                    res = False
                    break
            else:
                if abs(g) > self.con_tol:
                    #print('== type inequality constraint infeasible')
                    #print(con,g)
                    res = False
                    break
        
        return res
        #return np.all(self.eval_cons(self.X) <= self.con_tol)
    
    

    def type_of_vars(self):
        '''
        :return: "Variable", "IntegerVariable", "BinaryVariable",
        IndexVariable or "DiscreteVariable"
        '''
        var = self.vars[0]
        if isinstance(var, BinaryVariable):
            return "BinaryVariable"
        elif isinstance(var, IndexVariable):
            return "IndexVariable"
        elif isinstance(var, DiscreteVariable):
            return "DiscreteVariable"
        elif isinstance(var, IntegerVariable):
            return "IntegerVariable"
        elif isinstance(var, Variable):
            return "Variable"
        else:
            raise Exception("Error, type of type_of_vars is none of above")

    def clear_vars(self):
        """
            Deletes all variables
        """
        self.__vars.clear()
        
    def del_var(self,var):
        """
            Deletes variable 'var'

        """
        self.__vars.remove(var)
        
    def var_values(self):
        """
        Returns all current variable values
        -------

        """
        return np.array([var.value for var in self.vars])
    
    def random_values(self):
        """
        Generates random variable values within the range of allowable values
        --------
        
        Returns
        -------
        None.

        """
        xrand = np.zeros(self.nvars())
        rng = np.random.default_rng()

        for i, var in enumerate(self.vars):
            if isinstance(var,IntegerVariable) or isinstance(var,IndexVariable):
                xrand[i] = rng.integers(var.lb,var.ub+1)
            elif isinstance(var, DiscreteVariable):
                possible_values = var.values
                random_value = (var.ub - var.lb) * rng.random() + var.lb
                xrand[i] = min(possible_values, key=lambda x:abs(x-random_value))
            else:
                xrand[i] = (var.ub-var.lb)*rng.random() + var.lb
        return xrand

    def print_vars(self):
        """ Prints the current variable values """
        
        print("Current variables:")
        for var in self.vars:
            print("{0:s} = {1:5.4f}".format(var.name,var.value))

    def non_linear_constraint(self, *args, **kwargs):
        """
        Creates new NonLinearConstraint and adds it to cons -list

        :param args:
        :param kwargs:
        :return: NonLinearConstraint
        """
        con = NonLinearConstraint(*args, **kwargs)
        con.problem = self

        self.cons.append(con)
        return con

    def linear_constraint(self, *args, **kwargs):
        """
        Creates new LinearConstraint and adds it to cons -list
        :param args:
        :param kwargs:
        :return: LinearConstraint
        """

        con = LinearConstraint(*args, **kwargs)
        con.problem = self
        self.cons.append(con)
        return con

    #@lru_cache(CACHE_BOUND)
    def linearize(self, x, only_potential=False):
        """
        Linearizes problem around x

        Ax < B

        :param x:
        :param only_potential: True, if only potentially active constraints are included
        in the linearization
        :return:
        """
        #start = time.time()
        x = np.asarray(x)

        self.substitute_variables(x)

        """ Evaluate objective function at x """
        fx = self.obj(x)
        """ Evaluate nonlinear constraint functions at x """
        #b = self.eval_nonlin_cons(x)
        
        con_vals = self.eval_cons(x)
                
        """ Get number of nonlinear constraints """
        if only_potential:
            act_cons = [con for con in self.cons if con.potential]
            b = np.array([val for val, con in zip(con_vals,self.cons) if con.potential])
        else:
            act_cons = self.cons
            b = con_vals
            
        m = len(act_cons)
        #m = self.nnonlincons()
        """ Number of variables """
        n = self.nvars()
        """ Jacobian
            Rows of A are the gradients (transposed) of constraints
        """
        A = np.zeros((m, n))
        df = np.zeros(n)
        
        if self.grad is not None:
            df = self.grad(x)
                
        for i, var in enumerate(self.vars):
            if isinstance(var, BinaryVariable):
                h = 1
            elif isinstance(var, IndexVariable):
                h = 1
                prev_val = var.value
                if var.idx + h > var.ub:
                    h = -1
                var.substitute(var.idx + h)
            else:
                x[i] = var.value
                """ Get current value of variable """
                prev_val = var.value
                """ Step length """
                    
                h = max(1e-6 * abs(prev_val), 1e-8)
                var.substitute(prev_val + h)
                
            """ Get variable values """
            xh = [var.value for var in self.vars]
            #print(x,xh)
            if not isinstance(var,BinaryVariable) and np.linalg.norm(x-xh) == 0:
                print("Linearize: error with variable substitution - too small step size.")

            """ Evaluate objective function at x + hi*ei """
            if self.grad is None:
                if not isinstance(var, BinaryVariable):                 
                    f_val = self.obj(xh)
                    #print(f_val,fx)         
                    df[i] = (f_val - fx) / h
                    
            """ Evaluate constraint functions at x + hi*ei """
            #a = np.array([con(xh) for con in act_cons])
            dg = np.zeros(m)
            for j, con in enumerate(act_cons):
                if isinstance(con,LinearConstraint):
                    dg[j] = con.a[i]
                else:
                    if not isinstance(var,BinaryVariable):
                        dg[j] = (con(xh)-b[j])/h
                #a = self.eval_nonlin_cons(xh)
                #A[:, i] = (a - b) / h
            A[:, i] = dg
                
            """ Substitute the original value x[i] to current variable """
            var.substitute(prev_val)            
            

        #B = A @ x.T - b
        B = A.dot(x) - b
        # B = -b

        # Add linear constraints' values
        """
        for con in self.cons:
            if isinstance(con, LinearConstraint):
                A = np.vstack((A, con.a))
                if isinstance(con.b, Variable):
                    B = np.hstack((B, con.b.value))
                else:
                    B = np.hstack((B, con.b))
        """
        #end = time.time()
        # print("Linearization took: ", end - start, " s")

        return A, B, df, fx
    
    def conlin(self,x):
        """ Generates convex-linear approximation of the problem
            around point 'x'
            
            Assumption: all variables are positive
            
            Note: Linear function remain linear
        """
        
        Pconlin = OptimizationProblem(name="conlin")
        
        for i, var in enumerate(self.vars):
            """ Create variables for the conlin problem """
            Pconlin.add(Variable(name=var.name,lb=var.lb,ub=var.ub,value=x[i]))
        
        # Approximate objective
        if self.grad is None:
            dfx = NumGrad(self.obj,x)
        else:
            dfx = self.grad(x)
            
        Pconlin.obj = conlin_fun(x,self.obj(x),dfx)
        
        # Approximate constraints
        for con in self.cons:
            if isinstance(con,LinearConstraint):
                # Linear constraints are simply added
                Pconlin.add(con)
            else:
                # For nonlinear constraints, CONLIN approximation is made
                conlin_con = copy.copy(con)
                gx = con.con(x)
                dgx = con.df(x)
                conlin_con.con = conlin_fun(x,gx,dgx)
                
                Pconlin.add(conlin_con)
        
        return Pconlin
    
    def mma(self,x,U,L,mu=0.8):
        """
        Generates MMA approximation of the problem at point 'x' using the asymptotes 'U' and 'L'


        Parameters
        ----------
        x : numpy array
            Approximation point.
        U : numpy array
            Upper asymptote.
        L : numpy array
            Lower asymptote.

        Returns
        -------
        Optimization problem with nonlinear function replaced by their MMA approximations.

        """
        
        Pmma = OptimizationProblem(name="MMA")
        
        for i, var in enumerate(self.vars):
            """ Create variables for the mma problem """
            Pmma.add(Variable(name=var.name,lb=max(var.lb,L[i]+mu*(x[i]-L[i])),ub=min(var.ub,U[i]-mu*(U[i]-x[i])),value=x[i]))
        
        # Approximate objective
        if self.grad is None:
            dfx = NumGrad(self.obj,x)
        else:
            dfx = self.grad(x)
            
        Pmma.obj = mma_fun(x,self.obj(x),dfx,L,U)
        
        # Approximate constraints
        for con in self.cons:
            if isinstance(con,LinearConstraint):
                # Linear constraints are simply added
                Pmma.add(con)
            else:
                # For nonlinear constraints, MMA approximation is made
                mma_con = copy.copy(con)
                gx = con.con(x)
                dgx = con.df(x)
                mma_con.con = mma_fun(x,gx,dgx,L,U)
                
                Pmma.add(mma_con)
        
        return Pmma
        
    def nnonlincons(self):
        """ Returns number of non linear constraints"""
        non_lin_cons = [con for con in self.cons if
                        isinstance(con, NonLinearConstraint)]
        return len(non_lin_cons)

    def eval_nonlin_cons(self, X):
        """
        Calculates non linear constraints and returns their values as a numpy array
        :param X:
        :return:
        """
        g = []
        non_lin_cons = [con for con in self.cons if
                        isinstance(con, NonLinearConstraint)]
        start = time.time()
        for con in non_lin_cons:
            g.append(con(X))

        return np.asarray(g)

    #def fea(self, load_id=LoadIDs.ULS):
    def fea(self, load_id=None):
        """
        Runs finite element analysis on structure
        """
        if load_id is None:
            load_id = self.structure.load_ids[0]
        
        self.structure.calculate(load_id)
        self.fea_done = True
        self.num_fem_analyses += 1

    def __call__(self, x=None, prec=2, ncons=5):
        """ Call method evaluates the objective function and all constraints
            at x and returns their values
        """
        if x is None:
            x = [var.value for var in self.vars]
        self.substitute_variables(x)
        fx = self.obj(x)
        print("** {0} **".format(self.name))
        print(f'Solution is feasible: {self.feasible()}')

        #vals = [var.value for var in self.vars]
        #print(f"Optimal values: {vals}")

        print(f'X: {[round(val, prec) for val in x]}')

        if len(self.cons):
            g = self.eval_cons(x)
            print(F"Max constraint: {self.cons[np.argmax(g)].name}: {max(g):.{prec}f}")
        if len(x) < 100:
            print("Variables:")
            print("----------")
            for var in self.vars:
                if isinstance(var, IndexVariable):
                    print(f"{var.name} = {var.values[var.value]}")
                else:
                    print(f"{var.name} = {var.value:.{prec}f}")

            print("----------\n")

        print(f"Objective function = {fx:.{prec}f}\n")

        if len(self.cons):
            ncons = min(len(self.cons), ncons)
            print(f"{ncons} Maximum Constraints:")
            print("----------")

            idx = np.argpartition(g, -ncons)[-ncons:]
            for i in idx:
                gi = g[i]
                con = self.cons[i]
                print(f"{con.name}: {gi:.{prec}f} {con.type} 0")
            # for con in np.asarray(self.cons)[idx]:
            #     gi = con(x)
            #     print(f"{con.name}: {gi:.{prec}f} {con.type} 0")


    def eval_cons(self, x):
        """ Constraint evaluation

        """
        g = []

        for con in self.cons:
            g.append(con(x))

        return np.asarray(g)

    def eval_eq_con(self, x):
        """ Evaluates equality constraints """

        geq = []

        for con in self.cons:
            if con.type == "=":
                geq.append(con(x))

    def eval_ineq_con(self, x):
        """ Evaluate inequality constriants
            All constraints are evaluated in the form g(x) >= 0
        """

        g = []

        for con in self.cons:
            if con.type == ">":
                g.append(con(x))
            elif con.type == "<":
                g.append(-con(x))

    def nvars(self):
        """ Return number of variables """
        return len(self.vars)

    def nlincons(self):
        """ Number of nonlinear inequality constraints """

        return len(self.vars)

    def nlineqcons(self):
        """ Number of nonlinear equality constraints """

        return len(self.vars)

    def add_variable(self, var):
        """ Adds new variable to the problem """
        self.vars.append(var)

    def add_variables(self, variables):
        """ Adds new variables to the problem """
        for var in variables:
            self.vars.append(var)

    def add_constraint(self, con):
        """ Adds new constraint to the problem """
        self.cons.append(con)

    def add_constraints(self, constraints):
        """ Adds new constraints to the problem """
        for con in constraints:
            self.cons.append(con)

    def add_structure(self, structure):
        """ Add structure to the problem """
        self.structure = structure

    def var_bounds(self):
        """ Returns arrays lb and ub for variable bounds """
        lb = []
        ub = []

        for v in self.vars:
            lb.append(v.lb)
            ub.append(v.ub)

        return lb, ub

    # THIS METHOD IS NOT NEEDED!! USE substitute_variables INSTEAD!!
    def substitute_variable(self, i, xval):
        """ Substitutes the value 'xval' for variable i """
        var = self.vars[i]
        var.value = xval

        if var.target is not None:
            for obj in var.target['objects']:
                if var.target['property'] == 'A':
                    obj.A = xval
                elif var.target['property'] == 'IY':
                    obj.I[0] = xval
                elif var.target['property'] == 'IZ':
                    obj.I[1] = xval
                elif var.target['property'] == 'WELY':
                    obj.Wel[0] = xval
                elif var.target['property'] == 'WELZ':
                    obj.Wel[1] = xval
                elif var.target['property'] == 'WPLY':
                    obj.Wpl[0] = xval
                elif var.target['property'] == 'WPLZ':
                    obj.Wpl[1] = xval
                elif var.target['property'] == 'H':
                    obj.h = xval
                elif var.target['property'] == 'TW':
                    obj.tw = xval
                elif var.target['property'] == 'BF':
                    obj.b = xval
                elif var.target['property'] == 'TF':
                    obj.tf = xval
                elif var.target['property'] == 'PROFILE':
                    obj.profile = self.profiles(xval)

    def substitute_variables(self, xvals):
        """ Substitute variable values from xval to structure """

        """ Make sure all variable values are within lower and upper bounds """
        xvals = np.array([max(var.lb, min(x, var.ub)) for x, var in zip(xvals, self.vars)])
        #print(xvals)
        #print(self.var_values())
        #if np.linalg.norm(self.var_values()-xvals) > 1e-9:
        # Save starting point
        if not np.any(self.X):
            self.x0 = xvals.copy()

        # if np.any(self.X != xvals):
        self.X = xvals.copy()
        self.fea_done = False
        for x, var in zip(xvals, self.vars):
            var.substitute(x)
        #else:
            #print("No variable substitution.")
            #if self.nvars() < 10:
            #    print(self.var_values(),xvals)
            #
            # for i in range(len(xvals)):
            #     self.substitute_variable(i, xvals[i])

    # OBSOLETE! USE SOLVERS FROM solvers MODULE!
    def solve(self, solver="slsqp", **kwargs):
        """ Solve the optimization problem
            input:
                solver .. name of the solver
                **kwargs .. list of keyword argument
        """
        lb, ub = self.var_bounds()

        if solver == "slsqp":
            jac = "2-point"
            bnds = sciop.Bounds(lb, ub)

            if self.grad is not None:
                jac = self.grad

            if self.hess is not None:
                hess = self.hess

            """ Crerate constraints
            Constraints for COBYLA, SLSQP are defined as a 
            list of dictionaries. Each dictionary with fields:

                type: str

                Constraint type: ‘eq’ for equality, ‘ineq’ for inequality.

                fun: callable The function defining the constraint.

                jac: callable, optional. The Jacobian of fun (only for SLSQP).

                args: sequence, optional

                    Extra arguments to be passed to the function and Jacobian.

            """
            # cons = [{"type":"ineq", "fun":self.eval_ineq_con}]

            cons = []
            for con in self.cons:
                if con.type == "=":
                    new_con = {"type": "eq", "fun": con}
                elif con.type == "<":
                    new_con = {"type": "ineq", "fun": con.neg_call}
                else:
                    new_con = {"type": "ineq", "fun": con}

                cons.append(new_con)

            # print(cons)

            res = sciop.minimize(self.obj, kwargs["x0"], method=solver, \
                                 jac=jac, bounds=bnds, constraints=cons)
        elif solver == "trust-constr":
            jac = "2-point"
            hess = "2-point"
            bnds = sciop.Bounds(lb, ub)

            if self.grad is not None:
                jac = self.grad

            if self.hess is not None:
                hess = self.hess

            res = sciop.minimize(self.obj, kwargs["x0"], method=solver, \
                                 jac=jac, hess=hess, bounds=bnds)
        else:
            print("Solver unidentified.")
            res = None

        return res

    def discrete_neighborhood(self, x, k, method='nearest'):
        """ Creates a discrete neighborhood around point 'x' 
            input:
                k -- parameter used by the 'method'
                method .. procedure for defining the neigborhood
        """
        neighborhood = []
        
        for i, var in enumerate(self.vars):
            if method == 'nearest':
                if isinstance(var,DiscreteVariable):
                    values = np.array(var.values)
                    """ Sort values in ascending order of distance from 'x' 
                        np.argsort returns the indices of 'values'
                    """
                    sorted_ndx = np.argsort(abs(values-x[i]))
                    """ Get k actual values closest to x and transform
                        the Numpy array to list.
                        
                        Update the list of values for this variable
                    """
                    neighborhood.append(values[sorted_ndx[:k]].tolist())
                elif isinstance(var,IntegerVariable):
                    """ For integer variables, make a list of
                        integers between lb and ub
                        
                        The sorting is the same as for discrete variables,
                        but here only lb and ub need to be updated.
                    """
                    values = np.arange(var.lb,var.ub+1)
                    sorted_ndx = np.argsort(abs(values-x[i]))
                    new_values = values[sorted_ndx[:k]].tolist()
                    neighborhood.append(new_values)
                    #var.lb = min(new_values)
                    #var.ub = max(new_values)
                
        
        """
            Create iterable nList that contains the neighborhood points
            as a list of tuples
        """
        nList = product(*neighborhood)
        return nList
    
    def exact_penalty_fun(self,R=1.0):
        """ Exact penalty function 
            NOTE: Returns only the sum of constraint violations!
        """
        
        def penalty_fun(x):
            
            F = 0.0
            
            for con in self.cons:
                val = con(x)
                if con.type == '=':
                    """ For equality constraints, the penalty term
                        is simply the value of the constraint squared
                    """
                    F += abs(val)
                elif con.type == '<':
                    """ For less than or equal constraints, the penalty
                        term is max(0,val)
                    """
                    F += R*max(0,val)
                elif con.type == '>':
                    """ For geq type constraints, the penalty term is
                         min(0,val)
                    """
                    F += R*min(0,val)
            
            """ Variable bounds are treated separately """
            for (i,var) in enumerate(self.vars):                
                if var.lb > XLB:                    
                    #F += R*max(0,-1-x[i]/abs(var.lb))**2
                    F += R*max(0,-x[i]+var.lb)
                
                if var.ub < XUB:
                    #F += R*max(0,x[i]/abs(var.ub)-1)**2
                    F += R*max(0,x[i]-var.ub)
            
            return F
        
        return penalty_fun
            

    def exterior_penalty_fun(self,R=10,p=2):
        """ Return exterior quadratic penalty function """
        
        def penalty_fun(x):
            
            F = self.obj(x)
            
            for con in self.cons:
                val = con(x)
                if con.type == '=':
                    """ For equality constraints, the penalty term
                        is simply the value of the constraint squared
                    """
                    F += val**p
                elif con.type == '<':
                    """ For less than or equal constraints, the penalty
                        term is square of max(0,val)
                    """
                    F += R*max(0,val)**p
                elif con.type == '>':
                    """ For geq type constraints, the penalty term is
                        square of min(0,val)
                    """
                    F += R*min(0,val)**2
            
            
            for (i,var) in enumerate(self.vars):                
                if var.lb > XLB:                    
                    #F += R*max(0,-1-x[i]/abs(var.lb))**2
                    F += R*max(0,-x[i]+var.lb)**2
                
                if var.ub < XUB:
                    #F += R*max(0,x[i]/abs(var.ub)-1)**2
                    F += R*max(0,x[i]-var.ub)**2
            
            return F
        
        return penalty_fun

    def barrier_fun(self,R,barrier_type='log'):
        """ Returns the barrier function
            barrier_type: 'log' or 'inverse'
        """
        
        def bar_fun(x):
            
            F = self.obj(x)
            
            for con in self.cons:
                val = con(x)
                if con.type == '=':
                    """ For equality constraints, the barrier term
                        is not defined
                    """
                    
                elif con.type == '<':
                    """ For less than or equal constraints, the penalty
                        term is -1/R*log(-val)
                    """
                    if barrier_type == 'log':            
                        F -= 1/R*np.log(-val)
                    elif barrier_type == 'inverse':
                        F -= 1/R/val
                elif con.type == '>':
                    """ For geq type constraints, the penalty term is
                        square of min(0,val)
                    """
                    if barrier_type == 'log':            
                        F -= 1/R*np.log(val)
                    elif barrier_type == 'inverse':
                        F += 1/R/val
            
            for (i,var) in enumerate(self.vars):                
                if var.lb > XLB:                    
                    #F += R*max(0,-1-x[i]/abs(var.lb))**2
                    F -= 1/R/(-x[i]+var.lb)
                
                if var.ub < XUB:
                    #F += R*max(0,x[i]/abs(var.ub)-1)**2
                    F -= 1/R/(x[i]-var.ub)        
            
            return F
        
            
        
        return bar_fun

    def augmented_lagrangian(self,Rineq=1.0,Req=1.0):
        """ Constructs the augmented Lagrangian function for the
            problem.
            
            mult .. list of Lagrange multipliers
            R .. penalty factors
            
            NOTE! variable bounds are treated as regular inequalities
        """
        
        def aul_fun(x):
        
            F = self.obj(x)
            
            for con in self.cons:
                val = con(x)
                mult = con.mult
                
                if con.type == '=':
                    """ Equality constraints
                    """
                    F += -mult*val + Req*val**2
                    
                elif con.type == '<':
                    """ Less than or equal constraints
                    """                    
                    F += Rineq*max(0.5*mult/Rineq+val,0)**2
                    
                elif con.type == '>':
                    """ Less than or equal constraints
                        NOTE! CHECK THIS!!!!
                    """                    
                    F += Rineq*max(0.5*mult/Rineq-val,0)**2
          
            for (i,var) in enumerate(self.vars):                     
                if var.lb > XLB:
                    val = -x[i] + var.lb                            
                    F += Rineq*max(0.5*var.lb_mult/Rineq + val,0)**2
                
                if var.ub < XUB:                    
                    val = x[i] - var.ub
                    F += Rineq*max(0.5*var.ub_mult/Rineq + val,0)**2
          
            return F
        
        return aul_fun
    
    def update_lag_mult(self,x,Rineq=1.0,Req=1.0):
        """ Updates Lagrange multipliers according to the update rules
            associated with the augmented Lagrangian method.
        """
        
        """ First, iterate over regular constraints """
        for con in self.cons:
            val = con(x)
            mult = con.mult
            
            if con.type == '=':
                """ Equality constraints
                """
                con.mult -= 2*Req*val
                    
            elif con.type == '<':
                """ Less than or equal constraints
                """                           
                con.mult = max(mult + 2*Rineq*val,0)
                    
            elif con.type == '>':
                """ Less than or equal constraints
                NOTE! CHECK THIS!!!!
                """                    
                con.mult =  max(mult - 2*Rineq*val,0)
         
        """ Iterate over variable bound constraints """
        for (i,var) in enumerate(self.vars):                     
            if var.lb > XLB:
                val = -x[i] + var.lb                            
                var.lb_mult = max(var.lb_mult + 2*Rineq*val,0)
                                
            if var.ub < XUB:                    
                val = x[i] - var.ub
                var.ub_mult = max(var.lb_mult + 2*Rineq*val,0)
            




if __name__ == '__main__':

    b = BinaryVariable(name="oma",section=1,target=3)    

    """    
    prop = OptimizationProblem(name="Testi")

    def obj(x):
        a, b = x
        return a ** 2 + b + 10

    objective = ObjectiveFunction(name="Testi", obj_fun=obj)

    var1 = Variable("Var1", lb=0, ub=100)
    var2 = Variable("Var2", lb=0, ub=100)

    prop.add(objective)
    prop.add(var1)
    prop.add(var2)

    print(prop([5, 5]))

    # dvars = []
    # dvars.append(Variable("Height", 150, 1000))
    # dvars.append(Variable("Flange width", 150, 500))
    # dvars.append(Variable("Flange thickness", 5, 40))
    # dvars.append(Variable("Web thickness", 5, 40))
    #
    """
    """    
    def obj_fun(x):
        return x[0]**2 + 2*x[1]**3 -4*x[2]
    
    x1 = DiscreteVariable(name='x1',values = [4.5, 5.4, 7.23, 11.2, 13.4])
    x2 = IntegerVariable(name='x2',lb=1,ub =10)
    x3 = DiscreteVariable(name='x3',values = [3.77,6.987,11.223])
    x = [x1,x2,x3]
    x0 = [6.66,4.445,5.0]
    
    p = OptimizationProblem(name="I-Beam Weight Minimization",variables=x)
    
    p.obj = obj_fun
    
    N = p.discrete_neighborhood(x0,k=4)
    
    fvals = []
    for xi in N:
        print(xi)
        fvals.append((p.obj(xi)))
    # create grid of neighborhood points
    #M = np.meshgrid(*N)
    
    #s = M[0][0].shape
    """

    """
    for a in np.nditer(M):
        print(len(a),list(a))
    """
    
    """
    for i in range(s[0]):
        for j in range(s[1]):
            for k in 
            for I,m in enumerate(M):
                x[I] = m(k,i,j)
    """         
    
    # p = OptimizationProblem(name="I-Beam Weight Minimization",
    # variables=dvars)
    """
    problem = OptimizationProblem()
    vars = []
    for i in range(10):
        var = Variable("Name " + str(i),
                       lb=0,
                       ub=1)
        vars.append(var)

    problem.vars = vars
    var0 = problem.vars[0]
    var0.lock()
    var3 = problem.vars[3]
    var3.lock(0.2)
    print(len(problem.vars))
    var0.unlock()
    var3.unlock()
    print(len(problem.vars))
    """
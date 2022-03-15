# -*- coding: utf-8 -*-
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

import time
import math

import numpy as np
import scipy.optimize as sciop
from itertools import product
from collections.abc import Iterable
from functools import lru_cache
from copy import deepcopy
import copy
try:
    from metku.loadIDs import LoadIDs
except:
    from loadIDs import LoadIDs

CACHE_BOUND = 2**8
INT_TOL = 1e-4
""" Tolerance for discreteness violation of a discrete variable """
DISC_TOL = 1e-3
""" Tolerance for comparing variable vectors """
X_TOL = 1e-10

""" Artificial variable upper and lower bounds """
XUB = 1e8
XLB = -1e8


class Variable:
    """ Class for optimization variables
    """

    def __init__(self, name, lb=XLB, ub=XUB, id=None, target=None, profiles=None, value=None, target_fun=None):
        """
        Parameters:
            -----------
                :param name: string stating the name of the variable
                :param value: current value of the variable
                :param lb: lower bound
                :param ub: upper bound
                :param ub: upper bound
                :param target: dict that provides information about the part of
                the structure that the variable is changing
                :param profiles: list of available profiles (optional)

                target = {"property": string, "objects": list}
                    property .. unique identifier that shows which property of
                                the structure is affected by the variable
                                
                                this can also be a list of strings
                    objects .. list of objects (e.g. frame members) that are
                                affected by the variable

                if target is 'None', then the variable does not change anything
                in the structure
                :param target_fun: (optional) dict that provides functions for evaluating a functional relationship
                                    between the variable and its target. For cases, where the variable substitution
                                    is not a direct attribute insertion. For example, when the cross-sectional area
                                    of a member is used as a variable and it also affects the second moment of area
                                    by some function I = I(A)

                Allowed values for 'property':
                    'AREA' .. cross-sectional area of a member
                    'IY' .. second moment of area of a member (major axis)
                    'IZ' .. second moment of area of a member (minor axis)
                    'WELY' .. elastic section modulus of a member (major axis)
                    'WELZ' .. elastic section modulus of a member (minor axis)
                    'WPLY' .. plastic section modulus of a member (major axis)
                    'WPLZ' .. plastic section modulus of a member (minor axis)
                    'PROFILE' .. choosing a profile from a given catalogue
                    'BF' .. width of flange
                    'TF' .. thickness of flange
                    'H' .. height of the profile
                    'TW' .. thickness of the web
        """
        
        self.name = name
        self.value = value
        self.lb = lb
        self.ub = ub
        
        self.target = target
        self.profiles = profiles
        self.locked = False
        self.branch_priority = 0
        self.id = id            # MIKÄ TÄMÄ ON?
        self.scaling = 1.0
        self.target_fun = target_fun

        # Lagrange multiplier for lower bound
        self.lb_mult = 1.0
        # Lagrange multiplier for upper bound
        self.ub_mult = 1.0

    def __repr__(self):

        #return self.name + ": [" + str(self.lb) + "," + str(self.ub) + "]"
        if self.locked:
            fixed = " (fixed to {0:g})".format(self.value)
        else:
            fixed = ""

        return str(self.lb) + " <= " + self.name + " <= " + str(self.ub) + fixed

    def lock(self, val=None):
        """
        Locks the variable so it's value can't be changed

        Parameters:
        -----------
        :param val: value to give to variable before locking
        """
        if val is not None:
            self.substitute(val)

        self.locked = True

    def unlock(self):
        """
        Unlocks the variable so it's value can be changed
        """
        self.locked = False
        
    def scale(self,scaling_factor):
        """ Scale variable by the factor 'scaling_factor' 
            This alters the variable bounds. In 'substitute' method, scaling
            is reflected such that the variable value is multiplied by
            'scaling_factor'
        """
        
        self.scaling = scaling_factor
        self.lb = self.lb/scaling_factor
        self.ub = self.ub/scaling_factor

    def substitute(self, new_value):
        """
        Changes variable value and modifies target
        :param new_value:
        :return:
        """
        if not self.locked:
            """ Substitute a new value for the variable """
            self.value = new_value
            """ Modify target object(s) if any """            
            if not self.target is None:
                for obj in self.target['objects']:
                    prop = self.target['property']
                    """ If the variable affects several properties at a time
                        'prop' is a list of property names
                    """
                    if isinstance(prop, list):
                        for p in prop:
                            """ If the varible has a 'target_fun',
                                evaluate it for right properties
                            """
                            if self.target_fun is not None:
                                """ If 'p' is not a property that has a functional dependency,
                                    then do regular substitution
                                """
                                try:
                                    new_val = self.target_fun[p](new_value*self.scaling)
                                except KeyError:
                                    new_val = new_value*self.scaling
                                    
                                obj.__setattr__(p, new_val)
                            else:
                                obj.__setattr__(p, new_value)

                    elif isinstance(self, IndexVariable):
                        obj.__setattr__(prop, new_value)
                    elif isinstance(self, BinaryVariable):
                        """ If a binary variable has a target, assume that target
                            is the section to be substituted. The section is
                            given in the 'section' attribute of the variable.
                        """
                        if new_value == 1:
                            if prop == 'profile':
                                obj.__setattr__(prop,self.section)
                        else:
                            """ If the target of a binary variable is a member,
                                the binary variable is a member existence variable
                                used in topology optimization. In that case,
                                mark the member as passive.
                            """
                            if prop == 'member':
                                obj.__setattr__('active',False)
                    else:
                        obj.__setattr__(prop, new_value*self.scaling)


    def discrete_violation(self,x):
        """ Identify the violation of value 'x' from allowable (Discrete)
            value. Returns 0 for continous variables
        :return:

        """
        return 0



class IntegerVariable(Variable):
    """ Class for integer variables

        Integer variables can take integer values from the interval lb, ub
    """

    def __init__(self, name="", lb=0, ub=1e5, **kwargs):
        """ Constructor

            Arguments:
                name .. string stating the name of the variable
                lb .. lower bound
                ub .. upper bound
        """

        super().__init__(name, lb, ub, **kwargs)


    def is_allowed_value(self,x):
        """ Check if 'x' is an integer within the range of the variable """

        if x >= self.lb and x <= self.ub:
            """ If x is within INT_TOL of its rounded integer value,
                it is deemed integer
            """
            if math.isclose(x,round(x,0),abs_tol=INT_TOL):
                return True
            else:
                return False
        else:
            return False

    def discrete_violation(self,x):
        """ Identify the violation of value 'x' from allowable (Discrete)
            value. Returns 0 for continous variables
        """

        violation = abs(x-round(x,0))
        if violation <= INT_TOL:
            violation = 0

        return violation


class BinaryVariable(IntegerVariable):
    """ Class for binary variables

        Binary variables can only have values 0 or 1

    """

    def __init__(self, name="", section=None, target=None):
        """ Constructor

            Arguments:
                name .. string stating the name of the variable
        """
        super(IntegerVariable, self).__init__(name, 0, 1,target=target)

        self.section = section       

class IndexVariable(IntegerVariable):
    """
    Class for index variables.
    
    An index variable is usually used to select a profile for
    a member.
    
    """
    def __init__(self, name="", values=None, target=None, **kwargs):
        """ Constructor:
            :param values: list of possible 
        """
        self.values = values

        if values is not None:
            lb = 0
            ub = len(values) - 1
        else:
            lb = None
            ub = None

        super().__init__(name, lb, ub, target=target, **kwargs)

    def substitute(self, new_value):
        """ Substitute the profile corresponding to the integer new_value
        """
        try:
            # Take into account that new_value can be a float
            new_value = int(new_value)
            super().substitute(self.values[new_value])
        except:
            raise ValueError(
                f"Input {new_value} is erroneous "
                "IndexVariable's value must be either"
                " index or a value from the given list!\n"
            f"Values: {self.values}")

        self.value = new_value

    # @property
    # def idx(self):
    #     if self.value in self.values:
    #         return self.values.index(self.value)
    #     return self.value


class DiscreteVariable(Variable):
    """ Class for general discrete variables """

    def __init__(self, name="", values=None, target=None, **kwargs):
        """ Constructor

            Parameters:
            -----------
                :param name: string stating the name of the variable
                :param values: list of allowable discrete values

                Variables:
                ----------
                :ivar values: list of allowable discrete values
                :ivar name: name of the variable (string)
        """

        self.values = values

        if values is not None:
            lb = min(values)
            ub = max(values)

        else:
            lb = None
            ub = None

        super().__init__(name, lb, ub, target=target, **kwargs)

    def is_allowed_value(self,x):
        """ Check if 'x' is one of the values listed in
            self.values
        """

        if x >= self.lb and x <= self.ub:
            """ If x is within INT_TOL of any of the allowed values,
                it is deemed allowed
            """
            if min([abs(val-x) for val in self.values]) < INT_TOL:
                return True

                """
                for val in self.values:
                    if abs(x-val) < INT_TOL:            
                        return True        
                """
            return False
        else:
            return False

    def discrete_violation(self,x):
        """ Identify the violation of value 'x' from allowable (Discrete)
            value.
        """

        violation = min([abs(val-x) for val in self.values])

        if violation <= DISC_TOL:
            violation = 0

        return violation


    def smaller_discrete_value(self,x):
        """ Returns discrete value smaller than 'x' and closest to it  """

        return max(list(filter(lambda val: (val-x<0),self.values)))

            #diff = [val-x for val in self.values]

    def larger_discrete_value(self,x):
        """ Returns discrete value smaller than 'x' and closest to it  """

        return min(list(filter(lambda val: (val-x>0),self.values)))



class Constraint:
    """ General class for constraints """

    def __init__(self, name="", con_type="<",
                 problem=None, fea_required=False,
                 vars=None):
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
        """
        self.name = name
        self.type = con_type
        self.problem = problem
        self.fea_required = fea_required
        if not isinstance(vars, Iterable):
            vars = [vars]
        self.vars = vars
        
        # Lagrange multiplier
        self.mult = 1.0

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
        if self.fea_required and not self.problem.fea_done:            
            self.problem.fea()

    def __repr__(self):
        return f"{type(self).__name__}: {self.name}"

    def neg_call(self, x):
        """ Evaluate -g(x) """
        return -self(x)


class LinearConstraint(Constraint):
    """ Class for linear constraints """

    def __init__(self, a=None, b=None, con_type="<", name="", problem=None, **kwargs):
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

        return self.a.dot(np.array(X)) - b


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

        return self.con(fixed_vals)
    
    def df(self,x):
        """ Evaluate gradient at x """
        
        if self.grad is None:
            # Evaluate gradient numerically
            df = NumGrad(self.con,x)
        else:
            df = self.grad(x)
        
        return df


class ObjectiveFunction:

    def __init__(self, name, obj_fun, obj_type="MIN", problem=None, fea_required=False):

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
                    Because X was revesed, the pop method provides the correct
                    value here.
                """
                fixed_vals[i] = X.pop()
            if not X:
                """ If variable vars[i] is fixed, then use the corresponding
                    value from 'fixed_vals'
                """
                break

        if self.obj_type == "MIN":
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
    def __init__(self,name,c,obj_type="MIN",problem=None):
        """ Constructor:
            input:
                c .. constant array for performing the evaluation
                     c @ x
        """
        obj_fun = self.evaluate

        self.c = np.array(c)

        super().__init__(name,obj_fun,obj_type,problem)


    def evaluate(self,x):
        """ Evaluate function value """
        return self.c.dot(np.array(x))


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


    def add(self, this):
        """
        Adds given object to the problem
        :param this: object to be added
        
        :type: Variable, Constraint, ObjectiveFunction
        """
        # VARIABLES
        if isinstance(this, Variable):
            if this.id is None:
                this.id = len(self.vars)
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

    @property
    def feasible(self):
        """
        Returns problem's feasibility at X
        :return: True/False
        """
        
        res = True

           
        for con in self.cons:
            g = con(self.X)            
            if con.type == '<':
                if g > self.con_tol:
                    res = False
                    break
            elif con.type == '>':
                if g < -self.con_tol:
                    res = False
                    break
            else:
                if abs(g) > self.con_tol:
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

    @lru_cache(CACHE_BOUND)
    def linearize(self, *x):
        """
        Linearizes problem around x

        Ax < B

        :param x:
        :return:
        """
        start = time.time()
        x = np.asarray(x)
        self.substitute_variables(x)

        """ Evaluate objective function at x """
        fx = self.obj(x)
        """ Evaluate nonlinear constraint functions at x """
        b = self.eval_nonlin_cons(x)

        """ Get number of nonlinear constraints """
        m = self.nnonlincons()
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
            if not isinstance(var, BinaryVariable):
                if isinstance(var, IndexVariable):
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

                """ Make finite element analysis """
                #if self.structure and not self.fea_done:
                #    self.fea()
                """ Get variable values """
                xh = [var.value for var in self.vars]
                #print(x,xh)
                if np.linalg.norm(x-xh) == 0:
                    print("Linearize: error with variable substitution - too small step size.")
                """ Evaluate objective function at x + hi*ei """
                if self.grad is None:
                    f_val = self.obj(xh)
                    df[i] = (f_val - fx) / h
                    
                """ Evaluate constraint functions at x + hi*ei """
                a = self.eval_nonlin_cons(xh)
                A[:, i] = (a - b) / h
                
                """ Substitute the original value x[i] to current variable """
                var.substitute(prev_val)


        #B = A @ x.T - b
        B = A.dot(x) - b
        # B = -b

        # Add linear constraints' values
        for con in self.cons:
            if isinstance(con, LinearConstraint):
                A = np.vstack((A, con.a))
                if isinstance(con.b, Variable):
                    B = np.hstack((B, con.b.value))
                else:
                    B = np.hstack((B, con.b))

        end = time.time()
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
            """ Create variables for the conlin problem """
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
                # For nonlinear constraints, CONLIN approximation is made
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
        print(f'Solution is feasible: {self.feasible}')

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
        """ Evaluates equality constraionts """

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
            

def NumGrad(fun, x):
    """ Gradient by finite difference """
    fx = fun(x)
    n = len(x)

    """ if only a single value for step length 'h' is given,
        transform it to a list
    """
    
    
    #if isinstance(h, float):
    #    h = h * np.ones(n)

    df = np.zeros(n)
    #    
    for i in range(len(x)):
        #print("Täällä")
        xh = deepcopy(x)
        h = max(1e-6 * abs(xh[i]), 1e-8)
        xh[i] += h
        #print(xh)
        
        fh = fun(xh)
        #print(fx, fh)

        df[i] = (fh - fx) / h

    return df

def conlin_fun(xk,fk,dfk,pos_vars=True):
    """ Generates convex-linear approximation for a function
        at point 'xk'
        
        :param xk: approximation point
        :param fk: value of the function at 'xk'
        :parm dfk: value of the gradient at 'xk'
    """
    if isinstance(xk,float) or isinstance(xk,int):
        xk = np.array([xk])
        dfk = np.array([dfk])
        
    n = len(xk)

    """ Initialize 'alpha' as zeros. """
    a = np.zeros(n)
    for i, Xk in enumerate(xk):
        """ If the derivative of f with respect to xk is non-negative,
            set a[i] = 1
        """
        if dfk[i] >= 0:
            a[i] = 1

    def generate_a(x):
        
        a = np.ones(n)
        for i, Xk in enumerate(xk):
            if x[i]*dfk[i] < 0:
                a[i] = Xk/x[i]
                
    def eval_conlin(x):
        fC = fk
        if isinstance(x,float) or isinstance(x,int):
            x = np.array([x])
        
        for i in range(n):
            if dfk[i] >= 0:
                a[i] = 1.0
            else:
                a[i] = xk[i]/x[i]
                
        
        fC += (a*dfk).dot(x-xk)
        
        return fC          
            
    return eval_conlin

def mma_fun(xk,fk,dfk,L,U):
    """ Generates the MMA approximation of a function """
    
    if isinstance(xk,float) or isinstance(xk,int):
        xk = np.array([xk])
        dfk = np.array([dfk])
        
    n = len(xk)

    """ Initialize 'p' and 'q' as zeros. """
    p = np.zeros(n)
    q = np.zeros(n)
    for i, (Xk,Lk,Uk) in enumerate(zip(xk,L,U)):
        # Determine the constants p and q
        if dfk[i] >= 0:
            p[i] = (Uk-Xk)**2*dfk[i]
        else:
            q[i] = -(Xk-Lk)**2*dfk[i]
    
    # Evaluate r(xk)
    rk = fk-sum(p/(U-xk)+q/(xk-L))
    
    def eval_mma(x):
        # Evaluate MMA approximation
        return rk + sum(p/(U-x)+q/(x-L))


    return eval_mma

def Linearize(fun, x, grad=None):
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
        h = 0.01 * x  # check if this is OK
        # h = np.ones_like(x) * 1e-2
        a = NumGrad(fun, x, h)
    else:
        a = grad(x)

    c = a.dot(x)

    return a, b, c


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
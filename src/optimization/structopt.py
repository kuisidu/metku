# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:00:48 2019

Constructing and solving structural optimization problems

Idea: all the data is stored in OptimizationProblem object.
As different optimization solvers are applied, the methods of
OptimizationProblem class will be used to write the problem data in
a form suitable for each solver.

@author: kmela
"""

import math
import numpy as np
import scipy.optimize as sciop
from itertools import product

INT_TOL = 1e-4


class Variable:
    """ Class for optimization variables
    """

    def __init__(self, name, lb, ub, target=None, profiles=None):
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
                    objects .. list of objects (e.g. frame members) that are
                                affected by the variable

                if target is 'None', then the variable does not change anything
                in the structure

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
        self.value = None
        self.lb = lb
        self.ub = ub
        self.target = target
        self.profiles = profiles

    # @property
    # def value(self):
    #     return self.target['objects'][0].__dict__[self.target['property']]
    #
    #
    # @value.setter
    # def value(self, new_val):
    #     for obj in self.target['objects']:
    #         obj.__dict__[self.target['property']] = new_val

    def substitute(self, new_value):
        """ Substitute a new value for the variable """
        self.value = new_value

        """ Modify target object(s) if any """
        if self.target is not None:
            for obj in self.target['objects']:
                if self.target['property'] == 'AREA' \
                        or self.target['property'] == 'A':
                    obj.A = new_value
                elif self.target['property'] == 'IY':
                    obj.I[0] = new_value
                elif self.target['property'] == 'IZ':
                    obj.I[1] = new_value
                elif self.target['property'] == 'WELY':
                    obj.Wel[0] = new_value
                elif self.target['property'] == 'WELZ':
                    obj.Wel[1] = new_value
                elif self.target['property'] == 'WPLY':
                    obj.Wpl[0] = new_value
                elif self.target['property'] == 'WPLZ':
                    obj.Wpl[1] = new_value
                elif self.target['property'] == 'H':
                    obj.h = new_value
                elif self.target['property'] == 'TW':
                    obj.tw = new_value
                elif self.target['property'] == 'TT':
                    obj.tt = new_value
                elif self.target['property'] == 'TB':
                    obj.tb = new_value
                elif self.target['property'] == 'BT':
                    obj.bt = new_value
                elif self.target['property'] == 'BB':
                    obj.bb = new_value
                elif self.target['property'] == 'BF':
                    if isinstance(obj.b, list):
                        for i in range(len(obj.b)):
                            obj.b[i] = new_value
                    else:
                        obj.b = new_value
                elif self.target['property'] == 'TF':
                    if isinstance(obj.tf, list):
                        for i in range(len(obj.b)):
                            obj.tf[i] = new_value
                    else:
                        obj.tf = new_value
                elif self.target['property'] == 'PROFILE':
                    """ susbstitute profile """
                    if isinstance(new_value, str):
                        obj.profile = new_value
                    elif isinstance(new_value, int):
                        obj.profile = self.profiles[new_value]
                    else:
                        raise ValueError('Variable type must be either str or '
                                         'int!')

                elif self.target['property'] == 'x':
                    obj.x = new_value

                elif self.target['property'] == 'y':
                    obj.y = new_value

                elif self.target['property'] == 'loc':
                    obj.loc = new_value

    def is_allowed_value(self,x):
        """ Returns True, if 'x' is an allowed value for the variable
            and False otherise.
            
            This method is implemented separately for different variable types.
        """
        
        # For continuous variable, and a value is allowable, if it is
        # Between the lower and upper bounds
        if x < self.lb or x > self.ub:
            return False
        else:
            return True
    
    def discrete_violation(self,x):
        """ Identify the violation of value 'x' from allowable (Discrete)
            value. Returns 0 for continous variables
        """
        
        return 0

class IntegerVariable(Variable):
    """ Class for integer variables

        Integer variables can take integer values from the interval lb, ub
    """

    def __init__(self, name="", lb=0, ub=1e5):
        """ Constructor

            Arguments:
                name .. string stating the name of the variable
                lb .. lower bound
                ub .. upper bound
        """

        Variable.__init__(self, name, lb, ub)

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
        
        return abs(x-round(x,0))


class BinaryVariable(IntegerVariable):
    """ Class for binary variables

        Binary variables can only have values 0 or 1

    """

    def __init__(self, name=""):
        """ Constructor

            Arguments:
                name .. string stating the name of the variable
        """

        IntegerVariable.__init__(self, name, 0, 1)


class IndexVariable(Variable):

    def __init__(self, name="", values=None, target=None):
        super().__init__(name=name, target=target)

        self.values = values

        if values is not None:
            lb = 0
            ub = len(values) - 1
        else:
            lb = None
            ub = None

        super().__init__(name, lb, ub, target=target)

    def substitute(self, new_value):

        if new_value not in self.values:
            try:
                new_value = self.values[new_value]
            except:
                raise ValueError(
                    f"Input {new_value} is erroneous"
                    "IndexVariable's value must be either"
                    " index or a value from the given list!")

        super().substitute(new_value)

    @property
    def idx(self):
        if self.value in self.values:
            return self.values.index(self.value)
        return self.value


class DiscreteVariable(Variable):
    """ Class for general discrete variables """

    def __init__(self, name="", values=None, target=None):
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

        super().__init__(name, lb, ub, target=target)

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
        
        return min([abs(val-x) for val in self.values])
    
    def smaller_discrete_value(self,x):
        """ Returns discrete value smaller than 'x' and closest to it  """
        
        return max(list(filter(lambda val: (val-x<0),self.values)))
        
        #diff = [val-x for val in self.values]
    
    def larger_discrete_value(self,x):
        """ Returns discrete value smaller than 'x' and closest to it  """
        
        return min(list(filter(lambda val: (val-x>0),self.values)))
    
class Constraint:
    """ General class for constraints """

    def __init__(self, name="", con_type="<", parent=None):
        self.name = name
        self.type = con_type
        self.parent = parent
        self.fea_required = False

    def __call__(self, x):
        """
        Runs finite element analysis on problem's structure if needed
        """
        if self.fea_required and not self.parent.fea_done:
            if np.any(self.parent.X != x):  # Replace with norm
                self.parent.substitute_variables(x)
            self.parent.fea()

    def neg_call(self, x):
        """ Evaluate -g(x) """
        return -self(x)


class LinearConstraint(Constraint):
    """ Class for linear constraints """

    def __init__(self, a=None, b=None, con_type="<", name="", parent=None):
        """ Constructor """

        self.a = a
        self.b = b
        Constraint.__init__(self, name, con_type, parent)

    def __call__(self, x):
        """ Evaluate constraint at x """
        if isinstance(self.b, Variable):
            b = self.b.value
        else:
            b = self.b
        return np.array(self.a).dot(np.array(x)) - b


class NonLinearConstraint(Constraint):
    """ Class for linear constraints """

    def __init__(self, con_fun, con_type="<", name="", parent=None):
        """ Constructor
            con_fun: function that returns the value of constraint function
                     g(x) <= 0 (or >= 0 or = 0)
        """

        self.con = con_fun
        Constraint.__init__(self, name, con_type, parent)

    def __call__(self, x):
        """ Evaluate constraint at x """
        super().__call__(x)
        return self.con(x)


class OptimizationProblem:
    """ Class for defining and operating with optimization problems """

    def __init__(self, name="", variables=[], constraints=[], objective=None, \
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
                :param structure: structure to be optimized (optional)
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
        """

        self.con_tol = 1e-4
        self.name = name
        self.vars = variables
        self.cons = constraints
        self.obj = objective
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

    @property
    def feasible(self):
        """
        Returns problem's feasibility
        :return: True/False
        """
        return np.all(self.eval_cons(self.X) <= self.con_tol)

    def linearize(self, x):
        """
        Linearizes problem around x

        Ax < B

        :param x:
        :return:
        """
        import time
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
                    h = max(0.01 * abs(prev_val), 1e-4)
                    var.substitute(prev_val + h)

                """ Make finite element analysis """
                fea_start = time.time()
                if self.structure:
                    self.fea()
                fea_end = time.time()
                """ Get variable values """
                xh = [var.value for var in self.vars]
                """ Evaluate objective function at x + hi*ei """
                f_val = self.obj(xh)
                """ Evaluate constraint functions at x + hi*ei """
                a = self.eval_nonlin_cons(xh)
                A[:, i] = (a - b) / h
                df[i] = (f_val - fx) / h
                """ Substitute the original value x[i] to current variable """
                var.substitute(prev_val)

        B = A @ x.T - b
        # B = -b

        # Add linear constraints' values
        for con in self.cons:
            if isinstance(con, LinearConstraint):
                A = np.vstack((A, con.a))
                B = np.hstack((B, con.b))

        end = time.time()
        # print("Linearization took: ", end - start, " s")

        return A, B, df, fx

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

        for con in non_lin_cons:
            g.append(con(X))

        return np.asarray(g)

    def fea(self):
        """
        Runs finite element analysis on structure
        """
        self.structure.calculate()
        self.fea_done = True
        self.num_fem_analyses += 1

    def __call__(self, x, prec=2):
        """ Call method evaluates the objective function and all constraints
            at x and returns their values
        """
        self.substitute_variables(x)
        fx = self.obj(x)
        print("** {0} **".format(self.name))
        print(f'Solution is feasible: {self.feasible}')
        print(f'X: {[round(val, prec) for val in x]}')
        if len(x) < 10:
            print("Variables:")
            print("----------")
            for i in range(len(x)):
                print(f"{self.vars[i].name} = {x[i]:.{prec}f}")

            print("----------\n")

        print(f"Objective function = {fx:.{prec}f}\n")

        print("Constraints:")
        print("----------")

        for con in self.cons:
            g = con(x)
            print(f"{con.name}: {g:.{prec}f} {con.type} 0")

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

        # Save starting point
        if not np.any(self.X):
            self.x0 = xvals.copy()

        self.X = xvals
        self.fea_done = False
        for x, var in zip(xvals, self.vars):
            var.substitute(x)

        #
        # for i in range(len(xvals)):
        #     self.substitute_variable(i, xvals[i])

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
    
    def discrete_neighborhood(self,x,k,method='nearest'):
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



def NumGrad(fun, h, x):
    """ Gradient by finite difference """
    fx = fun(x)
    n = len(x)

    """ if only a single value for step length 'h' is given,
        transform it to a list
    """

    if isinstance(h, float):
        h = h * np.ones(n)

    df = np.zeros(n)

    np.eye(3, 1, -2)
    #
    for i in range(len(x)):
        xh = x + h[i] * np.eye(1, n, i)
        print(xh)
        fh = fun(xh)
        #     print(fx, fh)

        df[i] = (fh - fx) / h[i]

    return df


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
        a = NumGrad(fun, h, x)
    else:
        a = grad(x)

    c = a.dot(x)

    return a, b, c


if __name__ == '__main__':
    # dvars = []
    # dvars.append(Variable("Height", 150, 1000))
    # dvars.append(Variable("Flange width", 150, 500))
    # dvars.append(Variable("Flange thickness", 5, 40))
    # dvars.append(Variable("Web thickness", 5, 40))
    #
    
    
    
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
    
    
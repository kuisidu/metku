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

import numpy as np
import scipy.optimize as sciop

class Variable:
    """ Class for optimization variables     
    """
    
    def __init__(self,name,lb,ub,target=None,profiles=None):
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
                    property .. unique identifier that shows which property of the
                                structure is affected by the variable
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
        
    def substitute(self,new_value):
        """ Substitute a new value for the variable """
        
        self.value = new_value
        
        """ Modify target object(s) if any """
        if self.target is not None:
            for obj in self.target['objects']:
                if self.target['property'] == 'AREA':
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
                elif self.target['property'] == 'BF':
                    if isinstance(obj.b,list):
                        for i in range(len(obj.b)):
                            obj.b[i] = new_value
                    else:                        
                        obj.b = new_value
                elif self.target['property'] == 'TF':    
                    if isinstance(obj.tf,list):
                        for i in range(len(obj.b)):
                            obj.tf[i] = new_value
                    else:                        
                        obj.tf = new_value
                elif self.target['property'] == 'PROFILE':
                    """ susbstitute profile """
                    try:
                        obj.profile = self.profiles[new_value]
                    except ValueError:
                        print('Index variable value must be integer!')
                    
        
class IntegerVariable(Variable):
    """ Class for integer variables 
    
        Integer variables can take integer values from the interval lb, ub
    """
    
    def __init__(self,name="",lb=0,ub=1e5):
        """ Constructor
            
            Arguments:
                name .. string stating the name of the variable
                lb .. lower bound
                ub .. upper bound
        """

        Variable.__init__(self,name,lb,ub)
        
class BinaryVariable(IntegerVariable):
    """ Class for binary variables 
        
        Binary variables can only have values 0 or 1
    
    """
    
    def __init__(self,name=""):
        """ Constructor
            
            Arguments:
                name .. string stating the name of the variable
        """
        
        IntegerVariable.__init__(self,name,0,1)
        
class DiscreteVariable(Variable):
    """ Class for general discrete variables """
    
    def __init__(self,name="",values=None):
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
        
        if values != None:
            lb = min(values)
            ub = max(values)
        else:
            lb = None
            ub = None
        
        Variable.__init__(self,name,lb,ub)
        
class Constraint:
    """ General class for constraints """
    
    def __init__(self,name="",con_type="<"):
                
        self.name = name
        self.type = con_type

class LinearConstraint(Constraint):
    """ Class for linear constraints """

    def __init__(self,a=None,b=None,con_type="<",name=""):
        """ Constructor """
        
        self.a = a
        self.b = b
        Constraint.__init__(self,name,con_type)
        
    def __call__(self,x):
        """ Evaluate constraint at x """    
        return np.array(self.a).dot(np.array(x))-self.b
    
    def neg_call(self,x):
        """ Evaluate -g(x) """
        return -self(x)

class NonLinearConstraint(Constraint):
    """ Class for linear constraints """

    def __init__(self,con_fun,con_type="<",name=""):
        """ Constructor 
            con_fun: function that returns the value of constraint function
                     g(x) <= 0 (or >= 0 or = 0)
        """
        
        self.con = con_fun
        Constraint.__init__(self,name,con_type)
    
    def __call__(self,x):
        """ Evaluate constraint at x """        
        return self.con(x)
    
    def neg_call(self,x):
        """ Return -g(x) """
        return -self.con(x)
    
class OptimizationProblem:
    """ Class for defining and operating with optimization problems """
    
    def __init__(self,name="",variables=[],constraints=[],objective=None,\
                 gradient=None,hess=None, structure=None, profiles=None):
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

        self.name = name        
        self.vars = variables
        self.cons = constraints
        self.obj = objective
        self.grad = gradient
        self.hess = hess
        self.structure = structure
        self.profiles = profiles
        
    def __call__(self,x):
        """ Call method evaluates the objective function and all constraints
            at x and returns their values
        """
        
        fx = self.obj(x)
        print("** {0} **".format(self.name))
        if len(x) < 10:
            print("Variables:")
            print("----------")
            for i in range(len(x)):
                print("{0} = {1:4.2f}".format(self.vars[i].name,x[i]))
            
            print("----------\n")    
            
        print("Objective function = {0:4.2f}\n".format(fx))
        
        print("Constraints:")
        print("----------")
        
        for con in self.cons:
            g = con(x)            
            print("{0}: {1:4.2f} {2} 0".format(con.name,g,con.type))
        
    def eval_con(self,x):
        """ Constraint evaluation
            
        """
        g = []
        
        for con in self.cons:
            g.append(con(x))
        
        return g        
    
    def eval_eq_con(self,x):
        """ Evaluates equality constraionts """
        
        geq = []
        
        for con in self.cons:
            if con.type == "=":
                geq.append(con(x))
    
    def eval_ineq_con(self,x):
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

    def add_variables(self,variables):
        """ Adds new variables to the problem """
        
        self.vars.append(variables)
        
    def add_constraints(self,constraints):
        """ Adds new constraints to the problem """
        self.cons.append(constraints)
    
    def add_structure(self,structure):
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
    
    def substitute_variable(self,i,xval):
        """ Substitutes the value 'xval' for variable i """
        var = self.vars[i]
        var.value = xval
        
        if var.target is not None:
            for obj in var.target['objects']:
                if var.target['property'] == 'AREA':
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
                    
    
    
    def substitute_variables(self,xvals):
        """ Substitute variable values from xval to structure """
        
        for i in range(len(xvals)):
            self.substitute_variable(i,xvals[i])
            
    
    def solve(self,solver="slsqp",**kwargs):
        """ Solve the optimization problem 
            input:
                solver .. name of the solver
                **kwargs .. list of keyword argument
        """
        lb, ub = self.var_bounds()
                
        
        if solver == "slsqp":
            jac = "2-point"            
            bnds = sciop.Bounds(lb,ub)        
            
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
            #cons = [{"type":"ineq", "fun":self.eval_ineq_con}]
            
            
            cons = []
            for con in self.cons:
                if con.type == "=":
                    new_con = {"type":"eq", "fun":con}
                elif con.type == "<":
                    new_con = {"type":"ineq", "fun":con.neg_call}
                else:
                    new_con = {"type":"ineq", "fun":con}
                    
                cons.append(new_con)
            
            #print(cons)
            
            res = sciop.minimize(self.obj,kwargs["x0"],method=solver,\
                           jac=jac,bounds=bnds,constraints=cons)
        elif solver == "trust-constr":
            jac = "2-point"        
            hess = "2-point"
            bnds = sciop.Bounds(lb,ub)            
            
            if self.grad is not None:
                jac = self.grad
            
            if self.hess is not None:
                hess = self.hess          

            res = sciop.minimize(self.obj,kwargs["x0"],method=solver,\
                           jac=jac,hess=hess,bounds=bnds)
        else:
            print("Solver unidentified.")
            res = None
        
        return res


def NumGrad(fun,h,x):
    """ Gradient by finite difference """
    fx = fun(x)
    n = len(x)
   
    df = np.zeros(n)
    
    for i in range(len(x)):
        dx = np.zeros(n)
        dx[i] = h[i]
        print(x+dx)
        fh = fun(x+dx)        
        df[i] = (fh-fx)/h[i]

    return df

def Linearize(fun,x,grad=None):
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
        h = 0.01*x # check if this is OK
        a = NumGrad(fun,h,x)
    else:
        a = grad(x)
    
    c = a.dot(x)
    
    return a, b, c



if __name__ == '__main__':
    
    dvars = []
    dvars.append(Variable("Height",150,1000))
    dvars.append(Variable("Flange width",150,500))
    dvars.append(Variable("Flange thickness",5,40))
    dvars.append(Variable("Web thickness",5,40))
    
    p = OptimizationProblem(name="I-Beam Weight Minimization",variables=dvars)
    
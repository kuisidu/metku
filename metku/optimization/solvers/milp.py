# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Tue 30 Mar 2021

Solving mixed-integer linear programming problems.

@author: kmela
"""

# TODO! Add a MILP solver from ortools

from scipy.optimize import linprog
import numpy as np
import time

# import gurobipy as grb

from metku.optimization.solvers.optsolver import OptSolver
import metku.optimization.structopt as sopt

#from ortools.linear_solver import pywraplp

class MILP(OptSolver):

    def __init__(self,algorithm="gurobi"):
        super().__init__()
        self.algo = algorithm
        

    def solve(self,problem,verb=False):
        """
        
        """
    
        import gurobipy as grb
                
        lp = grb.Model("milp")
        
                
        x = []
        for var in problem.vars:
            """
                create variable, including the objective function
                coefficient and name
            """
            if isinstance(var,sopt.BinaryVariable):
                x.append(lp.addVar(var.lb,var.ub,vtype=grb.GRB.BINARY,name=var.name))
            else:
                x.append(lp.addVar(var.lb,var.ub,vtype=grb.GRB.CONTINUOUS,name=var.name))
                
        #nvars = len(x)
                
        if problem.obj.obj_type == "MIN":
            obj_sense = grb.GRB.MINIMIZE
        else:
            obj_sense = grb.GRB.MAXIMIZE
        
        lp.setObjective(grb.LinExpr(problem.obj.c,x),sense=obj_sense)
        
        
        for con in problem.cons:
            if isinstance(con,sopt.LinearConstraint):
                if con.type == '<':
                    con_sense = grb.GRB.LESS_EQUAL
                            #lp.addLConstr(LinExpr(con.a, x), grb.GRB.LESS_EQUAL, con.b)
                            #lp.addConstr(grb.quicksum([con.a[j]*x[j] for j in range(nvars)]) <= con.b[i])
                elif con.type == '>':
                    con_sense = grb.GRB.GREATER_EQUAL
                            #lp.addLConstr(LinExpr(con.a, x), grb.GRB.LESS_EQUAL, con.b)
                            #lp.addConstr(grb.quicksum([con.a[j]*x[j] for j in range(nvars)]) >= con.b[i])
                else:
                    con_sense = grb.GRB.EQUAL
                            #lp.addConstr(grb.quicksum([con.a[j]*x[j] for j in range(nvars)]) = con.b[i])
                
                lp.addLConstr(grb.LinExpr(con.a, x), con_sense, con.b)
                

        lp.update()                
       

        #for var in lp.getVars():
        #    print(var.lb,var.ub)
        
        
        if not verb:
            lp.setParam("LogToConsole",0)
                #qp.update()
        lp.optimize()
                
        if lp.Status == grb.GRB.OPTIMAL:
            self.feasible = True
        elif lp.Status == grb.GRB.INFEASIBLE:
            self.feasible = False
                
        X = []
        for v in lp.getVars():
            #print('%s %g' % (v.varName, v.x))                    
            X.append(v.x)
        
        # bounds = [x*self.move_limits for x in self.X]
        #
        # res = linprog(df, A, B, bounds=bounds)

        self.X = X
        #return X
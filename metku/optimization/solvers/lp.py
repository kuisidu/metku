# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Fri Nov 01 2019

Solving linear programming problems.

@author: kmela
"""

from scipy.optimize import linprog
import numpy as np
import time

# import gurobipy as grb

from metku.optimization.solvers.optsolver import OptSolver
import metku.optimization.structopt as sopt

#from ortools.linear_solver import pywraplp

class LP(OptSolver):

    def __init__(self,algorithm="gurobi"):
        super().__init__()
        self.algo = algorithm
        

    def solve(self,problem,verb=False):
        """
        
        """
    
        if self.algo == "gurobi":
            import gurobipy as grb
            
            try:
                lp = grb.Model("lp")
            except:
                self.algo = 'ortools'
                from ortools.linear_solver import pywraplp
                lp = pywraplp.Solver.CreateSolver('GLOP')
                
        else:
            from ortools.linear_solver import pywraplp
            lp = pywraplp.Solver.CreateSolver('GLOP')
                
        x = []
        for var in problem.vars:
            """
                create variable, including the objective function
                coefficient and name
            """
            if self.algo == "gurobi":
                x.append(lp.addVar(var.lb,var.ub,vtype=grb.GRB.CONTINUOUS,name=var.name))
            else:
                x.append(lp.NumVar(var.lb,var.ub,var.name))
                
        #nvars = len(x)
                
        if problem.obj.obj_type == "MIN":
            if self.algo == "gurobi":
                obj_sense = grb.GRB.MINIMIZE
            else:
                lp.Minimize(lp.Sum(C * X for C, X in zip(problem.obj.c,x)))
        else:
            if self.algo == "gurobi":
                obj_sense = grb.GRB.MAXIMIZE
            else:
                lp.Maximize(lp.Sum(C * X for C, X in zip(problem.obj.c,x)))
        
        if self.algo == "gurobi":
            lp.setObjective(grb.LinExpr(problem.obj.c,x),sense=obj_sense)
            
        # Generate constraints
        for con in problem.cons:
            if isinstance(con,sopt.LinearConstraint):
                if self.algo == "gurobi":
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
                else:
                    lhs = lp.Sum(A * X for A, X in zip(con.a,x))
            
                    if con.type == '<':
                        lp.Add(lhs <= con.b)
                    elif con.type == '>':
                        lp.Add(lhs >= con.b)
                    else:
                        lp.Add(lhs == con.b)
        
        if self.algo == "gurobi":
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
        else:
            status = lp.Solve()
            
            X = []
            if status == lp.OPTIMAL:
                for v in lp.variables():
                    X.append(v.solution_value())    

        self.X = X
        #return X

def Arora_Ex62():
    x1 = sopt.Variable('x1',0,1e4)
    x2 = sopt.Variable('x2',0,1e4)
    p = sopt.OptimizationProblem(name="Arora, Example 6.2", variables=[x1,x2])
    
    c = [-400,-600]
    
    obj = sopt.LinearObjective("Objective",c,problem=p)
    p.add(obj)
    
    cons = []
    cons.append(sopt.LinearConstraint([1,1],16,name="Con1",problem=p))
    cons.append(sopt.LinearConstraint([1/28,1/14],1,name="Con2",problem=p))
    cons.append(sopt.LinearConstraint([1/14,1/24],1,name="Con3",problem=p))
    
    for con in cons:
        p.add(con)
    
    solver = LP("gurobi")
    solver.solve(p)
    p(solver.X)
    
def Arora_Ex63():
    x1 = sopt.Variable('x1',0,1e4)
    x2 = sopt.Variable('x2',0,1e4)
    p = sopt.OptimizationProblem(name="Arora, Example 6.3", variables=[x1,x2])
    
    c = [4,5]
    
    obj = sopt.LinearObjective("Objective",c,obj_type="MAX",problem=p)
    p.add(obj)
    
    cons = []
    cons.append(sopt.LinearConstraint([-1,1],4,name="Con1",problem=p))
    cons.append(sopt.LinearConstraint([1,1],6,name="Con2",problem=p))

    
    for con in cons:
        p.add(con)
    
    #solver = LP("gurobi")
    solver = LP("ortools")
    solver.solve(p)
    p(solver.X)

    return solver, p

def Arora_Ex67():
    x = []
    for i in range(4):
        x.append(sopt.Variable('x'+str(i+1),0,1e4))
            
    p = sopt.OptimizationProblem(name="Arora, Example 6.7", variables=x)
    
    c = [4,5,0,0]
    
    obj = sopt.LinearObjective("Objective",c,obj_type="MAX",problem=p)
    p.add(obj)
    
    cons = []
    cons.append(sopt.LinearConstraint([-1,1,1,0],4,name="Con1",con_type='=',problem=p))
    cons.append(sopt.LinearConstraint([1,1,0,1],6,name="Con2",con_type='=',problem=p))

    
    for con in cons:
        p.add(con)
    
    #solver = LP("gurobi")
    solver = LP("ortools")
    solver.solve(p)
    p(solver.X)
    
    return solver, p

if __name__ == '__main__':
    
    solver, p = Arora_Ex63()
    

    
    
    
    


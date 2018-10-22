# Imported libraries
from tables_and_tuples import *
from classes import EndPlate
from timeit import default_timer as timer
import matplotlib.pyplot as plt

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import differential_evolution


class Optimizer:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, connection, beam_id):
        self.beam_id = beam_id
        self.cn = connection
        self.res = None

        self.vars = []
        self.x0 = np.array([])
        self.bnds = []
        self.var_name = []
        self.constraints = []

        self.fun0 = 0.0

    def add_design_variable(self, var, var_init, var_bounds, var_name=""):
        self.vars.append(var)
        self.x0 = np.append(self.x0, [var_init])
        self.bnds.append(var_bounds)
        self.var_name.append(var_name)

    def add_constraint(self, type, fun, name, var_update='global'):
        self.constraints = list(self.constraints)
        if var_update == 'global':
            self.constraints.append({'type': type, 'fun': lambda x: self.global_constraint_(fun, x), 'name': name})
        elif var_update == 'local':
            self.constraints.append({'type': type, 'fun': lambda x: self.local_constraint_(fun, x), 'name': name})
        self.constraints = tuple(self.constraints)

    def update_design_variables_(self, x):
        for i in range(0, len(self.vars)):
            self.vars[i](x[i])
            print("x[{0}] = {1}, {2}".format(i, x[i], self.var_name[i]))
        self.cn.update()

    def global_constraint_(self, fun, x):
        self.update_design_variables_(x)
        return fun()

    def local_constraint_(self, fun, x):
        for i in range(0, len(self.vars)):
            self.vars[i](x[i])
        return fun()

    def add_objective(self, fun):
        self.objective = fun

    def objective_(self, x):
        self.update_design_variables_(x)
        self.cn.info(values=False, cost_info=False, warnings_and_errors=False, figures=True, side=True, front=True)
        return self.objective()

    def solve(self, method='SLSQP'):
        print("Optimization started!")
        print('Number of variables:   {0}'.format(len(self.vars)))
        print('Number of constraints: {0} (Not including variable bounds)'.format(len(self.constraints)))
        self.fun0 = self.objective_(self.x0)
        self.method = method
        if self.method == 'SLSQP':
            self.SLSQP_optimizer()
        elif self.method == 'Diff_evolution':
            self.evolutionary_optimizer()
        else:
            print("Optimization method '{0}' not added! Usable optmization methods= 'SLSQP' or 'Diff_evolution'"
                  .format(method))

    def SLSQP_optimizer(self):
        # 'SLSQP'=Sequential Least SQuares Programming
        self.res = minimize(self.objective_,
                            self.x0,
                            method='SLSQP',
                            bounds=self.bnds,
                            constraints=self.constraints,
                            options={'disp': True, 'maxiter': 500, 'eps': 0.05, 'ftol': 1.0e-4})
        for const in self.constraints:
            if const['fun'](self.res.x) < -1.0e-3:
                self.res.message = 'No feasible results found.'
        print("\n" + self.res.message)

    def evolutionary_optimizer(self):
        # 'diff_evolution'=Differential evolution
        print(self.NOT_IMPLEMENTED + "Evolutionary optimization")
        self.res = differential_evolution(self.objective_, self.bnds, disp=True)
        print("\n" + self.res.message)

    def info(self):
        print("\nDesign variable! current! inital! bounds")
        for n in range(0, len(self.vars)):
            print("         {0:6}!  {1:{fm}}! {2:{fm}}! {3}"
                  .format(self.var_name[n], self.res.x[n], self.x0[n], self.bnds[n], fm="6.2f"))
        print("\nInitial objective function value: {0:{fm}}".format(self.fun0, fm=".2f"))
        print("Final objective function value:   {0:{fm}}".format(self.res.fun, fm=".2f"))



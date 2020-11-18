import time

import numpy as np
from scipy.optimize import minimize

try:
    from metku.optimization.solvers.optsolver import OptSolver
except:
    from optimization.solvers.optsolver import OptSolver


class SLSQP(OptSolver):

    def __init__(self):
        super().__init__()

    def calc_constraints(self, x=None):
        """
        Calculates constraint values and saves them to numpy array

        Parameters:
        ------------
        :param x: (Optional, default: self.X) Point where constraint values are calculated

        """
        constr_vals = []
        if x is None:
            x = self.X
        # Calculate constraints
        self.constr_vals = self.problem.eval_cons(x)

    def _create_eqcons(self):
        """ Creates a list of equality constraints

        :return: constraints
        """
        eqcons = []

        for con in self.problem.cons:
            if con.type == "=":
                eqcons.append(con)

        return eqcons

    def _create_ieqcons(self):
        """ Creates a list of inequality constraints

        :return: list of equality constraint functions
        """
        ieqcons = []

        for con in self.problem.cons:
            if con.type == "<" or con.type == "<=":
                ieqcons.append(con.neg_call)
            elif con.type == ">" or con.type == ">=":
                ieqcons.append(con)

        return ieqcons

    def _create_bounds(self):
        """ Creates a list of optimization variables' bounds

        :return: list of boundaries (lb, ub)
        """
        bounds = []
        for var in self.problem.vars:
            bound = (var.lb, var.ub)
            bounds.append(bound)

        return bounds

    def solve(self, problem, maxiter=100, x0=None, verb=False):
        """
        Solves given problem

        :param problem: problem to be solved
        :param maxiter: number of maximum iterations
        :param x0: initial guess

        :return: fopt, xopt
        """
        self.verb = verb
        self.problem = problem
        self.maxiter = maxiter
        bounds = self._create_bounds()
        eqcons = self._create_eqcons()
        ieqcons = self._create_ieqcons()
        
        print(bounds)
        
        # Create initial guess if one isn't provided
        if x0 is None:
            x0 = []
            for var in problem.vars:
                x0.append(var.ub)  # * np.random.uniform())

        """
            'eps' is the step length in numerical differentiation
        """
        options = {'maxiter': maxiter,
                   'disp': verb,
                   'iprint': 1,
                   'eps': 1e-6,
                   'ftol': 1e-5}

        constraints = []

        for eqcon in eqcons:
            con = {'type': 'eq', 'fun': eqcon}
            constraints.append(con)

        for ineqcon in ieqcons:
            con = {'type': 'ineq', 'fun': ineqcon}
            constraints.append(con)

        self.start = time.time()
        out = minimize(self.problem.obj,
                       x0,
                       method='slsqp',
                       tol=1e-5,
                       bounds=bounds,
                       constraints=constraints,
                       options=options,
                       #callback=self.callback,
                       )

        print(out)

        #print(self.calc_constraints(out.x))
        if out.success == True:
            self.feasible = True
        else:
            self.feasible = False
        self.result = out
        self.best_x = out.x
        self.best_f = out.fun
        self.X = out.x
        return self.best_f, self.best_x

    def callback(self, xk, iterations=[]):
        i = len(iterations)
        iterations.append(1)
        if self.verb:
            fval = self.problem.obj(xk)
            self.calc_constraints(xk)
            end = time.time()
            print(f'\r Iteration {i + 1} / {self.maxiter}: Obj: {fval} Feasible: {self.problem.feasible} ' \
                f'max g: {max(self.constr_vals):.4f}')
            print(f"Iteration took: {end - self.start :.2f} s")
            self.start = end


class COBYLA(SLSQP):

    def solve(self, problem, maxiter=100, x0=[]):
        """
        Solves given problem

        :param problem: problem to be solved
        :param maxiter: number of maximum iterations
        :param x0: initial guess

        :return: fopt, xopt
        """

        self.problem = problem
        eqcons = self._create_eqcons()
        ieqcons = self._create_ieqcons()

        # Create initial guess if one isn't provided
        if not len(x0):
            x0 = []
            for var in problem.vars:
                x0.append(var.ub)

        options = {'maxiter': maxiter,
                   'disp': True,
                   'rhobeg': 1e2}

        constraints = []

        for eqcon in eqcons:
            con = {'type': 'eq', 'fun': eqcon}
            constraints.append(con)

        for ineqcon in ieqcons:
            con = {'type': 'ineq', 'fun': ineqcon}
            constraints.append(con)

        out = minimize(self.problem.obj,
                       x0,
                       method='COBYLA',
                       constraints=constraints,
                       options=options)
        # print(out)
        return out


if __name__ == '__main__':
    from metku.optimization.benchmarks import *

    problem = TenBarTruss('continuous')
    solver = SLSQP()
    xopt, fopt = solver.solve(problem, maxiter=10)

    problem(xopt)

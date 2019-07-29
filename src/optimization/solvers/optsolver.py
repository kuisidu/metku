import time

import numpy as np
from scipy.optimize import basinhopping

from optimization.structopt import OptimizationProblem


class OptSolver:
    """
        Base class for optimization problem solvers
    """

    def __init__(self):
        self.constr_vals = np.array([])
        self.X = np.array([])
        self.problem = None

    def calc_constraints(self, x=[]):
        """
        Calculates constraint values and saves them to numpy array

        Parameters:
        ------------
        :param x: (Optional, default: self.X) Point where constraint values are calculated

        """
        constr_vals = []
        if not len(x):
            x = self.X
        # Calculate constraints
        for i in range(len(self.problem.cons)):
            constr_vals.append(self.problem.cons[i](x))
        self.constr_vals = np.asarray(constr_vals)
        return np.asarray(constr_vals)

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

    def solve(self, maxiter=100, maxtime=30):
        pass

class DiscreteVNS(OptSolver):
    """
        Variable Neighborhood Search (VNS) solver for discrete optimization problems
    """

    def __init__(self, step_length=1):
        super().__init__()
        self.step_length = int(step_length)
        self.X = None
        self.problem = None
        self.constr_vals = None
        self.fval = 10e3
        self.best_fval = 10e3

    def shake(self):
        action = np.random.randint(-self.step_length, self.step_length + 1,
                                   len(self.X))
        x = np.clip(self.X.copy() + action, 0, self.problem.vars[0].ub)
        return x

    def first_improvement(self, x):
        actions = []
        r = -1
        while r <= 0:
            action = np.random.randint(-self.step_length, self.step_length + 1,
                                       len(x))
            if list(action) not in actions:
                actions.append(list(action))
                s, r, d, _ = self.step(action, x)

        if d:
            print(f"DONE! Obj.val: {self.best_fval}")
            print(x)

        return x

    def neighbourhood_change(self, x):
        self.X = x

    def step(self, action, X=[]):

        if not len(X):
            X = self.X.copy()

        prev_cons_vals = self.constr_vals.copy()
        X += action

        for i in range(len(X)):
            X[i] = np.clip(X[i], self.problem.vars[i].lb,
                           self.problem.vars[i].ub)

        self.problem.substitute_variables(X)

        # Calculate objective
        fval = self.problem.obj(X)
        dfval = self.best_fval - fval
        self.fval = fval

        self.calc_constraints(X)

        prev_val = sum(np.clip(prev_cons_vals, 0, 100))
        cur_val = sum(np.clip(self.constr_vals, 0, 100))

        d_cons_vals = prev_val - cur_val

        if dfval > 0 and d_cons_vals >= 0 and np.all(self.constr_vals <= 0):
            self.best_fval = fval
            reward = 1
        else:
            reward = -1

        done = np.all(
            self.X == np.asarray([41, 0, 38, 31, 0, 0, 27, 38, 37, 0]))

        return X, reward, done, 'INFO'

    def worker(self):
        r = -1
        actions = []
        while r <= 0:

            action = np.random.randint(-self.step_length, self.step_length + 1,
                                       len(self.X))
            if list(action) not in actions:
                actions.append(list(action))
                s, r, d, _ = self.step(action)

        self.X = s
        if d:
            print("X: ", list(self.X))
            print(f"Weight: {self.problem.obj(self.X):2f}")

    def solve(self, problem, x0=[], maxiter=100, maxtime=60, subset_size=-1):
        """
        Solves given problem

        :param problem: problem to be solved
        :param maxiter: maximum number of iterations

        :type problem: OptimizationProblem
        :type maxiter: int

        :return: fopt, xopt

        """
        self.problem = problem
        # Srating point
        if not len(x0):
            self.X = np.zeros(len(problem.vars), dtype=int)
            self.X += int(problem.vars[0].ub)
        else:
            self.X = x0
        if subset_size == -1:
            subset_size = len(self.X)
        self.constr_vals = np.zeros(len(self.problem.cons))
        self.calc_constraints()

        X_vals = []

        for i in range(maxiter):
            r = -1
            actions = []
            d = False
            while r <= 0 and not d:
                if time.process_time() >= maxtime or d:
                    i = maxiter
                    break

                if len(actions) or i == 0:
                    action = np.random.randint(-self.step_length,
                                               self.step_length + 1,
                                               len(self.X))

                    # Choose random subset to stay still
                    subset = np.random.choice(len(self.X),
                                              len(self.X) - subset_size,
                                              replace=False)
                    action[subset] = 0

                if list(action) not in actions:
                    actions.append(list(action))
                    s, r, d, _ = self.step(action)

            if (time.process_time() < maxtime or i < maxiter) and not d:
                X_vals.append(s)
                self.X = s
                print("New Best! ", self.best_fval)
                self.best_X = s.copy()
                print(self.X)

        print("TIME: ", time.process_time(), " s")
        print("X: ", list(self.best_X))
        print(f"Weight: {problem.obj(self.best_X):2f}")

        return self.best_X.copy(), self.best_fval, X_vals


class Basinhopping(OptSolver):

    def __init__(self, step_length=1):
        self.step_length = step_length

    def step(self, x):
        print(x)
        self.X += x
        self.calc_constraints()

    def solve(self, problem, x0=[]):
        self.problem = problem
        # Srating point
        self.X = np.zeros(len(problem.vars), dtype=int)
        # self.X += problem.vars[0].ub
        self.constr_vals = np.zeros(len(self.problem.cons))
        self.calc_constraints()

        if not len(x0):
            x0 = self.X.copy()

        return basinhopping(problem.obj,
                            x0,
                            stepsize=self.step_length,
                            take_step=self.step)

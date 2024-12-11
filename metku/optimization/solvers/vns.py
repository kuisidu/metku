# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
import numpy as np
import time
import matplotlib.pyplot as plt
from ortools.linear_solver import pywraplp

from metku.optimization.solvers import *
from metku.optimization.structopt import DiscreteVariable, IndexVariable, IntegerVariable


class VNS(OptSolver):
    
    def __init__(self, step_length=1, stochastic=False, pop_size=5, maxiter=3, first_improvement=True,
                 solver='GA', solver_params=None, max_expansion=10):
        super().__init__()
        self.step_length = step_length
        self.initial_step_length = step_length
        self.stochastic = stochastic
        self.temp_X = np.array([])
        self.fopt = 1e100
        self.x0 = None
        self.prev_x = []
        self.pop_size = pop_size
        self.maxiter = maxiter
        self.space_expansion = 0
        self.max_expansion = max_expansion
        self.first_improvement = first_improvement
        self.solver = solver.upper()
        self.solver_params = solver_params
        
    def take_action(self, recursive=False):
        """
        Defines action to take
        :return:
        """


        bounds = [(var.lb, var.ub) for var in self.problem.vars]

        for var in self.problem.vars:
            if isinstance(var, DiscreteVariable):
                idx = np.argmin((np.asarray(var.values) - var.value) ** 2)
                delta = (var.ub - var.lb) / (self.step_length * len(var.values))
                lb, ub = var.lb - delta, var.ub + delta

                if idx == 0:
                    ub = max(ub, var.values[1] - var.value)
                elif idx == len(var.values) - 1:
                    lb = max(lb, var.value - var.values[idx - 1])
                else:
                    lb = max(lb, var.value - var.values[idx - 1])
                    ub = max(ub, var.values[idx + 1] - var.value)
                var.lb = lb
                var.ub = ub
                if not recursive:
                    self.prev_x.append(var.value)

            elif isinstance(var, IndexVariable):
                var.lb = max(0, var.value - self.step_length)
                var.ub = min(var.ub, var.value + self.step_length)
                if not recursive:
                    self.prev_x.append(var.value)

            else:
                delta = (var.ub - var.lb) / 10
                var.lb = max(var.lb, var.value - delta * self.step_length)
                var.ub = min(var.ub, var.value + delta * self.step_length)

                if not recursive:
                    self.prev_x.append(var.value)

        if self.space_expansion >= 10:
            return np.zeros_like(self.best_x)


        if self.x0 is None:
            self.x0 = [var.value for var in self.problem.vars]

        if type(self.solver).__name__ == 'GA':
            self.solver.counter = 0

        fopt, xopt = self.solver.solve(self.problem, maxiter=self.maxiter, x0=self.x0, plot=self.plot)

        if fopt < self.fopt and np.all(xopt != self.prev_x):
            self.fopt = fopt
            self.x0 = xopt
            self.solver.best_f, self.solver.prev_best = fopt, fopt
            self.space_expansion = 0
            for var, bound in zip(self.problem.vars, bounds):
                var.lb, var.ub = bound
            return np.asarray(xopt) - np.asarray(self.prev_x)

        else:
            if self.space_expansion >= self.max_expansion:
                print("EXPANDED OVER {} TIMES!".format(self.max_expansion))
                return np.zeros_like(self.best_x)
            print("EXPANDING SEARCH SPACE")
            self.space_expansion += 1
            for var, bound in zip(self.problem.vars, bounds):
                var.lb, var.ub = bound
            self.step_length += 1

            return self.take_action(recursive=True)

    def step(self, action):
        """
        Takes step
        """
        self.prev_x = []
        # Reset step_length
        self.step_length = self.initial_step_length
        # Take step
        self.X = np.asarray(self.X) + np.asarray(action)


        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                                self.problem.vars[i].ub)

        return self.X.copy(), 1, False, 'INFO'


    def solve(self, *args, **kwargs):

        try:
            self.plot = kwargs['plot']
        except:
            self.plot = False

        if self.solver == 'GA':
            if self.solver_params is None:
                def mutate(individual, prob=0.05, stepsize=1, multiplier=0.1):
                    """
                    Mutates individual
                    :param individual:
                    :param prob:
                    :return: mutated individual
                    """
                    for i in range(len(individual)):
                        if np.random.rand() < prob:
                            var = self.problem.vars[i]
                            if isinstance(var,
                                          (IndexVariable, IntegerVariable)):
                                step = min(var.ub - var.lb, stepsize)
                                val = np.random.randint(0, step)
                            else:
                                val = np.random.uniform(var.lb, var.ub)
                            if np.random.rand() < 0.5:
                                val *= -1
                            individual[i] += val

                    return individual
                self.solver_params = dict(
                    pop_size=self.pop_size,
                    mutation_fun=mutate,
                    mutation_kwargs={'prob': 0.05,
                                     "stepsize": 100,
                                     "multiplier": self.step_length}
                )
                self.solver_params["first_improvement"] = self.first_improvement

            self.solver = GA(**self.solver_params)


        elif self.solver == "SLP":
            self.solver = SLP(**self.solver_params)

        elif self.solver == "MISLP":
            self.solver = MISLP(**self.solver_params)


        return super().solve(*args, **kwargs)


class DiscreteVNS(OptSolver):
    """
        Variable Neighborhood Search (VNS) solver for discrete optimization problems
    """

    def __init__(self, step_length=1, subset_size=-1):
        self.step_length = int(step_length)
        self.X = None
        self.problem = None
        self.constr_vals = None
        self.fval = 10e3
        self.best_fval = 10e3
        self.subset_size = subset_size


    def shake(self):
        action = np.random.randint(-self.step_length, self.step_length + 1, len(self.X))
        x = np.clip(self.X.copy() + action, 0, self.problem.vars[0].ub)
        return x

    def first_improvement(self, x):
        actions = []
        r = -1
        while r <= 0:
            action = np.random.randint(-self.step_length, self.step_length + 1, len(x))
            if list(action) not in actions:
                actions.append(list(action))
                s, r, d, _ = self.step(action, x)

        if d:
            print(f"DONE! Obj.val: {self.best_fval}")
            print(x)

        return x


    def neighbourhood_change(self, x):
        self.X = x


    def take_action(self):
        """
        Defines action to take
        :return: action
        """
        if self.subset_size <= 0:
            self.subset_size = len(self.X)
        # Take random action in current neighborhood
        action = np.random.randint(-self.step_length,
                                   self.step_length + 1,
                                   len(self.X))

        # Choose random subset to stay still
        subset = np.random.choice(len(self.X),
                                  len(self.X) - self.subset_size,
                                  replace=False)
        action[subset] = 0

        return action

    def step(self, action, X=[]):
        """
        Takes a step
        :param action: direction of the step
        :param X: starting point of the step
        :return: state, reward, done, info
        """

        if not len(X):
            X = self.X.copy()
        prev_cons_vals = self.constr_vals.copy()
        X += action
        X = np.clip(X, 0, self.problem.vars[0].ub)

        for x, var in zip(X, self.problem.vars):
            var.substitute(int(x))

        # Calculate objective
        fval = self.problem.obj(X)
        dfval = self.best_fval - fval
        self.fval = fval

        self.calc_constraints(X)

        prev_val = sum(np.clip(prev_cons_vals, 0, 100)**2)
        cur_val = sum(np.clip(self.constr_vals, 0, 100)**2)

        d_cons_vals = prev_val - cur_val

        if dfval > 0 and d_cons_vals >= 0 and np.all(self.constr_vals <= 0):
            self.best_fval = fval
            reward = 1
        #
        # elif d_cons_vals > 0 and np.all(self.constr_vals <= 0):
        #     print(d_cons_vals)
        #     reward = 1

        # if self.best_fval > self.problem.obj(X):
        #     reward = 1
        
        else:
            reward = -1

        #done = np.all(self.constr_vals <= 0)
        done = False

        return X, reward, done, 'INFO'


    def solve(self, problem, x0=None, maxiter=None, maxtime=None, verb=None, plot=False):
        """
        Solves given problem

        :param problem: problem to be solved
        :param maxiter: maximum number of iterations

        :type problem: OptimizationProblem
        :type maxiter: int

        :return: fopt, xopt

        """
        if maxtime is not None:
            maxiter = 10000000
        elif maxiter is not None:
            maxtime = 10000000
        else:
            maxiter = 100
            maxtime = 100
        
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.ion()
            plt.show()
        
        
        self.problem = problem
        # Srating point
        if x0 is None:
            self.X = np.zeros(len(problem.vars), dtype=int)
            self.X += problem.vars[0].ub
        else:
            self.X = x0
        problem.substitute_variables(self.X)
        self.constr_vals = np.zeros(len(self.problem.cons))
        self.calc_constraints()

        start_time = time.time()
        
        for i in range(maxiter):
            r = -1
            actions = []
            while r <= 0:

                if time.time() - start_time >= maxtime:
                    print("MAXTIME")
                    i = maxiter
                    return self.X, problem.obj(self.X)
                if r != 1:
                    action = np.random.randint(-self.step_length,
                                               self.step_length + 1,
                                               len(self.X))
                if list(action) not in actions:
                    actions.append(list(action))
                    s, r, d, _ = self.step(action)
            if plot:
                self.update_plot(fig, ax)
            self.X = s
            if verb:
                print("New Best! ", problem.obj(s))
                print(self.X)

        print("TIME: ", time.time() - start_time, " s")
        print("X: ", list(self.X))
        print(f"Weight: {problem.obj(self.X):2f}")
        
        return problem.obj(self.X), self.X

        
if __name__ == '__main__':
    from metku.optimization.benchmarks import *
    
    problem = FifteenBarTruss('continuous')
    solver = VNS(step_length=1, stochastic=False, stoch_prob=[0.05, 0.9, 0.05])
    x0 = [var.ub for var in problem.vars]
    fopt, xopt = solver.solve(problem, maxiter=300, x0=x0, verb=True)
    problem(xopt)
    problem.structure.plot()

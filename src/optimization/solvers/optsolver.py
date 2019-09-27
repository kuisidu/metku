import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy.optimize import basinhopping

try:
    from src.optimization.structopt import OptimizationProblem
except:
    from optimization.structopt import OptimizationProblem



class OptSolver:
    """
        Base class for optimization problem solvers
    """

    def __init__(self):
        self.constr_vals = np.array([])
        self.X = np.array([])
        self.problem = None
        self.fvals = []
        self.xvals = []
        self.best_f = np.inf
        self.best_x = None

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

    def take_action(self):
        pass


    def step(self, action):
        """
        Takes a step
        :param action:
        :return:
        """

        self.step_length -= self.step_length * self.step_factor
        self.X += action
        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                           self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'



    def solve(self, problem, x0=None, maxiter=-1, maxtime=-1, log=False,
              min_diff=1e-5, verb=False):
        """
        Solves given problem

        :param problem:
        :param maxiter:
        :param maxtime:
        :return:
        """
        # If maxiter or maxtime aren't defined use 100 as maxiter
        if maxiter <= 0 and maxtime <= 0:
            maxiter = 100
        # If only maxtime is defined, use large maxiter val
        elif maxiter == -1:
            maxiter = 1000000000
        # If only maxiter is defined, use large maxtime val
        elif maxtime == -1:
            maxtime = 1000000000


        # Assign problem
        self.problem = problem
        # If initial starting point is/isn't defined
        if x0:
            self.X = x0
            problem.substitute_variables(x0)
        else:
            self.X = self.random_feasible_point()
            problem.substitute_variables(self.X)

        # Assign done
        done = False
        # Start iteration
        t_total = 0
        for i in range(maxiter):
            # Check time
            start = time.time()
            t_0 = time.process_time()
            if t_total >= maxtime or done:
                break
            # Save previous state
            prev_state = self.X.copy()
            # Define action to take
            action = self.take_action()
            # Take step
            state, reward, done, info = self.step(action)
            # If new state is almost same as previous
            if (np.linalg.norm(prev_state - state)/np.linalg.norm(prev_state) <= min_diff or \
                np.all(prev_state == state)) and problem.feasible:
                done = True
                break
            # Change current state
            self.X = state
            # Substitute new variables
            problem.substitute_variables(state)
            # Calculate constraints
            self.calc_constraints(self.X)
            #print(state)
            # Save best values

            fval = problem.obj(state)
            t_1 = time.process_time()
            t_total += t_1-t_0
            if verb:
                print(
                    f'\r Iteration {i + 1} / {maxiter}: Obj: {fval} Feasible: {problem.feasible} '\
                    f'max g: {max(self.constr_vals):.4f}')
            if fval < self.best_f and problem.feasible:
                self.best_f = fval
                self.best_x = state.copy()
                if verb:
                    print(self.best_x)
                    print(f"New best!: {fval:.2f} {[round(s, 2) for s in state]}")
            # Log objective vals per iteration
            if log:
                problem.num_iters += 1
                problem.fvals.append(problem.obj(self.X))
                problem.states.append(list(state))
                problem.gvals.append(list(self.constr_vals).copy())
                self.fvals.append(problem.obj(self.X))
                self.xvals.append(self.X)
            end = time.time()
            print(f"Iteration took: {end - start :.2f} s")

        return self.best_f, self.best_x

    def random_feasible_point(self):
        """
        Creates a random feasible starting point for optimization

        :return: random feasible starting point
        """
        print("Creating random feasible starting point!")

        X = [0] * self.problem.nvars()
        for i, var in enumerate(self.problem.vars):
            if self.problem.prob_type == 'continuous':
                X[i] = var.lb + (var.ub - var.lb) * np.random.rand()
            else:
                X[i] = np.random.randint(var.lb, var.ub)
        while np.any(self.calc_constraints(X) > -0.15):

            for i, var in enumerate(self.problem.vars):
                if self.problem.prob_type == 'continuous':
                    X[i] = var.lb + (var.ub - var.lb) * np.random.rand()
                else:
                    X[i] = np.random.randint(var.lb, var.ub)

            self.problem.substitute_variables(X)
            self.problem.fea()

        print("Starting point created!")

        return np.asarray(X)
        

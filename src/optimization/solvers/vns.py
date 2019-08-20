
import numpy as np
from .optsolver import OptSolver


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
            var.substitute(x)

        # Calculate objective
        fval = self.problem.obj(self.problem.vars)
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

        done = np.all(self.constr_vals <= 0)

        return X, reward, done, 'INFO'


    def solve(self, problem, maxiter=100, maxtime=100):
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
        self.X = np.zeros(len(problem.vars), dtype=int)
        self.X += problem.vars[0].ub
        self.constr_vals = np.zeros(len(self.problem.cons))
        self.calc_constraints()


        for i in range(maxiter):
            r = -1
            actions = []
            while r <= 0:
                if time.process_time() >= maxtime:
                    i = maxiter
                    break
                if r != 1:
                    action = np.random.randint(-self.step_length,
                                               self.step_length + 1,
                                               len(self.X))
                if list(action) not in actions:
                    actions.append(list(action))
                    s, r, d, _ = self.step(action)

            self.X = s
            print("New Best! ", self.best_fval)
            print(self.X)

        print("TIME: ", time.process_time(), " s")
        print("X: ", list(self.X))
        print(f"Weight: {problem.obj(self.X):2f}")
import numpy as np
import time
import matplotlib.pyplot as plt
from ortools.linear_solver import pywraplp


try:
    from src.optimization.solvers import *
    from src.optimization.structopt import DiscreteVariable, IndexVariable, IntegerVariable
except:
    from optimization.solvers import *
    from optimization.structopt import DiscreteVariable, IndexVariable



class VNS(OptSolver):
    
    def __init__(self, step_length=1, stochastic=False, stoch_vals=[-1, 0, 1],
                 stoch_prob=[0.05, 0.9, 0.05]):
        super().__init__()
        self.step_length = step_length
        self.initial_step_length = step_length
        self.stochastic = stochastic
        self.stoch_vals = stoch_vals
        self.stoch_prob = stoch_prob
        self.temp_X = np.array([])
        self.fopt = 1e100
        self.prev_x = []
        if sum(stoch_prob) != 1:
            raise ValueError(f"Sum of stoch_prob must be 1.0 !")
        
    def take_action(self, recursive=False):
        """
        Defines action to take
        :return:
        """


        bounds = [(var.lb, var.ub) for var in self.problem.vars]

        for var in self.problem.vars:
            if isinstance(var, DiscreteVariable):
                pass
            elif isinstance(var, IndexVariable):

                var.lb = max(0, var.idx - self.step_length)
                var.ub = min(var.ub, var.idx + self.step_length)
                if not recursive:
                    self.prev_x.append(var.idx)

            else:

                delta = (var.ub - var.lb) / 50

                var.lb = max(var.lb, var.value - delta * self.step_length)
                var.ub = min(var.ub, var.value + delta * self.step_length)

                if not recursive:
                    self.prev_x.append(var.value)


        fopt, xopt = self.solver.solve(self.problem, maxiter=3, plot=True)
        if fopt < self.fopt:
            self.fopt = fopt
            for var, bound in zip(self.problem.vars, bounds):
                var.lb, var.ub = bound
            return np.asarray(xopt) - np.asarray(self.prev_x)

        else:
            print("EXPANDING SEARCH SPACE")
            for var, bound in zip(self.problem.vars, bounds):
                var.lb, var.ub = bound
            self.step_length += self.initial_step_length

            return self.take_action(recursive=True)

    def step(self, action):
        """
        Takes step
        """
        print(action)
        self.prev_x = []
        # Reset step_length
        self.step_length = self.initial_step_length
        # Take step
        self.X += action

        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                                self.problem.vars[i].ub)

        return self.X.copy(), 1, False, 'INFO'


    def solve(self, *args, **kwargs):

        self.plot = kwargs['plot']

        def mutate(individual, prob=0.01, stepsize=1, multiplier=0.1):
            """
            Mutates individual
            :param individual:
            :param prob:
            :return: mutated individual
            """
            for i in range(len(individual)):
                var = self.problem.vars[i]
                if isinstance(var, (IndexVariable, IntegerVariable)):
                    step = min(var.ub - var.lb, stepsize)
                    val = np.random.randint(0, step)

                else:
                    val = np.random.uniform(var.lb, var.ub)

                if np.random.rand() < prob:
                    if np.random.rand() < 0.5:
                        val *= -1
                    individual[i] += val

            return individual


        self.solver = GA(pop_size=5,
                    mut_fun=mutate,
                    mutation_kwargs={'prob': 0.5, "stepsize": 10,
                                     "multiplier": self.step_length})
        super().solve(*args, **kwargs)


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

        prev_val = sum(np.clip(prev_cons_vals, 0, 100))
        cur_val = sum(np.clip(self.constr_vals, 0, 100))

        d_cons_vals = prev_val - cur_val

        if dfval > 0 and d_cons_vals >= 0 and np.all(self.constr_vals <= 0):
            self.best_fval = fval
            reward = 1
            
        elif d_cons_vals > 0 and np.all(self.constr_vals <= 0):
            print(d_cons_vals)
            reward = 1
        
        else:
            reward = -1

        done = np.all(self.constr_vals <= 0)

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
                print("New Best! ", self.best_fval)
                print(self.X)

        print("TIME: ", time.time() - start_time, " s")
        print("X: ", list(self.X))
        print(f"Weight: {problem.obj(self.X):2f}")
        
        return self.X, problem.obj(self.X)

        
if __name__ == '__main__':
    from src.optimization.benchmarks import *
    
    problem = FifteenBarTruss('continuous')
    solver = VNS(step_length=1, stochastic=False, stoch_prob=[0.05, 0.9, 0.05])
    x0 = [var.ub for var in problem.vars]
    fopt, xopt = solver.solve(problem, maxiter=300, x0=x0, verb=True)
    problem(xopt)
    problem.structure.plot()

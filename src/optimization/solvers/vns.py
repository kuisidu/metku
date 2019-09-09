import numpy as np
from src.optimization.solvers import OptSolver
from ortools.linear_solver import pywraplp
from src.optimization.structopt import DiscreteVariable


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
        if sum(stoch_prob) != 1:
            raise ValueError(f"Sum of stoch_prob must be 1.0 !")
        
    def take_action(self, linearized=[], X=[]):
        """
        Defines action to take
        """
        if len(linearized):
            A, B, df, fx = linearized
        else:
            if not len(X):
                X = self.X
            # Linearize problem
            A, B, df, fx = self.problem.linearize(X)
        
        if self.problem.prob_type == 'discrete':
            # Create CBC solver
            solver = pywraplp.Solver('MIP',
                               pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
        else:
            solver = pywraplp.Solver('LP',
                                 pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)
        # Number of constraints
        n = A.shape[0]
        # Number of variables
        m = A.shape[1]

        # Create continuous variables
        x = {}
        for i, var in enumerate(self.problem.vars):
            if isinstance(var, DiscreteVariable):
                idx = var.values.index(var.value)
                idx_lb = max(var.lb, idx - self.step_length)
                idx_ub = min(var.ub, idx + self.step_length)
                lb = var.values[idx_lb]
                ub = var.values[idx_ub]
            else:
                # This can be made into instance variable
                steps = 100
                step_length = (var.ub - var.lb) / steps * self.step_length
                lb, ub = var.value + np.array([-step_length, step_length])
                
            x[i] = solver.NumVar(lb,
                                    ub,
                                     str(var.target['property']) + str(i))
        # Linear Constraints
        # Ax <= B
        for i in range(n):
            solver.Add(solver.Sum([A[i, j] * x[j]
                                           for j in range(m)]) <= B[i])

        
        # Create binary variables
        discrete_vars = [var for var in
                         self.problem.vars if isinstance(var, DiscreteVariable)]
        y = {}
        for i, var in enumerate(discrete_vars):
            for j in range(var.ub):
                y[i, j] = solver.BoolVar(f'y{i},{j}')

        # Create binary constraints
        for i, var in enumerate(discrete_vars):
            # Binary constraint
            # sum(bin_vars) == 1
            solver.Add(solver.Sum([y[i, j]
                                           for j in range(var.ub)]) == 1)
            # Binary property constraint
            # Ai == sum(A[i][bin_idx])
            solver.Add(solver.Sum([y[i, j] * var.values[j]
                                           for j in range(var.ub)]) == x[i])

        # Objective
        solver.Minimize(solver.Sum(df[i] * x[i] for i in range(m)))

        # Solve
        sol = solver.Solve()

        # If solution is infeasible
        if sol == 2:
            print("Expanding search space!")
            self.step_length += 1
            return self.take_action(linearized=[A, B, df, fx])

        else:
            X = [j for i, var in enumerate(discrete_vars)
                 for j in range(var.ub) if y[i, j].solution_value() == 1.0]
    
            for i in range(len(X), len(x)):
                X.append(x[i].solution_value())
    
            X = np.asarray(X)

            return X - self.X
    
    
    def step(self, action):
        """
        Takes step
        """
        # Reset step_length
        self.step_length = self.initial_step_length
        # Stochastic push
        if self.stochastic:
            p = np.random.choice(self.stoch_vals,
                                 len(action),
                                 p=self.stoch_prob)
            action += p

        # Take step
        self.X += action

        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                                self.problem.vars[i].ub)

        return self.X.copy(), 1, False, 'INFO'


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
        

        
if __name__ == '__main__':
    from src.optimization.benchmarks import *
    
    problem = FiftyTwoBarTruss('discrete')
    solver = VNS(step_length=1, stochastic=False, stoch_prob=[0.1, 0.8, 0.1])
    x0 = [var.ub for var in problem.vars]
    fopt, xopt = solver.solve(problem, maxiter=100, x0=x0, verb=True
                              )
    problem(xopt)
    problem.structure.plot()

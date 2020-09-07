
from ortools.linear_solver import pywraplp

try:
    from metku.optimization.solvers.optsolver import OptSolver
    from metku.optimization.structopt import *
except:
    from optimization.solvers.optsolver import OptSolver
    from optimization.structopt import *


class MISLP(OptSolver):

    def __init__(self, move_limits=(0.1, 0.1), gamma=1e-3, C=2e5, beta=100):
        super().__init__()
        self.move_limits = np.asarray(move_limits)
        self.gamma = gamma
        self.C = C
        self.beta = beta

    def take_action(self):
        """
        Defines action to take
        :return: action
        """

        # Linearize problem
        A, B, df, fx = self.problem.linearize(*self.X.copy())
        # Create CBC solver
        CBC_solver = pywraplp.Solver('MISLP',
                           pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

        # Number of constraints
        n = A.shape[0]
        # Number of variables
        m = A.shape[1]

        # Create continuous variables
        x = {}
        # Save discrete variables' indices to map binary variables to
        # these relaxed variables
        discrete_ids = []
        for i, var in enumerate(self.problem.vars):
            delta = var.ub - var.lb
            lb, ub = delta * self.move_limits
            if isinstance(var, DiscreteVariable):
                discrete_ids.append(i)
                idx = np.argmin((np.asarray(var.values) - var.value)**2)

                if idx == 0:
                    ub = max(ub, var.values[1] - var.value)
                elif idx == len(var.values) -1:
                    lb = max(lb, var.value - var.values[idx - 1])
                else:
                    lb = max(lb, var.value - var.values[idx - 1])
                    ub = max(ub, var.values[idx + 1] - var.value)
            lb = max(var.value - lb, var.lb)
            ub = min(var.value + ub, var.ub)
            name = str(var.target['property']) + str(i)
            x[i] = CBC_solver.NumVar(float(lb),
                                     float(ub),
                                     name)

        beta = CBC_solver.NumVar(0, self.beta, 'Beta')
        # Linear Constraints
        # Ax <= B
        for i in range(n):
            CBC_solver.Add(CBC_solver.Sum([A[i, j] * x[j]
                                           for j in range(m)]) <= B[i] + beta)

        # Create binary variables
        discrete_vars = [var for var in
                         self.problem.vars if isinstance(var, DiscreteVariable)]
        y = {}
        for i, var in enumerate(discrete_vars):
            for j in range(len(var.values)):
                if (var.id, j) not in y.keys():
                    y[var.id, j] = CBC_solver.BoolVar(f'y{var.id},{j}')

        # Create binary constraints
        for i, var in enumerate(discrete_vars):
            # Binary constraint
            # sum(bin_vars) == 1
            CBC_solver.Add(CBC_solver.Sum([y[var.id, j]
                                           for j in range(len(var.values))]) == 1)
            # Binary property constraint
            # Ai == sum(A[i][bin_idx])
            idx = discrete_ids[i]
            CBC_solver.Add(CBC_solver.Sum([y[var.id, j] * var.values[j]
                                           for j in range(len(var.values))]) == x[idx])

        # Objective
        CBC_solver.Minimize(CBC_solver.Sum(df[i] * x[i] for i in range(m)) + self.C * beta)

        # Solve
        sol = CBC_solver.Solve()
        if sol == 2:
            print("SOLUTION INFEASIBLE")
            return np.zeros_like(self.X)

        self.beta = beta.solution_value()

        # X = [j for i, var in enumerate(discrete_vars)
        #      for j in range(len(var.values)) if y[var.id, j].solution_value() == 1.0]
        #
        #
        # print("MISLP: ", [a.solution_value() for a in x.values()])
        # for i in range(len(X), len(x)):
        #     X.append(x[i].solution_value())
        #
        # print(X)
        # # for i, var in enumerate(discrete_vars):
        # #     print(x[i].solution_value())
        # #     print(var.profiles[X[i]])
        #
        # # for i, var in enumerate(self.problem.vars):
        # #     for j in range(var.ub):
        # #         print(i, j, y[i,j].solution_value())
        #
        # X = np.asarray(X)
        #
        # print(X)
        # for i in discrete_ids:
        #     print(int(X[i]))
        #     X[i] = discrete_vars[i].values[int(X[i])]
        X = np.asarray([round(var.solution_value(), 3) for var in x.values()])
        return X - self.X

    def step(self, action):

        self.move_limits *= (1 - self.gamma)
        self.beta *= (1 - self.gamma)
        self.X += action
        # for i in range(len(self.X)):
        #     self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
        #                         self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'

if __name__ == '__main__':
    from metku.optimization.benchmarks import *
    import matplotlib.pyplot as plt

    ten_bar_DLM = [41, 0, 38, 31, 0, 0, 27, 38, 37, 0]
    # [40  0 38 34  0  0 27 38 38  0] 2497.67
    problem = TenBarTruss(prob_type='discrete')
    solver = MISLP(move_limits=[0.3, 0.3], gamma=1e-2)
    x0 = [var.ub for var in problem.vars]
    solver.solve(problem, x0=x0, maxiter=300, log=True, verb=True)
    problem(solver.X)

    # plt.plot(solver.fvals)
    # plt.show()
    problem.structure.plot()
    #problem(ten_bar_DLM)
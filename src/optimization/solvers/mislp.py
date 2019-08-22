from src.optimization.solvers.optsolver import OptSolver
from src.optimization.structopt import *
from ortools.linear_solver import pywraplp



class MISLP(OptSolver):

    def __init__(self, move_limits=[0.9, 1.1], gamma=1e-3, step_length=3):
        super().__init__()
        self.move_limits = np.asarray(move_limits)
        self.gamma = gamma
        self.step_length = step_length

    def take_action(self):
        """
        Defines action to take
        :return: action
        """
        # Linearize problem
        A, B, df, fx = self.problem.linearize(self.X.copy())
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
            if isinstance(var, DiscreteVariable):
                discrete_ids.append(i)
            lb, ub = var.value * self.move_limits
            x[i] = CBC_solver.NumVar(lb,
                                     ub,
                                     str(var.target['property']) + str(i))
        # Linear Constraints
        # Ax <= B
        for i in range(n):
            CBC_solver.Add(CBC_solver.Sum([A[i, j] * x[j]
                                           for j in range(m)]) <= B[i])

        # Create binary variables
        discrete_vars = [var for var in
                         self.problem.vars if isinstance(var, DiscreteVariable)]
        y = {}
        for i, var in enumerate(discrete_vars):
            for j in range(var.ub):
                y[i, j] = CBC_solver.BoolVar(f'y{i},{j}')

        # Create binary constraints
        for i, var in enumerate(discrete_vars):
            # Binary constraint
            # sum(bin_vars) == 1
            CBC_solver.Add(CBC_solver.Sum([y[i, j]
                                           for j in range(var.ub)]) == 1)
            # Binary property constraint
            # Ai == sum(A[i][bin_idx])
            idx = discrete_ids[i]
            CBC_solver.Add(CBC_solver.Sum([y[i, j] * var.profiles[j]
                                           for j in range(var.ub)]) == x[idx])

        # Objective
        CBC_solver.Minimize(CBC_solver.Sum(df[i] * x[i] for i in range(m)))

        # Solve
        sol = CBC_solver.Solve()

        X = [j for i, var in enumerate(discrete_vars)
             for j in range(var.ub) if y[i, j].solution_value() == 1.0]


        for i in range(len(X), len(x)):
            X.append(x[i].solution_value())

        # for i, var in enumerate(discrete_vars):
        #     print(x[i].solution_value())
        #     print(var.profiles[X[i]])

        # for i, var in enumerate(self.problem.vars):
        #     for j in range(var.ub):
        #         print(i, j, y[i,j].solution_value())

        X = np.asarray(X)

        return X - self.X

    def step(self, action):

        self.move_limits += (1 - self.move_limits) * self.gamma
        self.X += action
        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                                self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'

    def solve(self, *args, **kwargs):
        #if args[0].prob_type != 'discrete':
        #    raise ValueError("MISLP only works for discrete problems")
        super().solve(*args, **kwargs)

if __name__ == '__main__':
    from src.optimization.benchmarks import *
    import matplotlib.pyplot as plt

    ten_bar_DLM = [41, 0, 38, 31, 0, 0, 27, 38, 37, 0]
    # [40  0 38 34  0  0 27 38 38  0] 2497.67
    problem = TenBarTruss(prob_type='discrete')
    solver = MISLP(move_limits=[0.5, 2], gamma=1e-3)
    x0 = [var.ub for var in problem.vars]

    solver.solve(problem, x0=x0, maxiter=20, log=True)
    problem(solver.X)

    # plt.plot(solver.fvals)
    # plt.show()
    problem.structure.plot()
    #problem(ten_bar_DLM)
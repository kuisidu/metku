from src.optimization.solvers import OptSolver
import numpy as np


class TwoPhase:


    def __init__(self, solver1, solver2, limits=[0.9, 5]):
        super().__init__()

        if not isinstance(solver1, OptSolver)\
                and isinstance(solver2, OptSolver):
            raise ValueError("Solvers must be 'OptSolver' class solvers!")

        self.solver1 = solver1
        self.solver2 = solver2
        self.limits = np.asarray(limits)


    def discretize(self):
        """
        Discretizes relaxed problem
        :return:
        """
        problem = self.solver1.problem
        X = problem.X



    def solve(self, problem, **kwargs):
        """
        Solves given problem
        :param problem:
        :param kwargs:
        :return:
        """
        if problem.prob_type != 'continuous':
            problem.__init__(prob_type='continuous')
        print(problem.prob_type)
        print(f"Solving first phase using {type(self.solver1).__name__}")
        self.solver1.solve(problem, **kwargs)
        self.discretize()
        self.solver2.solve(problem, **kwargs)




if __name__ == '__main__':
    from src.optimization.solvers import *
    from src.optimization.benchmarks import *

    solver1 = SLP()
    solver2 = MISLP()

    solver = TwoPhase(solver1, solver2)
    problem = TenBarTruss(prob_type='discrete')
    solver.solve(problem, maxiter=10)
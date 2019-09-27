import numpy as np

try:
    from src.optimization.solvers import OptSolver
    from src.optimization.structopt import DiscreteVariable
except:
    from optimization.solvers import OptSolver
    from optimization.structopt import DiscreteVariable
    



class TwoPhase:


    def __init__(self, solver1, solver2, limits=[0.9, 5]):
        super().__init__()

        if not isinstance(solver1, OptSolver)\
                and isinstance(solver2, OptSolver):
            raise ValueError("Solvers must be 'OptSolver' class solvers!")

        self.solver1 = solver1
        self.solver2 = solver2
        self.problem = None
        self.limits = np.asarray(limits)


    def discretize(self):
        """
        Discretizes relaxed problem
        :return:
        """
        best_x = self.solver1.X
        self.problem.__init__(prob_type='discrete')
        for i, var in enumerate(self.problem.vars):
            if isinstance(var, DiscreteVariable):
                pass
            values = np.asarray(var.values)
            idx = np.argmin((values - best_x[i])**2)
            var.substitute(var.values[idx])
            ll, ul = self.limits
            lb = max(0, idx - ll)
            ub = min(len(values), idx + ul)
            var.values = var.values[ll:ub]
            
            
    def solve(self, problem, **kwargs):
        """
        Solves given problem
        :param problem:
        :param kwargs:
        :return:
        """
        self.problem = problem
        if problem.prob_type != 'continuous':
            problem.__init__(prob_type='continuous')
        print(f"Solving first phase using {type(self.solver1).__name__}")
        self.solver1.solve(problem, **kwargs)
        self.discretize()
        print(f"Discretizing first phase solution using {type(self.solver2).__name__}")
        x0 = [var.value for var in self.problem.vars]
        kwargs['x0'] = x0
        xopt, fopt = self.solver2.solve(problem, **kwargs)
            
        return xopt, fopt


if __name__ == '__main__':
    from src.optimization.solvers import *
    from src.optimization.benchmarks import *

    solver1 = SLP(move_limits=[0.75, 1.5])
    solver2 = MISLP(move_limits=[0.5, 5], gamma=1e-2)
    #solver2 = DiscreteVNS(step_length=3, subset_size=5)

    solver = TwoPhase(solver1, solver2, limits=[2, 5])
    problem = FiftyTwoBarTruss(prob_type='discrete')
    x0 = [var.ub for var in problem.vars]
    xopt, fopt = solver.solve(problem, x0=x0, maxiter=150, verb=True)
    print(xopt, fopt)
    problem(xopt)
import threading
from itertools import product

from src.optimization.result_exporter import ResultExporter


class GridSearch:

    def __init__(self, problem, solver, solver_params, problem_params=None):

        self.problem = problem
        self.solver = solver

        self.solver_grid = []
        for p in [solver_params]:
            items = sorted(p.items())
            keys, values = zip(*items)
            for v in product(*values):
                params = dict(zip(keys, v))
                self.solver_grid.append(params)

        self.problem_grid = []
        if problem_params:
            for p in [problem_params]:
                items = sorted(p.items())
                keys, values = zip(*items)
                for v in product(*values):
                    params = dict(zip(keys, v))
                    self.problem_grid.append(params)

    def solve(self, **params):
        """
        :param params:
        :return:
        """
        problem = self.problem('discrete')
        solver = self.solver(**params)

        solver.solve(problem, **self.solve_kwargs)
        best_x = solver.best_x
        solver.__init__(**params)
        solver.best_x = best_x
        exporter = ResultExporter(problem, solver)
        exporter.to_csv()

    def run(self, **kwargs):
        """
        Runs grid search
        :return:
        """
        self.solve_kwargs = kwargs
        solver_options = {
            'move_limits': [[0.75, 1.5], [0.5, 2]],
            'gamma': [1e-2, 1e-3, 1e-4]
        }

        threads = []
        for i in range(len(self.solver_grid)):
            options = self.solver_grid[i]
            # x = threading.Thread(target=self.solve,
            #                      kwargs=options)
            # threads.append(x)
            # x.start()
            # print("Starting thread: ", i)
            self.solve(**options)


if __name__ == '__main__':
    from src.optimization.benchmarks import *
    from src.optimization.solvers import *

    solver_options = {
        'step_length': [1, 2, 3],
        'stochastic': [True, False]
    }

    # MUISTA TARKISTAA TEHTÄVÄN TYYPPI SOLVE -METODISSA
    problem = TenBarTruss
    solver = VNS
    gs = GridSearch(problem, solver, solver_options)
    x0 = [var.ub for var in problem().vars]
    # MUISTA LAITTAA log=True JOTTA TULOKSET TALLENTUVAT!!
    gs.run(maxiter=5, x0=x0, log=True)

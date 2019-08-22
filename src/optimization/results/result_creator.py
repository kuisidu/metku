from src.optimization.benchmarks import *
from src.optimization.solvers import *
from src.optimization.result_exporter import ResultExporter


move_limits = [
    [0.1, 5],
    [0.3, 3],
    [0.6, 2],
    [0.75, 1.50]
]

gammas = [
    1e-1,
    1e-2,
    1e-3,
]

def initialize(move_limit, gamma):
    problem = TenBarTruss(prob_type='continuous')
    solver = SLP(move_limits=move_limit, gamma=gamma)
    return problem, solver

for move_limit in move_limits:
    for gamma in gammas:
        # Initialize problem and solver
        problem, solver = initialize(move_limit, gamma)
        x0 = [var.ub for var in problem.vars]
        solver.solve(problem, x0=x0, maxiter=1000, log=True)
        solver.move_limits = move_limit
        # Export results to csv
        exporter = ResultExporter(problem, solver)
        exporter.to_csv()
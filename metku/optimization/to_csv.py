import csv
import os

from datetime import datetime

def to_csv(problem, solver, name=""):
    """
    Saves problem's results as a csv file

    :param problem: optimized problem
    :type problem: OptimizationProblem
    :return:
    """
    if not name:
        time = datetime.now().strftime('_%Y_%m_%d')
        name = type(solver).__name__ + '_' + problem.name + time + '.csv'

    # Optimal variable values
    xopt = [round(var.value,2) for var in problem.vars]

    # Optimal result
    fopt = problem.obj(xopt)

    # Constraints' values
    con_vals = problem.eval_cons(xopt)

    # Initial point, x0
    x0 = [round(var.value,2) for var in problem.x0]

    # No. Iterations
    iters = problem.num_iters

    # No. FEM evaluations
    fem_analyses = problem.num_fem_analyses

    # Avg. subproblem iterations / time

    # States
    states = problem.states

    # Function values
    fvals = problem.fvals


    results = {'f*': fopt,
               'x*': xopt,
               'x0': x0,
               'Iterations': iters,
               'FEM Analyses': fem_analyses
               }

    write_header = not os.path.exists(name)

    with open(name, 'a', newline='') as csvfile:
        fieldnames = list(results.keys())

        writer = csv.DictWriter(csvfile, fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(results)


    #print(name, fopt, xopt, x0, iters, fem_analyses)

    print(f'{name} file created to {os.getcwd()}')


if __name__ == '__main__':
    from metku.optimization.benchmarks import *
    from metku.optimization.solvers import *
    solver = SLP(move_limits=[0.5, 1.5])
    maxiters = [10, 20, 30, 40, 50]
    for maxiter in maxiters:
        problem = ThreeBarTruss('continuous')
        solver.solve(problem, maxiter=maxiter, log=True)
        to_csv(problem, solver)



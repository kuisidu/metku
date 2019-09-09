from src.optimization.benchmarks import *
from src.optimization.solvers import *
from src.optimization.grid_search import GridSearch

solvers = [MISLP, VNS]

benchmarks = [TenBarTruss, FifteenBarTruss, FiftyTwoBarTruss]

VNS_options = {
    'step_length': [1, 2],
    'stochastic': [True],
    'stoch_prob': [[0.05, 0.9, 0.05],
                   [0.1, 0.9, 0],
                   [0.15, 0.7, 0.15]],  # This equals to stochastic=False
    'stoch_vals': [[-1, 0, 1]]
}


MISLP_options = {
    'move_limits': [[0.5, 2], [0.6, 5]],
    'gamma': [1e-2, 5e-3, 1e-3]
}

VNS_grid = []
for problem in benchmarks:
    gs = GridSearch(problem, VNS, VNS_options)
    VNS_grid.extend(gs.solver_grid)

MISLP_grid = []
for problem in benchmarks:
    gs = GridSearch(problem, MISLP, MISLP_options)
    MISLP_grid.extend(gs.solver_grid)

def run_optimization(args):
    benchmark, solver, opt = args
    for key, value in opt.items():
        opt[key] = [value]
    #print(benchmark.__name__, " ", options)
    # MUISTA TARKISTAA TEHTÄVÄN TYYPPI GIRDSEARCHIN SOLVE -METODISSA
    gs = GridSearch(benchmark, solver, opt)
    x0 = [var.ub for var in benchmark().vars]
    # MUISTA LAITTAA log=True JOTTA TULOKSET TALLENTUVAT!!
    gs.run(maxiter=200, x0=x0, log=True)
    #print(gs.solver_grid)

pairs = []

for solver in solvers:
    if solver.__name__ == 'VNS':
        #options = VNS_options
        options = VNS_grid
    elif solver.__name__ == 'MISLP':
        #options = MISLP_options
        options = MISLP_grid
    else:
        print(f"Solver {solver.__name__} options missing!")
        break

    for benchmark in benchmarks:
        for opt in options:
            pairs.append([benchmark, solver, opt])
        # print(benchmark.__name__, " ", options)
        # # MUISTA TARKISTAA TEHTÄVÄN TYYPPI GIRDSEARCHIN SOLVE -METODISSA
        # gs = GridSearch(benchmark, solver, options)
        # x0 = [var.ub for var in benchmark().vars]
        # # MUISTA LAITTAA log=True JOTTA TULOKSET TALLENTUVAT!!
        # gs.run(maxiter=500, x0=x0, log=True)

if __name__ == '__main__':
    # from multiprocessing import Pool
    # with Pool() as p:
    #     p.map(run_optimization, pairs)
    for pair in pairs:
        run_optimization(pair)

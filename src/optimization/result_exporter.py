import csv
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MaxNLocator

from src.optimization.structopt import DiscreteVariable


class ResultExporter:

    def __init__(self, problem, solver):

        self.problem = problem
        self.solver = solver

    def iteration_gif(self, name="", interval=50):
        """
        Creates a gif animation of the objective function at each iteration
        :param  name: name for the created file
        :param interval: ms between each frame
        """

        if not name:
            time = datetime.now().strftime('_%Y_%m_%d')
            name = type(self.solver).__name__ + '_' \
                   + self.problem.name + time + '.gif'
        elif '.gif' not in name:
            name += '.gif'

        fig, ax = plt.subplots()
        xdata, ydata = [], []
        fvals, = plt.plot([], [])

        def update(i):
            label = f'Iteration {i}'
            xdata.append(i)
            ydata.append(self.problem.fvals[i])

            fvals.set_data(xdata, ydata)
            ax.clear()
            ax.plot(xdata, ydata)
            ax.set_xlabel(label)
            plt.xlim(0, len(self.problem.fvals))
            plt.ylim(0, max(self.problem.fvals) + 100)
            return fvals, ax

        def init():
            xdata.clear()
            ydata.clear()
            ax.clear()

        anim = FuncAnimation(fig,
                             update,
                             frames=len(self.problem.fvals),
                             init_func=init,
                             interval=50,
                             repeat_delay=200)

        # plt.show()
        # anim.save(name)

    def iteration_jpg(self, name="", fopt=None):
        """
        Creates a jpg image of the objective function at each iteration
        """

        if not name:
            time = datetime.now().strftime('_%Y_%m_%d')
            name = type(self.solver).__name__ + '_' \
                   + self.problem.name + time
        fig, ax = plt.subplots()

        num_iters = len(self.problem.fvals)
        if fopt:
            ax.plot(np.arange(num_iters + 1),
                    [fopt] * (num_iters + 1),
                    linestyle='--',
                    linewidth=2)
        ax.plot(self.problem.fvals)
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Objective value")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        if num_iters <= 20:
            plt.xticks(np.arange(0, num_iters, 1))
        elif num_iters <= 50:
            plt.xticks(np.arange(0, num_iters, 2))
        elif num_iters <= 100:
            plt.xticks(np.arange(0, num_iters, 5))

        plt.savefig(name + '.jpg')
        plt.show()

    def to_csv(self, name=""):
        """
        Saves self.problem's results as a csv file

        :param self.problem: optimized self.problem
        :type self.problem: OptimizationProblem
        :return:
        """
        if not name:
            time = datetime.now().strftime('_%Y_%m_%d')
            name = type(
                self.solver).__name__ + '_' + self.problem.name + time + '.csv'

        # Optimal variable values
        xopt = self.solver.best_x

        # Optimal result
        if np.any(xopt):
            fopt = self.problem.obj(xopt)

            # Constraints' values
            con_vals = list(self.problem.eval_cons(xopt))
        else:
            fopt = None
            con_vals = None

        # Initial point, x0
        x0 = [round(var, 2) for var in self.problem.x0]

        # No. Iterations
        iters = self.problem.num_iters

        # No. FEM evaluations
        fem_analyses = self.problem.num_fem_analyses

        # Avg. subself.problem iterations / time

        # States
        states = self.problem.states

        # Function values
        fvals = self.problem.fvals

        # Constraint values
        gvals = self.problem.gvals

        # Upper and lower boundaries
        bounds = [(var.lb, var.ub) for var in self.problem.vars]

        # Targets
        targets = [var.target['property'] for var in self.problem.vars]

        # Objects
        objects = [var.target['objects'] for var in self.problem.vars]
        for i, obj in enumerate(objects):
            for j, mem in enumerate(obj):
                try:
                    objects[i][j] = mem.mem_id
                except:
                    objects[i][j] = mem.nid

        # Profiles and values for each discrete variable
        disc_vars = [var for var in self.problem.vars
                     if isinstance(var, DiscreteVariable)]

        profiles = [var.profiles for var in disc_vars]
        values = [var.values for var in disc_vars]

        # Solver parameters
        params = list(self.solver.__dict__.items())
        excluded_params = ['X',
                           'problem',
                           'fvals',
                           'constr_vals',
                           'xvals']
        params = [param for param in params if param[0] not in excluded_params]

        results = {
            'Feasible': self.problem.feasible,
            'f*': fopt,
            'x*': xopt,
            'g*': con_vals,
            'x0': x0,
            'Iterations': iters,
            'FEM Analyses': fem_analyses,
            'States': states,
            'Objective Values': fvals,
            'Constraint Values': gvals,
            'Variable boundaries': bounds,
            'Variable targets': targets,
            'Variable objects': objects,
            'Discrete variable profiles': profiles,
            'Discrete variable values': values
        }

        for param in params:
            results[param[0]] = param[1]

        write_header = not os.path.exists(name)

        with open(name, 'a', newline='') as csvfile:
            fieldnames = list(results.keys())

            writer = csv.DictWriter(csvfile, fieldnames, dialect='excel')
            if write_header:
                writer.writeheader()
            writer.writerow(results)

        # print(name, fopt, xopt, x0, iters, fem_analyses)

        print(f'{name} file created to {os.getcwd()}')

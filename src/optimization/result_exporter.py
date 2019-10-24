import csv
import os
import xlsxwriter
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MaxNLocator

try:
    from src.frame2d.frame2d import *
    from src.framefem.framefem import *
    from src.sections.steel import *
    
    from src.optimization.structopt import DiscreteVariable
    from src.optimization.solvers.trust_region import TrustRegionConstr
except:
    from frame2d.frame2d import *
    from framefem.framefem import *
    from sections.steel import *
    
    from optimization.structopt import DiscreteVariable
    from optimization.solvers.trust_region import TrustRegionConstr


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

        plt.savefig(name + '.png')
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
        con_vals = [f'{val:f}' for val in con_vals]

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
                if isinstance(mem, FrameMember):
                    objects[i][j] = mem.mem_id
                elif isinstance(mem, FEMNode):
                    objects[i][j] = mem.nid
                
                # Tähän joku tallennustapa poikkileikkaukselle
                # memberillä ja nodella tallennetaan niiden id
                # mikä kuvaisi hyvin poikkileikkausta?
                elif isinstance(mem, SteelSection):
                    # objects[i][j] = ???
                    pass
                
                else:
                    objects[i][j] = None

        # Profiles and values for each discrete variable
        disc_vars = [var for var in self.problem.vars
                     if isinstance(var, DiscreteVariable)]

        profiles = [var.profiles for var in disc_vars]
        values = [var.values for var in disc_vars]

        # Solver parameters
        params = list(self.solver.__dict__.items())
        problem_params = list(self.problem.__dict__.items())
        excluded_params = ['X',
                           'problem',
                           'fvals',
                           'constr_vals',
                           'xvals',
                           'name', 'vars', 'cons', 'obj', 'grad', 'hess',
                           'structure', 'profiles', 'fea_done',
                           'num_fem_analyses',
                           'states', 'gvals']
        params = [param for param in params if param[0] not in excluded_params]
        problem_params = [param for param in problem_params if param[0] not in
                          excluded_params]

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

        for param in problem_params:
            results[param[0]] = param[1]

        write_header = not os.path.exists(name)

        # with open(name, 'a', newline='') as csvfile:
        #     fieldnames = list(results.keys())
        #
        #     writer = csv.DictWriter(csvfile, fieldnames, dialect='excel')
        #     if write_header:
        #         writer.writeheader()
        #     writer.writerow(results)

        name2 = type(
                self.solver).__name__ + '_' + self.problem.name + time + \
                '_results.csv'

        write_header = not os.path.exists(name2)

        with open(name2, 'a', newline='') as csvfile:
            fieldnames = ["Feasible", "Iterations", "f*", "L", "Lpi", "Fy",
                          "Fx", "lcr", "buckling_z", "LT_buckling", "symmetry",
                          "top_flange_class", "bottom_flange_class",
                          "web_class"]
            h, tw, bf, tf = xopt
            results2 = {}
            for fieldname in fieldnames:
                results2[fieldname] = results[fieldname]

            results2["h"] = h
            results2["tw"] = tw
            results2["bf"] = bf
            results2["tf"] = tf
            fieldnames.extend(["h", "tw", "bf", "tf"])

            all_cons = {
                # 'Buckling_z 0': '',
                'Com_compression_bending_y 0': '',
                # 'Buckling_y 0': '',
                'Com_compression_bending_z 0': '',
                'LT_buckling 0': '',
                'Web in class 3': '',
                'Top flange in class 3': '',
                'Web in class 2': '',
                'Top flange in class 2': '',
                'Web height at least 50': ''
            }

            for con, g in zip(self.problem.cons, con_vals):
                all_cons[con.name] = g
            fieldnames.extend(list(all_cons.keys()))
            results2.update(all_cons)

            writer = csv.DictWriter(csvfile, fieldnames, dialect='excel')
            if write_header:
                writer.writeheader()
            writer.writerow(results2)

        # print(name, fopt, xopt, x0, iters, fem_analyses)

        print(f'{name2} file created to {os.getcwd()}')


    def csv_to_excel(self):
        """
        Saves self.problem's results as an excel file

        :param self.problem: optimized self.problem
        :type self.problem: OptimizationProblem
        :return:
        """

        time = datetime.now().strftime('_%Y_%m_%d')

        csvfile = type(
            self.solver).__name__ + '_' + self.problem.name + time + \
                '_results.csv'

        workbook = xlsxwriter.Workbook(csvfile[:-4] + '.xlsx', {
            'tmpdir': 'C:/Users/Victoria/Google Drive/Koulu/Diplomityö/Tuloksia',
            'strings_to_numbers': True,
            'strings_to_formulas': True,
            'strings_to_urls': True,
            'nan_inf_to_errors': True})

        worksheet = workbook.add_worksheet()
        worksheet.freeze_panes(1, 0)  # Freeze the first row.

        with open(csvfile, 'rt', encoding='utf8') as f:
            reader = csv.reader(f)
            for r, row in enumerate(reader):
                for c, col in enumerate(row):
                    worksheet.write(r, c, col)

        # Add a bold format to use to highlight cells.
        bold = workbook.add_format({'bold': True})

        # Add a number format for cells with results
        # (weight/cross section measures).
        result_format = workbook.add_format({'num_format': '0.00'})
        worksheet.set_column('C:C', None, result_format)
        worksheet.set_column('O:R', None, result_format)

        # Add a number format for cells with constraints.
        con_format = workbook.add_format({'num_format': '0.000'})
        worksheet.set_column('S:Z', None, con_format)

        # Add a number format for cells with inputs/source.
        input_format = workbook.add_format({'num_format': '0'})
        worksheet.set_column('D:G', None, input_format)

        # Add a colour format for cells.
        purple_color_format = workbook.add_format({'bg_color': '#CCCCFF'})
        yellow_color_format = workbook.add_format({'bg_color': '#FFFF99'})
        red_color_format = workbook.add_format({'bg_color': '#FF8080'})
        green_color_format = workbook.add_format({'bg_color': '#CCFFCC'})
        grey_color_format = workbook.add_format({'bg_color': '#C0C0C0'})

        worksheet.conditional_format('A2:A145', {
            'type': 'text',
            'criteria': 'containing',
            'value': 'False',
            'format': red_color_format
        })
        worksheet.conditional_format('A2:A145', {
            'type': 'text',
            'criteria': 'containing',
            'value': 'True',
            'format': green_color_format
        })
        worksheet.conditional_format('C2:C145', {
            'type': 'no_blanks',
            'format': purple_color_format
        })
        worksheet.conditional_format('O2:R145', {
            'type': 'no_blanks',
            'format': purple_color_format
        })
        worksheet.conditional_format('S2:Z145', {
            'type': 'blanks',
            'format': grey_color_format
        })
        worksheet.conditional_format('S2:Z145', {
            'type': 'cell',
            'criteria': '>',
            'value': 0,
            'format': red_color_format
        })
        worksheet.conditional_format('S2:Z145', {
            'type': 'cell',
            'criteria': 'between',
            'minimum': -0.01,
            'maximum': 0,
            'format': yellow_color_format
        })
        worksheet.conditional_format('S2:Z145', {
            'type': 'cell',
            'criteria': '<',
            'value': -0.01,
            'format': green_color_format
        })

        """
        # Adjust the column width.
        worksheet.set_column(1, 1, 15)
        """

        workbook.close()

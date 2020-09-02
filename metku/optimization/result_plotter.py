import ast
import os
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
import pandas as pd

from metku.frame2d.frame2d import *
from matplotlib.animation import FuncAnimation

class ResultPlotter:

    def __init__(self, filename):

        self.filename = filename
        self.df = pd.read_csv(filename, sep=',')
        numeric_columns = ['f*', 'Iterations', 'FEM Analyses']
        self.df[numeric_columns] = self.df[numeric_columns].apply(
            pd.to_numeric)

    def iteration_gif(self, idx=0, scaler=0, interval=100, save=False):
        """
        Creates a gif animation of the objective function at each iteration
        :param  df: pandas dataframe
        :param interval: ms between each frame
        """
        fvals = self.df['Objective Values'].values[idx]
        fvals = ast.literal_eval(fvals)

        gvals = self.df['Constraint Values'].values[idx]
        gvals = ast.literal_eval(gvals)
        gvals = [max(gval) for gval in gvals]

        states = self.df['States'][idx]
        states = ast.literal_eval(states)

        fig, (ax3, ax1, ax2) = plt.subplots(3, 1,
                                            constrained_layout=True)
        title = ""
        cols = self.df.columns[17:]
        vals = self.df[cols].values[idx]
        for i in range(len(cols)):
            title += cols[i] + ": " +str(vals[i])+ "\n "
        fig.suptitle(title)

        if not scaler:
            scaler = int(len(fvals)/5)
        best_fval = [np.inf]

        def update(i):
            gval = gvals[i]
            num_vals = len(fvals)
            j = max(0, i - scaler)
            fval = fvals[i]
            if fval < best_fval[0] and gval < 0:
                best_fval[0] = fval
            # Objective value
            # Line chart
            ax1.clear()

            ax1.plot(np.arange(j, i), fvals[j:i])
            ax1_label = f'Iteration {i}/{num_vals} Best: {best_fval[0]:.2f}'
            ax1.set_xlabel(ax1_label)
            ax1.set_xlim(j, min(num_vals, i + scaler))
            ax1.set_ylim(min(fvals[j:]) - 100, max(fvals[j:]) + 100)

            # Constraint value
            # Bar chart
            if gval > 0:
                ax2.bar(i, gval, color='r')
            else:
                ax2.bar(i, gval, color='g')
            ax2.set_xlim(j, min(num_vals, i + scaler))
            ax2.set_ylim(min(gvals[j:]) - 0.05, max(gvals[j:]) + 0.05)
            ax2.set_xlabel("Max Constraint val: " + str(round(gval, 4)))

            # States
            # Matrix plot
            # TODO: What if problem is continuous?
            state = states[i]
            ax3.clear()


            one_hot = np.eye(64)[state]
            if gval >= 0:
                cmap = "Reds_r"
            else:
                cmap = "Greens_r"
            ax3.set_yticks(np.arange(0, len(state)))
            ax3.set_yticks([y - 0.5 for y in ax3.get_yticks()][1:])
            ax3.matshow(one_hot, cmap=cmap)
            ax3.grid(False)


            return ax1, ax2, ax3

        def init():
            ax1.clear()
            ax2.clear()
            ax3.clear()

        anim = FuncAnimation(fig,
                             update,
                             frames=len(fvals),
                             init_func=init,
                             interval=interval,
                             repeat_delay=500)
        if save:
            filename = self.filename[:-4] + "_idx_" + str(idx) + ".mp4"
            anim.save(filename)
            print(f"Animation saved as: {filename} to {os.getcwd()}")
        else:
            plt.show()

    def structure_gif(self, problem, idx=0, save=False):
        """
        Animates how states change during iteration
        :return:
        """
        states = self.df['States'][idx]
        states = ast.literal_eval(states)
        gvals = self.df['Constraint Values'].values[idx]
        gvals = ast.literal_eval(gvals)
        gvals = [max(gval) for gval in gvals]
        fvals = self.df['Objective Values'].values[idx]
        fvals = ast.literal_eval(fvals)
        best_f = [np.inf]
        best_x = [states[0]]
        best_g = [gvals[0]]

        fig, ax = plt.subplots()

        def update(i):
            label = f'Iteration {i}'
            state = states[i]
            gval = gvals[i]
            fval = fvals[i]
            if gval > 0:
                color = 'r'

            else:
                color = 'g'
                if fval < best_f[0]:
                    best_f[0] = fval
                    best_x[0] = state
                    best_g[0] = gval

            problem.substitute_variables(state)
            ax.clear()
            ax.set_xlim(-500, problem.structure.L + 500)
            if problem.name == "FifteenBarTruss":
                ax.set_ylim(-1000, 4000)
            if i == len(states) - 1:
                problem.substitute_variables(best_x[0])
                ax.set_xlabel(f"Objective: {best_f[0]:.2f} "
                              f" Max g: {best_g[0]:.4f}")
                color = 'g'

                for var in problem.vars:
                    for obj in var.target['objects']:
                        if isinstance(obj, FrameMember):
                            (x1, y1), (x2, y2) = obj.coordinates
                            ind = var.profiles.index(var.value)
                            ind = max(0.2, ind / 5)
                            ax.plot([x1, x2], [y1, y2],
                                    linewidth=ind, color=color)
            else:
                ax.set_xlabel(f"Iteration {i} / {len(fvals)} "
                              f"Objective: {fvals[i]:.2f} "
                              f" Max g: {gval:.4f}")

                for var in problem.vars:
                    for obj in var.target['objects'] :
                        if isinstance(obj, FrameMember):
                            (x1, y1), (x2, y2) = obj.coordinates
                            idx = var.profiles.index(var.value)
                            idx = max(0.2, idx/5)
                            ax.plot([x1, x2], [y1, y2],
                                    linewidth=idx, color=color)





            return ax

        anim = FuncAnimation(fig,
                             update,
                             frames=len(states),
                             interval=50,
                             repeat_delay=5000)


        if save:
            filename = self.filename[:-4] + "_idx_" + str(idx) + ".mp4"
            anim.save(filename)
            print(f"Animation saved as: {filename} to {os.getcwd()}")
        else:
            plt.show()


    def optimal_results(self):
        """
        Plots objective values as barchart
        :return:
        """
        # Values from csv
        obj_vals = self.df['f*'].values
        min_val = min(obj_vals)
        max_val = max(obj_vals)
        mean = np.mean(obj_vals)
        std = np.std(np.nan_to_num(obj_vals))
        feasible = self.df['Feasible'].values
        feasible_vals = obj_vals[feasible]
        if np.any(feasible_vals):
            min_feasible = min(feasible_vals)
        else:
            min_feasible = np.inf
        legend_values = [[], []]
        X = []
        for i in range(len(obj_vals)):
            if self.df['Feasible'][i]:
                color = 'tab:green'
                if np.isclose(obj_vals[i], min_feasible):
                    color = 'tab:blue'
            else:
                color = 'tab:red'
            plt.bar(i, obj_vals[i], color=color)
            X.append(i)
        plt.ylim(min_val - std , max_val + std)
        ticks = np.arange(np.round(min_val - std, -2) , max_val + std, 100)
        #plt.yticks(ticks)
        plt.xticks(np.arange(0, len(obj_vals)))
        plt.ylabel("Objective Value")
        plt.show()


    def table(self):
        """
        Plots a table of results
        :return:
        """
        columns = ['Feasible',
                   'f*',
                   'x*',
                   'x0',
                   'Iterations',
                   'FEM Analyses']

        columns.extend(self.df.columns[17:])

        print(self.df[columns].to_string())

        # rows = ['gamma']
        # gammas = self.df['gamma'].values
        # move_limits = self.df['move_limits'].values
        #
        # rows = []
        # for i in range(len(gammas)):
        #     gamma = gammas[i]
        #     move_limit = move_limits[i]
        #     rows.append(f'{i}: {move_limit} {gamma}')

        #
        # fig = go.Figure(data=[go.Table(
        #     header=dict(values=columns),
        #     cells=dict(values=self.df[columns].values.T))])
        # fig.show()

        # plt.table(
        #     cellText=self.df[columns].values,
        #     rowLabels=rows,
        #     colLabels=columns,
        #     loc='center'
        # )
        # plt.show()


if __name__ == '__main__':
    filename = "results/VNS_TenBarTruss_2019_08_28.csv"
    #filename = "benchmarks/MISLP_FiftyTwoBarTruss_2019_08_21.csv"

    plotter = ResultPlotter(filename)
    # plotter.optimal_results()
    # plotter.iteration_gif(idx=3, scaler=10, save=False)
    # plotter.table()
    # from metku.optimization.benchmarks import *
    # problem = TenBarTruss()
    # plotter.structure_gif(problem, idx=14, save=True)
    print(plotter.df[['Feasible', 'f*']].to_latex())
    """
    HYVIÄ ESIMERKKEJÄ
    
    VNS Oskillaatio kahden arvon välissä
    # filename = "results/VNS_FiftyTwoBarTruss_2019_08_23.csv"
    # plotter = ResultPlotter(filename)
    # plotter.iteration_gif(idx=8, scaler=25)
    
    VNS "kuolee": juuttuu linearisoidun alitehtävän ratkaisuun,
    joka on oikeasti epäkäypä
    # filename = "results/VNS_FiftyTwoBarTruss_2019_08_23.csv"
    # plotter = ResultPlotter(filename)
    # plotter.iteration_gif(idx=2, scaler=10)
    
    
    """
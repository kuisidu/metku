import ast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation

class ResultPlotter:

    def __init__(self, filename):

        self.df = pd.read_csv(filename, sep=',')
        numeric_columns = ['f*', 'Iterations', 'FEM Analyses']
        self.df[numeric_columns] = self.df[numeric_columns].apply(
            pd.to_numeric)

    def iteration_gif(self, idx=0, scaler=0, interval=50):
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
                                           gridspec_kw={'hspace': 1})
        if not scaler:
            scaler = int(len(fvals)/5)

        def update(i):

            gval = gvals[i]
            num_vals = len(fvals)
            j = max(0, i - scaler)

            # Objective value
            # Line chart
            ax1.clear()
            ax1.plot(np.arange(j, i), fvals[j:i])
            ax1_label = f'Iteration {i}/{num_vals}'
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
                             interval=50,
                             repeat_delay=500)
        plt.show()

    def state_gif(self, idx=0):
        """
        Animates how states change during iteration
        :return:
        """
        states = self.df['States'][idx]
        states = ast.literal_eval(states)

        fig, ax = plt.subplots()

        def update(i):
            label = f'Iteration {i}'
            state = states[i]

            return ax

        anim = FuncAnimation(fig,
                             update,
                             frames=len(states),
                             interval=50,
                             repeat_delay=200)
        plt.show()


    def optimal_results(self):
        """
        Plots objective values as barchart
        :return:
        """
        # Values from csv
        obj_vals = self.df['f*'].values
        max_limits = self.df['move_limits'].values
        gammas = self.df['gamma'].values
        min_val = min(obj_vals)
        max_val = max(obj_vals)
        mean = np.mean(obj_vals)
        std = np.std(obj_vals)
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
                   'Iterations',
                   'FEM Analyses']
        rows = ['gamma']
        gammas = self.df['gamma'].values
        move_limits = self.df['move_limits'].values

        rows = []
        for i in range(len(gammas)):
            gamma = gammas[i]
            move_limit = move_limits[i]
            rows.append(f'{i}: {move_limit} {gamma}')



        plt.table(
            cellText=self.df[columns].values,
            rowLabels=rows,
            colLabels=columns,
            loc='center'
        )
        plt.show()


if __name__ == '__main__':
    filename = "MISLP_ThreeBarTruss_2019_08_22.csv"
    #filename = "benchmarks/MISLP_FiftyTwoBarTruss_2019_08_21.csv"

    plotter = ResultPlotter(filename)
    plotter.optimal_results()
    #plotter.iteration_gif(idx=1, scaler=10)
    plotter.table()

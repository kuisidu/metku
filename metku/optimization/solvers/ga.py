""" Genetic algorithm 
    Solver for optimization problems

"""

try:
    from metku.optimization.solvers import OptSolver
    from metku.optimization.structopt import Variable, OptimizationProblem, \
        IndexVariable, DiscreteVariable, IntegerVariable
except:
    from optimization.solvers import OptSolver
    from optimization.structopt import Variable, OptimizationProblem, \
        IndexVariable, DiscreteVariable, IntegerVariable

import random

import matplotlib.pyplot as plt
import numpy as np
from deap import base, creator, tools

import time

class GA(OptSolver):

    def __init__(self, pop_size=10, mut_rate=0.05, cx_rate=0.8, penalty=None,
                 mut_fun=tools.mutShuffleIndexes,
                 mutation_kwargs={'indpb': 0.05},
                 cx_fun=tools.cxTwoPoint,
                 cx_kwargs={},
                 first_improvement=False,
                 best_f=1e100):
        """ Constructor
            :param pop_size: population size
            :mut_rate: mutation rate
            :cx_rate: crossover rate
            :penalty: penalty function
            :mut_fun: mutation function
            :cx_fun: cross over function
            :cx_kwargs: arguments for the cross over function
            :first_improvement: if True, the iteration is stopped when a first better solution
            than the initial solution has been found. If False, the iteration is continued
            until no improvement has been made in the last 10 rounds.
            :best_f: initial best known value of the objective function
        
        """
        super().__init__()

        """ Toolbox contains basic evolutionary operators """
        self.toolbox = base.Toolbox()
        self.pop_size = pop_size
        self.mut_rate = mut_rate
        self.cx_rate = cx_rate
        self.mut_fun = mut_fun
        self.mutation_kwargs = mutation_kwargs
        self.cx_kwargs = cx_kwargs
        self.cx_fun = cx_fun
        self.plot = False
        self.best_f = best_f
        self.prev_best = best_f
        self.first_improvement = first_improvement
        self.best_ind = None
        if penalty is None:
            def penalty(constr_vals):
                pos_vals = np.clip(constr_vals, 0, 100)
                return 2e6 * sum(pos_vals ** 2)

            self.penalty = penalty
        else:
            self.penalty = penalty


        self.counter = 0

    def take_action(self):
        pass

    def step(self, action):
        pass

    def eval_individual(self, ind):
        """
        Evaluates one individual

        :param ind: 'individual', i.e. design variable values
        :return: objective value of optimization problem
        """
        self.problem.substitute_variables(ind)
        X = [round(var.value, 3) for var in self.problem.vars]
        
        """ Evaluate objective function """
        obj_val = self.problem.obj(X)
        """ Evaluate constraints """
        constr_vals = self.problem.eval_cons(X)
        """ Apply penalty on constraints """
        penalty_val = self.penalty(constr_vals)
        
        """ Fitness function equals objective function value + penalty """
        val = obj_val + penalty_val
        if self.plot:
            self.update_plot(self.fig, self.ax)


        # print("GAPS: ", [j.g1 for j in self.problem.structure.joints.values()])

        """ If a new best solution has been found, store it """
        if val < self.best_f:
            for i, x in enumerate(X):
                ind[i] = x
            self.best_ind = ind
            self.prev_best = self.best_f
            best_x = [var.value for var in self.problem.all_vars]
            self.best_x = best_x.copy()
            self.best_f = val
            self.counter = 0

        # NOTE! Need to return at least two values, hence the comma
        return val,  # <- DO NOT REMOVE THIS COMMA

    def create_fitness(self):

        """ Create Fitness class called FitnessMin 
            The weight -1.0 corresponds to minimization problem
        """
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        """ Create class Individual that inherits the list class and
            has the FitnessMin as its fitness attribute.
        """
        creator.create("Individual", list, fitness=creator.FitnessMin)

    def attribute_generator(self, counter_list=[0]):
        """
        Generates one gene attribute
        :param counter_list: list used for counting how
                             many times function has been called
        :return:
        """
        idx = counter_list[0] % len(self.problem.vars)
        counter_list[0] = counter_list[0] + 1

        var = self.problem.vars[idx]

        """ Creates random integer value for index or integer variables """
        if isinstance(var, (IndexVariable, IntegerVariable)):
            val = np.random.randint(var.lb, var.ub)

        else:
            val = np.random.uniform(var.lb, var.ub)

        return val

    def register(self):
        """ Registers a set of functions to be used in the algorithm 
            the arguments for the register method are:
                1. alias for the function to be registered
                2. the actual method to registered
                3. arguments to be passed to the aliased method
        
        """
        self.toolbox.register("attribute", self.attribute_generator)

        self.toolbox.register("individual", tools.initRepeat,
                              creator.Individual,
                              self.toolbox.attribute, len(self.problem.vars))

        # Evaluation function
        self.toolbox.register("evaluate", self.eval_individual)
        # Crossover operator
        self.toolbox.register("mate", self.cx_fun, **self.cx_kwargs)
        # Mutation
        self.toolbox.register("mutate", self.mut_fun, **self.mutation_kwargs)

        # operator for selecting individuals for breeding the next
        # generation: each individual of the current generation
        # is replaced by the 'fittest' (best) of three individuals
        # drawn randomly from the current generation.
        self.toolbox.register("select", tools.selBest, k=2)

        # Define the population to be a list of individuals
        self.toolbox.register("population", tools.initRepeat,
                              list, self.toolbox.individual)

    def new_population(self, offspring):
        """
        Creates a new population from offsprings by shuffling them
        :param offspring: list of offsprings
        :return: new population
        """
        idx = np.random.randint(len(offspring)-1)
        self.problem.substitute_variables(self.best_x)
        values = [var.value for var in self.problem.vars]
        for i, val in enumerate(offspring[idx]):
            offspring[idx][i] = values[i]
        multiplier = self.pop_size // len(offspring) + 1
        new_population = list(offspring) * multiplier
        np.random.shuffle(new_population[:self.pop_size])

        if self.best_ind is not None:
            new_population.append(self.best_ind)

        return new_population


    def solve(self, problem, x0=None, maxiter=100, maxtime=-1, log=False,
              min_diff=1e-5, verb=False, plot=False):


        
        self.toolbox = base.Toolbox() # THIS COMMAND IS ALREADY IN THE CONSTRUCTOR!
        self.problem = problem
        """ Creates the fitness function 
            but the create_fitness method does not return anything!
        """
        self.create_fitness()
        self.register()

        if plot and not self.plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.ion()
            plt.show()
            self.plot = True
            self.fig = fig
            self.ax = ax

        """ Create population """
        pop = self.toolbox.population(n=self.pop_size)
        
        """ Evaluate fitness """
        fitnesses = list(map(self.toolbox.evaluate, pop))

        # CXPB  is the probability with which two individuals
        #       are crossed
        #
        # MUTPB is the probability for mutating an individual
        CXPB, MUTPB = self.cx_rate, self.mut_rate

        # Main loop
        it = 0
        prev_val = None
        while it < maxiter:
            it += 1
            problem.num_iters += 1
            print(f"Iteration {it}/{maxiter}")
            # Select the next generation individuals
            offspring = self.toolbox.select(pop)
            new_pop = self.new_population(offspring)

            # Clone the selected individuals
            offspring = list(map(self.toolbox.clone, new_pop))

            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                # cross two individuals with probability CXPB
                if random.random() < CXPB:
                    self.toolbox.mate(child1, child2)

                    # fitness values of the children
                    # must be recalculated later
                    del child1.fitness.values
                    del child2.fitness.values
                for mutant in offspring[:-1]:
                    # mutate an individual with probability MUTPB
                    if random.random() < MUTPB:
                        self.toolbox.mutate(mutant)
                        del mutant.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if
                           not ind.fitness.valid]
            fitnesses = map(self.toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            # The population is entirely replaced by the offspring
            pop[:] = offspring
            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness.values[0] for ind in pop]
            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x * x for x in fits)
            std = abs(sum2 / length - mean ** 2) ** 0.5

            if verb:
                print("VERB: ", self.best_x)
                print(f"Min fitness: {min(fits)} \n"
                      f"Obj: {problem.obj(self.best_x)} \n"
                      f"Max con: {max(problem.eval_cons(self.best_x)):.3f}  "
                      f"{problem.cons[np.argmax(problem.eval_cons(self.best_x))].name}")
            problem.substitute_variables(self.best_x)

            if plot:
                self.update_plot(self.fig, self.ax)

            if self.first_improvement and self.best_f < self.prev_best:
                return self.best_f, self.best_x

            if self.counter > (10 - 3):
                print(f"The result has not improved in the last 10 generations."
                      f"\n Returning current best values.")
                return self.best_f, self.best_x

            elif it > 0 and prev_val == min(fits):
                self.counter += 1

            else:
                prev_val = min(fits)

            if log:
                problem.num_iters += 1
                problem.fvals.append(self.best_f)
                problem.states.append(list(pop[np.argmin(fits)]))
                problem.gvals.append(max(self.constr_vals))
                self.fvals.append(self.best_f)
                self.xvals.append(self.best_x)

        # best_idx = np.argmin(fits)
        # xopt = pop[best_idx]
        # print("XOPT: ", xopt)
        # fopt = fits[best_idx]

        return self.best_f, self.best_x


if __name__ == "__main__":
    from metku.optimization.benchmarks import *

    problem = TenBarTruss('index')
    solver = GA(pop_size=100)
    solver.solve(problem, maxiter=50)

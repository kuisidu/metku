
try:
    from src.optimization.solvers import OptSolver
    from src.optimization.structopt import Variable, OptimizationProblem,\
        IndexVariable, DiscreteVariable, IntegerVariable
except:
    from optimization.solvers import OptSolver
    from optimization.structopt import Variable, OptimizationProblem,\
        IndexVariable, DiscreteVariable, IntegerVariable

import numpy as np
import random

from deap import base
from deap import creator
from deap import tools


class GA(OptSolver):

    def __init__(self, pop_size=10, mut_rate=0.05, penalty=None, mutate=None,
                 mutation_kwargs={'indpb': 0.05}):
        super().__init__()

        self.toolbox = base.Toolbox()
        self.pop_size = pop_size
        self.mut_rate=mut_rate
        self.mutation_kwargs = mutation_kwargs
        if penalty is None:
            def penalty(constr_vals):

                pos_vals = np.clip(constr_vals, 0, 100)
                return 2e5 * sum(pos_vals**2)

            self.penalty = penalty
        else:
            self.penalty = penalty

        if mutate is None:
            self.mutate = tools.mutShuffleIndexes
        else:
            self.mutate = mutate



    def take_action(self):
        pass

    def step(self, action):
        pass


    def eval_individual(self, X):
        """
        Evaluates one individual

        :param X: one optimization problem solution
        :return: objective value of optimization problem
        """
        self.problem.substitute_variables(X)
        obj_val = self.problem.obj(X)
        constr_vals = self.problem.eval_cons(X)
        penalty_val = self.penalty(constr_vals)
        # NOTE! Need to return at least two values, hence the comma
        return obj_val + penalty_val,


    def create_fitness(self):

        creator.create("FitnessMin", base.Fitness, weights=(-1,0))
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

        if isinstance(var, (IndexVariable, IntegerVariable)):
            val = np.random.randint(var.lb, var.ub)

        else:
            val = np.random.uniform(var.lb, var.ub)

        return val



    def mutate(self, individual, prob=0.01, stepsize=1, multiplier=0.1):
        """
        Mutates individual
        :param individual:
        :param prob:
        :return:
        """

        for i in range(len(individual)):
            var = self.problem.vars[i]
            if isinstance(var, (IndexVariable, IntegerVariable)):
                val = np.random.randint(0, stepsize)


            else:
                val = np.random.uniform(0, multiplier * var.value)

            if np.random.rand() < prob:
                if np.random.rand() < 0.5:
                    val *= -1
                individual[i] += val

        return individual


    def register(self):

        self.toolbox.register("attribute", self.attribute_generator)

        self.toolbox.register("individual", tools.initRepeat, creator.Individual,
                              self.toolbox.attribute, len(self.problem.vars))


        # Evaluation function
        self.toolbox.register("evaluate", self.eval_individual)
        # Crossover operator
        self.toolbox.register("mate", tools.cxTwoPoint)
        # Mutation
        self.toolbox.register("mutate", self.mutate, **self.mutation_kwargs)

        # operator for selecting individuals for breeding the next
        # generation: each individual of the current generation
        # is replaced by the 'fittest' (best) of three individuals
        # drawn randomly from the current generation.
        self.toolbox.register("select", tools.selTournament, tournsize=3)

        # Define the population to be a list of individuals
        self.toolbox.register("population", tools.initRepeat,
                              list, self.toolbox.individual)


    def solve(self, problem, x0=None, maxiter=100, maxtime=-1, log=False,
              min_diff=1e-5, verb=False):

        self.problem = problem
        self.create_fitness()
        self.register()

        pop = self.toolbox.population(n=self.pop_size)
        fitnesses = list(map(self.toolbox.evaluate, pop))

        # CXPB  is the probability with which two individuals
        #       are crossed
        #
        # MUTPB is the probability for mutating an individual
        CXPB, MUTPB = 0.5, 0.05

        # Main loop
        it = 0
        while it < maxiter:
            it += 1
            print(f"Iteration {it}/{maxiter}")
            # Select the next generation individuals
            offspring = self.toolbox.select(pop, len(pop))
            # Clone the selected individuals
            offspring = list(map(self.toolbox.clone, offspring))
            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                # cross two individuals with probability CXPB
                if random.random() < CXPB:
                    self.toolbox.mate(child1, child2)

                    # fitness values of the children
                    # must be recalculated later
                    del child1.fitness.values
                    del child2.fitness.values

                for mutant in offspring:
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

            print("  Min %s" % min(fits))
            print("  Max %s" % max(fits))
            print("  Avg %s" % mean)
            print("  Std %s" % std)

        best_idx = np.argmin(fits)
        best_x = pop[best_idx]
        print(f"Best: {best_x} obj: {self.problem.obj(best_x)}")
        self.problem(best_x)
        self.problem.structure.plot()


if __name__ == "__main__":

    from src.optimization.benchmarks import *
    problem = TenBarTruss('index')
    solver = GA(pop_size=100)
    solver.solve(problem, maxiter=50)


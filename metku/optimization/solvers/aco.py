# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
'''
Ant colony optimization for truss structures.
Author: Olli Luukkonen

Limitations:
- Supports only truss structures.
- Works only if variable type is discrete. Does NOT work for continuous or index variables.
- Results are not perfect, for example pso.py gives better results than this algorithm.
'''


from metku.optimization.solvers import OptSolver
from metku.frame2d import simple_truss
from metku.optimization.structopt import Variable, OptimizationProblem, \
    IndexVariable, DiscreteVariable, IntegerVariable
from metku.optimization.solvers.optsolver import OptSolver, G_TOL

import numpy as np
import warnings

warnings.filterwarnings("ignore")

class AntColonyOptimizer(OptSolver):
    def __init__(self, pop_size=100, rho=1, alpha=1,
                                beta=0.5, xi=0.75, number_of_top_ants=4,
                                penalty_factor=10000):

        super().__init__()

        self.number_of_ants = pop_size

        # Factor for global pheromone update
        self.rho = rho

        self.penalty = None

        self.pheromone_matrix = None
        self.heuristic_matrix = None
        self.probability_matrix = None

        self.alpha = alpha
        self.beta = beta

        # Factor for local pheromone update
        self.xi = xi

        self.number_of_top_ants = number_of_top_ants

        # Factor which increases penalty value
        self.penalty_factor = penalty_factor

        self.design_variables = []
        self.best_fitness = np.inf
        self.best_x = None
        self.min_f = None

    def get_min_f(self):
        '''
        :return: Minimum objective function value
        '''
        min_f = 0
        for member in self.problem.structure.members.values():
            min_f += self.calculate_member_mass(member, min(self.design_variables))
        return min_f

    def calculate_member_mass(self, member, cross_section_area):
        return member.length * cross_section_area * member.rho

    def initialize_heuristic_matrix(self):
        self.heuristic_matrix = np.full((self.i, self.j), 0, dtype=np.float64)
        i = 0
        for member in self.problem.structure.members.values():
            for j in range(self.j):
                self.heuristic_matrix[i][j] = 1 / self.calculate_member_mass(member, self.design_variables[j])
            i += 1

    def update_probabilities(self):
        self.probability_matrix = (self.pheromone_matrix ** self.alpha) * (
                self.heuristic_matrix ** self.beta)

    def initialize(self):
        '''
        Initialize parameters and matrices
        '''
        self.design_variables = self.problem.vars[0].values
        number_of_members = len(self.problem.vars)

        self.min_f = self.get_min_f()

        self.i = number_of_members
        self.j = len(self.design_variables)

        initial_pheromone = 1 / self.get_min_f()
        self.pheromone_matrix = np.full((self.i, self.j), initial_pheromone, dtype=np.float64)
        self.initialize_heuristic_matrix()
        self.update_probabilities()

    def eval_ant(self, ant, x):
        """
        Evaluates one ant:
        Calculate objective function, penalty value and fitness
        :param ant: Ant class object
        :return: fitness
        """
        self.problem.substitute_variables(x)
        """ Evaluate objective function """
        ant.fval = self.problem.obj(x)
        """ Evaluate penalty function """

        ant.penalty = self.penalty(x) * self.penalty_factor

        """ Fitness function equals objective function value + penalty """
        val = ant.fval + ant.penalty

        return val

    def evaluate(self, ants):
        """
        Calculate fitness value for all ants
        :param ants: List of ant objects
        :return: best objective function, list of best design variables
        """
        best_x = []
        best_fitness = np.inf
        best_ant = None
        for ant in ants:
            x = ant.get_x()
            fitness = self.eval_ant(ant, x)
            if fitness < best_fitness:
                best_fitness = fitness
                best_x = ant.get_x()
                best_ant = ant
        return best_fitness, best_x, best_ant

    def local_pheromone_update(self, ant):
        x_index_and_value = ant.x_index_and_value
        members = ant.members
        for i in range(len(members)):
            j = x_index_and_value[i][0]
            self.pheromone_matrix[i][j] *= self.xi

    def global_pheromone_update_for_one_edge(self, fitness, rank, i, j):
        '''
        :param rank: if rank == 0, this is ant with best fitness value. Max value for rank
         is self.number_of_top_ants
        '''
        self.pheromone_matrix[i][j] += self.rho * (self.number_of_top_ants - rank) / fitness

    def global_pheromone_update(self, ants):
        '''
        Choose best ants based on fitness value and update their pheromone value.
        Number of best ants are self.number_of_top_ants
        '''
        # Sort ants based on their fitness
        sorted_ants = sorted(ants, key=lambda ant: ant.get_fitness())

        for it in range(self.number_of_top_ants):
            x_index_and_value = sorted_ants[it].x_index_and_value
            members = sorted_ants[it].members
            for i in range(len(members)):
                j = x_index_and_value[i][0]
                self.global_pheromone_update_for_one_edge(sorted_ants[it].get_fitness(),
                                                          it, i, j)

    def solve(self, problem, maxiter=100, nmax=10):
        self.penalty = problem.exact_penalty_fun()
        self.problem = problem
        self.initialize()

        # Early stop counter
        N: int = 0
        for t in range(maxiter):
            ants = []

            for k in range(self.number_of_ants):
                ants.append(Ant(self.problem.structure.members))
                ants[k].choose_design_variables(self.probability_matrix, self.design_variables)
                self.local_pheromone_update(ants[k])

            best_fitness, best_x, best_ant = self.evaluate(ants)

            self.global_pheromone_update(ants)
            self.update_probabilities()

            if best_fitness < self.best_fitness:
                self.best_fitness = best_fitness
                self.best_x = best_x
                N = 0
            else:
                N += 1

            if N > nmax:
                return self.best_fitness, self.best_x

        return self.best_fitness, self.best_x

class Ant:
    '''
    Class for storing a single ant
    '''

    def __init__(self, members):

        # list of pairs: pair[0]=index, pair[1] = x = design variables
        self.x_index_and_value = []
        self.members = members

        self.fval = None
        self.penalty = None

    def __choose_design_variable(self, i, probability_matrix, design_variables):
        '''
        Choose design variable (cross section area for truss) based on probability P
        '''
        numerator = probability_matrix[i]
        denominator = np.sum(numerator)
        probabilities = numerator / denominator
        j = np.random.choice(range(len(probabilities)), p=probabilities)
        self.x_index_and_value.append([j, design_variables[j]])

    def choose_design_variables(self, probability_matrix, design_variables):
        '''
        Choose design variable based on probability P for all members.
        '''
        for i in range(len(self.members)):
            self.__choose_design_variable(i, probability_matrix, design_variables)

    def get_x(self):
        '''
        Convert x_index_and_value to np array of x
        :return: design variables x
        '''
        x = map(lambda x: x[1], self.x_index_and_value)
        return np.fromiter(x, dtype=np.float64)

    def get_fitness(self):
        return self.fval + self.penalty

if __name__ == "__main__":
    from metku.optimization.benchmarks import simple_ten_bar_truss

    problem = simple_ten_bar_truss.SimpleTenBarTruss('discrete')

    solver = AntColonyOptimizer(pop_size=100, rho=1, alpha=1,
                                beta=0.5, xi=0.75, number_of_top_ants=4,
                                penalty_factor=10000)

    fopt, xopt = solver.solve(problem=problem, maxiter=100)
    problem(xopt)
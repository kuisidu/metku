""" Particle swarm optimization """

import operator
import random
import copy

import matplotlib.pyplot as plt
import numpy as np
from deap import base, creator, tools

from optimization.solvers.optsolver import OptSolver, G_TOL


class PSO(OptSolver):
    """ Particle swarm optimization solver """
    def __init__(self,pop_size=10, c = [2.0,2.0], inertia=0.9, crazy_prob=0.1, 
                 penalty=None,best_f = np.inf):
        """ Constructor
            :param pop_size: population size
            :c: [c1,c2] weighting factors
            :inertia: initial value for the inertia term
            :craziness: craziness probability            
        
        """
        
        super().__init__()

        self.pop_size = pop_size
        self.c = c
        self.inertia = inertia
        self._craziness = crazy_prob
        self.plot = False
        self.best_x = None
        self.best_f = best_f
        self.best_scv = np.inf
        self.prev_best = best_f
        self.best_part = None
        self.swarm = [0 for _ in range(pop_size)]
        self.vrel = 0.2
        self.vmin = None
        self.vmax = None
        self.vred = 0.95
        self.NMAX = 10
        self.KMAX = 5
        if penalty is None:
            def penalty(constr_vals):
                pos_vals = np.clip(constr_vals, 0, 100)
                return 2e6 * sum(pos_vals ** 2)

            self.penalty = penalty
        else:
            self.penalty = penalty

        self.counter = 0
        self.inertia_cnt = 0

    def eval_particle(self, part):
        """
        Evaluates one particle

        :param part: Particle class object
        """
        self.problem.substitute_variables(part.x)
        #X = [round(var.value, 5) for var in self.problem.vars]
        
        """ Evaluate objective function """
        part.fval = self.problem.obj(part.x)
        """ Evaluate penalty function """
        part.scv = self.penalty(part.x)
        
                    
        """ If the particle is currently infeasible, it is sufficient to
            check if the 'scv' value is smaller than the current value.
            
            NOTE! This part needs to be changed if penalty function approach
            is used.
        """
        if part.best_scv > G_TOL:
            #print("Particle is currently infeasible")
            if part.scv <= part.best_scv:
                #print("Improving...")
                part.update_best()
        else:
            """ If the particle is currently feasible, we need to check
                both the scv and val values
            """
            if part.scv <= G_TOL and part.fval <= part.best_f:
                part.update_best()
                
        if self.plot:
            self.update_plot(self.fig, self.ax)
        

    def update_best(self):
        """ Update current best solution 
            Assumption is that all particle of the swarm have already been
            evaluated.
        """
        
        self.counter += 1
        self.inertia_cnt += 1
        
        for part in self.swarm:
            """ If a new best solution has been found, store it """
            if self.best_scv > G_TOL:
                if part.scv <= G_TOL:
                    self.best_part = part
                    self.best_x = copy.deepcopy(part.x)            
                    self.best_f = copy.deepcopy(part.fval)
                    self.best_scv = copy.deepcopy(part.scv)
                    self.best_part = copy.deepcopy(part)
                    self.counter = 0
                    self.inertia_cnt = 0
            else:
                if part.scv <= G_TOL and part.fval < self.best_f:
                    self.best_part = part
                    self.best_x = copy.deepcopy(part.x)            
                    self.best_f = copy.deepcopy(part.fval)
                    self.best_scv = copy.deepcopy(part.scv)
                    self.best_part = copy.deepcopy(part)
                    self.counter = 0
                    self.inertia_cnt = 0
    
    def update_velocity(self,part,verb=False,part_nd=0):
        """ Apply the velocity update rule for particle 'part' """
        rng = np.random.default_rng()
        n = self.problem.nvars()
        
        # Generate multipliers for local and global velocities
        C1 = self.c[0]*rng.random()
        C2 = self.c[1]*rng.random()

        # Store values for the particle
        part.c1 = C1
        part.c2 = C2
        part.w = self.inertia        
        
        # Velocity update
        v0 = self.inertia*part.v + C1*(part.best-part.x) + C2*(self.best_x-part.x)
                        
        #print("v0 = ",v0)
        #print("vmin, vmax = ",self.vmin,self.vmax)
        
        # Keep velocity within velocity bounds
        vnew = np.minimum(np.maximum(v0,self.vmin),self.vmax)
        
        if verb:
            print('---- PARTICLE {0:2g} ----'.format(part_nd))
            print('x = [{0:5.4f},{1:5.4f}]'.format(*part.x))
            print('xbest = [{0:5.4f},{1:5.4f}]'.format(*part.best))
            print('gbest = [{0:5.4f},{1:5.4f}]'.format(*self.best_x))
            print('c1 = {0:5.4f}'.format(part.c1))
            print('c2 = {0:5.4f}'.format(part.c2))
            print('v0 = [{0:5.4f},{1:5.4f}]'.format(*v0))
            print('vnew = [{0:5.4f},{1:5.4f}]'.format(*vnew))
        
        #print("vnew = ",vnew)
        part.v = vnew
    
    def craziness(self,part):
        """ Apply craziness operator """
        rng = np.random.default_rng()
        
        if rng.random() <= self._craziness:
            # change the velocity randomly
            print("Craziness activated")
            part.v = (self.vmax-self.vmin)*rng.random(self.problem.nvars()) + self.vmin
            
    def elite_particle(self):
        """ Replace the worst particle by the best known solution """
        
        
        """ Detect any infeasible points 
            Worst particle is the one with the greatest infeasibility.
            
            NOTE! This has to be modified, if penalty function is used.
        """
        scv = [part.scv for part in self.swarm]
        
        """ Find the index of the particle with greatest constraint violation """
        nd_max = np.argmax(scv)
        if self.swarm[nd_max].scv > G_TOL:
            worst = self.swarm[nd_max]
        else:
            """ If all particles are feasible, find the one with the
                highest objective function value.
            """
            fval = [part.fval for part in self.swarm]
            worst = self.swarm[np.argmax(fval)]
        
        worst = copy.deepcopy(self.best_part)
        
        

    def initizalize_population(self):
        """ Creates initial population with random positions 
        
            NOTE! This has to be checked for discrete variables.
            The 'random_values()' method of the OptimizationProblem class
            should generate only allowable values for the discrete variables.
        """
        rng = np.random.default_rng()
        n = self.problem.nvars()
        
        xlb = np.array([var.lb for var in self.problem.vars])
        xub = np.array([var.ub for var in self.problem.vars])        
        
        self.vmax = self.vrel*(xub-xlb)
        self.vmin = -self.vmax
        
        for i in range(self.pop_size):
            """ Generate random values for the variables """   
            v0 = (self.vmax-self.vmin)*rng.random(n) + self.vmin
            self.swarm[i] = Particle(self.problem.random_values(),v0)
            self.eval_particle(self.swarm[i])
        
        self.update_best()
        
    def print_pop(self):
        
        print("--- Swarm best ---")
        print("fbest = {0:5.4f}".format(self.best_f))
        if self.best_scv <= G_TOL:
            feas = 'Feasible'
        else:
            feas = 'Infeasible'
            
        print("max sum of constraint violation = {0:5.4f} ({1:s})".format(self.best_scv,feas))
        
            
        print("--- Current Swarm ---")
        for i, part in enumerate(self.swarm):
            print(str(i) + ': ' + repr(part))

    def store_population(self):
        """ Stores the current population """
        
        """
        new_popu = []
        for part in self.swarm:
            new_popu.append(copy.deepcopy(part.x))
        """    
        new_popu = copy.deepcopy(self.swarm)
        self.xvals.append(new_popu)
        self.fvals.append(self.best_x)


    def solve(self, problem, maxiter=100, maxtime=-1, log=False,
              min_diff=1e-5, verb=False, plot=False):

        
        self.problem = problem
        xlb = np.array([var.lb for var in self.problem.vars])
        xub = np.array([var.ub for var in self.problem.vars])
        
        self.penalty = problem.exact_penalty_fun()
        """ Creates the fitness function 
            but the create_fitness method does not return anything!
        """
        self.initizalize_population()
        self.store_population()
        if verb:
            self.print_pop()
        """
        if plot and not self.plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.ion()
            plt.show()
            self.plot = True
            self.fig = fig
            self.ax = ax
        """
        
        # Main loop
        it = 0
        prev_val = None
        while it < maxiter:
            it += 1
            problem.num_iters += 1
            print(f"Iteration {it}/{maxiter}")            
            
            """ Perform iteration for particles """
            for i, part in enumerate(self.swarm):
                # Update velocity
                self.update_velocity(part,verb=verb,part_nd=i)
                # Perform craziness
                self.craziness(part)
                # Update design
                # NOTE! This must be modified for discrete variables!
                # (Employ rounding of discrete variables)
                part.move(xlb,xub)
                # evaluate particle at its new location
                self.eval_particle(part)
            
            self.update_best()
            self.elite_particle()
            
            self.store_population()
            if verb:
                self.print_pop()
            
            if plot:
                self.update_plot(self.fig, self.ax)

            if self.counter > self.NMAX:
                print(f"The result has not improved in the last 10 generations."
                      f"\n Returning current best values.")
                return self.best_f, self.best_x

            if self.inertia_cnt > self.KMAX:
                print("Reduce inertia and max velocity")
                self.inertia *= 0.9
                self.vmax *= self.vred
                self.vmin *= self.vred

            """
            if log:
                problem.num_iters += 1
                problem.fvals.append(self.best_f)
                problem.states.append(list(pop[np.argmin(fits)]))
                problem.gvals.append(max(self.constr_vals))
                self.fvals.append(self.best_f)
                self.xvals.append(self.best_x)
            """
        # best_idx = np.argmin(fits)
        # xopt = pop[best_idx]
        # print("XOPT: ", xopt)
        # fopt = fits[best_idx]

        return self.best_f, self.best_x

class Particle:
    """ Class for storing a single particle """

    def __init__(self, x, v=None):
        """ Constructor 
            :param x: value for design variables
            :param v: value for velocity
            :param best: best design the particle has visited (personal best)
            :param best_f: objective function value at 'best'
            :param best_scv: sum of constraint violations at 'best'
            :param scv: sum of constraint violations at current location
            :param fval: objective function value at current location
            :param c1, c2: multipliers for local and global velocity components
            :param w: inertia term
            
        """
        self.x = x
        if v is None:
            self.v = np.zeros(len(x))
        else:
            self.v = v
            
        self.best = copy.deepcopy(x)
        self.best_f = np.inf
        self.best_scv = np.inf
        # Sum of constraint violations
        self.scv = np.inf
        # Objective function value
        self.fval = np.inf
        
        self.c1 = None
        self.c2 = None
        self.w = None

    def move(self,xlb,xub):
        """
        Moves particle
        """
        xnew = self.x + self.v
        # Ensure that the particle does not exit the variable bounds
        # (Perhaps the velocity should be changed instead?)
        self.x = np.minimum(np.maximum(xnew,xlb),xub)
    
    def update_best(self):
        """ Updates the best value found so far """
        self.best_scv = copy.deepcopy(self.scv)
        self.best_f = copy.deepcopy(self.fval)
        self.best = copy.deepcopy(self.x)
    
    def __repr__(self):
        s = '\\x &= ['
        for i in range(len(self.x)):
            s += '{' + str(i) + ':5.4f}, '
        
        s = s[:-2]
        
        s += '] & \\vv &= ['
        s = s.format(*self.x)
        
        for i in range(len(self.v)):
            s += '{' + str(i) + ':5.4f}, '
        
        s = s[:-2]
        s += '] & '
        s = s.format(*self.v)
        s += 'SCV &= {0:5.3f} &'.format(self.scv)
        s += 'f &= {0:5.3f}\\\\'.format(self.fval)
        return s
        
        
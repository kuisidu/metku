# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
import time

import numpy as np
import matplotlib.pyplot as plt

from metku.optimization.constants import GRAD_TOL, ABS_F_TOL, REL_F_TOL, G_TOL, X_TOL

class OptSolver:
    """
        Base class for optimization problem solvers
    """

    def __init__(self):
        self.constr_vals = np.array([-1])
        self.X = np.array([], dtype=float) # Current iterate
        self.problem = None   # OptimizationProblem class object, problem to be solved
        self.fvals = []       # List of objective function values at different iterations
        self.xvals = []       # List of iterates
        self.best_f = np.inf  # Best known objective function value
        self.best_x = None    # Best found design
        self.feasible = False # Flag to determine if a feasible solution has been found
        self.message = None
        self.constr_vals = None

    def calc_constraints(self, x=[]):
        """
        Calculates constraint values and saves them to numpy array

        Parameters:
        ------------
        :param x: (Optional, default: self.X) Point where constraint values are calculated

        """
        constr_vals = []
        if len(x) == 0:
            x = self.X
        # Calculate constraints
        for i in range(len(self.problem.cons)):
            constr_vals.append(self.problem.cons[i](x))
        self.constr_vals = np.asarray(constr_vals)
        return np.asarray(constr_vals)
    
    def is_feasible(self, x=[]):
        
        if len(x) == 0:
            constr_vals = self.constr_vals
        else:
            constr_vals = self.calc_constraints(x)
        
        res = True
        for g, con in zip(constr_vals,self.problem.cons):
            if con.type == '<':
                if g > self.problem.con_tol:
                    res = False
                    break
            elif con.type == '>':
                if g < -self.problem.con_tol:
                    res = False
                    break
            else:
                if abs(g) > self.problem.con_tol:
                    res = False
                    break
        return res

    def _create_eqcons(self):
        """ Creates a list of equality constraints

        :return: constraints
        """
        eqcons = []

        for con in self.problem.cons:
            if con.type == "=":
                eqcons.append(con)

        return eqcons

    def _create_ieqcons(self):
        """ Creates a list of inequality constraints

        :return: list of equality constraint functions
        """
        ieqcons = []

        for con in self.problem.cons:
            if con.type == "<" or con.type == "<=":
                ieqcons.append(con.neg_call)
            elif con.type == ">" or con.type == ">=":
                ieqcons.append(con)

        return ieqcons

    def _create_bounds(self):
        """ Creates a list of optimization variables' bounds

        :return: list of boundaries (lb, ub)
        """
        bounds = []
        for var in self.problem.vars:
            bound = (var.lb, var.ub)
            bounds.append(bound)

        return bounds

    def stopping_criterion(self,criterion='x',**kwargs):
        """ Check various stopping criteria.
            input:
                criterion .. string that states which criterion
                            is checked. Possible values:
                                'x' .. compare two design variable vectors
                                'f' .. compare two objective function values
        """

        res = False

        if criterion == 'x':
            xtol = X_TOL
            for key, value in kwargs.items():
                if key == 'x':
                    x = value
                elif key == 'xtol':
                    xtol = value

            if np.linalg.norm(x[1]-x[0]) <  xtol*np.linalg.norm(x[0]):
                res = True
                self.message = "Stopping criterion fulfilled: relative change to x smaller than tolerance."
        elif criterion == 'f':
            abs_f_tol = ABS_F_TOL
            rel_f_tol = REL_F_TOL

            for key, value in kwargs.items():

                if key == 'f':
                    f = value
                elif key == 'abs_f_tol':
                    abs_f_tol = value
                elif key == 'rel_f_tol':
                    rel_f_tol = value

            # f[1] .. objective function at x(k+1)
            # f[0] .. objective function at x(k)
            if abs(f[1]-f[0]) <= abs_f_tol + rel_f_tol*abs(f[0]):                
                res = True
                self.message = "Stopping criterion fulfilled: relative change to f smaller than tolerance."

        return res

    def take_action(self):
        """ This method needs to be implemented for each solver separately """
        pass

    def update_plot(self, fig, ax):
        """
        Updates the plotted figure
        
        NOTE: This method seems redundant, maybe remove it? (KMe, 6.10.2022)
        
        :return:
        """
        ax.clear()
        for mem in self.problem.structure.members.values():
            (x0, y0), (x1, y1) = mem.coords()
            lw = mem.cross_section.A / 10000
            ax.plot((x0, x1), (y0, y1), 'k', linewidth=lw)
        #self.problem.structure.plot_loads()
        #self.problem.structure.plot_deflection(10, show=False)

        plt.title(f'Feasible: {self.problem.feasible()}')
        plt.axis('equal')
        fig.canvas.draw()
        plt.pause(0.0001)

    def step(self, action):
        """
        Takes a step
        :param action:
        :return:
        """

        self.step_length -= self.step_length * self.step_factor
        self.X += action
        for i in range(len(self.X)):
            self.X[i] = np.clip(self.X[i], self.problem.vars[i].lb,
                           self.problem.vars[i].ub)

        self.problem.substitute_variables(self.X)

        return self.X.copy(), 1, False, 'INFO'

    def solve(self, problem, x0=None, maxiter=-1, maxtime=-1, log=True,
              min_diff=1e-5, verb=False, plot=False):

        """
        Main iteration loop. Solves given problem

        :param problem:
        :param maxiter:
        :param maxtime:
        :return:
        """
        if maxiter < 0:
            maxiter = 100

        if maxtime < 0:
            maxtime = 1e9


        # Assign problem
        self.problem = problem
        # If initial starting point is/isn't defined        
        if x0 is not None:
            self.X = np.asarray(x0, dtype=float)
            problem.substitute_variables(x0)
        else:
            self.X = self.random_feasible_point()
            problem.substitute_variables(self.X)

        # Store initial point
        self.xvals.append(np.array(x0))
        self.fvals.append(problem.obj(x0))


        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.ion()
            plt.show()

        # Assign done
        done = False
        # Start iteration
        t_total = 0
        for i in range(maxiter):
            #if verb:
            #    print("*** Iteration {0:g} ***".format(i+1))
            if plot:
                self.update_plot(fig, ax)
            # Check time
            start = time.time()
            t_0 = time.process_time()
            if t_total >= maxtime or done:
                print(t_total,maxtime,done)
                print("T_total - Stopping criterion fulfilled.")
                break
            # Save previous state
            prev_state = self.X.copy()
            # Define action to take. For each method, this is the
            # procedure for updating the design variables.
            # The take_action method needs to be implemented separately for
            # each method.
            action = self.take_action()
            # Take step
            state, reward, done, info = self.step(action)
                    
            # if (np.linalg.norm(prev_state - state)/np.linalg.norm(prev_state) <= min_diff):
            #    print("Sufficiently small change in variable values in consecutive iterations.")
            #    print(state,prev_state)
            #    done = True
            #     #break

            #if np.all(prev_state == state):
            #    print("No change in variable values.")
            #    done = True
                #break
            # Change current state
            self.X = state

            # Substitute new variables
            problem.substitute_variables(state)
            # Calculate constraints
            #print(self.X)
            #self.calc_constraints(self.X)
            state_feasible = self.is_feasible(self.X)
            # print('Check state feasibility',state_feasible)
            #print(state)
            # Save best values

            fval = problem.obj(state)
            t_1 = time.process_time()
            t_total += t_1-t_0

            if fval < self.best_f and state_feasible:
                self.best_f = fval
                self.best_x = state.copy()
                # if verb:
                #     #print(self.best_x)
                #     if len(state) <= 10:
                #         print(f"New best!: {fval:.3f} {[round(s, 3) for s in state]}")
                #     else:
                #         print(f"New best!: {fval:.3f}")

            # Log objective vals per iteration
            if log:
                #if verb:
                #    print("logging X = ",self.X)
                problem.num_iters += 1
                problem.fvals.append(problem.obj(self.X))
                problem.states.append(list(state))
                problem.gvals.append(list(self.constr_vals).copy())
                self.fvals.append(problem.obj(self.X))
                self.xvals.append(self.X)
            end = time.time()

            if verb:
                print(
                    f'\r Iteration {i + 1} / {maxiter}: Obj: {fval:.4f} Feasible: {state_feasible} '\
                    f'max g: {max(self.constr_vals):.4f} best value: {self.best_f:.2f} '\
                    f'Iteration took: {end - start :.2f} s {" "*100}',end="")

            # CHECK STOPPING CRITERIA
            # If new state is almost same as previous AND the design is feasible,
            # the iteration can be stopped.
            # TODO! Add Feasibility Check
            if self.is_feasible():
                if self.stopping_criterion('x', x=[prev_state, state], xtol=X_TOL):
                    print(self.message)
                    done = True
    
                if self.stopping_criterion('f', f=[problem.obj(prev_state),problem.obj(state)], 
                                           abs_f_tol=ABS_F_TOL, rel_f_tol=REL_F_TOL):
                    print(self.message)                
                    done = True

            if done:
                break

        if self.best_x is None:
            self.best_x = self.X
            self.best_f = problem.obj(self.X)

        return self.best_f, self.best_x    

    def random_feasible_point(self):
        """
        Creates a random feasible starting point for optimization

        :return: random feasible starting point
        """
        print("Creating random feasible starting point!")

        X = [0] * self.problem.nvars()
        for i, var in enumerate(self.problem.vars):
            if self.problem.prob_type == 'continuous':
                X[i] = var.lb + (var.ub - var.lb) * np.random.rand()
            else:
                X[i] = np.random.randint(var.lb, var.ub)
        while np.any(self.calc_constraints(X) > -0.15):

            for i, var in enumerate(self.problem.vars):
                if self.problem.prob_type == 'continuous':
                    X[i] = var.lb + (var.ub - var.lb) * np.random.rand()
                else:
                    X[i] = np.random.randint(var.lb, var.ub)

            self.problem.substitute_variables(X)
            self.problem.fea()

        print("Starting point created!")

        return np.asarray(X)

    def plot_iters(self):
        """
        Plots iteration history

        Returns
        -------
        None.

        """
        
        niter = len(self.fvals)
        
        if niter > 0:        
            fig, ax = plt.subplots(1)
            it = np.arange(niter)
            ax.plot(it,self.fvals)
            
            ax.set_xlabel('Iteration')
            ax.set_ylabel('Objective function')
            
            ax.grid(True,axis='both',linestyle='--')
        
        

import numpy as np
from scipy.optimize import minimize, NonlinearConstraint, Bounds
from scipy.optimize import LinearConstraint as LinCon

try:
    import src.optimization.structopt as sopt
    from src.optimization.solvers.optsolver import OptSolver
except:
    import optimization.structopt as sopt
    from optimization.solvers.optsolver import OptSolver


class TrustRegionConstr(OptSolver):
    """ Trust-region algorithm of scipy.optimize """

    def __init__(self):
        super().__init__()


    def calc_constraints(self, x=[]):
        """
        Calculates constraint values and saves them to numpy array

        Parameters:
        ------------
        :param x: (Optional, default: self.X) Point where constraint values are calculated

        """
        constr_vals = []
        if not len(x):
            x = self.X
        # Calculate constraints
        for i in range(len(self.problem.cons)):
            constr_vals.append(self.problem.cons[i](x))

        return np.asarray(constr_vals)

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
                ieqcons.append(con)
            elif con.type == ">" or con.type == ">=":
                ieqcons.append(con.neg_call)

        return ieqcons


    def _create_bounds(self):
        """ Creates a list of optimization variables' bounds

        :return: list of boundaries (lb, ub)
        """
        lb = []
        ub = []
        #bounds = []
        for var in self.problem.vars:
            lb.append(var.lb)
            ub.append(var.ub)
            #bound = (var.lb, var.ub)
            #bounds.append(bound)
        
        bounds = Bounds(lb, ub, keep_feasible=False)

        return bounds



    def solve(self, problem, maxiter=100, x0=[]):
        """
        Solves given problem

        :param problem: problem to be solved
        :param maxiter: number of maximum iterations
        :param x0: initial guess

        :return: fopt, xopt
        """

        self.problem = problem
        
        constraints = []
        
        for con in self.problem.cons:           
            if isinstance(con,sopt.NonLinearConstraint):
                lb = -np.inf
                ub = np.inf
                if con.type == "<" or con.type == "<=":                
                    ub = 0
                elif con.type == ">" or con.type == ">=":
                    lb = 0
                elif con.type == "=":
                    lb = 0
                    ub = 0
                    
                constraints.append(NonlinearConstraint(
                    con,
                    lb,
                    ub,
                    jac='2-point',
                    # finite_diff_rel_step=1e-6
                    ))
            elif isinstance(con,sopt.LinearConstraint):
                lb = -np.inf
                ub = np.inf
                if con.type == "<" or con.type == "<=":                
                    ub = con.b
                elif con.type == ">" or con.type == ">=":
                    lb = con.b
                elif con.type == "=":
                    lb = con.b
                    ub = con.b
                                    
                
                constraints.append(LinCon(con.a,lb,ub))
            
        
        bounds = self._create_bounds()
        """
        eqcons = self._create_eqcons()
        ieqcons = self._create_ieqcons()
        """
        
        bnd_cons = LinCon(A=np.eye(problem.nvars()), lb=bounds.lb, ub=bounds.ub)

        constraints.append(bnd_cons)

        # Create initial guess if one isn't provided
        if not len(x0):
            x0 = []
            for var in problem.vars:
                x0.append(var.ub)  # * np.random.uniform())

        print(bounds.lb,bounds.ub)
        
        #print(constraints[-3].ub)

        options = {'maxiter': maxiter,
                   'verbose': 2,
                   'xtol': 1e-8,
                   'gtol': 1e-8,
                   # 'finite_diff_rel_step': 1e-6
                   }

        """
        constraints = []
        
        for eqcon in eqcons:
            con = {'type': 'eq', 'fun': eqcon}
            constraints.append(con)

        for ineqcon in ieqcons:
            con = {'type': 'ineq', 'fun': ineqcon}
            constraints.append(con)
        """

        try:
            out = minimize(self.problem.obj,
                           x0,
                           method='trust-constr',
                           tol = 1e-10,
                           bounds=bounds,
                           constraints=constraints,
                           options=options,
                           # callback=self.callback
                           )

            print(out)

            print(self.calc_constraints(out.x))
            self.best_x = out.x
            self.best_f = out.fun
            self.X = out.x

            return out.fun, out.x, out.nit

        except:
            self.best_x = x0
            self.best_f = -np.inf
            return -np.inf, x0, 0

    def callback(self, xk, state):
        fval = self.problem.obj(xk)
        self.problem.fvals.append(fval)
        return print('CALLBACK', xk)


if __name__ == '__main__':
    from src.optimization.benchmarks import *

    problem = TenBarTruss('continuous')
    solver = TrustRegionConstr()
    xopt, fopt = solver.solve(problem, maxiter=10)

    problem(xopt)



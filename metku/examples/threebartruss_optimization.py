from metku.optimization.structopt import OptimizationProblem, Variable, \
    NonLinearConstraint, LinearConstraint, DiscreteVariable
from metku.frame2d.frame2d import FrameMember, Frame2D, SteelBeam, PointLoad, XYHingedSupport
import numpy as np

class ThreeBarTruss(OptimizationProblem):

    # Problem parameters
    P = 20
    L = 100
    sigma_U = 20
    sigma_L = -15

    def __init__(self):
        super().__init__(name='ThreeBarTruss')

        # Problem type needs to be addressed
        # otherwise some solvers might not work
        self.prob_type = 'continuous'
        self.create_variables()
        self.create_constraints()
        self.create_objective()

    def create_variables(self):
        """
        Creates variables
        """
        self.vars = []

        var1 = Variable(name="X1",
                        lb=0.01,
                        ub=5)

        var2 = Variable(name="X2",
                        lb=0.01,
                        ub=5)

        self.vars = [var1, var2]


    def create_objective(self):
        """
        Creates objective function and assigns it to problem
        """

        def objective(X):
            """
            Weight of the truss
            :param X: design variables
            :return: weight
            """
            x1, x2 = X
            L1 = np.sqrt(2) * self.L
            L2 = self.L

            return 2 * L1 * x1 + L2 * x2

        # Assign function to problem
        self.obj = objective

    def create_constraints(self):
        """
        Creates constraints and assigns them to problem
        """
        L1 = np.sqrt(2) * self.L
        L2 = self.L

        self.cons = []

        # Constraint 1
        def g1(x):
            x1, x2 = x
            return self.P * (x2 + np.sqrt(2) * x1) / (
                        2 * x1 * x2 + np.sqrt(2) * x1 ** 2) - self.sigma_U

        # Constraint 2
        def g2(x):
            x1, x2 = x
            return self.P * (np.sqrt(2) * x1) / (
                        2 * x1 * x2 + np.sqrt(2) * x1 ** 2) - self.sigma_U

        # Constraint 3
        def g3(x):
            x1, x2 = x
            return self.P * (x2) / (2 * x1 * x2 + np.sqrt(2) * x1 ** 2) + self.sigma_L

        # Create NonLinearConstraint objects
        con1 = NonLinearConstraint(con_fun=g1, name="Con1", con_type="<")
        con2 = NonLinearConstraint(con_fun=g2, name="Con2", con_type="<")
        con3 = NonLinearConstraint(con_fun=g3, name="Con3", con_type="<")

        # Assign created constraints as a list to the problem
        self.cons = [con1, con2, con3]

if __name__ == '__main__':
    # Import all solvers
    from metku.optimization.solvers import *
    # Create optimization problem
    problem = ThreeBarTruss()
    # Choose solver
    solver = SLP(move_limits=[0.5, 5])
    # Create random initial starting point
    x0 = [np.random.randint(var.lb, var.ub) for var in problem.vars]
    # Solve problem
    # NOTE! if you use SLSQP the order is xopt, fopt
    fopt, xopt = solver.solve(problem, x0=x0, maxiter=500)
    # Print results by calling the problem with the optimal parameters
    problem(xopt)




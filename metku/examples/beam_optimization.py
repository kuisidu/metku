# Import dependencies
from metku.optimization.structopt import OptimizationProblem, Variable, NonLinearConstraint, LinearConstraint
from metku.frame2d.frame2d import Frame2D, SteelBeam, XYHingedSupport, XHingedSupport, LineLoad

import numpy as np

class BeamOptimizationProblem(OptimizationProblem):
    """
    Optimization example for dual-symmetric steel beam section optimization

    """
    def __init__(self):
        super().__init__(name="BeamOptimizationProblem")

        # Problem type needs to be declared, otherwise some solvers might not work
        self.prob_type = 'continuous'
        self.create_structure()
        self.create_variables()
        self.create_constraints()
        self.create_objective()


    def create_structure(self):
        """
        Creates the structure to be optimized

        In this case the structure is just the one beam
        """
        frame = Frame2D()

        # Lets create 6000 mm long beam with welded I-section profile
        coordinates = [[0, 0], [6000, 0]]
        beam = SteelBeam(coordinates,
                         profile='WI 500-12-10X300-10X300',
                         material="S355")

        # Add the beam to frame
        frame.add(beam)

        # Create supports
        sup1 = XYHingedSupport(coordinates[0])
        sup2 = XYHingedSupport(coordinates[1])

        # Add supports to frame
        frame.add(sup1)
        frame.add(sup2)

        # Create and add lineload
        ll = LineLoad(beam, [-39.15, -39.15], 'y')
        frame.add(ll)

        # Generate FE-model for the frame
        frame.generate()

        # Assing this frame as the structure to be optimized
        self.structure = frame


    def create_variables(self):
        """
        Creates optimization design variables

        In this case the variables are cross-section's dimensions
        """
        # Initialize as an empty array
        self.vars = []

        # Get the beam
        beam = self.structure.members[0]

        # Get the SteelSection object from beam
        section = beam.cross_section

        # CREATING DESIGN VARIABLES
        # -------------------------
        # Flanges' breadth
        var_BF = Variable(name="BF",
                     lb=100,
                     ub=700,
                     target={"property": "BF",
                             'objects': [section]})
        # Flanges' thickness
        var_TF = Variable(name="TF",
                     lb=5,
                     ub=100,
                     target={"property": "TF",
                             'objects': [section]})
        # Section's height
        var_H = Variable(name="H",
                     lb=300,
                     ub=3300,
                     target={"property": "H",
                             'objects': [section]})
        # Web's thickness
        var_TW = Variable(name="TW",
                     lb=6,
                     ub=100,
                     target={"property": "TW",
                             'objects': [section]})

        # Append created variables to list
        # This is the order of the design variables
        self.vars = [var_BF, var_TF, var_H, var_TW]




    def create_confuns(self, beam):
        """
        Creates constraint functions

        :param beam: FrameMember object
        :return: constraint functions
        """
        def bending_moment_confun(x):
            # Index 0 is the bending resistance about the major axis
            g = beam.med / beam.MRd[0] - 1

            return g

        def shear_confun(x):
            g = beam.ved / beam.VRd - 1

            return g

        def deflection_confun(x):

            min_v = min([disp[1] for disp in beam.nodal_displacements.values()])
            max_v = max([disp[1] for disp in beam.nodal_displacements.values()])
            abs_max_v = max(abs(min_v), max_v)

            g = abs_max_v / (beam.length / 250) - 1

            return g

        return bending_moment_confun, shear_confun, deflection_confun

    def flange_class_con(self):
        """
        Creates flange class constraint
        :return: constraint
        """
        # Get the beam
        beam = self.structure.members[0]
        eps = np.sqrt(235 / beam.fy)

        # The order of the variables is defined in create variables -method
        # x = BF, TF, H, TW

        # NOTE: This neglects the weld
        a = np.array([0.5, -14*eps, 0, -0.5])

        con = LinearConstraint(a=a,
                               b=0,
                               con_type='<',
                               name="Flange class")
        return con


    def web_class_con(self):
        """
        Creates flange class constraint
        :return: constraint
        """
        # Get the beam
        beam = self.structure.members[0]
        eps = np.sqrt(235 / beam.fy)

        # The order of the variables is defined in create variables -method
        # BF, TF, H, TW = x

        # NOTE: This neglects the weld
        a = np.array([0, -1, 1, -124*eps])

        con = LinearConstraint(a=a,
                               b=0,
                               con_type='<',
                               name="Web class")
        return con

    def create_constraints(self):
        """
        Creates constraints

        """
        # Initialize as an empty array
        self.cons = []

        # Get the beam
        beam = self.structure.members[0]

        # Create constraint functions
        bending_moment_confun, shear_confun, deflection_confun = self.create_confuns(beam)

        # Create constraints
        bending_constr = NonLinearConstraint(con_fun=bending_moment_confun,
                                             name="Bending Constraint",
                                             parent=self)

        shear_constr = NonLinearConstraint(con_fun=shear_confun,
                                             name="Shear Constraint",
                                             parent=self)

        deflection_constr = NonLinearConstraint(con_fun=deflection_confun,
                                             name="Deflection Constraint",
                                             parent=self)
        # Deflection requires FEM analysis
        deflection_constr.fea_required = True

        # Cross-section class constraints
        web_con = self.web_class_con()
        flange_con = self.flange_class_con()

        # Assign constraints
        self.cons = [bending_constr,
                     shear_constr,
                     deflection_constr,
                     web_con,
                     flange_con]


    def create_objective(self):
        """
        Creates the objective function, function to be minimized

        In this case the objective is to minimize beam's weight
        """
        def obj(x):
            BF, TF, H, TW = x

            beam = self.structure.members[0]
            return beam.weight

        self.obj = obj


if __name__ == '__main__':

    from metku.optimization.solvers import *

    problem = BeamOptimizationProblem()

    x0 = [var.lb for var in problem.vars]
    solver = SLSQP()

    xopt, fopt = solver.solve(problem, x0=x0, maxiter=20)

    problem(xopt)
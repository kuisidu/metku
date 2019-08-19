from src.frame2d.frame2d import *
from src.optimization.structopt import *

AISC_PROFILE_AREAS = [71.613, 90.968, 126.451, 161.29, 198.064, 252.258,
                      285.161, 363.225, 388.386, 494.193, 506.451, 641.289,
                      645.16, 792.256, 816.773, 939.998, 1008.385, 1045.159,
                      1161.288, 1283.868, 1374.191, 1535.481, 1690.319,
                      1696.771, 1858.061, 1890.319, 1993.544, 2019.351,
                      2180.641, 2283.705, 2290.318, 2341.931, 2477.414,
                      2496.769, 2503.221, 2696.769, 2722.575, 2896.768,
                      2961.284, 3096.768, 3206.445, 3303.219, 3703.218,
                      4658.055, 5141.925, 5503.215, 5999.988, 6999.986,
                      7419.34, 8709.66, 8967.724, 9161.272, 9999.98, 10322.56,
                      10903.204, 12129.008, 12838.684, 14193.52, 14774.164,
                      15806.42, 17096.74, 18064.48, 19354.8, 21612.86]


class FiftyTwoBarTruss(OptimizationProblem):
    # Problem parameters
    E = 207e3  # MPa
    rho = 7860e-9  # kg/mm^3
    Fx = 100e3  # N
    Fy = 200e3  # N
    L = 2000  # mm
    H = 3000  # mm

    properties = {
        'L': str(L) + " mm",
        'Fx': f'{Fx / 1e3 :.2f} kN',
        'Fy': f'{Fy / 1e3 :.2f} kN',
        'E': f'{E / 1e3 :.2f} MPa',
        'rho': f'{rho * 1e9 :.0f} kg/m3'
    }

    # Constraint limits
    sigma_max = 180  # MPa

    def __init__(self, prob_type='discrete'):

        super().__init__("FiftyTwoBarTruss")
        self.prob_type = prob_type
        self.structure = self.create_structure()
        self.create_variables()
        self.create_objective()
        self.create_constraints()

    def constraint_generator(self, mem):
        """
        Creates non-linear constraint functions fro each member
        :param mem: member for which the constraints are created
        :return: created functions
        """

        def compression_fun(x):
            return -mem.ned / (mem.A * mem.fy) - 1

        def tension_fun(x):
            return mem.ned / (mem.A * mem.fy) - 1

        # def disp_fun(x):
        #     displacements = mem.nodal_displacements.values()
        #     max_vals = [max(l[0:2]) for l in displacements]
        #     min_vals = [min(l[0:2]) for l in displacements]
        #     max_val = max(max_vals)
        #     min_val = min(min_vals)
        #     abs_max = max(max_val, abs(min_val))
        #     return abs_max / self.delta_max - 1

        return compression_fun, tension_fun  # , disp_fun

    def create_constraints(self):
        """
        Creates constraints
        :return:
        """
        # Initialize constraints as an empty list
        self.cons = []
        for var in self.vars:
            for mem in var.target["objects"]:
                if isinstance(mem, FrameMember):
                    compression_fun, tension_fun = self.constraint_generator(
                        mem)

                    comp_con = NonLinearConstraint(con_fun=compression_fun,
                                                   name="Compression "
                                                        + str(mem.mem_id),
                                                   parent=self)
                    comp_con.fea_required = True

                    tension_con = NonLinearConstraint(con_fun=tension_fun,
                                                      name="Tension "
                                                           + str(mem.mem_id),
                                                      parent=self)
                    tension_con.fea_required = True

                    # disp_con = NonLinearConstraint(con_fun=disp_fun,
                    #                                name='Displacement '
                    #                                + str(mem_id),
                    #                                parent=self)
                    # disp_con.fea_required = True
                    self.cons.append(comp_con)
                    self.cons.append(tension_con)
                    # self.cons.append(disp_con)

    def create_structure(self):
        """
        Creates the structure
        :return: structure
        """

        frame = Frame2D(num_elements=1)
        nodes = {}
        # Starting coordinates for nodes
        x = 0
        y = 0
        for i in range(20):
            if i % 4:
                x += self.L

            elif i == 0:
                x = 0
                y = 0
            else:
                x = 0
                y += self.H

            nodes[i] = [x, y]

        # Create bars
        self.vert_groups = []
        self.diag_groups = []
        self.horiz_groups = []
        for j in range(4):
            # Vertical bars
            vert_group = []
            for k in range(4):
                # Node indices
                idx1 = j * 4 + k
                idx2 = idx1 + 4
                n1 = nodes[idx1]
                n2 = nodes[idx2]
                beam = SteelBeam([n1, n2])
                frame.add(beam)
                vert_group.append(beam)
            # Diagonal bars
            diag_group = []
            for l in range(3):
                idx1 = j * 4 + l
                idx2 = idx1 + 1
                idx3 = idx1 + 4
                idx4 = idx1 + 5
                n1 = nodes[idx1]
                n2 = nodes[idx2]
                n3 = nodes[idx3]
                n4 = nodes[idx4]
                beam1 = SteelBeam([n2, n3])
                beam2 = SteelBeam([n1, n4])
                frame.add(beam1)
                frame.add(beam2)
                diag_group.append(beam1)
                diag_group.append(beam2)

            # Horizontal bars
            horiz_group = []
            for m in range(3):
                idx1 = j * 4 + m
                idx2 = idx1 + 1
                idx3 = idx1 + 4
                idx4 = idx1 + 5
                n3 = nodes[idx3]
                n4 = nodes[idx4]
                beam = SteelBeam([n3, n4])
                frame.add(beam)
                horiz_group.append(beam)

            self.vert_groups.append(vert_group)
            self.diag_groups.append(diag_group)
            self.horiz_groups.append(horiz_group)

        # Supports & Loads
        for i in range(4):
            # Support
            nsup = nodes[i]
            frame.add(XYHingedSupport(nsup))
            # Point load
            npload = nodes[i + 16]
            pl1 = PointLoad(npload, [self.Fx, 0, 0])
            pl2 = PointLoad(npload, [0, self.Fy, 0])
            frame.add(pl1)
            frame.add(pl2)

        # Change Young's modulus
        for mem in frame.members.values():
            mem.E = self.E
            mem.fy = self.sigma_max

        frame.hinge_joints()
        # frame.add_self_weight()
        frame.generate()
        frame.f.draw()
        return frame

    def create_variables(self):
        """
        Creates variables used in optimization

        Appends each variable to 'vars' 'list from where they can be accessed
        """
        groups = []
        for i in range(4):
            groups.append(self.vert_groups[i])
            groups.append(self.diag_groups[i])
            groups.append(self.horiz_groups[i])

        for i, group in enumerate(groups):
            name = 'A' + str(i + 1)
            if self.prob_type == "discrete":
                var = DiscreteVariable(name,
                                       profiles=AISC_PROFILE_AREAS,
                                       target={"property": "A",
                                               "objects": group})
            elif self.prob_type == "continuous":
                var = Variable(name,
                               lb=AISC_PROFILE_AREAS[0],
                               ub=AISC_PROFILE_AREAS[-1],
                               target={"property": "A",
                                       "objects": group})
            else:
                raise TypeError("Problem type must be either 'discrete' "
                                "or 'continuous")

            self.vars.append(var)

    def create_objective(self):
        """
        Creates the objective function
        """

        def objective(x):
            weight = 0
            for mem in self.structure.members.values():
                weight += self.rho * mem.A * mem.length

            return weight

        self.obj = objective


if __name__ == '__main__':
    from src.optimization.solvers import *
    from src.optimization.solvers.optsolver import DiscreteVNS
    from src.optimization.to_csv import to_csv
    from src.optimization.result_exporter import ResultExporter
    import matplotlib.pyplot as plt

    problem = FiftyTwoBarTruss('continuous')
    solver = SLP(move_limits=[0.15, 3], gamma=1e-3)
    #solver = DiscreteVNS(step_length=3)
    #x0 = [43, 18,  5, 42, 15,  5, 29, 16,  7, 19, 18,  9]
    x0 = [var.ub for var in problem.vars]
    solver.solve(problem, x0=x0, maxiter=500, log=True)
    problem(solver.X, prec=6)
    exporter = ResultExporter(problem, solver)
    exporter.iteration_jpg(fopt=1897)
    exporter.to_csv()
    #to_csv(problem, solver)

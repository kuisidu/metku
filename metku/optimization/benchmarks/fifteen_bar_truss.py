# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela

from metku.frame2d.frame2d import *
from metku.optimization.structopt import *

FIFTEEN_BAR_AREAS_in2 = [0.111, 0.141, 0.174, 0.220, 0.270, 0.287, 0.347,
                         0.440, 0.539, 0.954, 1.081, 1.174, 1.333, 1.488,
                         1.764, 2.142, 2.697, 2.800, 3.131, 3.565, 3.813,
                         4.805, 5.952, 6.572, 7.192, 8.525, 9.300, 10.850,
                         13.300, 14.290, 17.170, 19.180]

FIFTEEN_BAR_AREAS_mm2 = [71.613, 90.968, 112.258, 141.935, 174.193, 185.161,
                         223.871, 283.87, 347.741, 615.483, 697.418, 757.418,
                         859.998, 959.998, 1138.062, 1381.933, 1739.997,
                         1806.448, 2019.996, 2299.995, 2459.995, 3099.994,
                         3839.992, 4239.992, 4639.991, 5499.989, 5999.988,
                         6999.986, 8580.628, 9219.336, 11077.397, 12374.169]



class FifteenBarTruss(OptimizationProblem):
    # Problem parameters
    L = 3048  # mm
    F = 44537 # N # 1e5 lb
    E = 68950  # MPa # 10e7 psi
    rho = 2720e-9  # kg/mm3

    properties = {
        'L': str(L) + " mm",
        'F': f'{F / 1e3 :.2f} kN',
        'E': f'{E / 1e3 :.2f} MPa',
        'rho': f'{rho * 1e9 :.0f} kg/m3'
    }
    # Constraint limits

    sigma_max = 172.37  # MPa # 25000 psi
    delta_max = 50.8  # mm

    def __init__(self, prob_type='discrete'):
        super().__init__(name='FifteenBarTruss')
        self.prob_type = prob_type
        self.structure = self.create_structure()
        self.create_variables(profiles=FIFTEEN_BAR_AREAS_mm2)
        self.create_constraints()
        self.create_objective()

    def create_objective(self):

        def objective(X):
            self.substitute_variables(X)
            weight = 0
            for mem in self.structure.members.values():
                weight += self.rho * mem.A * mem.length

            return weight

        self.obj = objective

    def create_variables(self, profiles=[0, 1e6]):
        """
        Creates variables used in optimization

        Appends each variable to 'vars' 'list from where they can be accessed
        """
        self.vars = []
        # Member area variables
        for i, mem in enumerate(self.structure.members.values()):
            name = 'A' + str(i + 1)
            if self.prob_type == "discrete":
                var = DiscreteVariable(name,
                                       values=profiles,
                                       target={"property": "A",
                                               "objects": [mem]})

            elif self.prob_type == "continuous":
                var = Variable(name,
                               lb=profiles[0],
                               ub=profiles[-1],
                               target={"property": "A",
                                       "objects": [mem]})
            else:

                raise TypeError("Problem type must be either 'dicrete' "
                                "or 'continuous")

            self.add(var)

        # Node location variables
        # Nodes 2, 3, 6, 7 can move in x- and y -directions
        # x2 == x6 , x3 == x7
        # Nodes 4, 8 can only move in y -direction
        # Nodes 1, 5 are stationary
        x_nodes = [self.structure.f.nodes[i] for i in [1, 2, 5, 6]]
        y_nodes = [self.structure.f.nodes[i] for i in [1, 2, 3, 5, 6, 7]]

        x2 = Variable('x2',
                      lb=254,
                      ub=3556,
                      target={"property": "x",
                              "objects": [x_nodes[0], x_nodes[2]]})

        x3 = Variable('x3',
                      lb=5588,
                      ub=6604,
                      target={"property": "x",
                              "objects": [x_nodes[1], x_nodes[3]]})

        y2 = Variable('y2',
                      lb=2540,
                      ub=3556,
                      target={"property": "y",
                              "objects": [x_nodes[0]]})

        y3 = Variable('y3',
                      lb=2540,
                      ub=3556,
                      target={"property": "y",
                              "objects": [y_nodes[1]]})
        y4 = Variable('y4',
                      lb=1525,
                      ub=2286,
                      target={"property": "y",
                              "objects": [y_nodes[2]]})
        y6 = Variable('y6',
                      lb=-508,
                      ub=508,
                      target={"property": "y",
                              "objects": [y_nodes[3]]})
        y7 = Variable('y7',
                      lb=-508,
                      ub=508,
                      target={"property": "y",
                              "objects": [y_nodes[4]]})
        y8 = Variable('y8',
                      lb=508,
                      ub=1524,
                      target={"property": "y",
                              "objects": [y_nodes[5]]})

        for var in [x2, x3, y2, y3, y4, y6, y7, y8]:
            self.add(var)

    def create_structure(self):

        frame = Frame2D(num_elements=1)
        # Nodes
        n1 = [0, self.L]
        n2 = [self.L, self.L]
        n3 = [2 * self.L, self.L]
        n4 = [3 * self.L, self.L]
        n5 = [0, 0]
        n6 = [self.L, 0]
        n7 = [2 * self.L, 0]
        n8 = [3 * self.L, 0]
        # Create bars
        bar1 = SteelBeam([n1, n2])
        bar2 = SteelBeam([n2, n3])
        bar3 = SteelBeam([n3, n4])
        bar4 = SteelBeam([n5, n6])
        bar5 = SteelBeam([n6, n7])
        bar6 = SteelBeam([n7, n8])
        bar7 = SteelBeam([n6, n2])
        bar8 = SteelBeam([n7, n3])
        bar9 = SteelBeam([n8, n4])
        bar10 = SteelBeam([n1, n6])
        bar11 = SteelBeam([n5, n2])
        bar12 = SteelBeam([n2, n7])
        bar13 = SteelBeam([n6, n3])
        bar14 = SteelBeam([n3, n8])
        bar15 = SteelBeam([n7, n4])
        # Add bars
        frame.add(bar1)
        frame.add(bar2)
        frame.add(bar3)
        frame.add(bar4)
        frame.add(bar5)
        frame.add(bar6)
        frame.add(bar7)
        frame.add(bar8)
        frame.add(bar9)
        frame.add(bar10)
        frame.add(bar11)
        frame.add(bar12)
        frame.add(bar13)
        frame.add(bar14)
        frame.add(bar15)
        # Supports
        frame.add(XYHingedSupport(n1))
        frame.add(XYHingedSupport(n5))
        # Point loads
        frame.add(PointLoad(n8, [0, -self.F, 0]))
        # Change material properties
        for mem in frame.members.values():
            mem.fy = self.sigma_max
            mem.E = self.E  # MPa
            mem.A = FIFTEEN_BAR_AREAS_mm2[0]
            mem.Sj1 = 0
            mem.Sj2 = 0
        frame.generate()
        frame.calculate()
        return frame


    def constraint_generator(self, mem):

        def tension_fun(X):
            return mem.ned / mem.NRd - 1

        def compression_fun(X):
            return -mem.ned / mem.NRd - 1

        def buckling_fun(X):
            sigma_cr = 100 * mem.E * mem.A / (8 * mem.length ** 2)
            sigma = -mem.ned / mem.A
            return sigma / sigma_cr - 1


        return tension_fun, compression_fun, buckling_fun

    def create_constraints(self):

        # Initialize constraints as an empty list
        self.cons = []

        i = 0
        for j, var in enumerate(self.vars):
            for mem in var.target["objects"]:
                if isinstance(mem, FrameMember):
                    i += 1

                    compression_fun, tension_fun, buckling_fun = self.constraint_generator(mem)

                    comp_con = NonLinearConstraint(con_fun=compression_fun,
                                                   name="Compression " + str(
                                                       i),
                                                   problem=self)
                    comp_con.fea_required = True

                    tension_con = NonLinearConstraint(con_fun=tension_fun,
                                                      name="Tension " + str(i),
                                                      problem=self)
                    tension_con.fea_required = True

                    # buckl_con = NonLinearConstraint(con_fun=buckling_fun,
                    #                                name='Buckling ' + str(
                    #                                    i),
                    #                                problem=self)
                    # buckl_con.fea_required = True

                    self.add(comp_con)
                    self.add(tension_con)
                    # self.add(buckl_con)



if __name__ == '__main__':
    from metku.optimization.solvers import *
    problem = FifteenBarTruss(prob_type='continuous')
    #solver = GA(popsize=10)
    solver = SLP(move_limits=[0.05, 0.05], gamma=1e-2)
    x0 = [var.ub for var in problem.vars]
    solver.solve(problem, x0=x0, maxiter=200, verb=True, plot=True)
    problem(solver.X)
    problem.structure.plot()

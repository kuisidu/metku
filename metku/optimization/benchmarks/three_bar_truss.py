# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-

from metku.optimization.structopt import *
from metku.frame2d.frame2d import *


THREE_BAR_AREAS_mm2 = [71.613, 90.968, 112.258, 141.935, 174.193, 185.161,
                         223.871,
                         283.87, 347.741, 615.483, 697.418, 757.418, 859.998,
                         959.998,
                         1138.062, 1381.933, 1739.997, 1806.448, 2019.996,
                         2299.995,
                         2459.995, 3099.994, 3839.992, 4239.992, 4639.991,
                         5499.989,
                         5999.988, 6999.986, 8580.628, 9219.336, 11077.397,
                         12374.169]


class ThreeBarTruss(OptimizationProblem):

    # Problem parameters
    L = 3048  # mm
    F = 444890  # N # 1e5 lb
    E = 68950  # MPa # 10e7 psi
    rho = 2768e-9  # kg/mm3

    properties = {
        'L': str(L) + " mm",
        'F': f'{F / 1e3 :.2f} kN',
        'E': f'{E / 1e3 :.2f} MPa',
        'rho': f'{rho * 1e9 :.0f} kg/m3'
    }

    # Constraint limits
    sigma_max = 172.37  # MPa # 25000 psi
    delta_max = 50.8  # mm

    def __init__(self, prob_type="discrete"):

        super().__init__(name='ThreeBarTruss')

        self.binary_vars = None
        self.prob_type = prob_type.lower()
        self.structure = self.create_structure()
        self.create_constraints()
        self.create_objective()

    def create_objective(self):

        def objective(X):
            self.substitute_variables(X)
            weight = 0
            for mem in self.structure.members.values():
                weight += self.rho * mem.A * mem.length  # mem.weight

            return weight

        self.obj = objective

    def create_structure(self):

        frame = Frame2D(num_elements=1)
        # Nodes
        n1 = [0, 0]
        n2 = [self.L, 0]
        n3 = [2 * self.L, 0]
        n4 = [self.L, -self.L]
        # Create bars
        bar1 = SteelBeam([n1, n4])
        bar2 = SteelBeam([n2, n4])
        bar3 = SteelBeam([n3, n4])
        # Add bars
        frame.add(bar1)
        frame.add(bar2)
        frame.add(bar3)
        # Supports
        frame.add(XYHingedSupport(n1))
        frame.add(XYHingedSupport(n2))
        frame.add(XYHingedSupport(n3))
        # Point loads
        frame.add(PointLoad(n4, [-self.F, -self.F, 0]))
        # Change material properties
        for mem in frame.members.values():
            mem.fy = self.sigma_max
            mem.E = self.E  # MPa
            mem.A = THREE_BAR_AREAS_mm2[0]
            mem.Sj1 = 0
            mem.Sj2 = 0
        frame.generate()
        frame.calculate()

        if self.prob_type == "discrete":
            # Variables
            A1 = DiscreteVariable('A1', profiles=THREE_BAR_AREAS_mm2,
                                  target={"property": "A",
                                          "objects": [bar1, bar3]})
            A2 = DiscreteVariable('A2', profiles=THREE_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar2]})
            vars = [A1, A2]

        elif self.prob_type == "continuous":
            A1 = Variable('A1', lb=THREE_BAR_AREAS_mm2[0],
                          ub=THREE_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar1, bar3]})
            A2 = Variable('A2', lb=THREE_BAR_AREAS_mm2[0],
                          ub=THREE_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar2]})

            vars = [A1, A2]


        elif self.prob_type == 'binary':
            cont_A1 = Variable('A1', lb=THREE_BAR_AREAS_mm2[0],
                          ub=THREE_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar1, bar3]})
            cont_A2 = Variable('A2', lb=THREE_BAR_AREAS_mm2[0],
                          ub=THREE_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar2]})


            binary_A1 = []
            for i in range(len(THREE_BAR_AREAS_mm2)):
                binary = BinaryVariable("Bin A1"+str(i))
                binary.target_val = THREE_BAR_AREAS_mm2[i]
                binary_A1.append(binary)

            binary_A2 = []
            for i in range(len(THREE_BAR_AREAS_mm2)):
                binary = BinaryVariable("Bin A2"+str(i))
                binary.target_val = THREE_BAR_AREAS_mm2[i]
                binary_A2.append(binary)

            vars = [cont_A1, cont_A2]
            vars.extend(binary_A1)
            vars.extend(binary_A2)
            self.binary_vars = [binary_A1, binary_A2]

        self.vars = vars

        return frame


    def binary_constraints(self):

        # All available areas
        A = np.zeros((len(self.binary_vars), (len(THREE_BAR_AREAS_mm2))))

        for i, binary_A in enumerate(self.binary_vars):
            # Binary constraint
            # sum(bin_vars) == 1
            binary_idx = [self.vars.index(bvar) for bvar in binary_A]
            a = np.zeros(len(self.vars))
            a[binary_idx] = 1
            bin_con = LinearConstraint(a=a, b=1, con_type="=",
                                       name="Binary Constraint "
                                            + str(i+1))
            self.cons.append(bin_con)

            # Binary Area constraint
            # Ai == sum(A[i][bin_idx])
            a[binary_idx] = [bvar.target_val for bvar in binary_A]
            b = self.vars[i]
            bin_A_con = LinearConstraint(a=a, b=b, con_type="=",
                                         name="Binary Area Constraint "
                                         + str(i + 1))
            self.cons.append(bin_A_con)


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
            if not isinstance(var, BinaryVariable):
                for mem in var.target["objects"]:
                    if isinstance(mem, FrameMember):
                        i += 1

                        compression_fun, tension_fun, buckling_fun = self.constraint_generator(mem)

                        comp_con = NonLinearConstraint(con_fun=compression_fun,
                                                       name="Compression " + str(i),
                                                       parent=self)
                        comp_con.fea_required = True

                        tension_con = NonLinearConstraint(con_fun=tension_fun,
                                                          name="Tension " + str(i),
                                                          parent=self)
                        tension_con.fea_required = True

                        buckl_con = NonLinearConstraint(con_fun=buckling_fun,
                                                       name='Buckling ' + str(i),
                                                       parent=self)
                        buckl_con.fea_required = True

                        self.cons.append(comp_con)
                        self.cons.append(tension_con)
                        self.cons.append(buckl_con)

        if self.binary_vars:
            self.binary_constraints()

if __name__ == '__main__':
    from metku.optimization.solvers import *
    import matplotlib.pyplot as plt
    problem = ThreeBarTruss(prob_type='continuous')
    # X = [809, 7591, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    #      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    #      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # problem(X)
    # problem.structure.plot()




    solver = SLP(move_limits=[0.95, 1.2])
    solver.problem = problem
    solver.random_feasible_point()
    x0 = [THREE_BAR_AREAS_mm2[-1], THREE_BAR_AREAS_mm2[-1]]
    solver.solve(problem,
                 maxiter=15,
                 maxtime=5,
                 log=True)
    problem(solver.X)
    xvals = np.asarray(solver.xvals)
    A1 = np.array(x0[0])
    A2 = np.array(x0[1])
    A1 = np.hstack((A1, xvals[:, 0]))
    A2 = np.hstack((A2, xvals[:, 1]))
    X, Y = np.meshgrid(A1, A2)
    Z1 = problem.structure.members[0].length * X * problem.rho
    Z2 = problem.structure.members[1].length * Y * problem.rho
    Z = 2*Z1 + Z2

    plt.plot(A1, A2)
    plt.scatter(A1, A2, c='k', marker='x')
    CS = plt.contour(X, Y, Z, levels=10, cmap='Greens_r')
    plt.clabel(CS, inline=True)

    plt.show()


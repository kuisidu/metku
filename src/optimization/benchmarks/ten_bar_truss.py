
from src.optimization.structopt import *
from src.frame2d.frame2d import *
import numpy as np
from src.optimization.solvers import DiscreteVNS

ten_bar_DLM = [41, 0, 38, 31, 0, 0, 27, 38, 37, 0]
ten_bar_PSO = [40, 0, 38, 29, 0, 0, 27, 39, 37, 1]

TEN_BAR_AREAS_in2 = [1.62, 1.80, 1.99, 2.13, 2.38, 2.62, 2.63, 2.88, 2.93,
                    3.09, 3.13, 3.38, 3.47, 3.55, 3.63, 3.84, 3.87, 3.88,
                    4.18, 4.22, 4.49, 4.59, 4.80, 4.97, 5.12, 5.74, 7.22,
                    7.97, 11.50, 13.50, 13.90, 14.20, 15.50, 16.00, 16.90,
                    18.80, 19.90, 22.00, 22.90, 26.50, 30.00, 33.50]


TEN_BAR_AREAS_mm2 = [1044.9, 1161.0, 1283.55, 1373.85, 1535.1, 1689.9,
                     1696.35, 1857.6, 1889.85, 1993.05, 2018.85, 2180.1,
                     2238.15, 2289.75, 2341.35, 2476.8, 2496.15, 2502.6,
                     2696.1, 2721.9, 2896.05, 2960.55, 3096.0, 3205.65,
                     3302.4, 3702.3, 4656.9, 5140.65, 7417.5, 8707.5, 8965.5,
                     9159.0, 9997.5, 10320.0, 10900.5, 12126.0, 12835.5,
                     14190.0, 14770.5, 17092.5, 19350.0, 21607.5]


class TenBarTruss(OptimizationProblem):

    # Problem parameters
    L = 9144    # mm
    F = 444890  # N
    E = 68950   # MPa
    rho = 2768e-9  # kg/mm3

    properties = {
        'L': str(L) + " mm",
        'F': f'{F / 1e3 :.2f} kN',
        'E': f'{E / 1e3 :.2f} MPa',
        'rho': f'{rho * 1e9 :.0f} kg/m3'
    }

    # Constraint limits
    sigma_max = 172.37  # MPa
    delta_max = 50.8    # mm


    def __init__(self, prob_type="discrete"):

        super().__init__(name='TenBarTruss')
        self.prob_type = prob_type
        self.structure = self._create_structure()
        self._create_vars(profiles=TEN_BAR_AREAS_mm2)
        self.create_constraints()
        self.create_objective()


    def _create_vars(self, profiles=[0, 1e6]):
        """
        Creates variables used in optimization

        Appends each variable to 'vars' 'list from where they can be accessed
        """
        for i, mem in enumerate(self.structure.members.values()):
            name = 'A' + str(i + 1)
            if self.prob_type == "discrete":
                var = DiscreteVariable(name,
                                    profiles=profiles,
                                    target={"property": "A",
                                            "objects": [mem]})

            elif self.prob_type == "continuous":
                var = Variable(name,
                               lb=profiles[0],
                               ub=profiles[-1],
                               target={"property": "A",
                                       "objects": [mem]})

            else:

                raise TypeError("Problem type must be either 'discrete' "
                                "or 'continuous")


            self.vars.append(var)


    def _create_structure(self):
        """
        Creates the ten bar truss structure

        :return: structure
        :rtype: Frame2D
        """

        frame = Frame2D(num_elements=1)
        # Members
        n1 = [2 * self.L, self.L]
        n2 = [2 * self.L, 0]
        n3 = [self.L, self.L]
        n4 = [self.L, 0]
        n5 = [0, 0]
        n6 = [0, self.L]
        bar1 = SteelBeam([n6, n3])
        bar2 = SteelBeam([n3, n1])
        bar3 = SteelBeam([n5, n4])
        bar4 = SteelBeam([n4, n2])
        bar5 = SteelBeam([n3, n4])
        bar6 = SteelBeam([n1, n2])
        bar7 = SteelBeam([n6, n4])
        bar8 = SteelBeam([n5, n3])
        bar9 = SteelBeam([n3, n2])
        bar10 = SteelBeam([n4, n1])
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
        # Supports
        frame.add(XYHingedSupport(n5))
        frame.add(XYHingedSupport(n6))
        # Point loads
        frame.add(PointLoad(n4, [0, -self.F, 0]))
        frame.add(PointLoad(n2, [0, -self.F, 0]))
        # Change material properties
        for mem in frame.members.values():
            mem.cross_section.fy = self.sigma_max
            mem.E = self.E  # MPa
            mem.A = TEN_BAR_AREAS_mm2[0]
            mem.Sj1 = 0
            mem.Sj2 = 0
        frame.generate()
        frame.calculate()

        return frame


    def constraint_generator(self, mem):

        def compression_fun(x):

            return -mem.ned / (mem.A * mem.fy) - 1

        def tension_fun(x):

            return mem.ned / (mem.A * mem.fy) - 1

        def disp_fun(x):
            displacements = mem.nodal_displacements.values()
            max_vals = [max(l[0:2]) for l in displacements]
            min_vals = [min(l[0:2]) for l in displacements]
            max_val = max(max_vals)
            min_val = min(min_vals)
            abs_max = max(max_val, abs(min_val))
            return abs_max / self.delta_max - 1

        return compression_fun, tension_fun, disp_fun

    def create_constraints(self):

        # Initialize constraints as an empty list
        self.cons = []
        i = 0
        for var in self.vars:
            for mem in var.target["objects"]:
                if isinstance(mem, FrameMember):
                    i += 1
                    compression_fun, tension_fun, disp_fun = self.constraint_generator(mem)


                    comp_con = NonLinearConstraint(con_fun=compression_fun,
                                                 name="Compression " + str(i),
                                                 parent=self)
                    comp_con.fea_required = True

                    tension_con = NonLinearConstraint(con_fun=tension_fun,
                                                     name="Tension " + str(i),
                                                     parent=self)
                    tension_con.fea_required = True

                    disp_con = NonLinearConstraint(con_fun=disp_fun,
                                                 name='Displacement ' + str(i),
                                                 parent=self)
                    disp_con.fea_required = True
                    self.cons.append(comp_con)
                    self.cons.append(tension_con)
                    self.cons.append(disp_con)

    def create_objective(self):

        def objective(X):
            if np.any(self.X != X):
                self.substitute_variables(X)
            weight = 0
            for mem in self.structure.members.values():
                weight += self.rho * mem.A * mem.length

            return weight

        self.obj = objective


if __name__ == '__main__':
    import time
    sub_sizes = np.arange(2, 10)
    step_lengths = np.arange(1, 5)

    # for step_length in step_lengths:
    #     for sub_size in sub_sizes:
    problem = TenBarTruss()
    problem.structure.plot()
    # solver = DiscreteVNS(step_length=3)
    # start = time.time()
    # x, f, steps = solver.solve(problem,
    #                          maxiter=1000,
    #                          maxtime=300,
    #                          subset_size=3)
    # end = time.time()
    # total_time = round(end - start, 2)
    # with open("tenbartruss_results.csv", 'a') as f:
    #     f.write(f'{x};'
    #             f'{round(f,2)};'
    #             f'{len(steps)};'
    #             f'{step_length};'
    #             f'{sub_size};'
    #             f'{total_time}')


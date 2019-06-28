from optimization.structopt import *
from frame2d.frame2d import *
from optimization.optsolver import VNS

ten_bar_x_star = [41, 0, 38, 31, 0, 0, 27, 38, 37, 0]

TEN_BAR_AREAS_mm2 = [1044.9, 1161.0, 1283.55, 1373.85, 1535.1, 1689.9,
                     1696.35, 1857.6, 1889.85, 1993.05, 2018.85, 2180.1,
                     2238.15, 2289.75, 2341.35, 2476.8, 2496.15, 2502.6,
                     2696.1, 2721.9, 2896.05, 2960.55, 3096.0, 3205.65,
                     3302.4, 3702.3, 4656.9, 5140.65, 7417.5, 8707.5, 8965.5,
                     9159.0, 9997.5, 10320.0, 10900.5, 12126.0, 12835.5,
                     14190.0, 14770.5, 17092.5, 19350.0, 21607.5]


class ThreeBarTruss(OptimizationProblem):

    def __init__(self, name="",variables=[],constraints=[],objective=None,
                 gradient=None,hess=None, structure=None, profiles=None):

        super().__init__('ThreeBarTruss',variables,constraints,objective,
                        gradient,hess, structure, profiles)


class TenBarTruss(OptimizationProblem):

    # Problem parameters
    L = 9144    # mm
    F = 444890  # N
    E = 68950   # MPa
    rho = 2768e3  # kg/mm3

    # Variables
    A1 = IntegerVariable('A1', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A2 = IntegerVariable('A2', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A3 = IntegerVariable('A3', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A4 = IntegerVariable('A4', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A5 = IntegerVariable('A5', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A6 = IntegerVariable('A6', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A7 = IntegerVariable('A7', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A8 = IntegerVariable('A8', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A9 = IntegerVariable('A9', 0, len(TEN_BAR_AREAS_mm2) - 1)
    A10 = IntegerVariable('A10', 0, len(TEN_BAR_AREAS_mm2) - 1)

    int_vars = [A1, A2, A3, A4, A5,
                A6, A7, A8, A9, A10]

    # Constraint limits
    sigma_max = 172.37  # MPa
    delta_max = 50.8    # mm


    def __init__(self):

        super().__init__(name='TenBarTruss',
                         structure=self.create_structure(),
                         variables=self.int_vars,
                         objective=lambda x: self.rho * sum(x))
        self.create_constraints()


    def create_structure(self):

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
        bar7 = SteelBeam([n5, n3])
        bar8 = SteelBeam([n6, n4])
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
            mem.cross_section.E = self.E  # MPa
            mem.cross_section.A = TEN_BAR_AREAS_mm2[0]
            mem.Sj1 = 0
            mem.Sj2 = 0

        frame.generate()
        frame.calculate()
        return frame


    def create_constraints(self):

        # Initialize constraints as an empty list
        self.cons = []
        # Create stress constraints
        for i, mem in enumerate(self.structure.members.values()):
            def con_fun(A, i=i):
                idx = list(self.structure.members.keys()).index(mem.mem_id)
                mem.cross_section.A = A[idx]
                return mem.NRd
            constr = NonLinearConstraint(con_fun, name='g' + str(i+1))
            self.cons.append(constr)

        # Create displacement constraints
        for j, mem in enumerate(self.structure.members.values()):
            def con_fun(j=j):
                displacements = mem.nodal_displacements.values()
                max_vals = [max(l[0:2]) for l in displacements]
                min_vals = [min(l[0:2]) for l in displacements]
                max_val = max(max_vals)
                min_val = min(min_vals)
                abs_max = max(max_val, abs(min_val))
                return abs_max - self.delta_max
            constr = NonLinearConstraint(con_fun, name='g' + str(11 + j))
            self.cons.append(constr)


if __name__ == '__main__':
    problem = TenBarTruss()
    solver = VNS()
    #f_opt, x_opt = solver.solve(problem)


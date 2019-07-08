from optimization.structopt import *
from frame2d.frame2d import *
from optimization.optsolver import DiscreteVNS

ten_bar_DLM = [41, 0, 38, 31, 0, 0, 27, 38, 37, 0]
ten_bar_PSO = [40, 0, 38, 29, 0, 0, 27, 39, 37, 1]

TEN_BAR_AREAS_in2 = [1.62, 1.80, 1.99, 2.13, 2.38, 2.62, 2.63, 2.88, 2.93,
                    3.09, 3.13, 3.38, 3.47, 3.55, 3.63, 3.84, 3.87, 3.88,
                    4.18, 4.22, 4.49, 4.59, 4.80, 4.97, 5.12, 5.74, 7.22,
                    7.97, 11.50, 13.50, 13.90, 14.20, 15.50, 16.00, 16.90,
                    18.80, 19.90, 22.00, 22.90, 26.50, 30.00, 33.50,]


TEN_BAR_AREAS_mm2 = [1044.9, 1161.0, 1283.55, 1373.85, 1535.1, 1689.9,
                     1696.35, 1857.6, 1889.85, 1993.05, 2018.85, 2180.1,
                     2238.15, 2289.75, 2341.35, 2476.8, 2496.15, 2502.6,
                     2696.1, 2721.9, 2896.05, 2960.55, 3096.0, 3205.65,
                     3302.4, 3702.3, 4656.9, 5140.65, 7417.5, 8707.5, 8965.5,
                     9159.0, 9997.5, 10320.0, 10900.5, 12126.0, 12835.5,
                     14190.0, 14770.5, 17092.5, 19350.0, 21607.5]


FIFTEEN_BAR_AREAS_in2 = [0.111, 0.141, 0.174, 0.220, 0.270, 0.287, 0.347, 0.440,
                         0.539, 0.954, 1.081, 1.174, 1.333, 1.488, 1.764, 2.142,
                         2.697, 2.800, 3.131, 3.565, 3.813, 4.805, 5.952, 6.572,
                         7.192, 8.525, 9.300, 10.850, 13.300, 14.290, 17.170, 19.180]

in2_to_mm2 = 645.16

FIFTEEN_BAR_AREAS_mm2 = [71.613, 90.968, 112.258, 141.935, 174.193, 185.161, 223.871,
                         283.87, 347.741, 615.483, 697.418, 757.418, 859.998, 959.998,
                         1138.062, 1381.933, 1739.997, 1806.448, 2019.996, 2299.995,
                         2459.995, 3099.994, 3839.992, 4239.992, 4639.991, 5499.989,
                         5999.988, 6999.986, 8580.628, 9219.336, 11077.397, 12374.169]




class ThreeBarTruss(OptimizationProblem):

    pass

class FifteenBarTruss(OptimizationProblem):
    # Problem parameters
    L = 3048  # mm
    F = 444890  # N # 1e5 lb
    E = 68950  # MPa # 10e7 psi
    rho = 2768e-9  # kg/mm3

    # Variables
    A1 = IntegerVariable('A1', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A2 = IntegerVariable('A2', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A3 = IntegerVariable('A3', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A4 = IntegerVariable('A4', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A5 = IntegerVariable('A5', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A6 = IntegerVariable('A6', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A7 = IntegerVariable('A7', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A8 = IntegerVariable('A8', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A9 = IntegerVariable('A9', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A10 = IntegerVariable('A10', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A11 = IntegerVariable('A11', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A12 = IntegerVariable('A12', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A13 = IntegerVariable('A13', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A14 = IntegerVariable('A14', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)
    A15 = IntegerVariable('A15', 0, len(FIFTEEN_BAR_AREAS_mm2) - 1)

    int_vars = [A1, A2, A3, A4, A5,
                A6, A7, A8, A9, A10,
                A11, A12, A13, A14, A15]

    # Constraint limits

    sigma_max = 172.37  # MPa # 25000 psi
    delta_max = 50.8  # mm


    def __init__(self):

        super().__init__(name='TenBarTruss',
                         structure=self.create_structure(),
                         variables=self.int_vars)
        self.create_constraints()
        self.create_objective()


    def create_objective(self):

        def objective(X):
            if isinstance(X[0], IntegerVariable):
                X = [x.value for x in X]

            weight = 0
            for x, mem in zip(X, self.structure.members.values()):
                x = int(x)
                mem.A = TEN_BAR_AREAS_mm2[x]
                weight += self.rho * mem.A * mem.length
            # Calculate stucture with new values
            self.structure.calculate()
            return weight

        self.obj = objective

    def create_structure(self):

        frame = Frame2D(num_elements=1)
        # Nodes
        n1 = [0, self.L]
        n2 = [self.L, self.L]
        n3 = [2 * self.L, self.L]
        n4 = [3*self.L, self.L]
        n5 = [0, 0]
        n6 = [self.L, 0]
        n7 = [2 * self.L, 0]
        n8 = [3*self.L, 0]
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
            mem.cross_section.fy = self.sigma_max
            mem.E = self.E  # MPa
            mem.A = FIFTEEN_BAR_AREAS_mm2[0]
            mem.Sj1 = 0
            mem.Sj2 = 0
        frame.generate()
        frame.calculate()
        return frame

    def create_constraints(self):

        # Initialize constraints as an empty list
        self.cons = []
        # Create stress constraints
        for i in range(len(self.structure.members)):
            def con_fun(A, i=i):
                mem = self.structure.members[i]
                mem.A = TEN_BAR_AREAS_mm2[int(A[i])]
                return mem.ned / mem.NRd - 1
            constr = NonLinearConstraint(con_fun, name='Stress ' + str(i+1))
            self.cons.append(constr)

        # Create buckling constraints
        for j in range(len(self.structure.members)):
            def con_fun(A, j=j):
                mem = self.structure.members[j]
                sigma_cr = 100 * mem.E * mem.A / (8 * mem.length**2)
                sigma = abs(mem.ned) / mem.A
                return  sigma/sigma_cr  - 1
            constr = NonLinearConstraint(con_fun, name='Buckling ' + str(11 + j))
            self.cons.append(constr)


class TenBarTruss(OptimizationProblem):

    # Problem parameters
    L = 9144    # mm
    F = 444890  # N
    E = 68950   # MPa
    rho = 2768e-9  # kg/mm3

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
                         variables=self.int_vars)
        self.create_constraints()
        self.create_objective()


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


    def create_constraints(self):

        # Initialize constraints as an empty list
        self.cons = []
        # Create stress constraints
        for i in range(len(self.structure.members)):
            def con_fun(A, i=i):
                mem = self.structure.members[i]
                mem.A = TEN_BAR_AREAS_mm2[int(A[i])]
                return mem.ned / mem.NRd - 1
            constr = NonLinearConstraint(con_fun, name='Stress ' + str(i+1))
            self.cons.append(constr)

        # Create displacement constraints
        for j in range(len(self.structure.members)):
            def con_fun(A, j=j):
                mem = self.structure.members[j]
                displacements = mem.nodal_displacements.values()
                max_vals = [max(l[0:2]) for l in displacements]
                min_vals = [min(l[0:2]) for l in displacements]
                max_val = max(max_vals)
                min_val = min(min_vals)
                abs_max = max(max_val, abs(min_val))
                return abs_max / self.delta_max - 1
            constr = NonLinearConstraint(con_fun, name='Displacement ' + str(11 + j))
            self.cons.append(constr)

    def create_objective(self):

        def objective(X):
            if isinstance(X[0], IntegerVariable):
                X = [x.value for x in X]

            weight = 0
            for x, mem in zip(X, self.structure.members.values()):
                x = int(x)
                mem.A = TEN_BAR_AREAS_mm2[x]
                weight += self.rho * mem.A * mem.length
            # Calculate stucture with new values
            self.structure.calculate()
            return weight

        self.obj = objective

if __name__ == '__main__':
    problem = TenBarTruss()
    problem = FifteenBarTruss()
    solver = DiscreteVNS(step_length=3)
    solver.solve(problem, maxiter=int(1e3), maxtime=300)


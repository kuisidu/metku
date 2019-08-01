from gekko import GEKKO
from frame2d.frame2d import *


TEN_BAR_AREAS_mm2 = [1044.9, 1161.0, 1283.55, 1373.85, 1535.1, 1689.9,
                     1696.35, 1857.6, 1889.85, 1993.05, 2018.85, 2180.1,
                     2238.15, 2289.75, 2341.35, 2476.8, 2496.15, 2502.6,
                     2696.1, 2721.9, 2896.05, 2960.55, 3096.0, 3205.65,
                     3302.4, 3702.3, 4656.9, 5140.65, 7417.5, 8707.5, 8965.5,
                     9159.0, 9997.5, 10320.0, 10900.5, 12126.0, 12835.5,
                     14190.0, 14770.5, 17092.5, 19350.0, 21607.5]


class TenBarTruss(GEKKO):

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


    def __init__(self):

        super().__init__(name='TenBarTruss')
        self.create_vars()
        self.create_constraints()
        self.create_objective()
        self.structure = self.create_structure()



    def create_vars(self):

        # Variables
        # Var(value, lb, ub, integer, name)
        A1 = self.Var(name='A1', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A2 = self.Var(name='A2', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A3 = self.Var(name='A3', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A4 = self.Var(name='A4', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A5 = self.Var(name='A5', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A6 = self.Var(name='A6', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A7 = self.Var(name='A7', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A8 = self.Var(name='A8', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A9 = self.Var(name='A9', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)
        A10 = self.Var(name='A10', lb=0, ub=len(TEN_BAR_AREAS_mm2) - 1, integer=True)


        self.vars = [A1, A2, A3, A4, A5,
                    A6, A7, A8, A9, A10]


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


        for i in range(len(self.structure.members)):
            self.Equation(self.vars[i])

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
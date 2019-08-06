
from src.optimization.structopt import *
from src.frame2d.frame2d import *

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

        self.prob_type = prob_type.lower()
        self.structure = self.create_structure()
        self.create_constraints()
        self.create_objective()

    def create_objective(self):

        def objective(X):
            if np.any(self.X != X):
                self.substitute_variables(X)
            weight = 0
            for x, mem in zip(X, self.structure.members.values()):
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
            mem.cross_section.fy = self.sigma_max
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

        self.vars = vars

        return frame


    def _create_stress_constraint(self, mem, i, j):
        """
        THIS DOES NOT WORK YET!
        GIVES DIFFERENT SOLUTIONS THAN THE FUNCTION CREATED
        IN THE LOOP!
        :param mem:
        :param i:
        :param j:
        :return:
        """
        def con_fun(x, i=i, j=j):

            if self.prob_type == 'discrete':
                A = THREE_BAR_AREAS_mm2[x[j]]
            else:
                A = x[j]

            return mem.ned / (A * mem.fy) - 1

        return con_fun

    def create_constraints(self):

        # Initialize constraints as an empty list
        self.cons = []

        i = 0
        for j, var in enumerate(self.vars):
            for mem in var.target["objects"]:
                if isinstance(mem, FrameMember):
                    i += 1
                    def stress_fun(x, i=i, j=j):

                        return mem.ned / (mem.A * mem.fy) - 1

                    def buckling_fun(x, i=i, j=j):
                        sigma_cr = 100 * mem.E * mem.A / (8 * mem.length ** 2)
                        sigma = -mem.ned / mem.A

                        return sigma / sigma_cr - 1

                    def disp_fun(A, i=i):
                        displacements = mem.nodal_displacements.values()
                        max_vals = [max(l[0:2]) for l in displacements]
                        min_vals = [min(l[0:2]) for l in displacements]
                        max_val = max(max_vals)
                        min_val = min(min_vals)
                        abs_max = max(max_val, abs(min_val))
                        return abs_max / self.delta_max - 1

                    stress_con = NonLinearConstraint(con_fun=stress_fun,
                                                 name="Stress " + str(i),
                                                 parent=self)
                    stress_con.fea_required = True

                    buckling_con = NonLinearConstraint(con_fun=buckling_fun,
                                                     name="Buckling " + str(i),
                                                     parent=self)
                    buckling_con.fea_required = True

                    disp_con = NonLinearConstraint(con_fun=disp_fun,
                                                 name='Displacement ' + str(i),
                                                 parent=self)
                    disp_con.fea_required = True

                    self.cons.append(stress_con)
                    self.cons.append(buckling_con)
                    self.cons.append(disp_con)




if __name__ == '__main__':
    from src.optimization.solvers import SLP, SLSQP
    problem = ThreeBarTruss(prob_type='continuous')
    solver = SLP()
    solver.solve(problem, maxiter=200)
    problem(solver.X)

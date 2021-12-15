
try:
    from metku.optimization.structopt import *
    from metku.frame2d.frame2d import *
    from metku.frame2d.simple_truss import *
except:
    from optimization.structopt import *
    from frame2d.frame2d import *
#    from optimization.solvers import DiscreteVNS

import numpy as np

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
# Problem parameters
L = 9144  # mm
h = 3048 # mm
F = 444890  # N
E = 68950  # MPa
rho = 2768e-9  # kg/mm3

# Constraint limits
sigma_max = 172.37  # MPa # 25000 psi
delta_max = 50.8  # mm

class SimpleFifteenBarTruss(OptimizationProblem):

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
        self.structure = self.create_structure()
        self.create_variables(profiles=FIFTEEN_BAR_AREAS_mm2)
        self.create_constraints()
        self.create_objective()


    def create_variables(self, profiles=[0, 1e6]):

        """
        Creates variables used in optimization

        Appends each variable to 'vars' 'list from where they can be accessed
        """
        self.vars = []
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

            elif self.prob_type == 'index':
                var = IndexVariable(name,
                                   values=profiles,
                                   target={"property": "A",
                                           "objects": [mem]})


            else:

                raise TypeError("Problem type must be either 'discrete' "
                                "or 'continuous")


            self.add(var)

    def create_structure(self):

        nodes = []
        nodes.append([0, h])
        nodes.append([L / 3, h])
        nodes.append([2 * L / 3, h])
        nodes.append([L, h])
        nodes.append([0, 0])
        nodes.append([L / 3, 0])
        nodes.append([2 * L / 3, 0])
        nodes.append([L, 0])

        mems = [[0, 1], [1, 2], [2, 3], [4, 5], [5, 6], [6, 7],
                [5, 1], [6, 2], [7, 3],
                [0, 5], [4, 1],
                [1, 6], [5, 2],
                [2, 7], [6, 3]]

        profs = []
        for i in range(len(mems)):
            profs.append('SHS 50x50x3.0')

        t = SimpleTruss(nodes, mems, profs)

        t.add(PointLoad(nodes[7], [0, -F, 0]))

        t.add(XYHingedSupport(nodes[0]))
        t.add(XYHingedSupport(nodes[4]))

        # Change material properties
        for mem in t.members.values():
            mem.material.fy = self.sigma_max
            mem.E = self.E  # MPa
            mem.A = FIFTEEN_BAR_AREAS_mm2[0]
            mem.rho = self.rho

        t.generate()
        t.calculate()
        return t

    def constraint_generator(self, mem):

        def compression_fun(x):
            return -mem.ned / (mem.A * mem.fy) - 1

        def tension_fun(x):
            return mem.ned / (mem.A * mem.fy) - 1

        def disp_fun(x):
            displacementsList = list(mem.nodal_displacements.values())
            displacementsMultipleArrays = displacementsList[0].values()

            max_vals = [max(array) for array in displacementsMultipleArrays]
            min_vals = [min(array) for array in displacementsMultipleArrays]

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
                                                   )
                    comp_con.fea_required = True

                    tension_con = NonLinearConstraint(con_fun=tension_fun,
                                                      name="Tension " + str(i))
                    tension_con.fea_required = True

                    disp_con = NonLinearConstraint(con_fun=disp_fun,
                                                   name='Displacement ' + str(i),
                                                   )
                    disp_con.fea_required = True
                    self.add(comp_con)
                    self.add(tension_con)
                    self.add(disp_con)

    def create_objective(self):

        def objective(X):
            self.substitute_variables(X)
            weight = 0
            for mem in self.structure.members.values():
                weight += mem.weight
            return weight

        obj = ObjectiveFunction(name="Weight",
                                obj_fun=objective)

        self.add(obj)

def convert_to_Excel(list):
    '''
    Save cross section areas in result.xlsx - file.
    :param list: list[0] == weight of truss, list[1] == cross section areas
    '''
    # Importing pandas as pd
    import pandas as pd
    df = pd.DataFrame()

    # Creating column
    df['Cross sections of best truss'] = list[1]

    # Converting to excel
    df.to_excel('result.xlsx', index=False)

if __name__ == '__main__':

    from metku.optimization.solvers import ga
    from metku.optimization.solvers import pso
    from metku.optimization.solvers import aco

    problem = SimpleFifteenBarTruss('discrete')
    solver = aco.AntColonyOptimizer(pop_size=100)

    fopt, xopt = solver.solve(problem, maxiter=100)
    print(fopt, xopt)
    problem(xopt)
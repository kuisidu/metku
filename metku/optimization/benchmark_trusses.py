# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
from frame2d.frame2d import *
from optimization.structopt import *
from sections.steel.catalogue import rhs_profiles, ipe_profiles
from truss2d import *
from optimization.solvers import *

RHS_PROFILES = list(rhs_profiles.keys())

# Sorted by area
RHS_PROFILES = ['RHS 40X40X2.0', 'RHS 40X40X2.5', 'RHS 50X50X2.0',
                'RHS 40X40X3.0', 'RHS 60X60X2.0', 'RHS 50X50X2.5',
                'RHS 70X70X2.0', 'RHS 50X50X3.0', 'RHS 60X60X2.5',
                'RHS 80X80X2.0', 'RHS 70X70X2.5', 'RHS 60X60X3.0',
                'RHS 50X50X4.0', 'RHS 80X80X2.5', 'RHS 70X70X3.0',
                'RHS 60X60X4.0', 'RHS 90X90X2.5', 'RHS 80X80X3.0',
                'RHS 100X100X2.5', 'RHS 70X70X4.0', 'RHS 90X90X3.0',
                'RHS 60X60X5.0', 'RHS 100X100X3.0', 'RHS 120X120X2.5',
                'RHS 80X80X4.0', 'RHS 70X70X5.0', 'RHS 90X90X4.0',
                'RHS 140X140X2.5', 'RHS 120X120X3.0', 'RHS 80X80X5.0',
                'RHS 150X150X2.5', 'RHS 100X100X4.0', 'RHS 140X140X3.0',
                'RHS 90X90X5.0', 'RHS 150X150X3.0', 'RHS 120X120X4.0',
                'RHS 100X100X5.0', 'RHS 160X160X3.0', 'RHS 140X140X4.0',
                'RHS 100X100X6.0', 'RHS 120X120X5.0', 'RHS 150X150X4.0',
                'RHS 160X160X4.0', 'RHS 140X140X5.0', 'RHS 120X120X6.0',
                'RHS 180X180X4.0', 'RHS 150X150X5.0', 'RHS 120X120X7.1',
                'RHS 160X160X5.0', 'RHS 200X200X4.0', 'RHS 140X140X6.0',
                'RHS 150X150X6.0', 'RHS 120X120X8.0', 'RHS 180X180X5.0',
                'RHS 160X160X6.0', 'RHS 120X120X8.8', 'RHS 180X180X5.6',
                'RHS 200X200X5.0', 'RHS 150X150X7.1', 'RHS 120X120X10.0',
                'RHS 180X180X6.0', 'RHS 160X160X7.1', 'RHS 220X220X5.0',
                'RHS 150X150X8.0', 'RHS 200X200X6.0', 'RHS 160X160X8.0',
                'RHS 150X150X8.8', 'RHS 180X180X7.1', 'RHS 250X250X5.0',
                'RHS 220X220X6.0', 'RHS 160X160X8.8', 'RHS 150X150X10.0',
                'RHS 180X180X8.0', 'RHS 200X200X7.1', 'RHS 160X160X10.0',
                'RHS 180X180X8.8', 'RHS 250X250X6.0', 'RHS 300X300X5.0',
                'RHS 220X220X7.1', 'RHS 200X200X8.0', 'RHS 260X260X6.0',
                'RHS 180X180X10.0', 'RHS 200X200X8.8', 'RHS 220X220X8.0',
                'RHS 250X250X7.1', 'RHS 300X300X6.0', 'RHS 260X260X7.1',
                'RHS 220X220X8.8', 'RHS 200X200X10.0', 'RHS 250X250X8.0',
                'RHS 180X180X12.5', 'RHS 260X260X8.0', 'RHS 220X220X10.0',
                'RHS 300X300X7.1', 'RHS 250X250X8.8', 'RHS 260X260X8.8',
                'RHS 200X200X12.5', 'RHS 300X300X8.0', 'RHS 250X250X10.0',
                'RHS 400X400X6.0', 'RHS 260X260X10.0', 'RHS 300X300X8.8',
                'RHS 400X400X7.1', 'RHS 250X250X12.5', 'RHS 300X300X10.0',
                'RHS 260X260X12.5', 'RHS 400X400X8.0', 'RHS 400X400X8.8',
                'RHS 300X300X12.5', 'RHS 400X400X10.0', 'RHS 400X400X12.5']

RHS_A = [293.6991118430775, 358.9048622548086, 373.6991118430775,
         420.8230016469244, 453.6991118430775, 458.9048622548086,
         533.6991118430775, 540.8230016469245, 558.9048622548087,
         613.6991118430775, 658.9048622548087, 660.8230016469245,
         694.79644737231, 758.9048622548087, 780.8230016469245,
         854.79644737231, 858.9048622548087, 900.8230016469245,
         958.9048622548087, 1014.79644737231, 1020.8230016469245,
         1035.6194490192345, 1140.8230016469245, 1158.9048622548087,
         1174.79644737231, 1235.6194490192345, 1334.79644737231,
         1358.9048622548087, 1380.8230016469245, 1435.6194490192345,
         1458.9048622548087, 1494.79644737231, 1620.8230016469245,
         1635.6194490192345, 1740.8230016469245, 1814.79644737231,
         1835.6194490192345, 1860.8230016469245, 2134.79644737231,
         2163.292006587698, 2235.6194490192347, 2294.79644737231,
         2454.79644737231, 2635.6194490192347, 2643.292006587698,
         2774.79644737231, 2835.6194490192347, 3033.2707426698457,
         3035.6194490192347, 3094.79644737231, 3123.292006587698,
         3363.292006587698, 3364.247719318987, 3435.6194490192347,
         3603.292006587698, 3648.3397403759745, 3825.8010368497276,
         3835.6194490192347, 3885.2707426698457, 4056.6370614359175,
         4083.292006587698, 4169.270742669845, 4235.619449019235,
         4324.247719318987, 4563.292006587698, 4644.247719318987,
         4704.3397403759745, 4737.270742669846, 4835.619449019235,
         5043.292006587698, 5056.3397403759745, 5256.6370614359175,
         5284.247719318987, 5305.270742669846, 5656.6370614359175,
         5760.339740375975, 5763.292006587698, 5835.619449019235,
         5873.270742669846, 5924.247719318987, 6003.292006587698,
         6456.6370614359175, 6464.339740375975, 6564.247719318987,
         6725.270742669845, 6963.292006587698, 7009.270742669845,
         7168.339740375975, 7256.6370614359175, 7524.247719318987,
         7704.369260617026, 7844.247719318987, 8056.6370614359175,
         8145.270742669845, 8224.339740375975, 8576.339740375975,
         8704.369260617026, 9124.247719318988, 9256.637061435917,
         9363.292006587697, 9656.637061435917, 9984.339740375975,
         10985.270742669845, 11204.369260617026, 11256.637061435917,
         11704.369260617026, 12324.247719318988, 13504.339740375975,
         13704.369260617026, 15256.637061435917, 18704.369260617026]

RHS_Iy = [69402.1348577099, 82151.28824369762, 141469.80388201258,
          93235.56709132508, 251422.42849846912, 169438.8058049558,
          407260.0087070795, 194671.3623630676, 303421.5664789544,
          616982.5445078439, 494099.5702656935, 351348.30771715636,
          237358.76890471048, 751472.8171651728, 575266.4031535913,
          435510.7428089776, 1085541.3071773928, 878425.6486723726,
          1506305.0403023532, 721202.5390818602, 1272826.0442734999,
          504944.22072287824, 1770467.5899569737, 2647918.2358904947,
          1110434.1577233584, 846291.9300855394, 1619205.598733472,
          4256312.403929599, 3123474.1315709595, 1314420.611899162,
          5260552.35261826, 2263516.8621122013, 5033445.27351433,
          1929330.2661637464, 6227292.569609535, 4022758.855975506,
          2711020.892879293, 7596381.015787086, 6516160.139313272,
          3114741.7978090816, 4854745.06366327, 8078170.514535079,
          9871720.712125503, 7905593.124251096, 5621572.923474502,
          14217440.574412191, 9821188.61322145, 6235230.294796982,
          12023565.074642764, 19681319.726173345, 9204262.450457461,
          11459054.114443017, 6768760.707404428, 17368660.914838284,
          14054810.378757961, 7198754.151445843, 19184945.76423401,
          24100880.64483765, 12896997.263231725, 7768082.442480799,
          20365216.708375998, 15874082.660310294, 32380224.26464086,
          14118333.707433986, 28327481.43931158, 17412349.479375735,
          15131229.417059045, 23133398.0658679, 48050096.98772789,
          38133604.5715647, 18695912.47963438, 16525294.695811596,
          25458618.181157082, 32322239.619959485, 20476695.81973211,
          27417265.565841436, 56720023.77241476, 84168837.64189216,
          43667807.32258503, 35662536.426802225, 64045426.04002355,
          30167993.626788545, 38495934.600123696, 48280104.21631118,
          65227020.405024536, 99636681.11375256, 73737921.17343803,
          52213519.582481146, 42510618.84613217, 72292048.7953192,
          34064339.66473276, 81780188.42692047, 57824571.47776296,
          115161339.61842687, 78353864.45865832, 88691837.39142165,
          48594155.981109016, 128006870.81298491, 87066739.32324764,
          241042340.82113203, 98646458.97788613, 139105018.9926629,
          280318593.3302435, 101613144.87508953, 155193656.12715802,
          115478848.04297818, 312692443.7957626, 340887002.087082,
          183481345.34484133, 382159878.71536344, 458765381.0116587]

IPE_PROFILES = list(ipe_profiles.keys())
ten_bar_DLM = [41, 0, 38, 31, 0, 0, 27, 38, 37, 0]
ten_bar_PSO = [40, 0, 38, 29, 0, 0, 27, 39, 37, 1]

TEN_BAR_AREAS_in2 = [1.62, 1.80, 1.99, 2.13, 2.38, 2.62, 2.63, 2.88, 2.93,
                     3.09, 3.13, 3.38, 3.47, 3.55, 3.63, 3.84, 3.87, 3.88,
                     4.18, 4.22, 4.49, 4.59, 4.80, 4.97, 5.12, 5.74, 7.22,
                     7.97, 11.50, 13.50, 13.90, 14.20, 15.50, 16.00, 16.90,
                     18.80, 19.90, 22.00, 22.90, 26.50, 30.00, 33.50, ]

TEN_BAR_AREAS_mm2 = [1044.9, 1161.0, 1283.55, 1373.85, 1535.1, 1689.9,
                     1696.35, 1857.6, 1889.85, 1993.05, 2018.85, 2180.1,
                     2238.15, 2289.75, 2341.35, 2476.8, 2496.15, 2502.6,
                     2696.1, 2721.9, 2896.05, 2960.55, 3096.0, 3205.65,
                     3302.4, 3702.3, 4656.9, 5140.65, 7417.5, 8707.5, 8965.5,
                     9159.0, 9997.5, 10320.0, 10900.5, 12126.0, 12835.5,
                     14190.0, 14770.5, 17092.5, 19350.0, 21607.5]

FIFTEEN_BAR_AREAS_in2 = [0.111, 0.141, 0.174, 0.220, 0.270, 0.287, 0.347,
                         0.440,
                         0.539, 0.954, 1.081, 1.174, 1.333, 1.488, 1.764,
                         2.142,
                         2.697, 2.800, 3.131, 3.565, 3.813, 4.805, 5.952,
                         6.572,
                         7.192, 8.525, 9.300, 10.850, 13.300, 14.290, 17.170,
                         19.180]

in2_to_mm2 = 645.16

FIFTEEN_BAR_AREAS_mm2 = [71.613, 90.968, 112.258, 141.935, 174.193, 185.161,
                         223.871,
                         283.87, 347.741, 615.483, 697.418, 757.418, 859.998,
                         959.998,
                         1138.062, 1381.933, 1739.997, 1806.448, 2019.996,
                         2299.995,
                         2459.995, 3099.994, 3839.992, 4239.992, 4639.991,
                         5499.989,
                         5999.988, 6999.986, 8580.628, 9219.336, 11077.397,
                         12374.169]

THREE_BAR_AREAS_mm2 = FIFTEEN_BAR_AREAS_mm2


class DiscreteFrame(OptimizationProblem):
    """
    General class for creating optimization problems that utilize the Frame2D -module
    
    Finite element analysis is done in the objective function and it's should always
    be calculated before constraints.

    """

    def __init__(self, name="Discrete Frame", profiles=RHS_PROFILES):
        super().__init__(name)

        self.profiles = profiles

    def add_structure(self, structure):
        self.structure = structure
        self.create_variables()
        self.create_objective()

    def create_variables(self):
        self.vars = []
        for i in range(len(self.structure.members)):
            X = IntegerVariable('X' + str(i + 1), 0, len(self.profiles) - 1)
            self.vars.append(X)

    def create_objective(self):

        def objective(X):
            if isinstance(X[0], IntegerVariable):
                X = [x.value for x in X]

            weight = 0
            for x, mem in zip(X, self.structure.members.values()):
                x = int(x)
                mem.profile = self.profiles[x]
                weight += mem.weight
            # Calculate stucture with new values
            self.structure.calculate()
            return weight

        self.obj = objective

    def _create_normal_force_constraints(self):
        """
        Creates normal force constraints
        """
        for i in range(len(self.structure.members)):
            def con_fun(X, i=i):
                mem = self.structure.members[i]
                return mem.ned / mem.NRd - 1

            constr = NonLinearConstraint(con_fun,
                                         name='Normal Force ' + str(i + 1))
            self.cons.append(constr)

    def _create_shear_force_constraints(self):
        """
        Creates shear force constraints
        """
        for i in range(len(self.structure.members)):
            def con_fun(X, i=i):
                mem = self.structure.members[i]
                return mem.ved / mem.VRd - 1

            constr = NonLinearConstraint(con_fun,
                                         name='Shear Force ' + str(i + 1))
            self.cons.append(constr)

    def _create_bending_moment_constraints(self):
        """
        Creates bending moment constraints about y-axis
        """
        for i in range(len(self.structure.members)):
            def con_fun(X, i=i):
                mem = self.structure.members[i]
                return mem.med / mem.MRd[0] - 1

            constr = NonLinearConstraint(con_fun,
                                         name='Bending Moment ' + str(i + 1))
            self.cons.append(constr)

    def create_constraints(self, N=True, V=True, M=True):

        # Initialize constraints as an empty list
        self.cons = []
        # Create normal force constraints
        if N:
            self._create_normal_force_constraints()
        if V:
            self._create_shear_force_constraints()
        if M:
            self._create_bending_moment_constraints()


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
                if self.prob_type == 'discrete':
                    A = THREE_BAR_AREAS_mm2[x]
                else:
                    A = x
                weight += self.rho * A * mem.length
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
                A = TEN_BAR_AREAS_mm2[x[j]]
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

                        if self.prob_type == 'discrete':
                            A = TEN_BAR_AREAS_mm2[x[j]]
                        else:
                            A = x[j]

                        return mem.ned / (A * mem.fy) - 1

                    def buckling_fun(x, i=i, j=j):

                        if self.prob_type == 'discrete':
                            A = TEN_BAR_AREAS_mm2[x[j]]
                        else:
                            A = x[j]

                        sigma_cr = 100 * mem.E * A / (8 * mem.length ** 2)
                        sigma = -mem.ned / A
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












        # # Create stress constraints
        # for i in range(len(self.structure.members)):
        #     def con_fun(x, i=i):
        #         mem = self.structure.members[i]
        #         mem.A = x[i]
        #         return abs(mem.ned / mem.NRd) - 1
        #
        #     constr = NonLinearConstraint(con_fun, name='Stress ' + str(i + 1),
        #                                  parent=self)
        #     constr.fea_required = True
        #     self.cons.append(constr)
        #
        # # Create displacement constraints
        # for j in range(len(self.structure.members)):
        #     def con_fun(A, j=j):
        #         mem = self.structure.members[j]
        #         displacements = mem.nodal_displacements.values()
        #         max_vals = [max(l[0:2]) for l in displacements]
        #         min_vals = [min(l[0:2]) for l in displacements]
        #         max_val = max(max_vals)
        #         min_val = min(min_vals)
        #         abs_max = max(max_val, abs(min_val))
        #         return abs_max / self.delta_max - 1
        #
        #     constr = NonLinearConstraint(con_fun,
        #                                  name='Displacement ' + str(1 + j),
        #                                  parent=self)
        #     constr.fea_required = True
        #     self.cons.append(constr)
        #
        # # Create buckling constraints
        # for j in range(len(self.structure.members)):
        #     def con_fun(x, j=j):
        #         mem = self.structure.members[j]
        #         sigma_cr = 100 * mem.E * x[j] / (8 * mem.length ** 2)
        #         sigma = abs(mem.ned) / mem.A
        #         return sigma / sigma_cr - 1
        #
        #     constr = NonLinearConstraint(con_fun,
        #                                  name='Buckling ' + str(1 + j),
        #                                  parent=self)
        #     constr.fea_required = True
        #     self.cons.append(constr)


class FifteenBarTruss(OptimizationProblem):
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

            constr = NonLinearConstraint(con_fun, name='Stress ' + str(i + 1))
            self.cons.append(constr)

        # Create buckling constraints
        for j in range(len(self.structure.members)):
            def con_fun(A, j=j):
                mem = self.structure.members[j]
                sigma_cr = 100 * mem.E * mem.A / (8 * mem.length ** 2)
                sigma = abs(mem.ned) / mem.A
                return sigma / sigma_cr - 1

            constr = NonLinearConstraint(con_fun,
                                         name='Buckling ' + str(11 + j))
            self.cons.append(constr)


class TenBarTruss(OptimizationProblem):
    # Problem parameters
    L = 9144  # mm
    F = 444890  # N
    E = 68950  # MPa
    rho = 2768e-9  # kg/mm3

    properties = {
        'L': str(L) + " mm",
        'F': f'{F / 1e3 :.2f} kN',
        'E': f'{E / 1e3 :.2f} MPa',
        'rho': f'{rho * 1e9 :.0f} kg/m3'
    }

    # Constraint limits
    sigma_max = 172.37  # MPa
    delta_max = 50.8  # mm

    def __init__(self, prob_type="discrete"):

        self.prob_type = prob_type
        super().__init__(name='TenBarTruss')
        self.structure = self.create_structure()
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

        if self.prob_type == "discrete":
            # Variables
            A1 = DiscreteVariable('A1', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar1]})
            A2 = DiscreteVariable('A2', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar2]})
            A3 = DiscreteVariable('A3', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar3]})
            A4 = DiscreteVariable('A4', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar4]})
            A5 = DiscreteVariable('A5', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar5]})
            A6 = DiscreteVariable('A6', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar6]})
            A7 = DiscreteVariable('A7', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar7]})
            A8 = DiscreteVariable('A8', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar8]})
            A9 = DiscreteVariable('A9', profiles=TEN_BAR_AREAS_mm2,
                                  target={"property": "A", "objects": [bar9]})
            A10 = DiscreteVariable('A10', profiles=TEN_BAR_AREAS_mm2,
                                   target={"property": "A",
                                           "objects": [bar10]})

            vars = [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10]


        elif self.prob_type == "continuous":
            A1 = Variable('A1', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar1]})
            A2 = Variable('A2', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar2]})
            A3 = Variable('A3', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar3]})
            A4 = Variable('A4', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar4]})
            A5 = Variable('A5', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar5]})
            A6 = Variable('A6', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar6]})
            A7 = Variable('A7', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar7]})
            A8 = Variable('A8', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar8]})
            A9 = Variable('A9', lb=TEN_BAR_AREAS_mm2[0],
                          ub=TEN_BAR_AREAS_mm2[-1],
                          target={"property": "A", "objects": [bar9]})
            A10 = Variable('A10', lb=TEN_BAR_AREAS_mm2[0],
                           ub=TEN_BAR_AREAS_mm2[-1],
                           target={"property": "A",
                                   "objects": [bar10]})

            vars = [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10]

        else:
            raise TypeError("Problem type must be either 'dicrete' "
                            "or 'continuous")

        self.vars = vars

        return frame

    def create_constraints(self):
        """
        Creates constraints for the problem
        """
        # Initialize constraints as an empty list
        self.cons = []
        # Create stress constraints
        for i in range(len(self.structure.members)):
            def con_fun(A, i=i):
                mem = self.structure.members[i]
                return mem.ned / mem.NRd - 1

            constr = NonLinearConstraint(con_fun, name='Stress ' + str(i + 1),
                                         parent=self)
            constr.fea_required = True
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

            constr = NonLinearConstraint(con_fun,
                                         name='Displacement ' + str(11 + j),
                                         parent=self)
            constr.fea_required = True
            self.cons.append(constr)

    def create_objective(self):

        def objective(X):
            weight = 0
            for x, mem in zip(X, self.structure.members.values()):
                weight += self.rho * mem.A * mem.length
            return weight

        self.obj = objective


def draw_structure(structure):
    from pyautocad import Autocad, APoint
    acad = Autocad(create_if_not_exists=True)
    scale = min(structure.L, structure.H)
    text_size = scale / 40
    node_size = scale / 100

    for mem in structure.members.values():
        [x1, y1], [x2, y2] = mem.coordinates
        p1 = APoint(x1, y1)
        p2 = APoint(x2, y2)
        acad.model.AddLine(p1, p2)
        # Nodes
        acad.model.AddCircle(p1, node_size)
        acad.model.AddCircle(p2, node_size)

        p1.x += mem.length / 4
        p1.y = mem.shape(p1.x)
        acad.model.AddText(f'{mem.mem_id + 1}', p1, text_size)

    for pl in structure.point_loads.values():
        x1, y1 = pl.coordinate
        fx, fy, m = np.asarray(pl.v)

        p1 = APoint(x1, y1)
        p2 = APoint(x1 + fx / 5e2, y1 + fy / 5e2)
        acad.model.AddLine(p1, p2)
        p2.x += scale / 100
        p2.y += scale / 100
        acad.model.AddText(f'F: {abs(fy) / 1000 :.2f} kN', p2, text_size)


def draw_deflected(structure):
    from pyautocad import Autocad, APoint
    acad = Autocad(create_if_not_exists=True)
    for mem in structure.members.values():
        [x1, y1], [x2, y2] = mem.coordinates
        p1 = APoint(x1, y1)
        p2 = APoint(x2, y2)
        delta1 = list(mem.nodal_displacements.values())[0]
        delta2 = list(mem.nodal_displacements.values())[-1]
        dx1, dy1, dr1 = delta1
        dx2, dy2, dr2 = delta2

        p1.x += dx1
        p1.y += dy1

        p2.x += dx2
        p2.y += dy2

        acad.model.AddLine(p1, p2)
        acad.model.AddCircle(p1, 100)
        acad.model.AddCircle(p2, 100)

    for pl in structure.point_loads.values():
        x1, y1 = pl.coordinate
        fx, fy, m = np.asarray(pl.v) / 5e2
        p1 = APoint(x1, y1)
        p2 = APoint(x1 + fx, y1 + fy)
        acad.model.AddLine(p1, p2)


def plot_values(problem, loc=[22, 15], X_vals=[]):
    OBJ_VALS = np.zeros((len(THREE_BAR_AREAS_mm2), len(THREE_BAR_AREAS_mm2)))
    CON_VALS = np.zeros((len(THREE_BAR_AREAS_mm2), len(THREE_BAR_AREAS_mm2)))
    for i in range(len(THREE_BAR_AREAS_mm2)):
        for j in range(len(THREE_BAR_AREAS_mm2)):
            X = [j, i]
            obj_val = problem.obj(X)
            con_val = max(problem.eval_con(X))
            if con_val > 0:
                con_val = -1
            else:
                con_val = 1

            OBJ_VALS[i, j] = obj_val
            CON_VALS[i, j] = con_val

    OBJ_CON = CON_VALS * OBJ_VALS

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib import cm

    N = 16
    p = 0.2
    top = cm.get_cmap('Reds', N)
    bottom = cm.get_cmap('Greens', N)

    newcolors = np.vstack((top(np.linspace(p, 1, N)),
                           bottom(np.linspace(p, 1, N))))
    newcmp = ListedColormap(newcolors, name='OrangeBlue')

    vmax = np.max(OBJ_VALS)
    vmin = -vmax
    plt.matshow(OBJ_VALS * CON_VALS, cmap=newcmp, vmin=vmin, vmax=vmax,
                interpolation='nearest', aspect='equal')
    if len(X_vals):
        for i, loc in enumerate(X_vals):
            if i > 0:
                plt.plot([prev_loc[0], loc[0]],
                         [prev_loc[1], loc[1]], 'k')
            prev_loc = loc
    else:
        plt.text(loc[0] - 0.5, loc[1] + 0.5, "X")
    plt.colorbar()
    plt.xticks(np.arange(0, len(THREE_BAR_AREAS_mm2)))
    plt.yticks(np.arange(0, len(THREE_BAR_AREAS_mm2)))
    plt.gca().set_xticks([x - 0.5 for x in plt.gca().get_xticks()][1:],
                         minor='true')
    plt.gca().set_yticks([y - 0.5 for y in plt.gca().get_yticks()][1:],
                         minor='true')
    plt.grid(which='minor')
    plt.show()

    # X = np.arange(0, len(THREE_BAR_AREAS_mm2))
    # Y = np.arange(0, len(THREE_BAR_AREAS_mm2))
    # Z = CON_VALS
    # CON = plt.contour(X, Y, Z, levels=[1-1e-6])
    # plt.clabel(CON, inline=True)
    # CS = plt.contour(X, Y, OBJ_VALS, levels=20, cmap='Greens_r')
    # plt.clabel(CS, inline=True)
    # plt.grid(True)


def plot_one_hot(areas=TEN_BAR_AREAS_mm2, point=ten_bar_DLM):
    one_hot = np.eye(len(areas))[point]
    plt.matshow(one_hot, cmap="Greys_r")
    plt.xticks(np.arange(0, len(areas)))
    plt.yticks(np.arange(0, len(point)))
    plt.gca().set_xticks([x - 0.5 for x in plt.gca().get_xticks()][1:],
                         minor='true')
    plt.gca().set_yticks([y - 0.5 for y in plt.gca().get_yticks()][1:],
                         minor='true')
    plt.grid(False)
    plt.show()

def plot_convex_hull():
    from scipy.spatial import ConvexHull, convex_hull_plot_2d
    points = []
    for A, Iy in zip(RHS_A, RHS_Iy):
        points.append([A, Iy])
    points = np.asarray(points)
    hull = ConvexHull(points)
    import matplotlib.pyplot as plt

    plt.plot(points[:, 0], points[:, 1], 'o')
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

    plt.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'r--', lw=2)
    plt.plot(points[hull.vertices[0], 0], points[hull.vertices[0], 1], 'ro')
    plt.show()

if __name__ == '__main__':
    problem = TenBarTruss(prob_type='continuous')
    solver = SLSQP()
    solver.solve(problem, maxiter=100)


    # solver = GA()
    # solver.solve(problem, maxiter=100 )
    # solver = GA(popsize=15)
    # solver.solve(problem, maxiter=50)
    # X, fopt, X_vals = solver.solve(problem,maxiter=5000, maxtime=3, subset_size=-1)
    # plot_values(problem, X_vals=X_vals)
    # problem([20, 13])
    # plot_values(problem)
    # solver.solve(problem, maxiter=5000)

    #plot_convex_hull()
    #plot_one_hot()
    # draw_structure(problem.structure)
    # draw_deflected(problem.structure)

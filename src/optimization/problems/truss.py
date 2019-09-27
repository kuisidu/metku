"""
IDEAS:
-------

General layout
--------------
- problem = TrussProblem(name='16-Bar-24m-KTruss',
                prob_type='design',
                var_type='continuous',
                variables=['B', 'T'],
                profiles=RHS_PROFILES,
                groups=[[mem1, mem2], [mem3, mem4, mem5]],
                objective='weight',
                structure=truss)

    Parameters:
    ------------
        - name: name of the problem
            default: "TrussProblem"

        - prob_type: 'design', 'geometry', 'topology'
            default: 'design'

        - var_type: 'continuous', 'discrete', 'binary'
            default: 'continuous'

        - variables: list of design variables
            default: ['A']

        - profiles: list of available profiles
            default: RHS_PROFILES

        - groups: list of groups, one variable is assigned to each group
            default: None

        - objective: 'cost', 'weight' or user given function
            default: 'weight'

        - structure: Truss2D/ Fram2D object
            default: None

Constraints
-----------
- design standard
- stability
- deflection
- stress

Other
-----
- structure knows whether it's symmetrical or not
    - it will automatically create grouping for symmetrical members


# Groups could be list of dicts

groups = [
    {"objects": [mem1, mem2, mem3],
    "lb": 0,
    "ub": 100,
    "property": "A",
    "profiles": RHS_PROFILES[10:30],
    },
    {"objects": [mem4, mem5],
    "lb": 10,
    "ub": 300,
    "property": "IZ",
    "profiles": RHS_PROFILES[:10],
    },
    {"objects": [mem6, mem7, mem8],
    "lb": 500,
    "ub": 10000,
    "property": "A",
    "profiles": RHS_PROFILES[100:130],
    }
]
"""

from src.optimization.structopt import *

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

RHS_IY = [69402.1348577099, 82151.28824369762, 141469.80388201258,
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

RHS_IZ = [69402.1348577099, 82151.28824369762, 141469.80388201258,
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


class TrussProblem(OptimizationProblem):

    def __init__(self,
                 name='TrussProblem',
                 prob_type='design',
                 var_type='continuous',
                 variables=['A'],
                 profiles=RHS_PROFILES,
                 groups=None,
                 objective='weight',
                 structure=None,
                 lb=0,
                 ub=1e5):

        super().__init__(name, structure=structure)

        self.prob_type = prob_type
        self.var_type = var_type
        self.variables = variables
        self.profiles = profiles
        self.groups = groups
        self.objective = objective
        self.lb = lb
        self.ub = ub

        # Assign objective
        if isinstance(self.objective, str):
            if self.objective == 'weight':
                self.obj = self.weight_fun
            elif self.objective == 'cost':
                pass
                # TODO self.obj = self.cost_fun
            else:
                raise ValueError("Objective must be either "
                                 "'weight' or 'cost' or a function")
        elif callable(self.objective):
            self.obj = self.objective
        else:
            raise ValueError("Objective must be either "
                             "'weight' or 'cost' or a function")

        # Create variables
        self.create_variables()

        # Create constraints
        self.create_constraints()

    @property
    def weight_fun(self):
        if self.structure is not None:
            def obj(x):
                self.substitute_variables(x)
                return self.structure.weight

            return obj
        else:
            raise ValueError("Structure is not defined!")

    def create_variables(self):
        """
        Creates design variables
        """
        self.vars = []

        if self.groups is None:
            self.groups = []
            for mem in self.structure.members.values():
                group = {
                    "objects": [mem],
                    "profiles": self.profiles,
                    "property": "A",
                    "lb": self.lb,
                    "ub": self.ub
                }
                self.groups.append(group)

        if self.var_type == 'discrete':
            for i, group in enumerate(self.groups):
                var = DiscreteVariable(
                    name=f"Var {i + 1}",
                    profiles=group['profiles'],
                    target={"property": group["property"],
                            "objects": group["objects"]}
                )
                self.vars.append(var)
        elif self.var_type == 'continuous':

            for i, group in enumerate(self.groups):
                var = Variable(
                    name=f"Var {i + 1}",
                    lb=group['lb'],
                    ub=group['ub'],
                    target={"property": group["property"],
                            "objects": group["objects"]}
                )
                self.vars.append(var)

        elif self.var_type == 'binary':
            pass
        else:
            raise ValueError("var_type must be either 'discrete',"
                             " 'continuous' or 'binary")

    def cross_section_constraints(self, sect, forces):
        """
        Creates cross-section constraints
        :return:
        """

        N, V, M = forces

        def compression(x):
            return -N / sect.NRd - 1

        def tension(x):
            return N / sect.NRd - 1

        def shear(x):
            return abs(V) / sect.VRd - 1

        def bending_moment(x):
            # Moment about y
            # TODO: Moment about z
            return abs(M) / sect.MRd[0] - 1

        return compression, tension, shear, bending_moment

    def stability_constraints(self, mem):
        """
        Creates stability constraint functions

        :return:
        """

        def buckling_y(x):
            self.substitute_variables(x)
            return -mem.ned / mem.NbRd[0] - 1

        def buckling_z(x):
            self.substitute_variables(x)
            return -mem.ned / mem.NbRd[1] - 1

        # def com_compression_bending_y(x):
        #     self.substitute_variables(x)
        #     return mem.steel_member.check_beamcolumn()[0] - 1
        #
        # def com_compression_bending_z(x):
        #     self.substitute_variables(x)
        #     return mem.steel_member.check_beamcolumn()[1] - 1

        return buckling_y, \
               buckling_z, \
            # com_compression_bending_y,\
        # com_compression_bending_z

    def deflection_constraints(self, mem):
        """
        Creates deflection constraint functions
        :param mem:
        :return:
        """
        self.delta_max = self.structure.L / 300

        def disp_fun(x):
            displacements = mem.nodal_displacements.values()
            max_vals = [max(l[0:2]) for l in displacements]
            min_vals = [min(l[0:2]) for l in displacements]
            max_val = max(max_vals)
            min_val = min(min_vals)
            abs_max = max(max_val, abs(min_val))
            return abs_max / self.delta_max - 1

        return disp_fun

    def create_constraints(self):
        """
        Cretes constraints and saves them to self.cons
        """
        # Initialize cons as an empty list
        self.cons = []
        for mem in self.structure.members.values():

            # buckling_y, buckling_z, com_compression_bending_y, \
            # com_compression_bending_z = \
            #     self.stability_constraints(mem)

            buckling_y, buckling_z, = self.stability_constraints(mem)

            # BUCKLING Y
            buckling_y_con = NonLinearConstraint(con_fun=buckling_y,
                                                 name="Buckling_y " +
                                                      str(mem.mem_id),
                                                 parent=self)
            buckling_y_con.fea_required = True
            self.cons.append(buckling_y_con)

            # BUCKLING Z
            buckling_z_con = NonLinearConstraint(con_fun=buckling_z,
                                                 name="Buckling_z " +
                                                      str(mem.mem_id),
                                                 parent=self)
            buckling_z_con.fea_required = True
            self.cons.append(buckling_z_con)

            # # BENDING + COMPRESSION Y
            # com_compression_bending_con_y = NonLinearConstraint(
            #     con_fun=com_compression_bending_y,
            #     name="Com_compression_bending_y " + str(mem.mem_id),
            #     parent=self)
            # com_compression_bending_con_y.fea_required = True
            # self.cons.append(com_compression_bending_con_y)
            #
            # # BENDING + COMPRESSION Z
            # com_compression_bending_con_z = NonLinearConstraint(
            #     con_fun=com_compression_bending_z,
            #     name="Com_compression_bending_z " + str(mem.mem_id),
            #     parent=self)
            # com_compression_bending_con_z.fea_required = True
            # self.cons.append(com_compression_bending_con_z)

            for i, elem in enumerate(mem.elements.values()):
                forces = [elem.axial_force[0], elem.shear_force[0],
                          elem.bending_moment[0]]
                compression, tension, shear, bending_moment = \
                    self.cross_section_constraints(mem, forces)

                compression_con = NonLinearConstraint(con_fun=compression,
                                                      name="Compression " +
                                                           str(mem.mem_id) +
                                                           str(i), parent=self)
                compression_con.fea_required = True

                tension_con = NonLinearConstraint(con_fun=tension,
                                                  name="Tension " +
                                                       str(mem.mem_id) +
                                                       str(i), parent=self)
                tension_con.fea_required = True

                shear_con = NonLinearConstraint(con_fun=shear,
                                                name="Shear " + str(mem.mem_id)
                                                     + str(i), parent=self)
                shear_con.fea_required = True

                bending_moment_con = NonLinearConstraint(
                    con_fun=bending_moment, name="Bending_moment " +
                                                 str(mem.mem_id) +
                                                 str(i), parent=self)
                bending_moment_con.fea_required = True

                self.cons.extend([compression_con, tension_con, shear_con,
                                  bending_moment_con])

                # Last element's end node
                if i == len(mem.elements) - 1:
                    forces = [elem.axial_force[1], elem.shear_force[1],
                              elem.bending_moment[1]]
                    compression, tension, shear, bending_moment = \
                        self.cross_section_constraints(mem, forces)

                    compression_con = NonLinearConstraint(con_fun=compression,
                                                          name="Compression " +
                                                               str(mem.mem_id)
                                                               + str(i + 1),
                                                          parent=self)
                    compression_con.fea_required = True
                    tension_con = NonLinearConstraint(con_fun=tension,
                                                      name="Tension " +
                                                           str(mem.mem_id) +
                                                           str(i + 1),
                                                      parent=self)
                    tension_con.fea_required = True
                    shear_con = NonLinearConstraint(con_fun=shear,
                                                    name="Shear " +
                                                         str(mem.mem_id) +
                                                         str(i + 1),
                                                    parent=self)
                    shear_con.fea_required = True

                    bending_moment_con = NonLinearConstraint(
                        con_fun=bending_moment, name="Bending_moment " +
                                                     str(mem.mem_id) +
                                                     str(i + 1), parent=self)
                    bending_moment_con.fea_required = True

                    self.cons.extend([compression_con, tension_con, shear_con,
                                      bending_moment_con])


class PlaneTrussProblem(OptimizationProblem):

    def __init__(self, **kwargs):
        super().__init__(name='PlaneTruss')
        self.prob_type = kwargs['prob_type']

        self.structure = self.create_structure(**kwargs)
        self.create_variables()
        self.create_constraints()
        self.create_objective()

    def create_constraints(self):
        """
        Cretes constraints and saves them to self.cons
        """
        # Initialize cons as an empty list
        self.cons = []
        for mem in self.structure.members.values():

            # buckling_y, buckling_z, com_compression_bending_y, \
            # com_compression_bending_z = \
            #     self.stability_constraints(mem)

            buckling_y, buckling_z, = self.stability_constraints(mem)

            # BUCKLING Y
            buckling_y_con = NonLinearConstraint(con_fun=buckling_y,
                                                 name="Buckling_y " +
                                                      str(mem.mem_id),
                                                 parent=self)
            buckling_y_con.fea_required = True
            self.cons.append(buckling_y_con)

            # BUCKLING Z
            buckling_z_con = NonLinearConstraint(con_fun=buckling_z,
                                                 name="Buckling_z " +
                                                      str(mem.mem_id),
                                                 parent=self)
            buckling_z_con.fea_required = True
            self.cons.append(buckling_z_con)

            # # BENDING + COMPRESSION Y
            # com_compression_bending_con_y = NonLinearConstraint(
            #     con_fun=com_compression_bending_y,
            #     name="Com_compression_bending_y " + str(mem.mem_id),
            #     parent=self)
            # com_compression_bending_con_y.fea_required = True
            # self.cons.append(com_compression_bending_con_y)
            #
            # # BENDING + COMPRESSION Z
            # com_compression_bending_con_z = NonLinearConstraint(
            #     con_fun=com_compression_bending_z,
            #     name="Com_compression_bending_z " + str(mem.mem_id),
            #     parent=self)
            # com_compression_bending_con_z.fea_required = True
            # self.cons.append(com_compression_bending_con_z)

            for i, elem in enumerate(mem.elements.values()):
                forces = [elem.axial_force[0], elem.shear_force[0],
                          elem.bending_moment[0]]
                compression, tension, shear, bending_moment = \
                    self.cross_section_constraints(mem, forces)

                compression_con = NonLinearConstraint(con_fun=compression,
                                                      name="Compression " +
                                                           str(mem.mem_id) +
                                                           str(i), parent=self)
                compression_con.fea_required = True

                tension_con = NonLinearConstraint(con_fun=tension,
                                                  name="Tension " +
                                                       str(mem.mem_id) +
                                                       str(i), parent=self)
                tension_con.fea_required = True

                shear_con = NonLinearConstraint(con_fun=shear,
                                                name="Shear " + str(mem.mem_id)
                                                     + str(i), parent=self)
                shear_con.fea_required = True

                bending_moment_con = NonLinearConstraint(
                    con_fun=bending_moment, name="Bending_moment " +
                                                 str(mem.mem_id) +
                                                 str(i), parent=self)
                bending_moment_con.fea_required = True

                self.cons.extend([compression_con, tension_con, shear_con,
                                  bending_moment_con])

                # Last element's end node
                if i == len(mem.elements) - 1:
                    forces = [elem.axial_force[1], elem.shear_force[1],
                              elem.bending_moment[1]]
                    compression, tension, shear, bending_moment = \
                        self.cross_section_constraints(mem, forces)

                    compression_con = NonLinearConstraint(con_fun=compression,
                                                          name="Compression " +
                                                               str(mem.mem_id)
                                                               + str(i + 1),
                                                          parent=self)
                    compression_con.fea_required = True
                    tension_con = NonLinearConstraint(con_fun=tension,
                                                      name="Tension " +
                                                           str(mem.mem_id) +
                                                           str(i + 1),
                                                      parent=self)
                    tension_con.fea_required = True
                    shear_con = NonLinearConstraint(con_fun=shear,
                                                    name="Shear " +
                                                         str(mem.mem_id) +
                                                         str(i + 1),
                                                    parent=self)
                    shear_con.fea_required = True

                    bending_moment_con = NonLinearConstraint(
                        con_fun=bending_moment, name="Bending_moment " +
                                                     str(mem.mem_id) +
                                                     str(i + 1), parent=self)
                    bending_moment_con.fea_required = True

                    self.cons.extend([compression_con, tension_con, shear_con,
                                      bending_moment_con])

    def create_objective(self):
        """
        Creates the objective function
        """

        def objective(x):
            self.substitute_variables(x)
            weight = 0
            for mem in self.structure.members.values():
                weight += mem.weight

            return weight

        self.obj = objective

    def create_structure(self, **kwargs):
        """
        Creates truss
        :return:
        """
        truss = PlaneTruss(**kwargs)

        # Loads
        for tc in truss.top_chords:
            truss.add(LineLoad(tc, kwargs['q'], kwargs['dir']))

        # Supports
        truss.add(XYHingedSupport([0, truss.H1]))
        truss.add(YHingedSupport([truss.L1, truss.H1]))

        truss.generate()
        truss.calculate()

        return truss

    def create_variables(self, profiles=RHS_PROFILES):
        """
        Creates optimization variables
        """
        self.vars = []
        chords = self.structure.bottom_chords.copy()
        chords.extend(self.structure.top_chords)
        TC = [mem for mem in self.structure.top_chords]
        BC = [mem for mem in self.structure.bottom_chords]
        TW = [mem for mem in self.structure.members.values() if
              mem not in chords and mem.ned >= 0]
        CW = [mem for mem in self.structure.members.values() if
              mem not in chords and mem.ned < 0]

        groups = [TC, BC, TW, CW]
        names = ['TopChords', 'BottomChords', 'TensionWebs', 'CompressionWebs']
        for i, group in enumerate(groups):
            name = names[i]
            if self.prob_type == "discrete":
                var = DiscreteVariable(name,
                                       profiles=profiles,
                                       # TODO: THIS IS A PARAMETER
                                       values=RHS_A,
                                       target={"property": "A",
                                               "objects": group})

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

        #
        # WEBS = [mem for mem in self.structure.members.values() if mem.mem_type == "W"]
        # nodes = [w.n1 for w in WEBS]
        # nodes = sorted(nodes, key=lambda n: n.x)
        #
        # for i, node in enumerate(nodes):
        #     name = "Node " + str(node.nid)
        #     if i == 0:
        #         var = Variable(name,
        #                        lb=0,
        #                        ub=nodes[i+1].x,
        #                        target={"property": "x",
        #                                "objects": [node]})
        #     elif i == len(nodes) -1:
        #         var = Variable(name,
        #                        lb=nodes[i-1].x,
        #                        ub=self.structure.L1,
        #                        target={"property": "x",
        #                                "objects": [node]})
        #     else:
        #         var = Variable(name,
        #                        lb=nodes[i - 1].x,
        #                        ub=nodes[i+1].x,
        #                        target={"property": "x",
        #                                "objects": [node]})
        #     self.vars.append(var)

    def cross_section_constraints(self, sect, forces):
        """
        Creates cross-section constraints
        :return:
        """

        N, V, M = forces

        def compression(x):
            return -N / sect.NRd - 1

        def tension(x):
            return N / sect.NRd - 1

        def shear(x):
            return abs(V) / sect.VRd - 1

        def bending_moment(x):
            # Moment about y
            # TODO: Moment about z
            return abs(M) / sect.MRd[0] - 1

        return compression, tension, shear, bending_moment

    def stability_constraints(self, mem):
        """
        Creates stability constraint functions

        :return:
        """

        def buckling_y(x):
            self.substitute_variables(x)
            return -mem.ned / mem.NbRd[0] - 1

        def buckling_z(x):
            self.substitute_variables(x)
            return -mem.ned / mem.NbRd[1] - 1

        # def com_compression_bending_y(x):
        #     self.substitute_variables(x)
        #     return mem.steel_member.check_beamcolumn()[0] - 1
        #
        # def com_compression_bending_z(x):
        #     self.substitute_variables(x)
        #     return mem.steel_member.check_beamcolumn()[1] - 1

        return buckling_y, \
               buckling_z, \
            # com_compression_bending_y,\
        # com_compression_bending_z

    def deflection_constraints(self, mem):
        """
        Creates deflection constraint functions
        :param mem:
        :return:
        """
        self.delta_max = self.structure.L / 300

        def disp_fun(x):
            displacements = mem.nodal_displacements.values()
            max_vals = [max(l[0:2]) for l in displacements]
            min_vals = [min(l[0:2]) for l in displacements]
            max_val = max(max_vals)
            min_val = min(min_vals)
            abs_max = max(max_val, abs(min_val))
            return abs_max / self.delta_max - 1

        return disp_fun


if __name__ == '__main__':
    from src.optimization.solvers import *
    from src.frame2d.frame2d import *
    import matplotlib.pyplot as plt

    frame = Frame2D(simple=[1,1, 5e3, 5e3], supports='fixed')
    for mem in frame.members.values():
        if mem.mtype == "beam":
            frame.add(LineLoad(mem, [-20, -20], 'y'))

    frame.generate()
    frame.calculate()

    problem = TrussProblem(
        name="FrameTest",
        structure=frame,
        var_type='discrete',
        profiles=RHS_A
    )

    solver = DiscreteVNS()
    x0 = [var.ub for var in problem.vars]
    problem(x0)
    solver.solve(problem, x0=x0, maxtime=30)


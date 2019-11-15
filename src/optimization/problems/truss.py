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

try:
    from src.optimization.structopt import *
except:
    from optimization.structopt import *

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

RHS_COMP = ['RHS 40X40X3.0', 'RHS 50X50X3.0', 'RHS 50X50X4.0', 'RHS 60X60X3.0',
            'RHS 60X60X4.0', 'RHS 60X60X5.0', 'RHS 70X70X3.0', 'RHS 70X70X4.0',
            'RHS 70X70X5.0', 'RHS 80X80X3.0', 'RHS 80X80X4.0', 'RHS 80X80X5.0',
            'RHS 90X90X3.0', 'RHS 90X90X4.0', 'RHS 90X90X5.0',
            'RHS 100X100X3.0', 'RHS 100X100X4.0', 'RHS 100X100X5.0',
            'RHS 100X100X6.0', 'RHS 120X120X4.0', 'RHS 120X120X5.0',
            'RHS 120X120X6.0', 'RHS 120X120X7.1', 'RHS 120X120X8.0',
            'RHS 120X120X8.8', 'RHS 120X120X10.0', 'RHS 140X140X4.0',
            'RHS 140X140X5.0', 'RHS 140X140X6.0', 'RHS 150X150X5.0',
            'RHS 150X150X6.0', 'RHS 150X150X7.1', 'RHS 150X150X8.0',
            'RHS 150X150X8.8', 'RHS 150X150X10.0', 'RHS 160X160X5.0',
            'RHS 160X160X6.0', 'RHS 160X160X8.0', 'RHS 160X160X7.1',
            'RHS 160X160X8.8', 'RHS 160X160X10.0', 'RHS 180X180X5.0',
            'RHS 180X180X5.6', 'RHS 180X180X6.0', 'RHS 180X180X7.1',
            'RHS 180X180X8.0', 'RHS 180X180X8.8', 'RHS 180X180X10.0',
            'RHS 180X180X12.5', 'RHS 200X200X6.0', 'RHS 200X200X7.1',
            'RHS 200X200X8.0', 'RHS 200X200X8.8', 'RHS 200X200X10.0',
            'RHS 200X200X12.5', 'RHS 220X220X6.0', 'RHS 220X220X7.1',
            'RHS 220X220X8.0', 'RHS 220X220X8.8', 'RHS 220X220X10.0',
            'RHS 250X250X8.0', 'RHS 250X250X10.0', 'RHS 250X250X7.1',
            'RHS 250X250X8.8', 'RHS 250X250X12.5', 'RHS 260X260X7.1',
            'RHS 260X260X8.0', 'RHS 260X260X8.8', 'RHS 260X260X10.0',
            'RHS 260X260X12.5', 'RHS 300X300X8.0', 'RHS 300X300X8.8',
            'RHS 300X300X10.0', 'RHS 300X300X12.5', 'RHS 400X400X12.5']

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
            def obj_fun(x):
                return self.structure.weight

            obj = ObjectiveFunction(name="Objective ",
                                    obj_fun=obj_fun)
            obj.problem = self

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
                    "name": "var name",
                    "objects": [mem],
                    "profiles": self.profiles,
                    "property": "A",
                    "properties": self.variables,
                    "var_type": 'continuous',
                    "lb": self.lb,
                    "ub": self.ub
                }
                self.groups.append(group)

        for i, group in enumerate(self.groups):

            if group["var_type"] == 'discrete':
                var = DiscreteVariable(
                    name=group['name'],
                    values=group['profiles'],
                    target={"property": group["property"],
                            "objects": group["objects"]}
                )
                self.add(var)

            elif group["var_type"] == 'continuous':
                # Multiple properties
                if "properties" in group.keys():
                    bounds = group["bounds"]
                    if len(bounds) != len(group["properties"]):
                        raise ValueError(
                            "There must be same number of bounds as properties!"
                            f"{group['bounds']} != {group['properties']}")
                    for bounds, prop, value in zip(group["bounds"],
                                                   group["properties"],
                                                   group["values"]):
                        lb, ub = bounds
                        var = Variable(
                            name=group['name'],
                            lb=lb,
                            ub=ub,
                            value=value,
                            target={"property": prop,
                                    "objects": group["objects"]}
                        )
                        self.add(var)

                # Single property
                else:
                    var = Variable(
                        name=group['name'],
                        lb=group['lb'],
                        ub=group['ub'],
                        value=group['value'],
                        target={"property": group["property"],
                                "objects": group["objects"]}
                    )
                    self.add(var)

            elif group["var_type"] == 'index':
                var = IndexVariable(
                    name=group['name'],
                    values=group['profiles'],
                    target={"property": group["property"],
                            "objects": group["objects"]}
                )
                self.add(var)

            elif group["var_type"] == 'binary':
                pass
            else:
                raise ValueError("var_type must be either 'discrete',"
                                 " 'continuous', 'index' or 'binary")


    def index_to_binary(self):
        """
        Creates binary variables and continuous variables
        :return:
        """
        idx_vars = [var for var in self.vars if isinstance(var, IndexVariable)]

        # Temp member instance to get needed values
        mem = SteelBeam([[0, 0], [1000, 0]])
        vars = {}
        for i, idx_var in enumerate(idx_vars):
            A_list = []
            Iy_list = []
            Iz_list = []
            Wply_list = []
            Wplz_list = []
            bin_vars = []
            for j, profile in enumerate(idx_var.values):
                # Get attributes
                mem.profile = profile
                A = mem.cross_section.A
                Iy, Iz = mem.cross_section.I
                Wply, Wplz = mem.cross_section.Wpl

                # Append attributes
                A_list.append(A)
                Iy_list.append(Iy)
                Iz_list.append(Iz)
                Wply_list.append(Wply)
                Wplz_list.append(Wplz)

                bin_var = BinaryVariable(name=f"BinVar{i}{j}")
                bin_vars.append(bin_var)

            # Create continuous variables
            # A
            A_var = Variable(name=f"Continuous A{i}",
                             lb=min(A_list),
                             ub=max(A_list),
                             target={
                                 "objects": idx_var.target["objects"],
                                 "property": "A"
                             })
            # Iy
            Iy_var = Variable(name=f"Continuous Iy{i}",
                              lb=min(Iy_list),
                              ub=max(Iy_list),
                              target={
                                  "objects": idx_var.target["objects"],
                                  "property": "Iy"
                              })
            # Iz
            Iz_var = Variable(name=f"Continuous Iz{i}",
                              lb=min(Iz_list),
                              ub=max(Iz_list),
                              target={
                                  "objects": idx_var.target["objects"],
                                  "property": "Iz"
                              })
            # Wply
            Wply_var = Variable(name=f"Continuous Wply{i}",
                                lb=min(Wply_list),
                                ub=max(Wply_list),
                                target={
                                    "objects": idx_var.target["objects"],
                                    "property": "Wply"
                                })
            # Wplz
            Wplz_var = Variable(name=f"Continuous Wplz{i}",
                                lb=min(Wplz_list),
                                ub=max(Wplz_list),
                                target={
                                    "objects": idx_var.target["objects"],
                                    "property": "Wplz"
                                })

            cont_vars = [A_var, Iy_var, Iz_var, Wply_var, Wplz_var]
            lists = [A_list, Iy_list, Iz_list, Wply_list, Wplz_list]

            vars[i] = {'bin_vars': bin_vars, 'cont_vars': cont_vars,
                       'lists': lists}
            self.vars.remove(idx_var)
            self.vars.extend(bin_vars)
            self.vars.extend(cont_vars)

        # Create linear constraints
        for i, vars_dict in enumerate(vars.values()):
            bin_vars = vars_dict['bin_vars']
            cont_vars = vars_dict['cont_vars']
            lists = vars_dict['lists']
            A_var, Iy_var, Iz_var, Wply_var, Wplz_var = cont_vars
            A_list, Iy_list, Iz_list, Wply_list, Wplz_list = lists
            # Binary constraint sum(bin_vars) == 1
            bin_idx = [self.vars.index(bvar) for bvar in bin_vars]
            a = np.zeros(len(self.vars))
            a[bin_idx] = 1
            bin_con = LinearConstraint(a=a, b=1, con_type="=",
                                       name="Binary Constraint "
                                            + str(i))
            # A con; sum(Â*bin_vars) == A
            a_A = np.zeros(len(self.vars))
            a_A[bin_idx] = A_list
            b = A_var
            A_con = LinearConstraint(a=a_A,
                                     b=b,
                                     name=f"A Constraint {i}")
            # Iy con; sum(Îy*bin_vars) == Iy
            a_Iy = np.zeros(len(self.vars))
            a_Iy[bin_idx] = Iy_list
            b = Iy_var
            Iy_con = LinearConstraint(a=a_Iy,
                                      b=b,
                                      name=f"Iy Constraint {i}")
            # Iz con; sum(Îz*bin_vars) == Iz
            a_Iz = np.zeros(len(self.vars))
            a_Iz[bin_idx] = Iz_list
            b = Iz_var
            Iz_con = LinearConstraint(a=a_Iz,
                                      b=b,
                                      name=f"Iz Constraint {i}")
            # Wply con; sum(Wply*bin_vars) == Wply
            a_Wply = np.zeros(len(self.vars))
            a_Wply[bin_idx] = Wply_list
            b = Wply_var
            Wply_con = LinearConstraint(a=a_Wply,
                                        b=b,
                                        name=f"Wply Constraint {i}")
            # Wplz con; sum(Wplz*bin_vars) == Wplz
            a_Wplz = np.zeros(len(self.vars))
            a_Wplz[bin_idx] = Wplz_list
            b = Wplz_var
            Wplz_con = LinearConstraint(a=a_Wplz,
                                        b=b,
                                        name=f"Wplz Constraint {i}")

            self.cons.extend([bin_con, A_con, Iy_con,
                              Iz_con, Wply_con, Wplz_con])

    def joint_geometry_constraints(self, joint):



        if joint.joint_type == 'Y':
            pass



        elif joint.joint_type == 'K':
            w1, w2 = joint.webs.values()

            def g_min(x):
                gap = w1.cross_section.T + w2.cross_section.T
                val = gap / joint.g1 - 1
                return val

            def g_min_beta(x):
                beta = joint.rhs_joint.beta()
                val = 0.5 * (1 - beta) * joint.chord.cross_section.B
                val -= joint.g1
                return val

            def g_max_beta(x):
                beta = joint.rhs_joint.beta()
                val = - 1.5 * (1 - beta) * joint.chord.cross_section.B
                val += joint.g1
                return val

            def b1_min(x):
                val = 0.35 * joint.chord.cross_section.B
                val -= w1.cross_section.B
                return val

            def b1_max(x):
                val = -0.85 * joint.chord.cross_section.B
                val += w1.cross_section.B
                return val

            def b2_min(x):
                val = 0.35 * joint.chord.cross_section.B
                val -= w2.cross_section.B
                return val

            def b2_max(x):
                val = -0.85 * joint.chord.cross_section.B
                val += w2.cross_section.B
                return val

            def theta1(x):
                val = 30 - np.degrees(w1.angle - joint.chord.angle)
                return val

            def theta2(x):
                val = 30 - np.degrees(w2.angle - joint.chord.angle)
                return val

            def e_pos(x):
                return joint.e / ((0.25 * joint.chord.cross_section.H)) - 1

            def e_neg(x):
                return -0.55 * joint.chord.cross_section.H / joint.e - 1

            con_funs = {f'Joint {joint.jid} g_min': g_min,
                        f'Joint {joint.jid} g_min_beta': g_min_beta,
                        f'Joint {joint.jid} g_max_beta': g_max_beta,
                        f'Joint {joint.jid} b1_min': b1_min,
                        f'Joint {joint.jid} b1_max': b1_max,
                        f'Joint {joint.jid} b2_min': b2_min,
                        f'Joint {joint.jid} b2_max': b2_max,
                        f'Joint {joint.jid}theta1': theta1,
                        f'Joint {joint.jid}theta2': theta2,
                        f'Joint {joint.jid}e_pos': e_pos,
                        f'Joint {joint.jid}e_neg': e_neg}

            return con_funs

        elif joint.joint_type == 'KT':
            pass

        return {}

    def joint_strength_constraints(self, joint):
        """
        Creates joint strength costraint functions

        :param joint: joint to be calculated
        :return: constraint functions in a dict
        """

        if joint.joint_type == 'Y':
            pass

        elif joint.joint_type == 'K':
            w1, w2 = joint.webs.values()

            def chord_face_failure_1(x):
                NRd1, NRd2 = joint.rhs_joint.chord_face_failure()
                return w1.ned / NRd1 - 1

            def chord_face_failure_2(x):
                NRd1, NRd2 = joint.rhs_joint.chord_face_failure()
                return w2.ned / NRd2 - 1

            def punching_shear(x):
                NRd = joint.rhs_joint.punching_shear()
                return joint.N0 / NRd - 1

            def chord_shear(x):
                pass

            def brace_failure(x):
                pass

        elif joint.joint_type == 'KT':
            pass

        return {}





    def cross_section_constraints(self, sect, elem):
        """
        Creates cross-section constraints
        :return:
        """

        def compression(x):
            N = elem.axial_force[0]
            return -N / sect.NRd - 1

        def tension(x):
            N = elem.axial_force[0]
            return N / sect.NRd - 1

        def shear(x):
            V = elem.shear_force[0]
            return abs(V) / sect.VRd - 1

        def bending_moment(x):
            # Moment about y
            # TODO: Moment about z
            M = elem.bending_moment[0]
            return abs(M) / sect.MRd[0] - 1

        return compression, tension, shear, bending_moment

    def stability_constraints(self, mem):
        """
        Creates stability constraint functions

        :return:
        """

        def buckling_y(x):
            return -mem.NEd / mem.NbRd[0] - 1

        def buckling_z(x):
            return -mem.NEd / mem.NbRd[1] - 1

        def com_compression_bending_y(x):
            return mem.check_beamcolumn()[0] - 1

        def com_compression_bending_z(x):
            return mem.check_beamcolumn()[1] - 1

        return buckling_y, \
               buckling_z, \
               com_compression_bending_y, \
               com_compression_bending_z

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
            for j, smem in enumerate(mem.steel_members):
                buckling_y, buckling_z, com_compression_bending_y, \
                com_compression_bending_z = \
                    self.stability_constraints(smem)

                # buckling_y, buckling_z, = self.stability_constraints(mem)

                # BUCKLING Y
                buckling_y_con = self.non_linear_constraint(con_fun=buckling_y,
                                                            name="Buckling_y " +
                                                                 str(
                                                                     mem.mem_id) +
                                                                 str(j))
                buckling_y_con.fea_required = True

                # BUCKLING Z
                buckling_z_con = self.non_linear_constraint(con_fun=buckling_z,
                                                            name="Buckling_z " +
                                                                 str(
                                                                     mem.mem_id) + str(
                                                                j),
                                                            )
                buckling_z_con.fea_required = True

            # # BENDING + COMPRESSION Y
            # com_compression_bending_con_y = self.non_linear_constraint()(
            #     con_fun=com_compression_bending_y,
            #     name="Com_compression_bending_y " + str(mem.mem_id),
            #     )
            # com_compression_bending_con_y.fea_required = True
            # self.cons.append(com_compression_bending_con_y)
            #
            # # BENDING + COMPRESSION Z
            # com_compression_bending_con_z = self.non_linear_constraint()(
            #     con_fun=com_compression_bending_z,
            #     name="Com_compression_bending_z " + str(mem.mem_id),
            #     )
            # com_compression_bending_con_z.fea_required = True
            # self.cons.append(com_compression_bending_con_z)

            # CROSS-SECTION STRENGTH
            for i, elem in enumerate(mem.elements.values()):
                compression, tension, shear, bending_moment = \
                    self.cross_section_constraints(mem, elem)

                compression_con = self.non_linear_constraint(
                    con_fun=compression,
                    name="Compression " +
                         str(mem.mem_id) +
                         str(i), )
                compression_con.fea_required = True

                tension_con = self.non_linear_constraint(con_fun=tension,
                                                         name="Tension " +
                                                              str(mem.mem_id) +
                                                              str(i), )
                tension_con.fea_required = True

                shear_con = self.non_linear_constraint(con_fun=shear,
                                                       name="Shear " + str(
                                                           mem.mem_id)
                                                            + str(i), )
                shear_con.fea_required = True

                bending_moment_con = self.non_linear_constraint(
                    con_fun=bending_moment, name="Bending_moment " +
                                                 str(mem.mem_id) +
                                                 str(i), )
                bending_moment_con.fea_required = True

                # Last element's end node
                if i == len(mem.elements) - 1 and mem.mtype != 'web':
                    compression, tension, shear, bending_moment = \
                        self.cross_section_constraints(mem, elem)

                    compression_con = self.non_linear_constraint(
                        con_fun=compression,
                        name="Compression " +
                             str(mem.mem_id)
                             + str(i + 1),
                    )
                    compression_con.fea_required = True
                    tension_con = self.non_linear_constraint(con_fun=tension,
                                                             name="Tension " +
                                                                  str(
                                                                      mem.mem_id) +
                                                                  str(i + 1),
                                                             )
                    tension_con.fea_required = True
                    shear_con = self.non_linear_constraint(con_fun=shear,
                                                           name="Shear " +
                                                                str(
                                                                    mem.mem_id) +
                                                                str(i + 1),
                                                           )
                    shear_con.fea_required = True

                    bending_moment_con = self.non_linear_constraint(
                        con_fun=bending_moment, name="Bending_moment " +
                                                     str(mem.mem_id) +
                                                     str(i + 1), )
                    bending_moment_con.fea_required = True

        # JOINT CONS
        for joint in self.structure.joints.values():

            con_funs = self.joint_geometry_constraints(joint)
            for name, con_fun in con_funs.items():
                self.non_linear_constraint(con_fun=con_fun,
                                           name=name)


class Joint:

    def __init__(self, loc, chord):
        self.chord = chord
        self.node = None
        self._loc = loc

    def __repr__(self):
        return type(self).__name__ + str(self.loc)

    @property
    def coord(self):
        return self.chord.to_global(self.loc)

    @property
    def loc(self):
        return self._loc

    @loc.setter
    def loc(self, val):
        val = min(max(0, val), 1)
        self._loc = val
        x, y = self.coord
        self.node.x = x
        self.node.y = y
        self.chord.calc_nodal_coordinates()

    def generate(self, fem_model):
        coord = [round(c, 3) for c in self.coord]
        idx = fem_model.nodal_coords.index(coord)
        self.node = fem_model.nodes[idx]

    def plot(self, color=None):
        pass


def create_planetruss(L=25000, H1=1500, H2=2000, n=16, dx=0):
    truss = Frame2D(num_elements=2)

    # Top Chords
    helper = SteelBeam([[0, H1], [L / 2, H2]])
    helper.mtype = "top_chord"
    helper2 = SteelBeam([[L, H1], [L / 2, H2]])
    helper2.mtype = "top_chord"
    # Bottom Chords
    bc_helper = SteelBeam([[dx, 0], [L / 2, 0]])
    bc_helper.mtype = "bottom_chord"
    bc_helper2 = SteelBeam([[L - dx, 0], [L / 2, 0]])
    bc_helper2.mtype = "bottom_chord"

    sep = 1 / (n / 4)
    tc_locs = np.arange(0, 1 + sep, sep)
    tc_coords = []
    for loc in tc_locs:
        tc_coords.append(list(helper.to_global(loc)))

    # Left side
    for i in range(int(n / 4)):
        c0 = list(tc_coords[i])
        c1 = list(tc_coords[i + 1])
        tc = SteelBeam([c0, c1])
        truss.add(tc)
        tc.mtype = "top_chord"

    # Right side
    for i in range(int(n / 4)):
        c0 = list(tc_coords[i])
        c0[0] = L - c0[0]
        c1 = list(tc_coords[i + 1])
        c1[0] = L - c1[0]
        tc = SteelBeam([c0, c1])
        truss.add(tc)
        tc_coords.append(c0)
        tc.mtype = "top_chord"
    # Bottom Chords
    bc_split = (L - 2 * dx) / (n / 2 - 1)
    bc_coords = []
    sep = 1 / (int((n - 2) / 2))
    bc_locs = np.arange(0, 1 + sep, sep)
    for i in range(int((n - 2) / 2)):
        c0 = [i * bc_split + dx, 0]
        c1 = [(i + 1) * bc_split + dx, 0]
        bc = SteelBeam([c0, c1])
        truss.add(bc)
        bc.mtype = "BC"
        bc_coords.append(c0)
    bc_coords.append(c1)
    # Webs
    tc_coords.sort()
    bc_coords.sort()

    for i, bcoord in enumerate(bc_coords):
        tc0 = tc_coords[i]
        tc1 = tc_coords[i + 1]

        w0 = SteelBeam([bcoord, tc0], Sj1=0, Sj2=0)
        w1 = SteelBeam([bcoord, tc1], Sj1=0, Sj2=0)

        truss.add(w0)
        truss.add(w1)
        w0.mtype = "W"
        w1.mtype = "W"

    # Joints
    # TC joints
    for loc in tc_locs:
        jid = len(truss.joints)
        truss.joints[jid] = Joint(loc, helper)
        jid = len(truss.joints)
        truss.joints[jid] = Joint(loc, helper2)

    # BC joints
    for loc in bc_locs:
        if loc < 0.5:
            jid = len(truss.joints)
            truss.joints[jid] = Joint(loc, bc_helper)
        else:
            jid = len(truss.joints)
            loc = 1 - loc
            truss.joints[jid] = Joint(loc, bc_helper2)

    # Supports
    truss.add(XYHingedSupport([0, H1]))
    truss.add(YHingedSupport([L, H1]))

    return truss


if __name__ == '__main__':
    try: 
        from src.optimization.solvers import *
        from src.frame2d.frame2d import *
        from src.truss2d import *
    except:
        from optimization.solvers import *
        from frame2d.frame2d import *
        from truss2d import *

    import matplotlib.pyplot as plt

    n = 20

    truss = Truss2D(simple=dict(
        H0=0,
        H1=1000,
        H2=2000,
        L1=12500,
        dx=500,
        n=n
    ))
    truss.add(XYHingedSupport([0, truss.H0 + truss.H1]))
    truss.add(YHingedSupport([truss.L1 + truss.L2, truss.H0 + truss.H3]))

    for mem in truss.members.values():
        mem.profile = "RHS 300X300X5"
        if mem.mtype == "top_chord":
            truss.add(LineLoad(mem, [-30, -30], 'y'))

    truss.generate()
    truss.calculate()

    top_chords = []
    bot_chords = []
    webs_tension = []
    webs_compression = []
    webs = []
    bot_nodes = []
    for mem in truss.members.values():
        mem.material = "S420"
        if mem.mtype == "top_chord":
            top_chords.append(mem)
        elif mem.mtype == "bottom_chord":
            bot_chords.append(mem)
            for node in mem.nodes.values():
                if node not in bot_nodes:
                    bot_nodes.append(node)
        else:
            webs.append(mem)
            if mem.ned > 0:
                webs_tension.append(mem)
            else:
                webs_compression.append(mem)

    TC_group = {
        "objects": top_chords,
        "name": "TC index",
        # "values": [500, 500, 5],
        # "bounds": [[50, 800], [50, 800], [1, 8]],
        # "properties": ["H", "B", "T"],
        # "var_type": 'continuous',
        "var_type": 'index',
        'property': 'profile',
        "profiles": RHS_COMP,
    }

    BC_group = {
        "objects": bot_chords,
        "name": "BC index",
        # "values": [500, 500, 5],
        # "bounds": [[50, 800], [50, 800], [1, 8]],
        # "properties": ["H", "B", "T"],
        # "var_type": 'continuous',
        "var_type": 'index',
        'property': 'profile',
        "profiles": RHS_COMP,
    }

    WT_group = {
        "objects": webs_tension,
        "name": "W Tension index",
        # "values": [500, 500, 5],
        # "bounds": [[50, 500], [50, 500], [1, 8]],
        # "properties": ["H", "B", "T"],
        "var_type": 'index',
        # "var_type": 'continuous',
        'property': 'profile',
        "profiles": RHS_COMP,
    }

    WC_group = {
        "objects": webs_compression,
        "name": "W Compression index",
        # "values": [500, 500, 5],
        # "bounds": [[50, 500], [50, 500], [1, 8]],
        # "properties": ["H", "B", "T"],
        "var_type": 'index',
        # "var_type": 'continuous',
        'property': 'profile',
        "profiles": RHS_COMP,
    }

    bnode_group = {
        "name": "BC Y-loc",
        "value": 0,
        "objects": bot_nodes,
        "lb": -500,
        "ub": 500,
        "var_type": 'continuous',
        'property': 'y',
    }

    groups = [TC_group, BC_group]  # , WT_group, WC_group]

    for i in range(int(len(webs) / 2)):
        w = webs[i * 2: i * 2 + 2]
        W_group = {
            "objects": [web for web in w],
            "name": "W " + str(i),
            # "values": [500, 500, 5],
            # "bounds": [[50, 500], [50, 500], [1, 8]],
            # "properties": ["H", "B", "T"],
            "var_type": 'index',
            # "var_type": 'continuous',
            'property': 'profile',
            "profiles": RHS_COMP,
        }
        groups.append(W_group)

    tc_joints = [j for j in truss.joints.values() if
                 j.chord.mtype == "top_chord" and 0.1 < j.loc < 0.9]
    bc_joints = [j for j in truss.joints.values() if
                 j.chord.mtype == "bottom_chord"]

    # Top Joints
    m = int(len(tc_joints) / 2)
    for i in range(m):
        joint = tc_joints[2 * i]
        sym_joint = tc_joints[2 * i + 1]
        j_group = {
            "value": joint.loc,
            "name": f"TopJoint {i} loc",
            "objects": [joint, sym_joint],
            "lb": max(0, joint.loc - 2 / n),
            "ub": min(1, joint.loc + 2 / n),
            "var_type": 'continuous',
            'property': 'loc'
        }
        groups.append(j_group)
        je_group = {
            "value": joint.e,
            "name": f"TopJoint {i} e",
            "objects": [joint, sym_joint],
            "lb": -50,
            "ub": 50,
            "var_type": 'continuous',
            'property': 'e'
        }
        groups.append(je_group)

    # Bottom joints
    n = int(len(bc_joints) / 2)
    for i in range(m):
        joint = bc_joints[2 * i]
        sym_joint = bc_joints[2 * i + 1]
        sym_joint.reverse = True

        j_group = {
            "value": joint.loc,
            "name": f"BottomJoint {i} loc",
            "objects": [joint, sym_joint],
            "lb": max(0, joint.loc - 0.1 / n),
            "ub": min(1, joint.loc + 0.1 / n),
            "var_type": 'continuous',
            'property': 'loc'
        }
        groups.append(j_group)
        je_group = {
            "value": joint.e,
            "name": f"BottomJoint {i} e",
            "objects": [joint, sym_joint],
            "lb": -50,
            "ub": 50,
            "var_type": 'continuous',
            'property': 'e'
        }
        groups.append(je_group)

    problem = TrussProblem(
        name="FrameTest",
        structure=truss,
        groups=groups
    )

    tjoint_vars = [var for var in problem.vars if "TopJoint" in var.name]
    bjoint_vars = [var for var in problem.vars if "BottomJoint" in var.name]

    # Joint loc constraints
    for i, jvar in enumerate(tjoint_vars):
        if i < len(tjoint_vars) - 1:
            a = np.zeros_like(problem.vars)
            idx = problem.vars.index(jvar)
            a[idx] = 1
            con = LinearConstraint(a=a, b=tjoint_vars[i + 1],
                                   name="Joint constraint")
            problem.add(con)

    # Joint loc constraints
    for i, jvar in enumerate(bjoint_vars):
        if i < len(bjoint_vars) - 1:
            a = np.zeros_like(problem.vars)
            idx = problem.vars.index(jvar)
            a[idx] = 1
            con = LinearConstraint(a=a, b=bjoint_vars[i + 1],
                                   name="Joint constraint")
            problem.add(con)

        else:
            a = np.zeros_like(problem.vars)
            idx = problem.vars.index(jvar)
            a[idx] = 1
            con = LinearConstraint(a=a, b=0.5,
                                   name="Joint constraint")
            problem.add(con)


    def mutate(individual, prob=0.01, stepsize=1, multiplier=0.1):
        """
        Mutates individual
        :param individual:
        :param prob:
        :return: mutated individual
        """
        for var in problem.vars:
            i = problem.all_vars.index(var)

            if isinstance(var, (IndexVariable, IntegerVariable)):
                val = np.random.randint(0, stepsize)

            else:
                val = np.random.uniform(0, multiplier * (var.ub - var.lb))

            if np.random.rand() < prob:
                if np.random.rand() < 0.5:
                    val *= -1
                individual[i] += val

        return individual


    # problem.index_to_binary()

    def cx_fun(A, B):

        C = np.arange(len(A))
        C = np.random.choice(C, 3, replace=False)
        for val in C:
            A[val], B[val] = B[val], A[val]
        return A, B


    solver = GA(pop_size=20,
                mut_fun=mutate,
                mutation_kwargs={'prob': 0.5, "stepsize": 5,
                                 "multiplier": 0.25})

    solver = VNS(step_length=3)
    x0 = [var.ub for var in problem.vars]
    fopt, xopt = solver.solve(problem, x0=x0, maxiter=50, verb=True, plot=True)


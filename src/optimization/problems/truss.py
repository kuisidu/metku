from src.optimization.structopt import *
from src.truss2d import *
from src.frame2d.frame2d import *

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

class PlaneTruss(OptimizationProblem):

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

            buckling_y, buckling_z, com_compression_bending_y, \
            com_compression_bending_z = \
                self.stability_constraints(mem)

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

            # BENDING + COMPRESSION Y
            com_compression_bending_con_y = NonLinearConstraint(
                con_fun=com_compression_bending_y,
                name="Com_compression_bending_y " + str(mem.mem_id),
                parent=self)
            com_compression_bending_con_y.fea_required = True
            self.cons.append(com_compression_bending_con_y)

            # BENDING + COMPRESSION Z
            com_compression_bending_con_z = NonLinearConstraint(
                con_fun=com_compression_bending_z,
                name="Com_compression_bending_z " + str(mem.mem_id),
                parent=self)
            com_compression_bending_con_z.fea_required = True
            self.cons.append(com_compression_bending_con_z)

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
        truss = Truss2D(simple=kwargs)

        sup1 = XYHingedSupport(truss.top_chords[0].coordinates[0])
        sup2 = YHingedSupport(truss.top_chords[-1].coordinates[1])

        # Top chord load
        for tc in truss.top_chords:
            truss.add(LineLoad(tc, kwargs['q'], kwargs['dir']))

        truss.add(sup1)
        truss.add(sup2)
        truss.generate()

        return truss


    def create_variables(self, profiles=RHS_PROFILES):
        """
        Creates optimization variables
        """
        self.vars = []
        for i, mem in enumerate(self.structure.members.values()):
            name = 'A' + str(i + 1)
            if self.prob_type == "discrete":
                var = DiscreteVariable(name,
                                       profiles=profiles,
                                       # TODO: THIS IS A PARAMETER
                                       values=RHS_A,
                                       target={"property": "PROFILE",
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

        def com_compression_bending_y(x):
            self.substitute_variables(x)
            return mem.steel_member.check_beamcolumn()[0] - 1

        def com_compression_bending_z(x):
            self.substitute_variables(x)
            return mem.steel_member.check_beamcolumn()[1] - 1

        return buckling_y,\
               buckling_z,\
               com_compression_bending_y,\
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

if __name__ == '__main__':

    from src.optimization.solvers import *

    kwargs = {'H1': 1000,
              'H2': 2000,
              'L1': 12500,
              'n': 16,
              'q': [-25, -25],
              'dir': 'y',
              'prob_type': 'discrete'}

    problem = PlaneTruss(**kwargs)
    x0 = [var.ub for var in problem.vars]
    print(x0)
    solver = VNS()
    solver.solve(problem, maxiter=100, x0=x0)


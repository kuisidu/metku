"""
@author:
"""

from src.frame2d.frame2d import *
from src.optimization.structopt import *


class WIColumn(OptimizationProblem):
    def __init__(self, L=5000, Fx=10e3, Fy=-200e3, Mz=0, top_flange_class=2,
                 bottom_flange_class=2, web_class=2):
        super().__init__("WIColumn")
        self.create_structure(L, Fx, Fy, Mz)
        self.create_variables()
        self.create_constraints()
        self.create_objective()

    def create_objective(self):
        def obj(x):
            return self.structure.weight

        self.obj = obj

    def create_structure(self, L, Fx, Fy, Mz):
        # Luo tyhjän kehän
        frame = Frame2D(num_elements=4)
        # Luo pilarin (koordinaatit, profile=vapaaehtoinen)
        col = SteelColumn([[0, 0], [0, L]], profile='WI 800-12-30X450-25X300')
        # Lisää pilarin kehälle
        # col.profile = 'WI 800-12-30X450-25X300'
        frame.add(col)
        # Lisää jäykän tuen pisteeseen (0,0)
        frame.add(FixedSupport([0, 0]))
        # Lisää pistekuorman pilarin yläpäähän
        frame.add(PointLoad([0, L], [Fx, Fy, Mz]))
        # Luo kehän fem -mallin
        frame.generate()
        # Laskee
        frame.calculate()
        # Piirtää kehän (pilarin) näytölle
        # frame.plot()

        self.structure = frame

    def create_variables(self):

        col = self.structure.members[0].cross_section

        var_h = Variable("h", 100, 800,
                         target={"property": "H", "objects": [col]})
        var_tt = Variable("tt", 5, 50,
                         target={"property": "TT", "objects": [col]})
        var_tb = Variable("tb", 5, 50,
                         target={"property": "TB", "objects": [col]})
        var_tw = Variable("tw", 5, 50,
                          target={"property": "TW", "objects": [col]})
        var_bt = Variable("bt", 100, 500,
                          target={"property": "BT", "objects": [col]})
        var_bb = Variable("bb", 100, 500,
                          target={"property": "BB", "objects": [col]})

        self.vars = [var_h, var_tt, var_tb, var_tw, var_bt, var_bb]

    # def section_class_constraint(self, mem):
    #
    #     def class_constraint(x):
    #         return mem.cross_section.C - 3
    #
    #     return class_constraint

    def WIColumnTopFlangeClassCon(self, top_flange_class):
        """
        Builds a constraint for cross-section class of a WI column
        0.5*bt - 0.5*tw - C0*e*tt <= sqrt(2)*aw
        """

        a = [0, 0, 0, 0, 0, 0]
        con_type = '<'

        e = self.structure.members.eps

        if top_flange_class == 1:
            C0 = 9
        elif top_flange_class == 2:
            C0 = 10
        elif top_flange_class > 2:
            C0 = 14

        for i in range(len(self.vars)):

            if self.vars[i].target.property == "BT":
                a[i] = 0.5
            elif self.vars[i].target.property == "TT":
                a[i] = - C0 * e
            elif self.vars[i].target.property == "TW":
                a[i] = 0.5

        b = math.sqrt(2) * self.structure.members.weld_throat

        # print(a)

        if top_flange_class > 3:
            """ if class 4 is required, the cf/tf ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type = '>'

        con_name = "Flange in class " + str(top_flange_class)
        con = self.structure.members.sop.LinearConstraint(a, b, con_type,
                                                          name=con_name)

        return con

    def WIColumnWebClassCon(self, web_class):
        """
        Builds a constraint for cross-section class of a WI column
        h - tt - tb - C1*e*tw <= 2*sqrt(2)*aw
        """

        a = [0, 0, 0, 0, 0, 0]
        con_type = '<'

        e = self.structure.members.eps

        # web is assumed to be in bending
        if web_class == 1:
            C1 = 72
        elif web_class == 2:
            C1 = 83
        elif web_class > 2:
            C1 = 124

        # web is assumed to be in compression
        # if web_class == 1:
        #     C1 = 33
        # elif web_class == 2:
        #     C1 = 38
        # elif web_class > 2:
        #     C1 = 42

        for i in range(len(self.vars)):

            if self.vars[i].target.property == "H":
                a[i] = 1
            elif self.vars[i].target.property == "TT":
                a[i] = - 1
            elif self.vars[i].target.property == "TB":
                a[i] = 1
            elif self.vars[i].target.property == "TW":
                a[i] = - C1 * e

        b = 2 * math.sqrt(2) * self.structure.members.weld_throat

        # print(a)

        if web_class > 3:
            """ if class 4 is required, the cw/tw ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type[1] = '>'

        con_name = "Web in class " + str(web_class)
        con = self.structure.members.sop.LinearConstraint(a, b, con_type,
                                                          name=con_name)

        return con

    def create_stability_constraints(self, mem, forces):

        N, V, M = forces

        def buckling(x):
            return -N / mem.NbRd - 1

        def com_compression_bending(x):
            pass

        def lt_buckling(x):
            pass

        return buckling, com_compression_bending, lt_buckling

    def create_section_constraints(self, sect, forces):

        N, V, M = forces

        def compression(x):
            return -N / sect.NRd - 1

        def tension(x):
            return N / sect.NRd - 1

        def shear(x):
            return abs(V) / sect.VRd - 1

        def bending_moment(x):
            return abs(M) / sect.MRd[0] - 1

        return compression, tension, shear, bending_moment

    def create_constraints(self):
        for mem in self.structure.members.values():
            for i, elem in enumerate(mem.elements.values()):
                forces = [elem.axial_force[0], elem.shear_force[0],
                          elem.bending_moment[0]]
                compression, tension, shear, bending_moment = \
                    self.create_section_constraints(mem, forces)

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

                if i == len(mem.elements) - 1:
                    forces = [elem.axial_force[1], elem.shear_force[1],
                              elem.bending_moment[1]]
                    compression, tension, shear, bending_moment = \
                        self.create_section_constraints(mem, forces)

                    compression_con = NonLinearConstraint(con_fun=compression,
                                                          name="Compression " +
                                                               str(mem.mem_id)
                                                               + str(i+1),
                                                          parent=self)
                    compression_con.fea_required = True
                    tension_con = NonLinearConstraint(con_fun=tension,
                                                      name="Tension " +
                                                           str(mem.mem_id) +
                                                           str(i+1),
                                                      parent=self)
                    tension_con.fea_required = True
                    shear_con = NonLinearConstraint(con_fun=shear,
                                                    name="Shear " +
                                                         str(mem.mem_id) +
                                                         str(i+1),
                                                    parent=self)
                    shear_con.fea_required = True

                    bending_moment_con = NonLinearConstraint(
                        con_fun=bending_moment, name="Bending_moment " +
                                                     str(mem.mem_id) +
                                                     str(i+1), parent=self)
                    bending_moment_con.fea_required = True

                    self.cons.extend([compression_con, tension_con, shear_con,
                                      bending_moment_con])

            class_constraint = self.section_class_constraint(mem)

            class_con = NonLinearConstraint(
                con_fun=class_constraint,
                name="Class_constraint " + str(mem.mem_id))
            self.cons.append(class_con)


if __name__ == "__main__":
    from src.optimization.solvers import *
    problem = WIColumn()
    x0 = [200, 5, 5, 5, 100, 100]
    solver = SLP(step_length=10)
    solver.solve(problem, maxiter=10000, maxtime=10, x0=x0)
    problem(solver.X)



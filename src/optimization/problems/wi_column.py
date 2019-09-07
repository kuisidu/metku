"""
@author: Victoria
"""

from src.frame2d.frame2d import *
from src.optimization.structopt import *


class WIColumn(OptimizationProblem):
    """
        Pilarin kuormitukset ja lähtötiedot

                Lpi -- pilarin pituus
                F -- pistekuormat
                    Fx .. vaakasuuntainen
                    Fy .. pystysuuntainen
                Q -- tasaiset kuormat
                    Qx .. vaakasuuntainen
                    Qy .. pystysuuntainen
                Mz .. pistemomentti
                Poikkileikkausluokat ylä- ja alalaipalle sekä uumalle
                    top_flange_class
                    bottom_flange_class
                    web_class
                symmetry .. kaksoissymmetrinen: dual tai monosymmetrinen: mono
                buckling_z -- heikomman suunnan nurjahdus
                    True .. rajoitusehto huomioidaan
                    False .. rajoitusehtoa ei huomioida
                LT_buckling -- kiepahdus
                    True .. rajoitusehto huomioidaan
                    False .. rajoitusehtoa ei huomioida
    """
    def __init__(self, Lpi=6000, Fx=800, Fy=-280e3, Qx=5.85, Qy=0, Mz=0, lcr=2,
                 top_flange_class=2, bottom_flange_class=2, web_class=2,
                 symmetry="dual", buckling_z=True, LT_buckling=False):
        super().__init__("WIColumn")
        self.cons.clear()
        self.LT_buckling = LT_buckling
        self.buckling_z = buckling_z
        self.symmetry = symmetry
        self.top_flange_class = top_flange_class
        self.bottom_flange_class = bottom_flange_class
        self.web_class = web_class
        self.create_structure(Lpi, Fx, Fy, Qx, Qy, Mz, lcr, LT_buckling)
        self.create_variables()
        self.create_constraints()
        self.create_objective()
        self.prob_type = 'continuous'

    def create_objective(self):
        def obj(x):
            return self.structure.weight

        self.obj = obj

    def create_structure(self, Lpi, Fx, Fy, Qx, Qy, Mz, lcr, LT_buckling):
        # Luo tyhjän kehän
        frame = Frame2D(num_elements=4)
        # Luo pilarin (koordinaatit, profile=vapaaehtoinen)
        col = SteelColumn([[0, 0], [0, Lpi]], LT_buckling,
                          profile='WI 500-12-10X300-10X300')
        # Lisätään nurjahduspituus (masto lcr[0]=2, muuten lcr[0]=0.7)
        col.steel_member.lcr[0] = lcr
        # Lisää pilarin kehälle
        # col.profile = 'WI 800-12-30X450-25X300'
        frame.add(col)
        # Lisää niveltuen pisteeseen (0,Lpi)
        if lcr == 0.7:
            frame.add(XHingedSupport([0, Lpi]))
        # Lisää jäykän tuen pisteeseen (0,0)
        frame.add(FixedSupport([0, 0]))
        # Lisää pistekuorman pilarin yläpäähän
        frame.add(PointLoad([0, Lpi], [Fx, Fy, Mz]))
        # Lisää tasaisen kuorman rakenneosalle
        frame.add(LineLoad(col, [Qx, Qx], "x"))
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
        var_tw = Variable("tw", 5, 50,
                          target={"property": "TW", "objects": [col]})
        if self.symmetry == "mono":
            var_tt = Variable("tt", 5, 50,
                             target={"property": "TT", "objects": [col]})
            var_tb = Variable("tb", 5, 50,
                             target={"property": "TB", "objects": [col]})
            var_bt = Variable("bt", 100, 500,
                              target={"property": "BT", "objects": [col]})
            var_bb = Variable("bb", 100, 500,
                              target={"property": "BB", "objects": [col]})

            self.vars = [var_h, var_tt, var_tb, var_tw, var_bt, var_bb]

        elif self.symmetry == "dual":
            var_tf = Variable("tf", 5, 50,
                              target={"property": "TF", "objects": [col]})
            var_bf = Variable("bf", 100, 500,
                              target={"property": "BF", "objects": [col]})

            self.vars = [var_h, var_tw, var_bf, var_tf]

        else:
            raise ValueError("Symmetry must be either dual or mono")

    # def section_class_constraint(self, mem):
    #
    #     def class_constraint(x):
    #         return mem.cross_section.C - 3
    #
    #     return class_constraint

    def WIColumnTopFlangeClassCon(self, mem):
        """
        Builds a constraint for cross-section class of a WI column
        0.5*bt - 0.5*tw - C0*e*tt <= sqrt(2)*aw
        """

        a = np.zeros_like(self.vars)
        con_type = '<'

        e = mem.cross_section.eps

        if self.top_flange_class == 1:
            C0 = 9
        elif self.top_flange_class == 2:
            C0 = 10
        elif self.top_flange_class > 2:
            C0 = 14

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "BT" or \
                    self.vars[i].target["property"] == "BF":
                a[i] = 0.5
            elif self.vars[i].target["property"] == "TT" or \
                    self.vars[i].target["property"] == "TF":
                a[i] = - C0 * e
            elif self.vars[i].target["property"] == "TW":
                a[i] = - 0.5

        b = math.sqrt(2) * mem.cross_section.weld_throat

        # print(a)

        if self.top_flange_class > 3:
            """ if class 4 is required, the cf/tf ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type = '>'

        con_name = "Flange in class " + str(self.top_flange_class)
        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def WIColumnBottomFlangeClassCon(self, mem):
        """
        Builds a constraint for cross-section class of a WI column
        0.5*bt - 0.5*tw - C0*e*tt <= sqrt(2)*aw
        """

        a = np.zeros_like(self.vars)
        con_type = '<'

        e = mem.cross_section.eps

        if self.bottom_flange_class == 1:
            C0 = 9
        elif self.bottom_flange_class == 2:
            C0 = 10
        elif self.bottom_flange_class > 2:
            C0 = 14

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "BT":
                a[i] = 0.5
            elif self.vars[i].target["property"] == "TT":
                a[i] = - C0 * e
            elif self.vars[i].target["property"] == "TW":
                a[i] = - 0.5

        b = math.sqrt(2) * mem.cross_section.weld_throat

        # print(a)

        if self.top_flange_class > 3:
            """ if class 4 is required, the cf/tf ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type = '>'

        con_name = "Flange in class " + str(self.bottom_flange_class)
        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def WIColumnWebClassCon(self, mem):
        """
        Builds a constraint for cross-section class of a WI column
        h - tt - tb - C1*e*tw <= 2*sqrt(2)*aw
        """

        a = np.zeros_like(self.vars)
        con_type = '<'

        e = mem.cross_section.eps

        # web is assumed to be in bending
        # if self.web_class == 1:
        #     C1 = 72
        # elif self.web_class == 2:
        #     C1 = 83
        # elif self.web_class > 2:
        #     C1 = 124

        # web is assumed to be in compression
        if self.web_class == 1:
            C1 = 33
        elif self.web_class == 2:
            C1 = 38
        elif self.web_class > 2:
            C1 = 42

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "H":
                a[i] = 1
            elif self.vars[i].target["property"] == "TT":
                a[i] = - 1
            elif self.vars[i].target["property"] == "TB":
                a[i] = - 1
            elif self.vars[i].target["property"] == "TW":
                a[i] = - C1 * e
            elif self.vars[i].target["property"] == "TF":
                a[i] = -2

        b = 2 * math.sqrt(2) * mem.cross_section.weld_throat

        # print(a)

        if self.web_class > 3:
            """ if class 4 is required, the cw/tw ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type = '>'

        con_name = "Web in class " + str(self.web_class)
        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def create_stability_constraints(self, mem):

        def buckling_y(x):
            return -mem.ned / mem.NbRd[0] - 1

        def buckling_z(x):
            return -mem.ned / mem.NbRd[1] - 1

        def com_compression_bending_y(x):
            return mem.steel_member.check_beamcolumn()[0] - 1

        def com_compression_bending_z(x):
            return mem.steel_member.check_beamcolumn()[1] - 1

        def lt_buckling(x):
            #  med = max(abs(np.min(np.array(mem.steel_member.myed))),
            #  np.max(np.array(mem.steel_member.myed)))
            return mem.med / mem.MbRd - 1

        return \
            buckling_y, buckling_z, com_compression_bending_y, \
            com_compression_bending_z, lt_buckling

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

            buckling_y, buckling_z, com_compression_bending_y, \
                com_compression_bending_z, lt_buckling = \
                self.create_stability_constraints(mem)

            buckling_y_con = NonLinearConstraint(con_fun=buckling_y,
                                                 name="Buckling_y " +
                                                      str(mem.mem_id),
                                                 parent=self)
            buckling_y_con.fea_required = True
            self.cons.append(buckling_y_con)

            if self.buckling_z:

                buckling_z_con = NonLinearConstraint(con_fun=buckling_z,
                                                     name="Buckling_z " +
                                                          str(mem.mem_id),
                                                     parent=self)
                buckling_z_con.fea_required = True
                self.cons.append(buckling_z_con)

            com_compression_bending_con_y = NonLinearConstraint(
                con_fun=com_compression_bending_y,
                name="Com_compression_bending_y " + str(mem.mem_id),
                parent=self)
            com_compression_bending_con_y.fea_required = True
            self.cons.append(com_compression_bending_con_y)

            com_compression_bending_con_z = NonLinearConstraint(
                con_fun=com_compression_bending_z,
                name="Com_compression_bending_z " + str(mem.mem_id),
                parent=self)
            com_compression_bending_con_z.fea_required = True
            self.cons.append(com_compression_bending_con_z)

            if self.LT_buckling:

                lt_buckling_con = NonLinearConstraint(con_fun=lt_buckling,
                                                      name="LT_buckling " +
                                                           str(mem.mem_id),
                                                      parent=self)
                lt_buckling_con.fea_required = True
                self.cons.append(lt_buckling_con)

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

            self.cons.append(self.WIColumnWebClassCon(mem))
            self.cons.append(self.WIColumnTopFlangeClassCon(mem))
            if self.symmetry == "mono":
                self.cons.append(self.WIColumnBottomFlangeClassCon(mem))


if __name__ == "__main__":
    from src.optimization.solvers import *
    problem = WIColumn()
    x0 = [800, 50, 500, 50]
    solver = SLP(move_limits=[0.9, 6])
    solver.solve(problem, maxiter=50000, maxtime=30, x0=x0)
    problem(solver.X, prec=5)
    from src.optimization.result_exporter import *
    name = "WIColumn_buckling_z:{0}_LT_buckling:{1}"\
        .format(problem.buckling_z, problem.LT_buckling)
    ResultExporter(problem, solver).to_csv()

    # solver = SLSQP()
    # solver.solve(problem, maxiter=50000, x0=x0)
    # problem(solver.X)

    # r0 = problem.solve("slsqp", x0=x0)
    # print(r0)

    # print(problem.structure.f.elements[0].bending_moment)
    # print(problem.structure.f.elements[0].axial_force)
    # print(problem.structure.f.loads[1].qval)


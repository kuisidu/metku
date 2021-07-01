"""

Minimum cost design of end-plate connections

@author: Kristo Mela
@date: 31.1.2020
"""

import numpy as np
import math
import time

MIN_HEIGHT = 100
MAX_HEIGHT = 500
step_height = 10
MIN_WIDTH = 100
MAX_WIDTH = 500
step_width = 10

HEIGHTS = np.arange(MIN_HEIGHT, MAX_HEIGHT+step_height, step_height)
WIDTHS = np.arange(MIN_WIDTH, MAX_WIDTH+step_width, step_width)
#THICKNESSES = [5, 6, 8, 10, 12, 14, 15, 16, 18, 20, 22, 25, 30, 35, 40, 50]
THICKNESSES = [4, 5, 6, 8, 10, 12, 14, 15, 16, 18, 20]
MAX_THICK = max(THICKNESSES)
MIN_THICK = min(THICKNESSES)

try:
    from metku.frame2d.frame2d import *
    import metku.frame2d.frame2d as f2d
    from metku.optimization.structopt import *
    # import metku.optimization.structopt as sopt
    from metku.optimization.solvers.trust_region import TrustRegionConstr
    from metku.optimization.solvers.bnb import BnB
    from metku.optimization.solvers.lp import LP
except:
    
    from optimization.structopt import *
    from optimization.solvers.trust_region import TrustRegionConstr    
        
class EndPlateJointOpt(OptimizationProblem):
    """
        
    """
    def __init__(self,end_plate_joint):
        super().__init__("End plate connection")

        self.prob_type = prob_type
        self.cons.clear()        
        self.clear_vars()
        self.structure = None        
        self.LT_buckling = LT_buckling
        self.buckling_z = buckling_z
        self.symmetry = symmetry
        self.top_flange_class = top_flange_class
        self.bottom_flange_class = bottom_flange_class
        self.web_class = web_class
        self.L = L
        self.Lpi = Lpi
        self.Fx = Fx
        self.Fy = Fy
        self.Qx = Qx
        self.lcr = lcr
        
        section_class = max(top_flange_class, bottom_flange_class)
        self.create_structure(Lpi, Fx, Fy, Qx, Qy, Mz, lcr, LT_buckling)
        self.create_variables()
        self.create_constraints(section_class)
        self.create_objective()
        self.prob_type = 'continuous'

    def create_objective(self):
        def obj_fun(x):
            self.substitute_variables(x)
            return self.structure.weight
        obj = ObjectiveFunction("weight", obj_fun,problem=self)
        self.add(obj)

    def create_structure(self, Lpi, Fx, Fy, Qx, Qy, Mz, lcr, LT_buckling):
        # Luo tyhjän kehän
        frame = Frame2D(num_elements=4)
        # Luo pilarin (koordinaatit, profile=vapaaehtoinen)
        col = SteelColumn([[0, 0], [0, Lpi]], LT_buckling,
                          profile='WI 500-12-10X300-10X300',
                          material='S355MC')
        # Lisätään pilarille materiaali
        col.material = "S355MC"
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
        frame.add(f2d.PointLoad([0, Lpi], [Fx, Fy, Mz]))
        # Lisää tasaisen kuorman rakenneosalle
        frame.add(f2d.LineLoad(col, [Qx, Qx], "x"))
        # Luo kehän fem -mallin
        frame.generate()
        # Laskee
        frame.calculate()
        # Piirtää kehän (pilarin) näytölle
        # frame.plot()

        self.structure = frame

    def create_variables(self):
        """
        Creates design variables
        """

        col = self.structure.members[0].cross_section

        if self.prob_type == "continuous":

            var_h = Variable("h", MIN_HEIGHT, MAX_HEIGHT,
                             target={"property": "h", "objects": [col]})
            var_tw = Variable("tw", MIN_THICK, MAX_THICK,
                              target={"property": "tw", "objects": [col]})
            if self.symmetry == "mono":
                var_tt = Variable("tt", MIN_THICK, MAX_THICK,
                                 target={"property": "tt", "objects": [col]})
                var_tb = Variable("tb", MIN_THICK, MAX_THICK,
                                 target={"property": "tb", "objects": [col]})
                var_bt = Variable("bt", MIN_WIDTH, MAX_WIDTH,
                                  target={"property": "bt", "objects": [col]})
                var_bb = Variable("bb", MIN_WIDTH, MAX_WIDTH,
                                  target={"property": "bb", "objects": [col]})

                vars = [var_h, var_tw, var_bt, var_tt, var_bb, var_tb]
                for var in vars:
                    self.add(var)

            elif self.symmetry == "dual":
                var_tf = Variable("tf", MIN_THICK, MAX_THICK,
                                  target={"property": "tf", "objects": [col]})
                var_bf = Variable("b", MIN_WIDTH, MAX_WIDTH,
                                  target={"property": "b", "objects": [col]})

                vars = [var_h, var_tw, var_bf, var_tf]
                for var in vars:
                    self.add(var)

            else:
                raise ValueError("Symmetry must be either dual or mono")

        elif self.prob_type == "discrete":
            var_h = DiscreteVariable(
                "h", values=HEIGHTS,
                target={"property": "h", "objects": [col]})
            var_tw = DiscreteVariable(
                "tw", values=THICKNESSES,
                target={"property": "tw", "objects": [col]})

            var_tw.branch_priority = 1
            if self.symmetry == "mono":
                var_tt = DiscreteVariable(
                    "tt", values=THICKNESSES,
                    target={"property": "tt", "objects": [col]})
                var_tb = DiscreteVariable(
                    "tb", values=THICKNESSES,
                    target={"property": "tb", "objects": [col]})
                var_bt = DiscreteVariable(
                    "bt", values=WIDTHS,
                    target={"property": "bt", "objects": [col]})
                var_bb = DiscreteVariable(
                    "bb", values=WIDTHS,
                    target={"property": "bb", "objects": [col]})

                vars = [var_h, var_tw, var_bt, var_tt, var_bb, var_tb]
                for var in vars:
                    self.add(var)

            elif self.symmetry == "dual":
                var_tf = DiscreteVariable(
                    "tf", values=THICKNESSES,
                    target={"property": "tf", "objects": [col]})
                
                var_tf.branch_priority = 1
                var_bf = DiscreteVariable(
                    "b", values=WIDTHS,
                    target={"property": "b", "objects": [col]})

                vars = [var_h, var_tw, var_bf, var_tf]
                for var in vars:
                    self.add(var)

            else:
                raise ValueError("Symmetry must be either dual or mono")

        elif self.prob_type == "index":
            var_h = IndexVariable(
                "h", values=HEIGHTS,
                target={"property": "h", "objects": [col]})
            var_tw = IndexVariable(
                "tw", values=THICKNESSES,
                target={"property": "tw", "objects": [col]})

            if self.symmetry == "mono":
                var_tt = IndexVariable(
                    "tt", values=THICKNESSES,
                    target={"property": "tt", "objects": [col]})
                var_tb = IndexVariable(
                    "tb", values=THICKNESSES,
                    target={"property": "tb", "objects": [col]})
                var_bt = IndexVariable(
                    "bt", values=WIDTHS,
                    target={"property": "bt", "objects": [col]})
                var_bb = IndexVariable(
                    "bb", values=WIDTHS,
                    target={"property": "bb", "objects": [col]})

                vars = [var_h, var_tw, var_bt, var_tt, var_bb, var_tb]
                for var in vars:
                    self.add(var)

            elif self.symmetry == "dual":
                var_tf = IndexVariable(
                    "tf", values=THICKNESSES,
                    target={"property": "tf", "objects": [col]})
                var_bf = IndexVariable(
                    "b", values=WIDTHS,
                    target={"property": "b", "objects": [col]})

                vars = [var_h, var_tw, var_bf, var_tf]
                for var in vars:
                    self.add(var)

            else:
                raise ValueError("Symmetry must be either dual or mono")

    def WIColumnTopFlangeClassCon(self, mem):
        """
        Builds a constraint for cross-section class of a WI column
        0.5*bt - 0.5*tw - C0*e*tt <= sqrt(2)*aw
        """

        #a = np.zeros_like(self.vars)
        a = np.zeros(self.nvars())
        con_type = '<'

        e = mem.cross_section.eps

        if self.top_flange_class == 1:
            C0 = 9
        elif self.top_flange_class == 2:
            C0 = 10
        elif self.top_flange_class > 2:
            C0 = 14

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "bt" or \
                    self.vars[i].target["property"] == "b":
                a[i] = 0.5
            elif self.vars[i].target["property"] == "tt" or \
                    self.vars[i].target["property"] == "tf":
                a[i] = - C0 * e
            elif self.vars[i].target["property"] == "tw":
                a[i] = - 0.5

        b = math.sqrt(2) * mem.cross_section.weld_throat

        # print(a)

        if self.top_flange_class > 3:
            """ if class 4 is required, the cf/tf ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type = '>'

        con_name = "Top flange in class " + str(self.top_flange_class)
        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def WIColumnBottomFlangeClassCon(self, mem):
        """
        Builds a constraint for cross-section class of a WI column
        0.5*bt - 0.5*tw - C0*e*tt <= sqrt(2)*aw
        """

        #a = np.zeros_like(self.vars)
        a = np.zeros(self.nvars())
        con_type = '<'

        e = mem.cross_section.eps

        if self.bottom_flange_class == 1:
            C0 = 9
        elif self.bottom_flange_class == 2:
            C0 = 10
        elif self.bottom_flange_class > 2:
            C0 = 14

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "bt":
                a[i] = 0.5
            elif self.vars[i].target["property"] == "tt":
                a[i] = - C0 * e
            elif self.vars[i].target["property"] == "tw":
                a[i] = - 0.5

        b = math.sqrt(2) * mem.cross_section.weld_throat

        # print(a)

        if self.top_flange_class > 3:
            """ if class 4 is required, the cf/tf ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type = '>'

        con_name = "Bottom flange in class " + str(self.bottom_flange_class)
        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def WIColumnWebClassCon(self, mem):
        """
        Builds a constraint for cross-section class of a WI column
        h - tt - tb - C1*e*tw <= 2*sqrt(2)*aw
        """

        #a = np.zeros_like(self.vars)
        a = np.zeros(self.nvars())
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

            if self.vars[i].target["property"] == "h":
                a[i] = 1
            elif self.vars[i].target["property"] == "tt":
                a[i] = - 1
            elif self.vars[i].target["property"] == "tb":
                a[i] = - 1
            elif self.vars[i].target["property"] == "tw":
                a[i] = - C1 * e
            elif self.vars[i].target["property"] == "tf":
                a[i] = -2

        b = 2 * math.sqrt(2) * mem.cross_section.weld_throat

        # print("PRINT", a)

        if self.web_class > 3:
            """ if class 4 is required, the cw/tw ratio needs to
                be greater than the class 3 limit. This  changes the
                direction of the constraint from < to >.
            """
            con_type = '>'

        con_name = "Web in class " + str(self.web_class)
        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def WIColumnWebHeightCon(self, h_min=50):
        """
        Constraint for ensuring that the web height is positive, at least 'h_min'
        h-tt-tb >= h_min, or
        -h + tt + tb <= -h_min
        """

        #a = np.zeros_like(self.vars)
        a = np.zeros(self.nvars())

        con_type = '<'

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "h":
                a[i] = -1
            elif self.vars[i].target["property"] == "tt":
                a[i] = + 1
            elif self.vars[i].target["property"] == "tb":
                a[i] = + 1            
            elif self.vars[i].target["property"] == "tf":
                a[i] = +2

        b = -h_min

        con_name = "Web height at least " + str(h_min)

        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def TopFlangeShearLagCon(self, mem):
        """
        Builds a constraint for shear lag
        0.5 * bt - 0.5 * tw <= Lpi/50 + sqrt(2)*aw
        """

        a = np.zeros(self.nvars())
        con_type = '<'

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "bt" or \
                    self.vars[i].target["property"] == "b":
                a[i] = +0.5
            elif self.vars[i].target["property"] == "tw":
                a[i] = -0.5

        b = self.Lpi / 50 + math.sqrt(2) * mem.cross_section.weld_throat

        con_name = "Top flange shear buckling"

        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    def BottomFlangeShearLagCon(self, mem):
        """
        Builds a constraint for shear lag
        0.5 * bb - 0.5 * tw <= Lpi/50 + sqrt(2)*aw
        """

        a = np.zeros(self.nvars())
        con_type = '<'

        for i in range(len(self.vars)):

            if self.vars[i].target["property"] == "bb":
                a[i] = +0.5
            elif self.vars[i].target["property"] == "tw":
                a[i] = -0.5

        b = self.Lpi / 50 + math.sqrt(2) * mem.cross_section.weld_throat

        con_name = "Bottom flange shear buckling"

        con = LinearConstraint(a, b, con_type, name=con_name)

        return con

    # def ShearBucklingCon(self, mem):
    #     """
    #     Builds a constraint for shear buckling resistance
    #     hw - (72 * e / ny) * tw <= 0
    #     """
    #
    #     a = np.zeros(self.nvars())
    #     con_type = '<'
    #
    #     e = mem.cross_section.eps
    #
    #     fy =
    #     if fy > 460:
    #         ny = 1
    #     else:
    #         ny = 1.2
    #
    #     for i in range(len(self.vars)):
    #
    #         if self.vars[i].target["property"] == "HW":
    #             a[i] = +1
    #         elif self.vars[i].target["property"] == "tw":
    #             a[i] = -0.5
    #
    #     b = self.Lpi / 50 + math.sqrt(2) * mem.cross_section.weld_throat
    #
    #     con_name = "Shear lag "
    #
    #     con = LinearConstraint(a, b, con_type, name=con_name)
    #
    #     return con

    def create_stability_constraints(self, mem, section_class):
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
            return mem.member.check_beamcolumn(section_class=section_class)[0] - 1

        def com_compression_bending_z(x):
            self.substitute_variables(x)
            return mem.member.check_beamcolumn(section_class=section_class)[1] - 1

        def lt_buckling(x):
            self.substitute_variables(x)
            #  med = max(abs(np.min(np.array(mem.steel_member.myed))),
            #  np.max(np.array(mem.steel_member.myed)))
            # print("mem.med = {0:4.2f}".format(mem.med*1e-6),
            #       "mem.MbRd = {0:4.2f}".format(mem.MbRd*1e-6))
            if x[0] < 0:
                print(x)
            return abs(mem.med) / mem.MbRd - 1

        return \
            buckling_y, buckling_z, com_compression_bending_y, \
            com_compression_bending_z, lt_buckling

    def create_section_constraints(self, sect, elem):
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
            M = elem.bending_moment[0]
            return abs(M) / sect.MRd[0] - 1

        return compression, tension, shear, bending_moment

    def create_constraints(self, section_class):
        """ Create constraints for each member 
            (for the single column the list is a single member)
        """
        for mem in self.structure.members.values():

            
            buckling_y, buckling_z, com_compression_bending_y, \
                com_compression_bending_z, lt_buckling = \
                self.create_stability_constraints(mem, section_class)

            buckling_y_con = NonLinearConstraint(con_fun=buckling_y,
                                                 name="Buckling_y " +
                                                      str(mem.mem_id),
                                                 )
            buckling_y_con.fea_required = True
            
            com_compression_bending_con_y = NonLinearConstraint(
                con_fun=com_compression_bending_y,
                name="Com_compression_bending_y " + str(mem.mem_id),
                )
            com_compression_bending_con_y.fea_required = True
            self.add(com_compression_bending_con_y)

            # self.add(buckling_y_con)
            
            if self.buckling_z:

                buckling_z_con = NonLinearConstraint(con_fun=buckling_z,
                                                     name="Buckling_z " +
                                                          str(mem.mem_id),
                                                     )
                buckling_z_con.fea_required = True

                # self.add(buckling_z_con)

                com_compression_bending_con_z = NonLinearConstraint(
                    con_fun=com_compression_bending_z,
                    name="Com_compression_bending_z " + str(mem.mem_id),
                    )
                com_compression_bending_con_z.fea_required = True
                self.add(com_compression_bending_con_z)
            
            if self.LT_buckling:

                lt_buckling_con = NonLinearConstraint(con_fun=lt_buckling,
                                                      name="LT_buckling " +
                                                           str(mem.mem_id),
                                                      )
                lt_buckling_con.fea_required = True
                self.add(lt_buckling_con)

            """
            for i, elem in enumerate(mem.elements.values()):
                forces = [elem.axial_force[0], elem.shear_force[0],
                          elem.bending_moment[0]]
                compression, tension, shear, bending_moment = \
                    self.create_section_constraints(mem, elem)

                compression_con = NonLinearConstraint(con_fun=compression,
                                                      name="Compression " +
                                                           str(mem.mem_id) +
                                                           str(i), )
                compression_con.fea_required = True

                tension_con = NonLinearConstraint(con_fun=tension,
                                                  name="Tension " +
                                                       str(mem.mem_id) +
                                                       str(i), )
                tension_con.fea_required = True

                shear_con = NonLinearConstraint(con_fun=shear,
                                                name="Shear " + str(mem.mem_id)
                                                     + str(i), )
                shear_con.fea_required = True

                bending_moment_con = NonLinearConstraint(
                    con_fun=bending_moment, name="Bending_moment " +
                                                 str(mem.mem_id) +
                                                 str(i), )
                bending_moment_con.fea_required = True

                self.cons.extend([compression_con, tension_con, shear_con,
                                  bending_moment_con])

                if i == len(mem.elements) - 1:
                    forces = [elem.axial_force[1], elem.shear_force[1],
                              elem.bending_moment[1]]
                    compression, tension, shear, bending_moment = \
                        self.create_section_constraints(mem, elem)

                    compression_con = NonLinearConstraint(con_fun=compression,
                                                          name="Compression " +
                                                               str(mem.mem_id)
                                                               + str(i+1),
                                                          )
                    compression_con.fea_required = True
                    tension_con = NonLinearConstraint(con_fun=tension,
                                                      name="Tension " +
                                                           str(mem.mem_id) +
                                                           str(i+1),
                                                      )
                    tension_con.fea_required = True
                    shear_con = NonLinearConstraint(con_fun=shear,
                                                    name="Shear " +
                                                         str(mem.mem_id) +
                                                         str(i+1),
                                                    )
                    shear_con.fea_required = True

                    bending_moment_con = NonLinearConstraint(
                        con_fun=bending_moment, name="Bending_moment " +
                                                     str(mem.mem_id) +
                                                     str(i+1), )
                    bending_moment_con.fea_required = True

                    self.cons.extend([compression_con, tension_con, shear_con,
                                      bending_moment_con])
            """

            self.add(self.WIColumnWebClassCon(mem))
            self.add(self.WIColumnTopFlangeClassCon(mem))
            self.add(self.TopFlangeShearLagCon(mem))
            if self.symmetry == "mono":
                self.add(self.WIColumnBottomFlangeClassCon(mem))
                self.add(self.BottomFlangeShearLagCon(mem))
                
            self.add(self.WIColumnWebHeightCon(h_min=50))


if __name__ == "__main__":
    from optimization.solvers import *
    from optimization.solvers.bnb import BnB
    from optimization.result_exporter import *

    problem = WIColumn(prob_type='discrete',
                       Lpi=8000,                       
                       top_flange_class=3,
                       bottom_flange_class=3,
                       web_class=3,
                       buckling_z=True,
                       LT_buckling=True)

    # x0 = [300, 8, 200, 10, 200, 10]
    # x0 = [300, 8, 200, 10]
    x0 = [200, 10, 170, 12]

    # x0 = [var.ub for var in problem.vars]

    # # TrustRegionConstr
    # solver = TrustRegionConstr()
    # f_best, x_best, nit = solver.solve(problem,
    #                                    maxiter=200,
    #                                    x0=x0)
    # # print(x_best)
    # problem.num_iters = nit
    # problem(x_best, prec=5)

    # BnB
    lb_solver = TrustRegionConstr()
    solver = BnB(problem, lb_solver)

    # return solver, problem

    # break

    solver.solve(problem, x0=x0, verb=1)

    print(solver.X)
    print(solver.best_x)
    print(solver.best_f)

    # SLP
    # solver = SLP(move_limits=[0.9, 6])
    # solver.solve(problem,
    #              maxiter=50000,
    #              maxtime=30,
    #              x0=x0)
    # problem(solver.X, prec=5)

    # SLSQP
    # solver = slsqp.SLSQP()
    # f_best, x_best = solver.solve(problem,
    #                               maxiter=100,
    #                               x0=x0)
    # problem(solver.best_x, prec=5)

    # MISLP
    # solver = MISLP(move_limits=[0.5, 5])
    # # problem(x0)
    # solver.solve(problem,
    #              maxiter=100,
    #              x0=x0,
    #              min_diff=1e-2,
    #              verb=True)
    # problem(solver.X, prec=5)

    # 2-vaihetekniikalla
    # solver1 = SLP(move_limits=[0.9, 4])
    # solver1 = TrustRegionConstr()
    # solver2 = MISLP(move_limits=[0.85, 1.5])
    # solver = TwoPhase(solver1, solver2, limits=[3, 3])
    # fopt, xopt = solver.solve(problem,
    #                           x0=x0,
    #                           maxiter=50,
    #                           # min_diff=1e-6,
    #                           # verb=True
    #                           )
    # problem(xopt)

    import matplotlib.pyplot as plt

    """
    fvals = problem.fvals
    X = np.arange(len(fvals))
    plt.plot(X, fvals)
    plt.show()
    """

    #ResultExporter(problem, solver).to_csv()
    #ResultExporter(problem, solver).csv_to_excel()

    seconds = time.process_time()
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    print("Process time:", "%d:%02d:%02d" % (h, m, s))

    # problem.structure.members[0].cross_section.draw()

    # G = {}
    # for con in problem.cons:
    #     name = con.name
    #     G[name] = 0

    # print(problem.structure.f.elements[0].bending_moment)
    # print(problem.structure.f.elements[0].axial_force)
    # print(problem.structure.f.loads[1].qval)
    # print("PINTA_ALA=", problem.structure.members[0].A)



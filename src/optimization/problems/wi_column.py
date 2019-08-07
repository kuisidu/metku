
from src.frame2d.frame2d import *
from src.optimization.structopt import *


class WIColumn(OptimizationProblem):
    def __init__(self, L=5000, Fx=10e3, Fy=-200e3, Mz=0):
        super().__init__("WIColumn")
        self.create_structure(L, Fx, Fy, Mz)

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

    def create_stability_constraints(self, mem):

        def buckling(x):
            pass

        def com_compression_bending(x):
            pass

        def lt_buckling(x):
            pass

    def create_section_constraints(self, sect, forces):

        N, V, M = forces

        def compression(x):
            return -N-sect.NRd

        def tension(x):
            return N-sect.NRd

        def shear(x):
            return abs(V)-sect.VRd

        def bending_moment(x):
            return abs(M)-sect.MRd

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
                shear_con = NonLinearConstraint(con_fun=shear, name="Shear " + str(mem.mem_id) + str(i), parent=self)
                shear_con.fea_required = True
                bending_moment_con = NonLinearConstraint(con_fun=bending_moment,
                                                name="Bending_moment " +
                                                     str(mem.mem_id) +
                                                     str(i), parent=self)
                bending_moment_con.fea_required = True

                self.cons.extend([compression_con, tension_con, shear_con, bending_moment_con])

                if i == len(mem.elements) - 1:
                    forces = [elem.axial_force[1], elem.shear_force[1],
                              elem.bending_moment[1]]
                    compression, tension, shear, bending_moment = \
                        self.create_section_constraints(mem, forces)

                    compression_con = NonLinearConstraint(con_fun=compression,
                                                          name="Compression " +
                                                               str(
                                                                   mem.mem_id) +
                                                               str(i+1),
                                                          parent=self)
                    compression_con.fea_required = True
                    tension_con = NonLinearConstraint(con_fun=tension,
                                                      name="Tension " +
                                                           str(mem.mem_id) +
                                                           str(i+1), parent=self)
                    tension_con.fea_required = True
                    shear_con = NonLinearConstraint(con_fun=shear,
                                                    name="Shear " + str(
                                                        mem.mem_id) + str(i+1),
                                                    parent=self)
                    shear_con.fea_required = True
                    bending_moment_con = NonLinearConstraint(
                        con_fun=bending_moment,
                        name="Bending_moment " +
                             str(mem.mem_id) +
                             str(i+1), parent=self)
                    bending_moment_con.fea_required = True

                    self.cons.extend([compression_con, tension_con, shear_con,
                                      bending_moment_con])

if __name__ == "__main__":
    from src.optimization.solvers import *
    problem = WIColumn()
    solver = SLP(step_length=10)
    solver.solve(problem, maxiter=50)
    problem(solver.X)



"""
@author: Victoria
"""

from src.frame2d.frame2d import *
from src.optimization.structopt import *
from src.optimization.problems.wi_column import WIColumn
from src.optimization.solvers import *
from src.optimization.result_exporter import *


class ColumnCalculation(WIColumn):

    q_wind = 0.00065  # N/mm2
    q_snow = 0.002  # N/mm2
    g_truss = 1  # N/mm
    g_roof = 0.0005  # N/mm2
    L_list = [24000, 30000, 36000]  # mm
    Lpi_list = [6000, 8000, 10000]  # mm
    c = 6000  # mm
    lcr_list = [2, 0.7]
    buckling_z = [True, False]
    LT_buckling = [True, False]
    cross_section_class_list = [2, 3]
    sym = "dual"

    Mz = 0
    Qy = 0
    Qx = q_wind * c

    #  x0 = [800, 50, 500, 50, 500, 50]
    #  x0 = [400, 20, 200, 20]
    x0 = [600, 30, 300, 30]
    #  x0 = [700, 40, 400, 40]
    #  x0 = [800, 50, 500, 50]

    # for LT in LT_buckling:
    for buck_z in buckling_z:
        for cross_section_class in cross_section_class_list:
            for L in L_list:
                Fy = -(1.15 * (g_truss + g_roof * c) + (
                        1.5 * q_snow * c)) * L / 2

                for Lpi in Lpi_list:
                    # phi_0 = 1 / 200
                    # m = 2
                    # alpha_m = math.sqrt(0.5 * (1 + 1 / m))
                    # alpha_h = max(2 / 3, min(2 / math.sqrt(Lpi), 1))
                    # phi = phi_0 * alpha_m * alpha_h
                    # Fx = phi * Fy

                    phi = en1993_1_1.sway_imperfection(Lpi, m=2)

                    Fx = -phi * Fy

                    for lcr in lcr_list:
                        problem = WIColumn(
                            Lpi, Fx, Fy, Qx, Qy, Mz, lcr=lcr,
                            top_flange_class=cross_section_class,
                            bottom_flange_class=cross_section_class,
                            web_class=cross_section_class, symmetry=sym,
                            buckling_z=buck_z, LT_buckling=False)

                        #  print(problem.nnonlincons())

                        # print("L={0}, Lpi={1}, Fx={2}, Fy={3}, Qx={4}, Qy={5},"
                        #       " Mz={6}, lcr={7}, cross_section_class={8}, "
                        #       "buckling_z={9}"
                        #       .format(L, Lpi, Fx, Fy, Qx, Qy, Mz, lcr,
                        #               cross_section_class, buck_z))
                        solver = slsqp.SLSQP()
                        f_best, x_best = solver.solve(problem, maxiter=100,
                                                      x0=x0)
                        #  print(x_best)
                        problem(solver.best_x, prec=5)

                        # solver = SLP(move_limits=[0.9, 6])
                        # solver.solve(problem, maxiter=500, maxtime=40, x0=x0)
                        # problem(solver.X, prec=5)

                        ResultExporter(problem, solver).to_csv()

                        breakpoint()

                        #  wi.cross_section.draw()



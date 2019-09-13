"""
@author: Victoria
"""

from src.frame2d.frame2d import *
from src.optimization.structopt import *
from src.optimization.problems.wi_column import WIColumn
from src.optimization.solvers import *


class ColumnCalculation(WIColumn):

    q_wind = 0.00065  # N/mm2
    q_snow = 0.002  # N/mm2
    g_truss = 1  # N/mm
    g_roof = 0.0005  # N/mm2
    L_list = [24000, 30000, 36000]  # mm
    Lpi_list = [6000, 8000, 10000]  # mm
    c = 6000  # mm
    lcr_list = [2, 2]
    buckling_z = [True, False]
    LT_buckling = [True, False]
    cross_section_class_list = [2, 3]

    Qy = 0
    Qx = q_wind * c

    for LT in LT_buckling:
        for buck_z in buckling_z:
            for cross_section_class in cross_section_class_list:
                for L in L_list:
                    Fy = (1.15 * (g_truss + g_roof * c)
                          + (1.5 * q_snow * c)) * L / 2

                    for Lpi in Lpi_list:
                        phi_0 = 1 / 200
                        m = 2
                        alpha_m = math.sqrt(0.5 * (1 + 1 / m))
                        alpha_h = max(2 / 3, min(2 / math.sqrt(Lpi), 1))
                        phi = phi_0 * alpha_m * alpha_h
                        Fx = phi * Fy

                        for lcr in lcr_list:
                            problem = WIColumn(
                                Fy, Fx, Qx, lcr, LT_buckling=LT,
                                buckling_z=buck_z,
                                top_flange_class=cross_section_class,
                                bottom_flange_class=cross_section_class,
                                web_class=cross_section_class)
                            x0 = [500, 20, 300, 20]
                            print("Debug")

                            solver = SLP(step_length=3)
                            solver.solve(problem, maxiter=10, maxtime=10,
                                         x0=x0)
                            #  problem(solver.X)
                            #  wi.cross_section.draw()



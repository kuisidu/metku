"""
@author: Victoria
"""

import time

try:
    from src.frame2d.frame2d import *
    from src.optimization.solvers.optsolver import OptSolver
    # from src.optimization.structopt import *
    import src.optimization.structopt as sopt
    from src.optimization.problems.wi_column import WIColumn
    from src.optimization.solvers import slsqp, slp, slqp, mislp, two_phase
    from src.optimization.solvers.trust_region import TrustRegionConstr
    from src.optimization.solvers.bnb import BnB
    from src.optimization.solvers.lp import LP
    from src.optimization.result_exporter import *
    from src.eurocodes.en1993 import en1993_1_1
except:
    from frame2d.frame2d import *
    from optimization.solvers.optsolver import OptSolver
    from optimization.structopt import *
    import optimization.structopt as sopt
    from optimization.problems.wi_column import WIColumn
    from optimization.solvers import slsqp, slp, slqp, mislp, two_phase
    from optimization.solvers.trust_region import TrustRegionConstr
    from optimization.solvers.bnb import BnB
    from optimization.solvers.lp import LP
    from optimization.result_exporter import *
    from eurocodes.en1993 import en1993_1_1
    from copy import deepcopy


class ColumnCalculation(WIColumn):
    # def ColumnCalculation():

    q_wind = 0.00065  # N/mm2
    q_snow = 0.002  # N/mm2
    g_truss = 1  # N/mm
    g_roof = 0.0005  # N/mm2

    L_list = [24000, 30000, 36000]  # mm
    # L_list = [30000]  # mm
    Lpi_list = [6000, 8000, 10000]  # mm
    # Lpi_list = [6000]  # mm

    c = 6000  # mm
    lcr_list = [2, 0.7]
    buckling_z = [True, False]
    LT_buckling = [True, False]
    cross_section_class_list = [3, 2]
    sym = "dual"  # dual or mono
    prob_type = "discrete"  # continuous or discrete

    Mz = 0
    Qy = 0
    Qx = 1.5 * q_wind * c

    #  x0 = [800, 50, 500, 50, 500, 50]
    #  x0 = [250, 8, 150, 10]
    #  x0 = [340, 6, 180, 10]
    x0 = [300, 8, 200, 10]

    #  x0 = [400, 20, 200, 20]
    #  x0 = [300, 8, 200, 10, 200, 10]
    #  x0 = [var.ub for var in problem.vars]

    for lcr in lcr_list:
        for cross_section_class in cross_section_class_list:
            for Lpi in Lpi_list:
                for L in L_list:
                    for LT in LT_buckling:
                        for buck_z in buckling_z:
                            Fy = -(1.15 * (g_truss + g_roof * c) + (
                                    1.5 * q_snow * c)) * L / 2

                            # phi_0 = 1 / 200
                            # m = 2
                            # alpha_m = math.sqrt(0.5 * (1 + 1 / m))
                            # alpha_h = max(2 / 3, min(2 / math.sqrt(Lpi), 1))
                            # phi = phi_0 * alpha_m * alpha_h

                            phi = en1993_1_1.sway_imperfection(Lpi, m=2)

                            Fx = phi * (-Fy)

                            problem = WIColumn(
                                L, Lpi, Fx, Fy, Qx, Qy, Mz, lcr=lcr,
                                top_flange_class=cross_section_class,
                                bottom_flange_class=cross_section_class,
                                web_class=cross_section_class, symmetry=sym,
                                buckling_z=buck_z, LT_buckling=LT,
                                prob_type=prob_type)

                            print("Constraints:", len(problem.cons))
                            print("Nonlinear Constraints:", problem.nnonlincons())
                            
                            # problem([340, 10, 150, 8])

                            # problem = WIColumn(
                            #     L=24000, Lpi=6000, Fx=782.89, Fy=-271200,
                            #     Qx=5.85, Qy=0, Mz=0, lcr=2,
                            #     top_flange_class=3,
                            #     bottom_flange_class=3,
                            #     web_class=3,
                            #     symmetry="dual",
                            #     buckling_z=True,
                            #     LT_buckling=True,
                            #     prob_type='continuous')

                            # print("L={0}, Lpi={1}, Fx={2}, Fy={3}, Qx={4}, "
                            #       "Qy={5}, Mz={6}, lcr={7}, "
                            #       "cross_section_class={8}, sym={9}, "
                            #       "buckling_z={10}, LT={11}, prob_type={12}"
                            #       .format(L, Lpi, Fx, Fy, Qx, Qy, Mz, lcr,
                            #               cross_section_class, sym, buck_z,
                            #               LT, prob_type))

                            #vnew = deepcopy(problem.vars[0])
                            #solver = TrustRegionConstr()
                            #f_best, x_best, nit = solver.solve(
                            #    problem, maxiter=200, x0=x0)

                            #problem.num_iters = nit
                            #problem(x_best, prec=5)

                            # BnB
                            lb_solver = TrustRegionConstr()
                            # lb_solver = slp.SLP()
                            solver = BnB(problem, lb_solver)

                            # return solver, problem

                            # break

                            solver.solve(problem, x0=x0, verb=2)

                            x0 = solver.best_x
                            print("x0=", x0)
                            print("Solver done.")
                            #print(solver.X)
                            print("Best found solution.",solver.best_x)
                            print("Best found objective function value:",solver.best_f)
                                                        
                            # problem(solver.best_x)

                            #return solver, problem

                            # TrustRegionConstr
                            # solver = TrustRegionConstr()
                            # f_best, x_best, nit = solver.solve(
                            #     problem,
                            #     maxiter=200,
                            #     x0=x0)
                            # problem.num_iters = nit
                            # problem(x_best, prec=5)

                            # SLP
                            # solver = slp.SLP(move_limits=[0.9, 6])
                            # solver.solve(problem,
                            #              maxiter=500,
                            #              maxtime=40,
                            #              x0=x0)
                            # problem(solver.X, prec=5)

                            # SLSQP
                            # solver = slsqp.SLSQP()
                            # f_best, x_best = solver.solve(problem,
                            #                               maxiter=100,
                            #                               x0=x0)
                            # problem(solver.best_x, prec=5)

                            # MISLP
                            # solver = mislp.MISLP(move_limits=[0.5, 5])
                            # # problem(x0)
                            # solver.solve(problem,
                            #              maxiter=200,
                            #              x0=x0,
                            #              min_diff=1e-2,
                            #              verb=True)
                            # problem(solver.X, prec=5)

                            # 2-vaihetekniikalla
                            # solver1 = SLP(move_limits=[0.9, 4])
                            # solver1 = TrustRegionConstr()
                            # solver2 = mislp.MISLP(move_limits=[0.85, 1.5])
                            # solver = two_phase.TwoPhase(
                            #     solver1, solver2, limits=[3, 3])
                            # fopt, xopt = solver.solve(problem,
                            #                           x0=x0,
                            #                           maxiter=200,
                            #                           # min_diff=1e-6,
                            #                           # verb=True
                            #                           )
                            # problem(xopt)
                            # ResultExporter(problem, solver2).to_csv()
                            # ResultExporter(problem, solver2).csv_to_excel()

                            # # ResultExporter muille kuin 2-vaihetekniikalle
                            ResultExporter(problem, solver).to_csv()
                            ResultExporter(problem, solver).csv_to_excel()

                            seconds = time.process_time()
                            m, s = divmod(seconds, 60)
                            h, m = divmod(m, 60)

                            print("")
                            print("Process time:", "%d:%02d:%02d" % (h, m, s))
                            print("----------------------------")
                            print("")

                            # breakpoint()

                            #  wi.cross_section.draw()


# if __name__ == '__main__':
#
#     solver, p = ColumnCalculation()
#
#     p.vars[1].lock(6)

"""
@author: Victoria
"""

from frame2d.frame2d import *
from optimization.structopt import *
from optimization.problems.wi_column import WIColumn
from optimization.solvers import slsqp, slp, slqp
from optimization.solvers.trust_region import TrustRegionConstr
#from optimization.result_exporter import *
from eurocodes.en1993 import en1993_1_1

q_wind = 0.00065  # N/mm2    
q_snow = 0.002  # N/mm2    
g_truss = 1  # N/mm    
g_roof = 0.0005  # N/mm2    
L_list = [24000]  # mm    
Lpi_list = [8000]  # mm    
c = 6000  # mm    
lcr_list = [0.7]    
buckling_z = [True]    
LT_buckling = [True]    
cross_section_class_list = [2]
sym = "dual"  # dual or mono
prob_type = "discrete"  # continuous or discrete
    
Mz = 0    
Qy = 0    
Qx = q_wind * c    
x0 = [300, 8, 200, 10]
    
# for LT in LT_buckling:
    
for buck_z in buckling_z:
    for cross_section_class in cross_section_class_list:
        for L in L_list:
            Fy = -(1.15 * (g_truss + g_roof * c) + (1.5 * q_snow * c)) * L / 2
                   
            for Lpi in Lpi_list:
                    
                phi = en1993_1_1.sway_imperfection(Lpi,m=2)
                    
                Fx = phi * Fy

                for lcr in lcr_list:
                    problem = WIColumn(
                            Lpi, Fx, Fy, Qx, Qy, Mz, lcr=lcr,
                            top_flange_class=cross_section_class,
                            bottom_flange_class=cross_section_class,
                            web_class=cross_section_class, symmetry="dual",
                            buckling_z=buck_z, LT_buckling=False)
                    print(problem.nnonlincons())
                    #breakpoint()

                    print("L={0}, Lpi={1}, Fx={2}, Fy={3}, Qx={4}, Qy={5},"
                          "Mz={6}, lcr={7}, cross_section_class={8}, "
                          "buckling_z={9}"
                          .format(L, Lpi, Fx, Fy, Qx, Qy, Mz, lcr,
                                  cross_section_class, buck_z))

                    
                    #solver = slsqp.SLSQP()
                    solver = TrustRegionConstr()
                    f_best, x_best, nit = solver.solve(problem, maxiter=200, x0=x0)
                    print(x_best)
                    #solver = slp.SLP(move_limits=[0.95, 1.05])
                    #solver.solve(problem, maxiter=100, maxtime=400, x0=x0,verb=True)
                    #solver = slqp.SLQP()
                    #solver.solve(problem, maxiter=100, maxtime=400, x0=x0,min_diff=1e-3,verb=False)
                    problem(x_best,prec=5)
                    #solver.solve(problem, maxiter=500, maxtime=40, x0=x0)
                    #r0 = problem.solve("slsqp",x0=x)                    
                    #problem(r0.x, prec=5)
                    
                    #ResultExporter(problem, solver).to_csv()
                    #breakpoint()




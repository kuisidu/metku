# Import dependencies
import numpy as np

MIN_HEIGHT = 100
MAX_HEIGHT = 500
step_height = 10
MIN_WIDTH = 100
MAX_WIDTH = 500
step_width = 10

HEIGHTS = list(np.arange(MIN_HEIGHT, MAX_HEIGHT+step_height, step_height))
WIDTHS = list(np.arange(MIN_WIDTH, MAX_WIDTH+step_width, step_width))
#THICKNESSES = [5, 6, 8, 10, 12, 14, 15, 16, 18, 20, 22, 25, 30, 35, 40, 50]
THICKNESSES = [4, 5, 6, 8, 10, 12, 14, 15, 16, 18, 20]
MAX_THICK = max(THICKNESSES)
MIN_THICK = min(THICKNESSES)

try:
    from src.frame2d.frame2d import *
    from src.truss2d import Truss2D
    from src.optimization.solvers import *
    from src.optimization.problems.structural_problem import StructuralProblem
except:
    from frame2d.frame2d import *
    from truss2d import Truss2D
    from optimization.solvers import *


def create_structure(L, H0, H1, H2, dx, n):

    simple_frame = [1, 1, H0, L]
    frame = Frame2D(simple=simple_frame, create_beams=False,
                    supports='fixed')  # Tuet, aina fixed?
    for col in frame.columns:
        col.profile = 'WI 500-12-10X300-10X300'

    simple_truss = {
        "L1": L/2,
        "H0": H0,
        "H1": H1,
        "H2": H2,
        "dx": dx,
        "n": n
    }

    truss = Truss2D(simple=simple_truss, fem_model=frame.f)
    frame.add(truss)

    columns = frame.columns

    # Kuormat
    for tc in truss.top_chords:
        frame.add(LineLoad(tc, [-22.6, -22.6], 'y'))
    for bc in truss.bottom_chords:
        frame.add(LineLoad(bc, [-0.7, -0.7], 'y'))
    # frame.add(PointLoad([0, 6500], [100e3, 0, 0]))
    frame.add(LineLoad(columns[0], [5.85, 5.85], 'x'))
    frame.add(LineLoad(columns[1], [0.2, 0.2], 'x'))
    # frame.add(LineLoad(columns), [])

    frame.generate()
    frame.f.draw()
    frame.plot()
    frame.calculate()

    return frame


def create_continuous_variable_groups(structure, col_bounds, col_values,
                                      tc_values, tc_bounds, bc_values,
                                      bc_bounds, web_values, web_bounds):

    # HUOM! objectit oltava poikkileikkausolioita
    # mem.cross_section

    groups = []
    col_sections = [col.cross_section for col in structure.columns]
    
    COL_group = {
        'name': 'Columns',
        'var_type': 'continuous',
        'values': col_values, # [300, 100, 10, 5],
        'bounds': col_bounds, # [[100, 800], [100, 300], [5, 50], [5, 50]],
        'properties': ['h', 'b', 'tf', 'tw'],
        'objects': col_sections
    }
    groups.append(COL_group)
    
    # Ristikon osille samaan tapaan
    truss = structure.truss[0]
    # top_chords = truss.top_chords

    # H1_group = {
    #     'name': 'H1',
    #     'var_type': 'continuous',
    #     'value': 1500,
    #     'lb': 500,
    #     'ub': 2000,
    #     'property': 'H1',
    #     'objects': [truss]
    # }
    # groups.append(H1_group)

    TC_group = {
        'name': 'TopChords',
        'var_type': 'continuous',
        'values': tc_values,  # [200, 200, 10],
        'bounds': tc_bounds,  # [[100, 300], [100, 300], [4, 12.5]],
        'properties': [['H', 'B'], 'T'],
        'objects': [tc.cross_section for tc in truss.top_chords],
        }
    groups.append(TC_group)

    BC_group = {
        'name': 'BottomChords',
        'var_type': 'continuous',
        'values': bc_values,  # [150, 150, 8],
        'bounds': bc_bounds,  # [[80, 200], [80, 200], [4, 10]],
        'properties': [['H', 'B'], 'T'],
        'objects': [bc.cross_section for bc in truss.bottom_chords]
    }
    groups.append(BC_group)

    WEB_group = {
        'name': 'Webs',
        'var_type': 'continuous',
        'values': web_values,  # [100, 100, 5],
        'bounds': web_bounds,  # [[50, 150], [50, 150], [3, 8]],
        'properties': [['H', 'B'], 'T'],
        'objects': [web.cross_section for web in truss.webs.values()]
    }
    groups.append(WEB_group)
    
    return groups


def create_discrete_variable_groups(structure):

    groups = []
    col_sections = [col.cross_section for col in structure.columns]

    COL_group_h = {
        'name': 'Columns h',
        'var_type': 'index',
        'value': 15,
        'values': HEIGHTS,
        'property': 'h',
        'objects': col_sections
    }
    groups.append(COL_group_h)

    COL_group_tw = {
        'name': 'Columns tw',
        'var_type': 'index',
        'value': 6,
        'values': THICKNESSES,
        'property': 'tw',
        'objects': col_sections
    }
    groups.append(COL_group_tw)

    COL_group_b = {
        'name': 'Columns b',
        'var_type': 'index',
        'value': 10,
        'values': WIDTHS,
        'property': 'b',
        'objects': col_sections
    }
    groups.append(COL_group_b)

    COL_group_tf = {
        'name': 'Columns tf',
        'var_type': 'index',
        'value': 6,
        'values': THICKNESSES,
        'property': 'tf',
        'objects': col_sections
    }
    groups.append(COL_group_tf)

    truss = structure.truss[0]

    # H1_group = {
    #     'name': 'H1',
    #     'var_type': 'continuous',
    #     'value': 1500,
    #     'lb': 500,
    #     'ub': 2000,
    #     'property': 'H1',
    #     'objects': [truss]
    # }
    # groups.append(H1_group)

    TC_group = {
        'name': 'TopChords',
        'var_type': 'index',
        'value': 10,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': truss.top_chords,
        }
    groups.append(TC_group)

    BC_group = {
        'name': 'BottomChords',
        'var_type': 'index',
        'value': 15,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': truss.bottom_chords
    }
    groups.append(BC_group)

    WEB_group = {
        'name': 'Webs',
        'var_type': 'index',
        'value': 5,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': list(truss.webs.values())
    }
    groups.append(WEB_group)

    return groups


def create_binary_discrete_variable_groups(structure):

    groups = []
    col_sections = [col.cross_section for col in structure.columns]

    COL_group_h = {
        'name': 'Columns h',
        'var_type': 'discrete',
        'value': 300,
        'values': HEIGHTS,
        'property': 'h',
        'objects': col_sections
    }
    groups.append(COL_group_h)

    COL_group_tw = {
        'name': 'Columns tw',
        'var_type': 'discrete',
        'value': 5,
        'values': THICKNESSES,
        'property': 'tw',
        'objects': col_sections
    }
    groups.append(COL_group_tw)

    COL_group_b = {
        'name': 'Columns b',
        'var_type': 'discrete',
        'value': 150,
        'values': WIDTHS,
        'property': 'b',
        'objects': col_sections
    }
    groups.append(COL_group_b)

    COL_group_tf = {
        'name': 'Columns tf',
        'var_type': 'discrete',
        'value': 6,
        'values': THICKNESSES,
        'property': 'tf',
        'objects': col_sections
    }
    groups.append(COL_group_tf)

    truss = structure.truss[0]

    # H1_group = {
    #     'name': 'H1',
    #     'var_type': 'continuous',
    #     'value': 1500,
    #     'lb': 500,
    #     'ub': 2000,
    #     'property': 'H1',
    #     'objects': [truss]
    # }
    # groups.append(H1_group)

    TC_h_group = {
        'name': 'TopChords',
        'var_type': 'discrete',
        'value': 180,
        'values': list(np.arange(100, 310, 10)),
        'property': ['H', 'B'],
        'objects': [tc.cross_section for tc in truss.top_chords],
        }
    groups.append(TC_h_group)
    TC_t_group = {
        'name': 'TopChords',
        'var_type': 'discrete',
        'value': 8,
        'values': list(np.arange(3, 13, 1)),
        'property': ['T'],
        'objects': [tc.cross_section for tc in truss.top_chords],
    }
    groups.append(TC_t_group)

    BC_h_group = {
        'name': 'BottomChords',
        'var_type': 'discrete',
        'value': 140,
        'values': list(np.arange(100, 310, 10)),
        'property': ['H', 'B'],
        'objects': [bc.cross_section for bc in truss.bottom_chords]
    }
    groups.append(BC_h_group)
    BC_t_group = {
        'name': 'TopChords',
        'var_type': 'discrete',
        'value': 8,
        'values': list(np.arange(3, 13, 1)),
        'property': ['T'],
        'objects': [bc.cross_section for bc in truss.bottom_chords],
    }
    groups.append(BC_t_group)

    WEB_h_group = {
        'name': 'Webs',
        'var_type': 'discrete',
        'value': 120,
        'values': list(np.arange(80, 200, 10)),
        'property': ['H', 'B'],
        'objects': [web.cross_section for web in truss.webs.values()]
    }
    groups.append(WEB_h_group)

    WEB_t_group = {
        'name': 'Webs',
        'var_type': 'discrete',
        'value': 5,
        'values': list(np.arange(3, 13, 1)),
        'property': ['T'],
        'objects': [web.cross_section for web in truss.webs.values()]
    }
    groups.append(WEB_t_group)

    return groups


def create_constraint_groups():
    pass
    

if __name__ == '__main__':

    structure = create_structure(L=24000,
                                 H0=8000,
                                 H1=1800,
                                 H2=2400,
                                 dx=0,
                                 n=14
                                 )

    col_bounds = [[100, 800], [100, 300], [5, 50], [5, 50]]
    col_values = [300, 100, 10, 5]
    tc_values = [200, 10]
    tc_bounds = [[100, 300], [4, 12.5]]
    bc_values = [150, 8]
    bc_bounds = [[80, 200], [4, 10]]
    web_values = [100, 5]
    web_bounds = [[50, 150], [3, 8]]

    cont_var_groups = create_continuous_variable_groups(structure=structure,
                                                        col_bounds=col_bounds,
                                                        col_values=col_values,
                                                        tc_bounds=tc_bounds,
                                                        tc_values=tc_values,
                                                        bc_bounds=bc_bounds,
                                                        bc_values=bc_values,
                                                        web_bounds=web_bounds,
                                                        web_values=web_values,
                                                        )

    disc_var_groups = create_discrete_variable_groups(structure=structure)

    binary_disc_var_groups = create_binary_discrete_variable_groups(
        structure=structure)

    problem = StructuralProblem(name="Example",
                                structure=structure,
                                var_groups=binary_disc_var_groups,
                                constraints={
                                    'buckling_y': True,
                                    'buckling_z': True,
                                    'compression_bending_y': True,
                                    'compression_bending_z': True,
                                    'compression': True,
                                    'tension': True,
                                    'shear': True,
                                    # 'deflection_y': structure.L / 200,
                                    #  'deflection_x': structure.H / 300,
                                    'web_class': True,
                                    'flange_class': True
                                })

    # GA
    # solver = GA(pop_size=50, mut_rate=0.15, cx_rate=0.9)
    # x0 = [var.value for var in problem.vars]
    # # problem(x0)
    # fopt, xopt = solver.solve(problem, x0=x0, maxiter=10, plot=False, verb=True)
    # problem(xopt, ncons=10)
    # structure.plot_deflection(100)
    # structure.to_robot("lopputulos")

    # VNS
    # solver = VNS(pop_size=10, maxiter=5, first_improvement=True, solver="GA",
    #              solver_params={'pop_size': 10,
    #                             'mut_rate': 0.15,
    #                             'cx_rate': 0.9,
    #                             'first_improvement': True})
    # x0 = [var.value for var in problem.vars]
    # solver.solve(problem, x0=x0, maxiter=10, verb=True)

    # SLP
    # solver = SLP(move_limits=[0.1, 0.1], beta=100)
    # x0 = [var.value for var in problem.vars]
    # solver.solve(problem,
    #              maxiter=50000,
    #              maxtime=30,
    #              x0=x0,
    #              verb=True)
    # problem(solver.X, prec=5)

    # MISLP
    solver = MISLP(move_limits=[0.1, 0.1], beta=100)
    # problem(x0)
    x0 = [var.value for var in problem.vars]
    solver.solve(problem,
                 maxiter=100,
                 x0=x0,
                 min_diff=1e-2,
                 verb=True)
    problem(solver.X, prec=5)






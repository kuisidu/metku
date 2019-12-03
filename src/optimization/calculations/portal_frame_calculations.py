# Import dependencies
import numpy as np

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
        col.profile = "HE 300 A"

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
    frame.to_robot("testi")

    return frame


def create_continuous_variable_groups(structure, col_bounds, col_values,
                                      tc_values, tc_bounds, bc_values,
                                      bc_bounds, web_values, web_bounds):

    # HUOM! objectit oltava poikkileikkausolioita
    # mem.cross_section

    groups = []
    
    COL_group = {
        'name': 'Columns',
        'var_type': 'continuous',
        'values': col_values, # [300, 100, 10, 5],
        'bounds': col_bounds, # [[100, 800], [100, 300], [5, 50], [5, 50]],
        'properties': ['h', 'b', 'tf', 'tw'],
        'objects': structure.columns
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
        'properties': ['H', 'B', 'T'],
        'objects': truss.top_chords,
        'constraints': {
            'buckling_y': True,
            'buckling_z': False,
            'deflection_y': truss.L / 300,
            'deflection_x': truss.H / 300
        }}
    groups.append(TC_group)

    BC_group = {
        'name': 'BottomChords',
        'var_type': 'continuous',
        'values': bc_values,  # [150, 150, 8],
        'bounds': bc_bounds,  # [[80, 200], [80, 200], [4, 10]],
        'properties': ['H', 'B', 'T'],
        'objects': truss.bottom_chords
    }
    groups.append(BC_group)

    WEB_group = {
        'name': 'Webs',
        'var_type': 'continuous',
        'values': web_values,  # [100, 100, 5],
        'bounds': web_bounds,  # [[50, 150], [50, 150], [3, 8]],
        'properties': ['H', 'B', 'T'],
        'objects': list(truss.webs.values())
    }
    groups.append(WEB_group)
    
    return groups


def create_discrete_variable_groups(structure):

    groups = []

    COL_group_h = {
        'name': 'Columns h',
        'var_type': 'index',
        'value': [300],
        'values': [MIN_HEIGHT, MAX_HEIGHT],
        'property': 'h',
        'objects': structure.columns
    }
    groups.append(COL_group_h)

    COL_group_tw = {
        'name': 'Columns tw',
        'var_type': 'index',
        'value': [6],
        'values': [MIN_THICK, MAX_THICK],
        'property': 'tw',
        'objects': structure.columns
    }
    groups.append(COL_group_tw)

    COL_group_b = {
        'name': 'Columns b',
        'var_type': 'index',
        'value': [200],
        'values': [MIN_WIDTH, MAX_WIDTH],
        'property': 'b',
        'objects': structure.columns
    }
    groups.append(COL_group_b)

    COL_group_tf = {
        'name': 'Columns tf',
        'var_type': 'index',
        'value': [10],
        'values': [MIN_THICK, MAX_THICK],
        'property': 'tf',
        'objects': structure.columns
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
        'constraints': {
            'buckling_y': True,
            'buckling_z': False,
            'deflection_y': structure.L / 300,
            'deflection_x': structure.H / 300
        }}
    groups.append(TC_group)

    BC_group = {
        'name': 'BottomChords',
        'var_type': 'index',
        'value': 10,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': truss.bottom_chords
    }
    groups.append(BC_group)

    WEB_group = {
        'name': 'Webs',
        'var_type': 'index',
        'value': 10,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': list(truss.webs.values())
    }
    groups.append(WEB_group)

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
    tc_values = [200, 200, 10]
    tc_bounds = [[100, 300], [100, 300], [4, 12.5]]
    bc_values = [150, 150, 8]
    bc_bounds = [[80, 200], [80, 200], [4, 10]]
    web_values = [100, 100, 5]
    web_bounds = [[50, 150], [50, 150], [3, 8]]

    # var_groups = create_continuous_variable_groups(structure=structure,
    #                                                col_bounds=col_bounds,
    #                                                col_values=col_values,
    #                                                tc_bounds=tc_bounds,
    #                                                tc_values=tc_values,
    #                                                bc_bounds=bc_bounds,
    #                                                bc_values=bc_values,
    #                                                web_bounds=web_bounds,
    #                                                web_values=web_values,
    #                                                )

    var_groups = create_discrete_variable_groups(structure=structure)

    problem = StructuralProblem(name="Example",
                                structure=structure,
                                var_groups=var_groups,
                                constraints={
                                    'buckling_y': True,
                                    'buckling_z': True,
                                    'compression_bending_y': True,
                                    'compression_bending_z': False,
                                    'compression': True,
                                    'tension': True,
                                    'shear': True,
                                    # 'deflection_y': structure.L / 200,
                                    # 'deflection_x': structure.H / 300
                                })
    solver = GA(pop_size=50, mut_rate=0.15)
    x0 = [1 for var in problem.vars]
    fopt, xopt = solver.solve(problem, x0=x0, maxiter=10, plot=True)
    problem(xopt, ncons=10)
    structure.plot_deflection(100)

    # 2-vaihetekniikalla
    # solver1 = SLP(move_limits=[0.9, 4])
    # solver1 = TrustRegionConstr()
    # solver2 = GA(pop_size=50, mut_rate=0.15)
    # solver = two_phase.TwoPhase(
    #     solver1, solver2, limits=[3, 3])
    # fopt, xopt = solver.solve(problem,
    #                           x0=x0,
    #                           maxiter=200,
    #                           # min_diff=1e-6,
    #                           # verb=True
    #                           )
    # problem(xopt)



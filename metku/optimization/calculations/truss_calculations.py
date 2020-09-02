# Import dependencies
import numpy as np

# I-profiilin mittojen ylä- ja alarajat sekä askelpituus
MIN_HEIGHT = 100
MAX_HEIGHT = 500
step_height = 10
MIN_WIDTH = 100
MAX_WIDTH = 500
step_width = 10

# I-profiilin mittojen listat
HEIGHTS = list(np.arange(MIN_HEIGHT, MAX_HEIGHT + step_height, step_height))
WIDTHS = list(np.arange(MIN_WIDTH, MAX_WIDTH + step_width, step_width))
# THICKNESSES = [5, 6, 8, 10, 12, 14, 15, 16, 18, 20, 22, 25, 30, 35, 40, 50]
THICKNESSES = [4, 5, 6, 8, 10, 12, 14, 15, 16, 18, 20]
MAX_THICK = max(THICKNESSES)
MIN_THICK = min(THICKNESSES)

# Vapaavälin ylä- ja alaraja sekä askelpituus ja lista mahdollisista
# vapaaväleistä
min_g = 20
max_g = 90
step_g = 1
g_values = list(np.arange(min_g, max_g, step_g))

try:
    from metku.frame2d.frame2d import *
    from metku.truss2d import Truss2D
    from metku.optimization.solvers import *
    from metku.optimization.solvers.bnb import BnB
    from metku.optimization.problems.structural_problem import StructuralProblem
    from metku.sections.steel.steel_section import SteelSection
    from metku.materials.steel_data import Steel
except:
    from frame2d.frame2d import *
    from truss2d import Truss2D
    from optimization.solvers import *
    from optimization.solvers.bnb import BnB
    from optimization.problems.structural_problem import StructuralProblem
    from materials.steel_data import Steel


def create_structure(L, H0, H1, H2, dx, n):
    """
    Luodaan rakenne, eli pilarit ja ristikko
    :param L: ristikon jänneväli
    :param H0: pilarin vapaa korkeus
    :param H1: ristikon korkeus reunalla
    :param H2: ristikon korkeus harjalla
    :param dx: epäkeskisyys pilarille
    :param n: uumasauvojen lukumäärä
    :return: frame, eli kehän rakenne
    """

    # Luodaan ristikko
    simple_truss = {
        "L1": L / 2,
        "H0": H0,
        "H1": H1,
        "H2": H2,
        "dx": dx,
        "n": n
    }

    truss = Truss2D(simple=simple_truss)

    # Määritetään ristikon sauvojen liitosten paikat
    dloc = 1 / (n / 2)

    tc = truss.top_chords[0]
    for i, joint in enumerate(tc.joints):
        if i:
            joint.loc = 2 * dloc * i + dloc
        else:
            joint.loc = dloc

    tc2 = truss.top_chords[1]
    joints = sorted(tc2.joints, key=lambda j: j.loc)
    for i, joint in enumerate(joints):
        if i:
            joint.loc = 2 * dloc * i

    # Määritetään ristikon osille materiaali, tällä vaihtuu kaikkien materiaali
    for mem in truss.members.values():
        mem.material = "S420"

    # print(frame.truss[0].webs[2].j1.e, frame.truss[0].webs[2].j1.g1)
    # frame.truss[0].webs[2].profile = "SHS 90x5"
    # print(frame.truss[0].webs[2].j1.e, frame.truss[0].webs[2].j1.g1)

    # Lisätään kehälle kuorma
    for tc in truss.top_chords:
        truss.add(LineLoad(tc, [-21.45, -21.45], 'y'))
    for bc in truss.bottom_chords:
        truss.add(LineLoad(bc, [-0.69, -0.69], 'y'))


    # Generoidaan kehä ja lasketaan voimasuureet
    truss.generate()
    # frame.f.draw()
    # frame.plot()
    truss.calculate()

    return truss


def create_continuous_variable_groups(structure,
                                      tc_values, tc_bounds, bc_values,
                                      bc_bounds, web_values, web_bounds):
    # HUOM! objectit oltava poikkileikkausolioita
    # mem.cross_section
    """
    Luodaan jatkuvien muuttujien joukko kehälle, muuttujina
    :param structure:
    :param col_bounds: pilarin muuttujien rajat
    :param col_values: pilarin muuttujien arvot
    :param tc_values: yläpaarteen muuttujien arvot
    :param tc_bounds: yläpaarteen muuttujien rajat
    :param bc_values: alapaarteen muuttujien arvot
    :param bc_bounds: alapaarteen muuttujien rajat
    :param web_values: uumasauvojen muuttujien arvot
    :param web_bounds: uumasauvojen muuttujien rajat
    :return: jatkuvat muuttujat
    """

    groups = []

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

    # Yläpaarteen muuttujat
    TC_group = {
        'name': 'TopChords',
        'var_type': 'continuous',
        'values': tc_values,  # [200, 200, 10],
        'bounds': tc_bounds,  # [[100, 300], [100, 300], [4, 12.5]],
        'properties': [['H', 'B'], 'T'],
        'objects': [tc.cross_section for tc in truss.top_chords],
    }
    groups.append(TC_group)

    # Alapaarteen muuttujat
    BC_group = {
        'name': 'BottomChords',
        'var_type': 'continuous',
        'values': bc_values,  # [150, 150, 8],
        'bounds': bc_bounds,  # [[80, 200], [80, 200], [4, 10]],
        'properties': [['H', 'B'], 'T'],
        'objects': [bc.cross_section for bc in truss.bottom_chords]
    }
    groups.append(BC_group)

    # Uumasauvojen muuttujat
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


def create_index_variable_groups(structure):
    """
    Luodaan kehän indeksimuuttujat
    :param structure:
    :return: indeksimuuttujat
    """
    groups = []

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

    # Yläpaarteen profiili valikoimasta
    TC_group = {
        'name': 'TopChords',
        'var_type': 'index',
        'value': 64,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': truss.top_chords,
    }
    groups.append(TC_group)

    # Alapaarteen profiili valikoimasta
    BC_group = {
        'name': 'BottomChords',
        'var_type': 'index',
        'value': 47,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': truss.bottom_chords
    }
    groups.append(BC_group)

    # Uumasauvojen profiili valikoimasta
    WEB_group = {
        'name': 'Webs',
        'var_type': 'index',
        'value': 33,
        'values': list(shs_profiles.keys()),
        'property': 'profile',
        'objects': list(truss.webs.values())
    }
    groups.append(WEB_group)

    return groups

def create_continuous_variable_g(structure):
    """
    Luodaan vapaavälin jatkuvat muuttujat
    :param structure:
    :return: vapaavälin jatkuvat muuttujat
    """
    groups = []

    truss = structure.truss[0]

    # Järjestetään uumasauvat pareiksi
    tc1 = truss.top_chords[0]
    tc_joints_1 = sorted(tc1.joints, key=lambda j: j.loc)
    tc2 = truss.top_chords[1]
    tc_joints_2 = sorted(tc2.joints, key=lambda j: j.loc, reverse=True)

    tc_joints = [[tc_joints_1[i], tc_joints_2[i]]
                 for i in range(0, len(tc_joints_1))]

    bc1 = truss.bottom_chords[0]
    bc_joints_list = sorted(bc1.joints, key=lambda j: j.loc)

    bc_joints = []

    ok = True

    while ok:
        joint_par = [bc_joints_list[0], bc_joints_list[-1]]
        bc_joints.append(joint_par)
        bc_joints_list.remove(bc_joints_list[0])
        bc_joints_list.remove(bc_joints_list[-1])
        if len(bc_joints_list) == 0:
            ok = False

    truss_joints = tc_joints + bc_joints

    for i, joint_pair in enumerate(truss_joints):
        groups.append({
            'name': 'Joint ' + str(i) + ' g',
            'var_type': 'continuous',
            'value': 30,
            'lb': 20,
            'ub': 90,
            'property': 'g1',
            'objects': joint_pair
            })

    return groups

def create_index_variable_g(structure):
    """
    Luodaa indeksimuuttujat vapaavälille
    :param structure:
    :return: vapaavälin indeksimuuttujat
    """
    groups = []

    truss = structure.truss[0]

    # Järjestetään uumasauvat pareiksi
    tc1 = truss.top_chords[0]
    tc_joints_1 = sorted(tc1.joints, key=lambda j: j.loc)

    tc2 = truss.top_chords[1]
    tc_joints_2 = sorted(tc2.joints, key=lambda j: j.loc, reverse=True)

    tc_joints = [[tc_joints_1[i], tc_joints_2[i]]
                 for i in range(0, len(tc_joints_1))]

    bc1 = truss.bottom_chords[0]
    bc_joints_list = sorted(bc1.joints, key=lambda j: j.loc)

    bc_joints = []

    ok = True

    while ok:
        joint_par = [bc_joints_list[0], bc_joints_list[-1]]
        bc_joints.append(joint_par)
        bc_joints_list.remove(bc_joints_list[0])
        bc_joints_list.remove(bc_joints_list[-1])
        if len(bc_joints_list) == 0:
            ok = False

    truss_joints = tc_joints + bc_joints
    # print(truss_joints)

    for i, joint_pair in enumerate(truss_joints):
        groups.append({
            'name': 'Joint ' + str(i) + ' g',
            'var_type': 'index',
            'value': 30,
            'values': g_values,
            'property': 'g1',
            'objects': joint_pair
            })

    return groups


def create_binary_discrete_variable_groups(structure):
    """
    Luodaan kehälle binäärimuuttujat
    :param structure:
    :return: kehän binäärimuuttujat
    """
    groups = []

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

    # Yläpaarteen leveys
    TC_h_group = {
        'name': 'TopChords',
        'var_type': 'discrete',
        'value': 180,
        'values': list(np.arange(100, 310, 10)),
        'property': ['H', 'B'],
        'objects': [tc.cross_section for tc in truss.top_chords],
    }
    groups.append(TC_h_group)
    # Yläpaarteen paksuus
    TC_t_group = {
        'name': 'TopChords',
        'var_type': 'discrete',
        'value': 10,
        'values': list(np.arange(3, 13, 1)),
        'property': ['T'],
        'objects': [tc.cross_section for tc in truss.top_chords],
    }
    groups.append(TC_t_group)

    # Alapaarteen leveys
    BC_h_group = {
        'name': 'BottomChords',
        'var_type': 'discrete',
        'value': 140,
        'values': list(np.arange(100, 310, 10)),
        'property': ['H', 'B'],
        'objects': [bc.cross_section for bc in truss.bottom_chords]
    }
    groups.append(BC_h_group)
    # Alapaarteen paksuus
    BC_t_group = {
        'name': 'TopChords',
        'var_type': 'discrete',
        'value': 10,
        'values': list(np.arange(3, 13, 1)),
        'property': ['T'],
        'objects': [bc.cross_section for bc in truss.bottom_chords],
    }
    groups.append(BC_t_group)

    # Uumasauvan leveys
    WEB_h_group = {
        'name': 'Webs',
        'var_type': 'discrete',
        'value': 120,
        'values': list(np.arange(80, 200, 10)),
        'property': ['H', 'B'],
        'objects': [web.cross_section for web in truss.webs.values()]
    }
    groups.append(WEB_h_group)
    # Uumasauvan paksuus
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
    """
    Luodaan rajoitusehtojen joukko
    :return:
    """
    pass


def discrete_neighbors(structure, n, params=('A', 'Iy'),
                       profiles=shs_profiles):
    """
    Luodaan diskreetti lähiympäristö
    :param structure:
    :param n: lähiympäristön koko
    :param params: muuttujat, joiden mukaan lähiympäristö valitaan
    :param profiles: profiilivalikoima
    :return: lähiympäristö
    """
    neighbors = []

    helper = SteelBeam([[0, 0], [0, 1000]])
    for mem in structure.members.values():
        vals = []
        for profile in profiles:
            helper.profile = profile
            d = 0
            for param in params:
                d += np.linalg.norm(mem.cross_section._getattribute_(
                    param) - helper.cross_section._getattribute_(param))
            vals.append(d)
        neighbors.append(list(np.argsort(vals)[:n]))

    return neighbors


if __name__ == '__main__':
    # Luodaan rakenne
    # 24 m ristikko
    # structure = create_structure(L=24000,
    #                              H0=8000,
    #                              H1=1800,
    #                              H2=2400,
    #                              dx=0,
    #                              n=14
    #                              )

    # 30 m ristikko
    structure = create_structure(L=30000,
                                 H0=8000,
                                 H1=2100,
                                 H2=2900,
                                 dx=0,
                                 n=14
                                 )

    # 36 m ristikko
    # structure = create_structure(L=36000,
    #                              H0=6000,
    #                              H1=2500,
    #                              H2=3500,
    #                              dx=0,
    #                              n=14
    #                              )

    # Jatkuvien muutttujien ylä- ja alarajat
    tc_values = [200, 10]
    tc_bounds = [[100, 300], [4, 12.5]]
    bc_values = [150, 8]
    bc_bounds = [[80, 200], [4, 10]]
    web_values = [100, 5]
    web_bounds = [[50, 150], [3, 8]]

    # Luodaan jatkuvien muuttujien joukko
    cont_var_groups = create_continuous_variable_groups(structure=structure,
                                                        tc_bounds=tc_bounds,
                                                        tc_values=tc_values,
                                                        bc_bounds=bc_bounds,
                                                        bc_values=bc_values,
                                                        web_bounds=web_bounds,
                                                        web_values=web_values,
                                                        )

    # Luodaan indeksimuuttujien joukko
    index_var_groups = create_index_variable_groups(structure=structure)

    # Luodaan diskreettien binäärimuuttujien joukko
    binary_disc_var_groups = create_binary_discrete_variable_groups(
        structure=structure)

    # Luodaan tehtävä/ongelma ja valitaan, mitä rajoitusehtoja otetaan mukaan
    problem = StructuralProblem(name="Example 1",
                                structure=structure,
                                var_groups=index_var_groups,
                                constraints={
                                    'joint_geometry_constraints': False,
                                    'joint_strength_constraints': True,
                                    'buckling_y': True,
                                    'buckling_z': True,
                                    'lt_buckling': True,
                                    'compression_bending_y': True,
                                    'compression_bending_z': True,
                                    'compression': True,
                                    'tension': True,
                                    'shear': True,
                                    # 'deflection_y': structure.L / 200 * 1.5,
                                    # 'deflection_x': structure.H / 300 * 1.5,
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
    # # structure.to_robot("lopputulos")

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

    # muutetaan muuttuja indeksimuuttujista binäärimuuttujiksi
    problem.index_to_binary()
    # for var in problem.vars:
    #     print(type(var))

    # MISLP
    solver = MISLP(move_limits=[0.05, 0.05], beta=100)
    # problem(x0)
    x0 = [var.value for var in problem.vars]
    print(x0)
    fopt, xopt = solver.solve(problem,
                              maxiter=50,
                              x0=x0,
                              min_diff=1e-2,
                              verb=True)
    problem(xopt, prec=5, ncons=100)

    # Luodaan tehtävä/ongelma toisen vaiheen optimointia varten
    index_problem = StructuralProblem(name="Example 2",
                                structure=structure,
                                var_groups=index_var_groups,
                                constraints={
                                    'joint_geometry_constraints': False,
                                    'joint_strength_constraints': True,
                                    'buckling_y': True,
                                    'buckling_z': True,
                                    'compression_bending_y': True,
                                    'compression_bending_z': True,
                                    'compression': True,
                                    'tension': True,
                                    'shear': True,
                                    'deflection_y': structure.L / 200 * 1.5,
                                    'deflection_x': (structure.H / 300) * 1.5,
                                    'web_class': True,
                                    'flange_class': True
                                })

    idx_vars = [var for var in index_problem.vars]
    ids = []
    for var in problem.vars:
        if 'Column' in var.name:
            idx = var.values.index(var.value)
            idx_vars[var.id].value = idx
        elif var.id not in ids:
            ids.append(var.id)
            mem = var.target['objects'][0]
            idx = list(shs_profiles.keys()).index(str(mem))
            idx_vars[var.id].value = idx

    index_problem.fea()
    index_problem([var.value for var in index_problem.vars])

    # VNS
    from deap.tools import mutShuffleIndexes
    solver = VNS(pop_size=50, maxiter=5, first_improvement=True, solver="GA",
                 solver_params={'pop_size': 50,
                                'mut_rate': 0.15,
                                'cx_rate': 0.9,
                                'first_improvement': True})
    x0 = [var.value for var in index_problem.vars]
    fopt, xopt = solver.solve(index_problem, x0=x0, maxiter=10, verb=True)
    index_problem(xopt, prec=5, ncons=100)

    index_variable_g_group = create_index_variable_g(structure=structure)

    continuous_variable_g_group = create_continuous_variable_g(structure=structure)

    g_problem = StructuralProblem(name="Example 3",
                                structure=structure,
                                var_groups=continuous_variable_g_group,
                                constraints={
                                    'joint_geometry_constraints': True,
                                    'joint_strength_constraints': False,
                                    'buckling_y': False,
                                    'buckling_z': False,
                                    'compression_bending_y': False,
                                    'compression_bending_z': False,
                                    'compression': False,
                                    'tension': False,
                                    'shear': False,
                                    # 'deflection_y': structure.L / 200,
                                    # 'deflection_x': (structure.H / 300) * 1.5,
                                    'web_class': False,
                                    'flange_class': False
                                })

    # TrustRegionConstr
    # solver = TrustRegionConstr()
    # x0 = [var.value for var in g_problem.vars]
    # f_best, x_best, nit = solver.solve(
    #     g_problem,
    #     maxiter=200,
    #     x0=x0)
    # g_problem.num_iters = nit
    # g_problem(x_best, prec=5)

    # SLP
    # solver = SLP(move_limits=[0.1, 0.1], beta=100)
    # x0 = [var.value for var in g_problem.vars]
    # solver.solve(g_problem,
    #              maxiter=500,
    #              # maxtime=3000,
    #              x0=x0,
    #              verb=True)
    # g_problem(solver.X, prec=5)

    # MISLP
    # solver = MISLP(move_limits=[0.05, 0.05], beta=100)
    # # problem(x0)
    # x0 = [var.value for var in g_problem.vars]
    # fopt, xopt = solver.solve(g_problem,
    #                           maxiter=50,
    #                           x0=x0,
    #                           min_diff=1e-2,
    #                           verb=True)
    # g_problem(xopt, prec=5)

    # VNS
    # solver = VNS(pop_size=50, maxiter=5, first_improvement=True, solver="GA",
    #              solver_params={'pop_size': 50,
    #                             'mut_rate': 0.15,
    #                             'cx_rate': 0.9,
    #                             'first_improvement': True})
    # x0 = [var.value for var in g_problem.vars]
    # fopt, xopt = solver.solve(g_problem, x0=x0, maxiter=10, verb=True)
    # g_problem(xopt, prec=5)

    # g_problem.structure.plot()

    # problem.structure.to_robot("lopputulos")

    print("")
    print("----PAINOT---------------------------------------------------")
    # Painot
    tc_weight = 0
    for tc in structure.truss[0].top_chords:
        tc_weight += tc.weight
    print("TOP CHORD WEIGHT =", tc_weight)

    bc_weight = 0
    for bc in structure.truss[0].bottom_chords:
        bc_weight += bc.weight
    print("BOTTOM CHORD WEIGHT =", bc_weight)

    web_weight = 0
    for web in structure.truss[0].webs.values():
        print("WEB WEIGHT =", web.weight)
        web_weight += web.weight
    print("WEBS WEIGHTS =", web_weight)

    # Ristikosta löytyy top_chords ja bottom_chords listoina, webs on
    # dicti(structure.truss[0])

    seconds = time.process_time()
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)

    print("")
    print("----------------------------")
    print("Process time:", "%d:%02d:%02d" % (h, m, s))
    print("----------------------------")
    print("")

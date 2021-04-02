from frame2d.simple_truss import ten_bar_truss
import optimization.structopt as sop
import frame2d.frame2d as f2d
from sections.timber.timber_section import TimberSection, HarjaPalkki, KaarevaPalkki
from materials.timber_data import T
import eurocodes.en1991.en1991 as en91
import numpy as np
from time import time
import datetime as dt

class Keha(sop.OptimizationProblem):

    def __init__(self):
        super().__init__()
        self.luo_rakenne()
        self.luo_muuttujat()
        self.luo_rajoitusehdot()
        self.luo_kohdefunktio()


    def luo_rakenne(self):
        fr = f2d.Frame2D()

        coord1 = [[0, 0], [0, 6000]]
        coord2 = [[0, 6000], [20000, 6000]]
        coord4 = [[20000, 6000], [20000, 0]]

        fr.add(f2d.TimberFrameMember(coord1, TimberSection(240, 360, material=T.GL30c), mtype='column', ldc='inst',
                                     num_elements=5, beta=0.85))
        fr.add(f2d.TimberFrameMember(coord2, HarjaPalkki(240, 600, 1300, material=T.GL30c), mtype='rigid-beam', ldc='mt',
                                     num_elements=20, edge_load='compression', lateral_support_z=[0.25, 0.5, 0.75]))
        fr.add(f2d.TimberFrameMember(coord4, TimberSection(240, 360, material=T.GL30c), mtype='column', ldc='inst',
                                     num_elements=5, beta=0.85))

        # fr.members[1].add_hinge(0)
        # fr.members[1].add_hinge(1)

        fr.add(f2d.FixedSupport([0, 0]))
        fr.add(f2d.FixedSupport([20000, 0]))

        qw = en91.WindLoad(0.5)
        qs = en91.SnowLoad(2.0)
        sw = en91.Load(1)
        sls = en91.SLSCombiner('B', kk_jako=6)
        sls.add([qw, 'x', [0]])
        sls.add([qs, 'y', [1]])
        sls.add([sw, 'y', [0, 1, 2]])
        sls_result = sls.get_result()
        uls = en91.ULSCombiner('B', kk_jako=6)
        uls.add([qw, 'x', [0]])
        uls.add([qs, 'y', [1]])
        uls.add([sw, 'y', [0, 1, 2]])
        uls_result = uls.get_result()
        acc = en91.ACCCombiner('B', kk_jako=6)
        acc.add([qw, 'x', [0]])
        acc.add([qs, 'y', [1]])
        acc.add([sw, 'y', [0, 1, 2]])
        acc_result = acc.get_result()

        for i, g in enumerate(sls_result[0][0]):
            if g == 0 or g is None:
                continue
            fr.add(f2d.LineLoad(fr.members[i], [g, g], 'x', load_id=1))

        for i, g in enumerate(sls_result[0][1]):
            if g == 0 or g is None:
                continue
            fr.add(f2d.LineLoad(fr.members[i], [-g, -g], 'y', load_id=1))

        for i, u in enumerate(uls_result[0]):
            if u == 0 or u is None:
                continue
            fr.add(f2d.LineLoad(fr.members[i], [u, u], 'x', load_id=2))

        for i, u in enumerate(uls_result[1]):
            if u == 0 or u is None:
                continue
            fr.add(f2d.LineLoad(fr.members[i], [-u, -u], 'y', load_id=2))

        for i, u in enumerate(acc_result[0]):
            if u == 0 or u is None:
                continue
            fr.add(f2d.LineLoad(fr.members[i], [u, u], 'x', load_id=3))

        for i, u in enumerate(acc_result[1]):
            if u == 0 or u is None:
                continue
            fr.add(f2d.LineLoad(fr.members[i], [-u, -u], 'y', load_id=3))

        #
        # #fr.add(f2d.LineLoad(fr.members[0], [5, 5], 'x'))
        # fr.add(f2d.LineLoad(fr.members[1], [-28.3, -28.3], 'y'))

        self.structure = fr
        fr.generate()
        fr.calculate('all')

    def luo_muuttujat(self):
        var_b = sop.Variable(name=f'Pilarit B', lb=100, ub=290,
                             target={'property': 'B', 'objects': [self.structure.members[0].cross_section, self.structure.members[2].cross_section]})
        var_h = sop.Variable(name=f'Pilarit H', lb=100, ub=990,
                             target={'property': 'H', 'objects': [self.structure.members[0].cross_section, self.structure.members[2].cross_section]})
        self.add(var_b)
        self.add(var_h)

        palk_b = sop.Variable(name=f'HarjaPalkki B', lb=100, ub=290,
                             target={'property': 'B', 'objects': [self.structure.members[1].cross_section]})
        palk_h = sop.Variable(name=f'HarjaPalkki H0', lb=100, ub=1500,
                             target={'property': 'H0', 'objects': [self.structure.members[1].cross_section]})
        palk_hap = sop.Variable(name=f'HarjaPalkki Hap', lb=500, ub=3000,
                             target={'property': 'Hap', 'objects': [self.structure.members[1].cross_section]})
        self.add(palk_b)
        self.add(palk_h)
        self.add(palk_hap)
        # al = sop.Variable(name=f'Palkin alpha', lb=20, ub=50, target={'property': 'Alpha', 'objects': [self.structure.members[1].cross_section]})
        # self.add(al)
        # be = sop.Variable(name=f'Palkin beta', lb=5, ub=20, target={'property': 'Beta', 'objects': [self.structure.members[1].cross_section]})
        # self.add(be)

    def mrt_rajoitusehto_generaattori(self, sauva):
        def UN(x):
            self.fea(2)
            return sauva.member.check_sections()[0] - 1

        def UV(x):
            return sauva.member.check_sections()[1] - 1

        def UM(x):
            return sauva.member.check_sections()[2] - 1

        def UBT(x):
            return sauva.member.check_sections()[4] - 1

        def UBC(x):
            return sauva.member.check_sections()[5] - 1

        def bucklingY(x):
            return sauva.member.check_buckling()[0] - 1

        def bucklingZ(x):
            return sauva.member.check_buckling()[1] - 1

        def buclingLT(x):
            return sauva.member.check_LT_buckling() - 1

        return UN, UV, UM, UBT, UBC, bucklingY, bucklingZ, buclingLT

    def krt_rajoitusehto_generaattori(self, sauva):
        def deflection(x):
            self.fea(1)
            maxdef = max([abs(max(md[1], key=abs)) for md in sauva.nodal_displacements.values()])
            return maxdef / (sauva.length / 200) - 1

        return deflection

    def acc_rajoitusehto_generaattori(self, sauva):
        def Palo_UN(x):
            sauva.R = 30
            self.fea(3)
            r = sauva.member.check_sections()[0] - 1
            sauva.R = 0
            return r

        return Palo_UN

    def luo_rajoitusehdot(self):
        self.add(sop.NonLinearConstraint(lambda x: x[0] - x[1], name=f'pilarit B < H'))
        self.add(sop.NonLinearConstraint(lambda x: x[3] - x[4], name=f'Palkki H0 < Hap'))
        self.add(sop.NonLinearConstraint(lambda x: x[2] - x[3], name=f'Palkki B < H0'))
        for i, sauva in enumerate(self.structure.members.values()):
            deflection = self.krt_rajoitusehto_generaattori(sauva)

            def_raj_ehto = sop.NonLinearConstraint(deflection, name=f'Taipuma sauvassa {i+1}')
            self.add(def_raj_ehto)

        for i, sauva in enumerate(self.structure.members.values()):
            UN, UV, UM, UBT, UBC, bucklingy, bucklingz, bucklinglt = self.mrt_rajoitusehto_generaattori(sauva)

            un_raj_ehto = sop.NonLinearConstraint(UN, name=f'UN sauvassa {i+1}', fea_required=True)
            self.add(un_raj_ehto)
            uv_raj_ehto = sop.NonLinearConstraint(UV, name=f'UV sauvassa {i+1}', fea_required=True)
            self.add(uv_raj_ehto)
            um_raj_ehto = sop.NonLinearConstraint(UM, name=f'UM sauvassa {i+1}', fea_required=True)
            self.add(um_raj_ehto)
            ubt_raj_ehto = sop.NonLinearConstraint(UBT, name=f'UBT sauvassa {i+1}', fea_required=True)
            self.add(ubt_raj_ehto)
            ubc_raj_ehto = sop.NonLinearConstraint(UBC, name=f'UBC sauvassa {i+1}', fea_required=True)
            self.add(ubc_raj_ehto)
            buck_raj_ehto = sop.NonLinearConstraint(bucklingy, name=f'buckling y sauvassa {i+1}', fea_required=True)
            self.add(buck_raj_ehto)
            buckz_raj_ehto = sop.NonLinearConstraint(bucklingz, name=f'buckling z sauvassa {i+1}', fea_required=True)
            self.add(buckz_raj_ehto)

            if sauva.mtype == 'rigid-beam':
                bucklt_raj_ehto = sop.NonLinearConstraint(bucklinglt, name=f'buckling lt sauvassa {i+1}', fea_required=True)
                self.add(bucklt_raj_ehto)
                self.add(sop.NonLinearConstraint(lambda x: sauva.member.check_perpendicular_compression() - 1, name=f'UNcp', fea_required=True))
                # self.add(sop.NonLinearConstraint(lambda x: sauva.member.check_apex_bending_stress() - 1, name='APXBS', fea_required=True))
                # self.add(sop.NonLinearConstraint(lambda x: sauva.member.check_apex_perpendicular_tension() - 1, name='APXPT', fea_required=True))
                # self.add(sop.NonLinearConstraint(lambda x: sauva.member.check_apex_shear_perpendicular_tension_combined() - 1, name='APXSPT', fea_required=True))

        for i, sauva in enumerate(self.structure.members.values()):
            Palo_UN = self.acc_rajoitusehto_generaattori(sauva)

            palo_un_raj_ehto = sop.NonLinearConstraint(Palo_UN, name=f'Palo UN sauvassa {i+1}')
            self.add(palo_un_raj_ehto)

    def luo_kohdefunktio(self):
        def kohde(x):
            return self.structure.weight
        obj = sop.ObjectiveFunction(name="Massa",
                                    obj_fun=kohde)
        self.add(obj)


if __name__ == '__main__':

    from scipy.optimize import minimize
    from optimization.solvers import *
    from log_to_file import log_to_file
    @log_to_file
    def keha():
        ongelma = Keha()

        start = time()
        #x0 = [600] * len(ongelma.vars)
        x0 = [240, 600, 240, 600, 1300]
        solver = PSO(pop_size=5)
        fopt, xopt = solver.solve(ongelma, maxiter=5, verb=False, log=True)

        ongelma(xopt)
        print('Laskentaan kulunut', round((time()-start), 2), 's')
        print(f'fea num {ongelma.num_fem_analyses}')
        print(f'states {ongelma.states}')
        print(f'xvals {solver.xvals}')

        ongelma.structure.calculate(2)
        ongelma.structure.members[0].print_utility_factor()
        ongelma.structure.members[1].print_utility_factor()
        ongelma.structure.members[2].print_utility_factor()

    keha()

    #ongelma.structure.plot()
    # ongelma.structure.plot_normal_force()
    # ongelma.structure.plot_deflection()
    # for i, sauva in enumerate(ongelma.structure.members.values()):
    #     print('sauvassa', i, 'on', sauva.ned // 1000, 'kN')

from frame2d.simple_truss import ten_bar_truss
import optimization.structopt as sop
import frame2d.frame2d as f2d
from sections.timber.timber_section import TimberSection
from materials.timber_data import T
import numpy as np
from time import time

class Keha(sop.OptimizationProblem):

    def __init__(self):
        super().__init__()
        self.luo_rakenne()
        self.luo_muuttujat()
        self.luo_rajoitusehdot()
        self.luo_kohdefunktio()


    def luo_rakenne(self):
        #ristikko = ten_bar_truss(L=3000, F=200e3)
        fr = f2d.Frame2D()

        coord1 = [[0, 0], [0, 6000]]
        coord2 = [[0, 6000], [16000, 6000]]
        coord4 = [[16000, 6000], [16000, 0]]

        fr.add(f2d.TimberFrameMember(coord1, TimberSection(240, 360, material=T.GL30c), mtype='column', ldc='inst',
                                     num_elements=5, beta=[[0.85], [0.85]]))
        fr.add(f2d.TimberFrameMember(coord2, TimberSection(240, 360, material=T.GL30c), mtype='rigid-beam', ldc='mt',
                                     num_elements=20, edge_load='compression'))
        fr.add(f2d.TimberFrameMember(coord4, TimberSection(240, 360, material=T.GL30c), mtype='column', ldc='inst',
                                     num_elements=5, beta=[[0.85], [0.85]]))

        fr.add(f2d.FixedSupport([0, 0]))
        fr.add(f2d.FixedSupport([16000, 0]))

        fr.add(f2d.LineLoad(fr.members[0], [5, 5], 'x'))
        fr.add(f2d.LineLoad(fr.members[1], [-15, -15], 'y'))

        self.structure = fr
        fr.generate()
        fr.calculate()

    def luo_muuttujat(self):
        var_b = sop.Variable(name=f'Pilarit B', lb=100, ub=1000,
                             target={'property': 'B', 'objects': [self.structure.members[0].cross_section, self.structure.members[2].cross_section]})
        var_h = sop.Variable(name=f'Pilarit H', lb=100, ub=10000,
                             target={'property': 'H', 'objects': [self.structure.members[0].cross_section, self.structure.members[2].cross_section]})
        self.add(var_b)
        self.add(var_h)

        palk_b = sop.Variable(name=f'Palkki B', lb=100, ub=1000,
                             target={'property': 'B', 'objects': [self.structure.members[1].cross_section]})
        palk_h = sop.Variable(name=f'Palkki H', lb=100, ub=10000,
                             target={'property': 'H', 'objects': [self.structure.members[1].cross_section]})
        self.add(palk_b)
        self.add(palk_h)


    def rajoitusehto_generaattori(self, sauva):
        # def NRd(x):
        #     A = sauva.cross_section.A
        #     fy = 355
        #     NEd = abs(sauva.ned)
        #     return NEd/(A*fy)-1
        #
        # def nurjahdus(x):
        #     if sauva.ned <= 0:
        #         return -sauva.ned / sauva.NbRd[0] - 1
        #     else:
        #         return -1
        def UN(x):
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


    def luo_rajoitusehdot(self):
        for i, sauva in enumerate(self.structure.members.values()):
            UN, UV, UM, UBT, UBC, bucklingy, bucklingz, bucklinglt = self.rajoitusehto_generaattori(sauva)

            un_raj_ehto = sop.NonLinearConstraint(UN, name=f'UN sauvassa {i+1}', fea_required=True)
            uv_raj_ehto = sop.NonLinearConstraint(UV, name=f'UV sauvassa {i+1}', fea_required=True)
            um_raj_ehto = sop.NonLinearConstraint(UM, name=f'UM sauvassa {i+1}', fea_required=True)
            ubt_raj_ehto = sop.NonLinearConstraint(UBT, name=f'UBT sauvassa {i+1}', fea_required=True)
            ubc_raj_ehto = sop.NonLinearConstraint(UBC, name=f'UBC sauvassa {i+1}', fea_required=True)
            buck_raj_ehto = sop.NonLinearConstraint(bucklingy, name=f'buckling y sauvassa {i+1}', fea_required=True)
            buckz_raj_ehto = sop.NonLinearConstraint(bucklingz, name=f'buckling z sauvassa {i+1}', fea_required=True)
            if sauva.mtype == 'rigid-beam':
                bucklt_raj_ehto = sop.NonLinearConstraint(bucklinglt, name=f'buckling lt sauvassa {i+1}', fea_required=True)
                self.add(bucklt_raj_ehto)
            self.add(un_raj_ehto)
            self.add(uv_raj_ehto)
            self.add(um_raj_ehto)
            self.add(ubt_raj_ehto)
            self.add(ubc_raj_ehto)
            self.add(buck_raj_ehto)
            self.add(buckz_raj_ehto)

            # NRd, nurjahdus = self.rajoitusehto_generaattori(sauva)
            # rajoitusehto = sop.NonLinearConstraint(name=f"Normaalivoimakestävyys sauvassa {i+1}",
            #                                        con_fun=NRd,
            #                                        fea_required=True)
            # self.add(rajoitusehto)
            #
            # nurjahdusehto = sop.NonLinearConstraint(name=f'Nurjahduskestävyys sauvassa {i+1}',
            #                                         con_fun=nurjahdus,
            #                                         fea_required=True)
            # self.add(nurjahdusehto)

    def luo_kohdefunktio(self):
        def kohde(x):
            return self.structure.weight
        obj = sop.ObjectiveFunction(name="Massa",
                                    obj_fun=kohde)
        self.add(obj)


if __name__ == '__main__':

    from scipy.optimize import minimize
    from optimization.solvers import *

    ongelma = Keha()

    start = time()
    #x0 = [600] * len(ongelma.vars)
    x0 = [300, 700, 300, 700]
    solver = SLP(move_limits=(0.01, 0.01))
    fopt, xopt = solver.solve(ongelma, x0=x0, maxiter=10, verb=False, log=True)

    ongelma(xopt)
    print('Laskentaan kulunut', round((time()-start), 2), 's')
    print(f'fea num {ongelma.num_fem_analyses}')
    print(f'states {ongelma.states}')
    print(f'xvals {solver.xvals}')

    ongelma.structure.members[0].print_utility_factor()
    ongelma.structure.members[1].print_utility_factor()
    ongelma.structure.members[2].print_utility_factor()

    #ongelma.structure.plot()
    # ongelma.structure.plot_normal_force()
    # ongelma.structure.plot_deflection()
    # for i, sauva in enumerate(ongelma.structure.members.values()):
    #     print('sauvassa', i, 'on', sauva.ned // 1000, 'kN')

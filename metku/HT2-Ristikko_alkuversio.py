try:
    from metku.frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, \
        Support, Hinge, XYHingedSupport, PREC
    from metku.framefem.elements import EBBeam, EBSemiRigidBeam
    # from fem.elements.eb_semi_rigid_beam import EBSemiRigidBeam
    from metku.sections.steel import HEA

    from metku.framefem import FrameFEM, BeamSection
    from metku.structures.steel.steel_member import SteelMember
    from metku.eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint, \
        RHSYJoint
except:
    from frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, \
        Support, Hinge, XYHingedSupport, PREC
    from framefem.elements import EBBeam, EBSemiRigidBeam
    # from fem.elements.eb_semi_rigid_beam import EBSemiRigidBeam
    from sections.steel import HEA
    from sections.steel.RHS import SHS

    from framefem import FrameFEM, BeamSection
    from structures.steel.steel_member import SteelMember
    from eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint, RHSYJoint

from frame2d.simple_truss import ten_bar_truss
import optimization.structopt as sop
from frame2d.simple_truss import SimpleTruss
import numpy as np


class Ristikko(sop.OptimizationProblem):

    def __init__(self):
        super().__init__()
        self.luo_rakenne()
        self.luo_muuttujat()
        self.luo_rajoitusehdot()
        self.luo_kohdefunktio()

    def luo_rakenne(self):
            """ (n % 4 == 0)-bar truss """

    L = 24000
    alpha = 1 / 20
    b = 60  # käyttäjä antaa asteissa (pitää olla välillä 20-80 astetta, iteroituu parhaimpaan mahdolliseen)
    beta = b * 2 * np.pi / 360  # muunnos asteista radiaaneiksi
    h_1 = L / 10  # käyttäjä syöttää jotakin (voidaan antaa vapasti)
    yp_vv = 4  # tehtävänannon mukainen vapaavälien lkm (voitaisiin antaa vapasti --> geometrinen optimointi)
    ap_vv = 7  # tehtävänannon mukainen vapaavälien lkm (voitaisiin antaa vapasti --> geometrinen optimointi)
    F = 21.5 / 9

    x1 = (L/2)/np.sqrt((L/2)**2+((L/2)*alpha)**2)           # yläpaarteen vasemman puolen suuntavektorin x komponentti
    y1 = (L/2*alpha)/np.sqrt((L/2)**2+((L/2)*alpha)**2)     # yläpaarteen vasemman puolen suuntavektorin y komponentti
    x2 = (L/2)/np.sqrt((L/2)**2+((L/2)*alpha)**2)           # yläpaarteen oikean puolen suuntavektorin x komponentti
    y2 = -(L/2*alpha)/np.sqrt((L/2)**2+((L/2)*alpha)**2)    # yläpaarteen oikean puolen suuntavektorin y komponentti

    # print(x1, y1, np.sqrt(x1 ** 2 + y1 ** 2))

    sv_v = [x1, y1]             # vasemman yläpaarteen suuntavektori
    sv_o = [x2, y2]             # oikean yläpaarteen suuntavektori

    vyp_p = np.sqrt((L/2)**2+((L/2)*alpha)**2)          # vasemman yläpaarteen pituus
    oyp_p = np.sqrt((L/2)**2+(-(L/2)*alpha)**2)         # oikean yläpaarteen pituus
    ap_p = L - 2*(h_1/beta)                             # alapaarteen pituus
    print(vyp_p, oyp_p, ap_p)

    vyp_vv = vyp_p / yp_vv       # vasemman yläpaarteen vapaavälin pituus
    oyp_vv = oyp_p / yp_vv       # oikean yläpaarteen vapaavälin pituus
    ap_vv = ap_p / ap_vv         # alapaarteen vapaavälin pituus
    print(round(vyp_vv, -2), round(ap_vv, -2))

    p9 = [L/2, L/2 * alpha]
    p7 = [L/2 - round(x1 * vyp_vv, 0), L/2 * alpha - round(y1 * vyp_vv, 0)]
    p5 = [L/2 - round(x1 * 2 * vyp_vv, 0), L/2 * alpha - round(y1 * 2 * vyp_vv, 0)]
    p3 = [L/2 - round(x1 * 3 * vyp_vv, 0), L/2 * alpha - round(y1 * 3 * vyp_vv, 0)]
    p1 = [L/2 - round(x1 * 4 * vyp_vv, 0), L/2 * alpha - round(y1 * 4 * vyp_vv, 0)]

    p11 = [L/2 + round(x2 * oyp_vv, 0), L/2 * alpha + round(y2 * oyp_vv, 0)]
    p13 = [L/2 + round(x2 * 2 * oyp_vv, 0), L/2 * alpha + round(y2 * 2 * oyp_vv, 0)]
    p15 = [L/2 + round(x2 * 3 * oyp_vv, 0), L/2 * alpha + round(y2 * 3 * oyp_vv, 0)]
    p17 = [L/2 + round(x2 * 4 * oyp_vv, 0), L/2 * alpha + round(y2 * 4 * oyp_vv, 0)]

    p2 = [round(h_1/beta, 0), -h_1]
    p4 = [p2[0] + round(ap_vv, 0), -h_1]
    p6 = [p4[0] + round(ap_vv, 0), -h_1]
    p8 = [p6[0] + round(ap_vv, 0), -h_1]

    p16 = [L - round(h_1/beta, 0), -h_1]
    p14 = [p16[0] - round(ap_vv, 0), -h_1]
    p12 = [p14[0] - round(ap_vv, 0), -h_1]
    p10 = [p12[0] - round(ap_vv, 0), -h_1]

    print(p1, p3, p5, p7, p9, p11, p13, p15, p17)
    print(p2, p4, p6, p8, p10, p12, p14, p16)


    nodes = []
    nodes.append(p1)
    nodes.append(p2)
    nodes.append(p3)
    nodes.append(p4)
    nodes.append(p5)
    nodes.append(p6)
    nodes.append(p7)
    nodes.append(p8)
    nodes.append(p9)
    nodes.append(p10)
    nodes.append(p11)
    nodes.append(p12)
    nodes.append(p13)
    nodes.append(p14)
    nodes.append(p15)
    nodes.append(p16)
    nodes.append(p17)

    sauvat = [[0, 2], [2, 4], [4, 6], [6, 8], [8, 10], [10, 12], [12, 14], [14, 16], [1, 3], [3, 5], [5, 7], [7, 9], [9, 11], [11, 13], [13, 15], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [12, 13], [13, 14], [14, 15], [15, 16]]

    profs = []
    for i in range(31):
        profs.append('SHS 50x50x3.0')

    t = SimpleTruss(nodes, sauvat, profs)

    t.add(PointLoad(nodes[0], [0, -F, 0]))
    t.add(PointLoad(nodes[2], [0, -F, 0]))
    t.add(PointLoad(nodes[4], [0, -F, 0]))
    t.add(PointLoad(nodes[6], [0, -F, 0]))
    t.add(PointLoad(nodes[8], [0, -F, 0]))
    t.add(PointLoad(nodes[10], [0, -F, 0]))
    t.add(PointLoad(nodes[12], [0, -F, 0]))
    t.add(PointLoad(nodes[14], [0, -F, 0]))
    t.add(PointLoad(nodes[16], [0, -F, 0]))
#    t.add(LineLoad)

    t.add(XYHingedSupport(nodes[0]))
    t.add(XYHingedSupport(nodes[16]))

    t.generate()
    t.calculate()

    t.plot()
        #self.structure = t
        #return t

    def luo_muuttujat(self):
        for i, sauva in enumerate(self.structure.members.values()):
            var = sop.Variable(name=f"sauva {i}",
                               lb=173,
                               ub=13704,
                               target={'property': 'A',
                                       'objects': [sauva.cross_section]})
            self.add(var)


    def rajoitusehto_generaattori(self, sauva):
        def NRd(x):
            A = sauva.cross_section.A
            fy = 355
            NEd = abs(sauva.ned)
            return NEd/(A*fy)-1

        def nurjahdus(x):
            #return -sauva.ned / sauva.NbRd[0] - 1
            if sauva.ned <= 0:
                sigma_cr = (np.pi ** 2 * sauva.cross_section.E * 0.252 * sauva.cross_section.A ** 1.142) / sauva.length ** 2
                sigma = -sauva.ned / sauva.cross_section.A
                return sigma / sigma_cr - 1
            else:
                return -1

        return NRd, nurjahdus


    def luo_rajoitusehdot(self):
        for i, sauva in enumerate(self.structure.members.values()):
            NRd, nurjahdus = self.rajoitusehto_generaattori(sauva)
            plk = sop.NonLinearConstraint(name=f"NRd {i}",
                                                con_fun=NRd,
                                                fea_required=True)
            self.add(plk)

            nurjahdusehto = sop.NonLinearConstraint(name=f'NbRd {i}',
                                                    con_fun=nurjahdus,
                                                    fea_required=True)
            self.add(nurjahdusehto)

    def luo_kohdefunktio(self):
        def kohde(x):
            return self.structure.weight
        obj = sop.ObjectiveFunction(name="Massa",
                                    obj_fun=kohde)
        self.add(obj)


if __name__ == '__main__':
    # t = three_bar_truss(L=3000,F1=-200e3,F2=-250e3)


    from metku.optimization.solvers import *

    ongelma = Ristikko()
    ongelma.structure.calculate()
#    ongelma.structure.plot_normal_force()

    x0 = [1000] * (19)
    solver = GA(pop_size=50)
    fopt, xopt = solver.solve(ongelma, x0=x0, maxiter=30, verb=False)
    print('xopt', xopt)
    ongelma(xopt)

# class Ristikko(sop.OptimizationProblem):
#
#     def __init__(self):
#         super().__init__()
#         self.luo_rakenne()
#         self.luo_muuttujat()
#         self.luo_rajoitusehdot()
#         self.luo_kohdefunktio()
#
#     def luo_rakenne(self):
#         truss = ristikko(L=3000, F=200e3)
#         self.structure = truss
#         truss.plot()
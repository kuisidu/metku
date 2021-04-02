"""
Timber member

@author Viktor Haimi
"""
from math import log, tan, pi, sqrt, isclose, sin, cos, radians
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from materials.timber_data import Timber, T
from eurocodes.en1995.en1995_1_1 import kmod, kdef, k_to_d, get_gammaM, k_m
from sections.timber.timber_section import TimberSection, PulpettiPalkki, HarjaPalkki, KaarevaHarjaPalkki, KaarevaPalkki
from framefem.elements import EBBeam, EBSemiRigidBeam, Rod
from sections.steel.ISection import HEA

PREC = 5

class TimberMember:
    def __init__(self, profile, material, length, ldc, sc, mtype='beam',
                 lcr=(0.9, 0.9), varnished=False, Sj1=np.inf, Sj2=np.inf, nodes=None, parent=None, lateral_support_y=None,
                 lateral_support_z=None, edge_load='compression', beta=None):
        """
        Initialization method for timber member
        @param material:    sauvan materiaali muodossa: C14..D70 / GL24c..GL32h
        @param H:           sauvan poikkileikkauksen korkeus
        @param B:           sauvan poikkileikkauksen leveys
        @param length:      sauvan pituus
        @param ldc:         kuorman aikaluokka
        @param sc:          rakenneosan käyttöluokka
        @param lcr:         nurjahduskerroin
        """
        self.material = material
        self.lcr = lcr
        self.ldc = ldc
        self.sc = sc
        self.varnished = varnished
        self.length = length

        self.type = mtype
        self.Sj1 = Sj1
        self.Sj2 = Sj2
        self.beta = beta

        self.ned = []
        self.myed = []
        self.mzed = []
        self.vyed = []
        self.vzed = []
        self.ted = []
        self.loc = []

        self.nodes = nodes
        self.parent = parent

        self.lateral_support_y = lateral_support_y
        self.lateral_support_z = lateral_support_z
        self.L_c = np.array([])

        self.edge_load = edge_load

        self.kmod = kmod(self.material.type, self.sc, self.ldc)

        self.section = profile

        self.R = 0

        self.k_L = 1.0
        if self.material.type == 'lvl':
            if self.section.Ned > 0 and self.length != 3000:
                self.k_L = min(1.1, (3000 / self.length) ** (self.material.s / 2))

        self.gammaM = get_gammaM(self.material.type)

    @property
    def R(self):
        return self.section.R

    @R.setter
    def R(self, val):
        self.section.R = val

    """
    Design values of material parameters are calculated from corresponding characteristic values
    """
    @property
    def fmd(self):
        return self.section.k_h * k_to_d(self.material.fmk, self.kmod, self.gammaM)

    @property
    def unmodified_fmd(self):
        return k_to_d(self.material.fmk, self.kmod, self.gammaM)

    @property
    def ft0d(self):
        if self.material.type == 'solid_timber' or self.material.type == 'glt':
            return self.section.k_h * k_to_d(self.material.ft0k, self.kmod, self.gammaM)
        elif self.material.type == 'lvl':
            return self.k_L * k_to_d(self.material.ft0k, self.kmod, self.gammaM)

    @property
    def unmodified_ft0d(self):
        return k_to_d(self.material.ft0k, self.kmod, self.gammaM)

    @property
    def ft90d(self):
        # TODO kysy samilta mikä on lvl materiaalin parametrien vastaavuus glt:hen nähden
        if self.material.type == 'lvl':
            raise Exception('lvl does not have ft90d, try ft90edged')
        return k_to_d(self.material.ft90k, self.kmod, self.gammaM)

    @property
    def fc0d(self):
        return k_to_d(self.material.fc0k, self.kmod, self.gammaM)

    @property
    def fc90d(self):
        if self.material.type == 'lvl':
            raise Exception('lvl does not have fc90d, try fc90edged or fc90flatd')
        return k_to_d(self.material.fc90k, self.kmod, self.gammaM)

    @property
    def fvd(self):
        return k_to_d(self.material.fvk, self.kmod, self.gammaM)

    @property
    def fm0flatd(self):
        if self.material.type == 'glt' or self.material.type == 'solid_timber':
            raise Exception(f'{self.material.type} does not have fm0flatd, try fmd')
        return k_to_d(self.material.fm0flatk, self.kmod, self.gammaM)

    @property
    def ft90edged(self):
        if self.material.type == 'glt' or self.material.type == 'solid_timber':
            raise Exception(f'{self.material.type} does not have ft90edged, try ft90d')
        return k_to_d(self.material.ft90edgek, self.kmod, self.gammaM)

    @property
    def fc90edged(self):
        if self.material.type == 'glt' or self.material.type == 'solid_timber':
            raise Exception(f'{self.material.type} does not have fc90edged, try fc90d')
        return k_to_d(self.material.fc90edgek, self.kmod, self.gammaM)

    @property
    def fc90flatd(self):
        if self.material.type == 'glt' or self.material.type == 'solid_timber':
            raise Exception(f'{self.material.type} does not have fc90flatd, try fc90d')
        return k_to_d(self.material.fc90flatk, self.kmod, self.gammaM)

    @property
    def fr0d(self):
        if self.material.type == 'glt' or self.material.type == 'solid_timber':
            raise Exception(f'{self.material.type} does not have fr0d')
        return k_to_d(self.material.fr0k, self.kmod, self.gammaM)

    @property
    def A(self):
        try:
            return self.section.A
        except:
            raise Exception('Use get_A')

    def get_A(self, at):
        return self.section.get_A(at)

    @property
    def Iy(self):
        return self.section.I[0]

    @Iy.setter
    def Iy(self, val):
        self.section.I[0] = val

    @property
    def Iz(self):
        return self.section.I[1]

    @Iz.setter
    def Iz(self, val):
        self.section.I[1] = val

    @property
    def E0mean(self):
        return self.material.E0mean

    @E0mean.setter
    def E0mean(self, val):
        self.material.E0mean = val

    # @property
    # def weight(self):
    #     """
    #     Function that calculates member's weight
    #     :return weight: float, weight of the member
    #     """
    #
    #     weight = self.A * self.length * self.material.rhomean
    #     return weight

    @property
    def h(self):
        """
        Property, returns member's cross-section's height
        """
        return self.section.get_H(0)

    @h.setter
    def h(self, val):
        """
        Sets cross-section's height to given value.
        """
        self.section.H = val

    @property
    def b(self):
        """
        Property, returns member's cross-section's breadth
        """
        return self.section.B

    @b.setter
    def b(self, val):
        """
        Sets cross-section's height to given value.
        """
        self.section.B = val

    def sigma_t0d(self, n):
        return self.ned[n] / self.section.get_A(self.loc[n])

    def sigma_c0d(self, n):
        return self.ned[n] / self.section.get_A(self.loc[n])

    def sigma_m_alpha_d(self):
        return max([self.myed[self.loc.index(i)] / self.section.get_Wel(i)[0] for i in self.loc], key=abs)

    def sigma_t90d(self):
        # (6.54)
        # print(f'kp {self.section.k_p}')
        # print(f'Mapd {self.M_apd()}')
        # print(f'Wap {self.W_ap()[0]}')

        return self.k_p() * self.M_apd() / self.W_ap()[0]

    def M_apd(self):
        try:
            return self.myed[self.loc.index(self.section.xap)]
        except:
            minusloc = None
            for loc in self.loc:
                if minusloc is None:
                    minusloc = loc
                    continue
                diff1 = self.section.xap - loc
                diff2 = self.section.xap - minusloc
                if diff2 > diff1 >= 0:
                    minusloc = loc
            plusloc = None
            for loc in self.loc[::-1]:
                if plusloc is None:
                    plusloc = loc
                    continue
                diff1 = loc - self.section.xap
                diff2 = plusloc - self.section.xap
                if diff2 > diff1 >= 0:
                    plusloc = loc
            return (self.myed[self.loc.index(minusloc)] + self.myed[self.loc.index(plusloc)]) / 2

    def W_ap(self):
        return self.section.get_Wel(self.section.xap)

    def k_r(self):
        if isinstance(self.section, (HarjaPalkki, PulpettiPalkki)):
            return 1
        elif isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki)):
            if self.section.r_in * 1000 / self.section.t >= 240:
                return 1
            else:
                return 0.76 + 0.001 * (self.section.r_in * 1000 / self.section.t)
        else:
            raise NotImplementedError

    def V(self):
        # V <= 2/3 * koko palkin tilavuus
        if isinstance(self.section, HarjaPalkki):
            return self.section.B * np.power(self.section.hap, 2)
        elif isinstance(self.section, KaarevaHarjaPalkki):

            l1 = self.section.h0 * cos(radians(self.alfa()))
            l2 = self.section.up_y1 - self.section.up_y2 / cos(radians(self.alfa()))
            l3 = self.section.up_x1 / cos(radians(self.alfa()))
            Vj = self.section.B * l3 * (l1 + l2) * 0.5
            # print(f'Vj: {Vj}')
            # print(f'2 Vj: {2* Vj}')
            # print(f'V_tot: {self.V_tot()}')
            # print(f'V_tot - 2 Vj: {self.V_tot() - 2 * Vj}')
            return self.V_tot() - 2 * Vj
            #
            # alpha = self.alfa() * 2 * pi / 360
            # return self.section.B * (tan(alpha) * (self.section.r_in + self.section.hap) ** 2 -
            #         sin(alpha) * cos(alpha) * (tan(alpha) * (self.section.r_in + self.section.hap)) ** 2 -
            #         pi * self.section.r_in ** 2 * alpha / (2 * pi)) / (1000 ** 3)
        elif isinstance(self.section, KaarevaPalkki):
            alpha = self.alfa() * 2 * pi / 360
            return self.section.B * self.section.hap * (self.section.r_in + self.section.hap / 2) * alpha
        else:
            raise NotImplementedError

    def V_tot(self):
        if isinstance(self.section, HarjaPalkki):
            return self.section.B * self.length * (self.section.hap - self.length * np.tan(self.alfa()) / 4)
        elif isinstance(self.section, (KaarevaPalkki, KaarevaHarjaPalkki)):
            return self.section.A_face() * self.section.B
        else:
            raise NotImplementedError

    def k_vol(self):
        if self.material.type == 'solid timber':
            return 1.0
        v = self.V()
        if v > (2 / 3) * self.V_tot():
            v = (2/3) * self.V_tot()
        k_vol = np.power(1e7 / v, 0.2)
        return k_vol

    def k_dis(self):
        if isinstance(self.section, (HarjaPalkki, PulpettiPalkki)):
            return 1.4
        elif isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki)):
            return 1.7

    def k_p(self):
        if isinstance(self.section, HarjaPalkki):
            return 0.2 * np.tan(self.alfa())
        elif isinstance(self.section, KaarevaPalkki):
            r = self.section.r_in + 0.5 * self.section.hap
            return 0.25 * (self.section.hap / r)
        elif isinstance(self.section, KaarevaHarjaPalkki):
            k5 = 0.2 * np.tan(self.alfa())
            k6 = 0.25 - 1.5 * np.tan(self.alfa()) + 2.6 * np.tan(self.alfa()) ** 2
            k7 = 2.1 * np.tan(self.alfa()) - 4 * np.tan(self.alfa()) ** 2
            r = self.section.r_in + 0.5 * self.section.hap
            return k5 + k6 * (self.section.hap / r) + k7 * (self.section.hap / r) ** 2

    def k_l(self):
        k1 = 1 + 1.4 * np.tan(self.alfa()) + 5.4 * np.tan(self.alfa()) ** 2
        if isinstance(self.section, HarjaPalkki):
            return k1
        elif isinstance(self.section, KaarevaPalkki):
            r = self.section.r_in + 0.5 * self.section.hap
            return 1 + 0.35 * (self.section.hap / r) + 0.6 * (self.section.hap / r) ** 2
        elif isinstance(self.section, KaarevaHarjaPalkki):
            r = self.section.r_in + 0.5 * self.section.hap
            k2 = 0.35 * - 8 * np.tan(self.alfa())
            k3 = 0.6 + 8.3 * np.tan(self.alfa()) - 7.8 * np.tan(self.alfa()) ** 2
            k4 = 6 * np.tan(self.alfa()) ** 2
            return k1 + k2 * (self.section.hap / r) + k3 * (
                    self.section.hap / r) ** 2 + k4 * (self.section.hap / r) ** 3

    def k_m_alpha(self):
        sigma_malphad = self.sigma_m_alpha_d()
        if not self.section.flipped:
            if sigma_malphad >= 0:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (4 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1/2)
            else:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (2 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1 / 2)
        else:
            if sigma_malphad <= 0:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (4 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1/2)
            else:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (2 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1 / 2)
        return kma

    def alfa(self):
        if isinstance(self.section, PulpettiPalkki):
            a = (self.section.hap - self.section.h0) / self.length
        elif isinstance(self.section, HarjaPalkki):
            a = (self.section.hap - self.section.h0) / (self.length / 2)
        elif isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki)):
            a = self.section.alpha
        #elif mahapalkki:
        #     a = (h_ap - h0)*2 / (self.length / 2)
        return a


    #
    #
    # print(Vkahapa(200, 10000, 800, 5))
    #
    #
    # def Vbumpa(h, b, r_in, alpha):
    #     return 4 * h * b * (alpha * math.pi / 360) * (r_in + h/2) / (1000 ** 3)
    #
    # print(Vbumpa(200, 600, 10000, 5))

    def sigma_c90d(self, node):
        """
        @param l:   sen sauvan kosketuspinnan pituus, joka rasitetaan syysuuntaa vastaan esim. pilarin päällä oleva palkki
        @param b:   palkin leveys, joka puristaa toista palkkia sen syysuuntaan nähden kohtisuorasti
        @param ned: kosketuspinnalla vallitseva normaalivoima puristavasta sauvasta
        @return:    sigma_c90d
        """
        '''
        Pitää 1) tarkistaa missä nodessa ollaan
        2) otetaan sen noden vanhemmat, katsotaan omasta memberistä molemmat elementit 
        3) katsotaan elementtien periusteelkla -- > elementtien normaalivoimiin --> saadaan seuraavien nodien normaalivoimat
        4) saadaan toisen vanhemman seuraava node --> tiedossa 4 noden tiedot --> ympäröivistä nodeista lasketaan geometrian perusteella halutun noden voimien komponentit
        '''

        if len(node.parents) == 1:
            return 0

        other_members = list(filter(lambda v: v.member != self if True else False, node.parents.copy()))

        if len(other_members) == 0:
            return 0

        def get_sigma(om):
            other_locs = om.member.loc
            other_idx = None
            for i in range(om.member.nsect()):
                if isclose(node.x, om.to_global(other_locs[i])[0], rel_tol=10) and isclose(node.y,
                                                                                           om.to_global(other_locs[i])[
                                                                                               1],
                                                                                           rel_tol=10):
                    other_idx = i

            own_idx = None
            for i in range(self.nsect()):
                if isclose(node.x, self.parent.to_global(self.loc[i])[0], rel_tol=10) and isclose(node.y,
                                                                                                  self.parent.to_global(
                                                                                                  self.loc[i])[1],
                                                                                                  rel_tol=10):
                    own_idx = i

            if om.member.ned[other_idx] >= 0:
                return 0

            if own_idx == 0 or own_idx == self.loc[-1]:
                # Nolla edustaa a mittaa. RIL 205-1-2009 kohta 6.1.5
                lc90ef = max(om.cross_section.H, om.cross_section.B) + min(30, 0, om.cross_section.H) + min(30, om.cross_section.H)
            else:
                lc90ef = max(om.cross_section.H, om.cross_section.B) + min(30, om.cross_section.H) + min(30, om.cross_section.H)

            # SFS-EN 1995-1-1 § 6.1.5
            Aef = lc90ef * min(om.cross_section.B, self.section.B)
            sigma = om.member.ned[other_idx] / Aef
            return sigma

        if len(other_members) == 1:
            om = other_members[0]
            return get_sigma(om)
        else:
            return [get_sigma(om) for om in other_members]

    def kc_perpendicular(self, l, a, l1) -> float:
        """
        Tukipainekerroin: RIL 205-1-2009,
        @param l:   sen sauvan kosketuspinnan pituus, joka rasitetaan syysuuntaa vastaan esim. pilarin päällä oleva palkki
        @param a:   kosketuspinnan ulkopuolinen, puristetun palkin päätyyn asti ulottuva mitta, tai seuraavaan puristuspintaan
        @param l1:  puristuspintojen lähin välimatka toisiinsa nähden
        @return:    tukipainekerroin
        """
        self.update_kc90(l1, l)
        lc90ef = l + min(30, a, l) + min(30, l, l1 / 2)
        return lc90ef / l * self.k_c90


    def update_kc90(self, l1, l):
        """
        Materiaalikohtainen, tukipainekertoimen määritykseen tarkoitettu kerroin: SFS-EN 1995-1-1 § 6.1.5

        Rajataan SFS-EN 1995-1-1 § 6.1.5 Kuva 6.2 b kohtaan

        @param l1:  kosketuspintojen reunojen välinen etäisyys
        @param l: Toisen palkin poikkileikkauksen korkeus
        @return:
        """
        self.k_c90 = 1.0
        if self.section.H * 2 <= l1:
            if self.material.type == 'solid_timber':
                self.k_c90 = 2.5
            elif self.material.type == 'glt' and l <= 400:
                self.k_c90 = 1.75


    def sigma_md(self, n):
        return [self.myed[n] / self.section.get_Wel(self.loc[n])[0],
                self.mzed[n] / self.section.get_Wel(self.loc[n])[1]]

    def get_max_sigma_per_segment(self):
        supp_y = self.lateral_support_y.copy()
        supp_z = self.lateral_support_z.copy()
        supp_y.append(0.0)
        supp_y.append(1.0)
        supp_z.append(0.0)
        supp_z.append(1.0)
        supp_y.sort()
        supp_z.sort()

        sig_md_list = [self.sigma_md(i) for i in range(self.nsect())]

        iy = []
        iz = []
        for i in range(len(supp_y) - 1):
            sig_md_in_seg = []
            locs = []
            for loc in self.loc:
                #print(f'Member {self}')
                #print(f'loc {loc}, supp i {supp_y[i]}, supp i+1 {supp_y[i+1]}')
                if supp_y[i] <= loc <= supp_y[i+1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iy.append(max(sig_md_in_seg, key=lambda x: abs(x[0])))

        for i in range(len(supp_z) - 1):
            sig_md_in_seg = []
            locs = []
            for loc in self.loc:
                if supp_z[i] <= loc <= supp_z[i+1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iz.append(max(sig_md_in_seg, key=lambda x: abs(x[0])))
        return [iy, iz]

    def tau_d(self, n):
        """
        Poikkileikkauksessa vaikuttava työntövoima muutetaan leikkausjännitykseski Jourawskin kaavaa hyväksikäyttäen
        @return: palauttaa suorakaidepoikkileikkauksen leikkausjännityksen arvon, ottaa huomioon k_cr kertoimen
        """
        k_cr = 1.0
        if self.material.type == 'solid_timber' or self.material.type == 'glt':
            if not self.varnished:
                k_cr = 0.67

        return [(3/2 * self.vyed[n]) / (self.section.get_H(self.loc[n]) * self.section.B * k_cr),
                (3/2 * self.vzed[n]) / (self.section.get_H(self.loc[n]) * self.section.B * k_cr)]


    def tau_tord(self, n) -> float:
        """
        Evaluation of torsional constant Iv (or J) and torsional resistance Wv for any arbitrarily selected cross section
        require rather sophisticated calculations. In structural wood design practice, cross sections are selected using
        manufacturer's catalogs and only sometimes random section dimensions are used. In this method only rectangular
        cross sections are permitted, therefore it's possible to estimate Iv or Wv values using simple mathematical
        expressions. Equations developed here by far are not unique, one can develop better ones, but for the time being
        they are giving sufficiently close enough results to compare with commercial programs.

        Torsional constant Iv can be estimated as stated below (change the base of logarithm to get better results):
        alpha = 0.08 * log(self.H / self.B, 2.5) + 0.14
        Iv = alpha * self.H * self.B ** 3
        but... to calculate torsional stress, torsional resistance Wv is used, which can be obtained from an analogous
        equation pair (here manipulating the formula one can also get better results for some particular cross section):
        beta = 0.05 * log(self.H / self.B, 2.6) + 0.2
        Wv = beta * self.H * self.B ** 2

        Alpha and beta are based on the theoretical tables given in: T. Salmi 2000, Lujuusopin perusteet, p. 258
        Functions were tested with logarithmic base 2..3 and middle value was selected, which underestimates close to
        square-shaped section properties and overestimates high and thin section properties, optimality is at 400x70
        """
        section_height_to_width_ratio = self.section.get_H(self.loc[n]) / self.section.B
        if section_height_to_width_ratio < 1:
            section_height_to_width_ratio = 1 / section_height_to_width_ratio
        if 0.95 < section_height_to_width_ratio < 1.05:
            beta = 0.208
        elif 1.45 < section_height_to_width_ratio < 1.55:
            beta = 0.231
        elif 1.95 < section_height_to_width_ratio < 2.05:
            beta = 0.246
        elif 2.45 < section_height_to_width_ratio < 2.55:
            beta = 0.258
        elif 2.95 < section_height_to_width_ratio < 3.05:
            beta = 0.267
        elif 3.95 < section_height_to_width_ratio < 4.05:
            beta = 0.282
        elif 4.95 < section_height_to_width_ratio < 5.05:
            beta = 0.299
        elif 5.95 < section_height_to_width_ratio < 6.05:
            beta = 0.307
        elif 7.95 < section_height_to_width_ratio < 8.05:
            beta = 0.313
        elif 9 < section_height_to_width_ratio < np.inf:
            beta = 0.333
        else:
            beta = 0.05 * log(section_height_to_width_ratio, 2.6) + 0.2
        Wv = beta * self.section.get_H(self.loc[n]) * self.section.B ** 2
        return self.ted[n] / Wv

    def I_tor(self, loc) -> float:
        """
        See help(tau_tord) for full description
        @return: torsional constant for rectangular cross sections
        """
        ratio = self.section.get_H(loc) / self.section.B
        if ratio < 1:
            ratio = 1 / ratio
        if 0.95 < ratio < 1.05:
            alpha = 0.141
        elif 1.45 < ratio < 1.55:
            alpha = 0.196
        elif 1.95 < ratio < 2.05:
            alpha = 0.229
        elif 2.45 < ratio < 2.55:
            alpha = 0.249
        elif 2.95 < ratio < 3.05:
            alpha = 0.263
        elif 3.95 < ratio < 4.05:
            alpha = 0.281
        elif 4.95 < ratio < 5.05:
            alpha = 0.299
        elif 5.95 < ratio < 6.05:
            alpha = 0.307
        elif 7.95 < ratio < 8.05:
            alpha = 0.313
        elif 9 < ratio < np.inf:
            alpha = 0.333
        else:
            alpha = 0.08 * log(ratio, 2.5) + 0.14
        return alpha * self.section.get_H(loc) * self.section.B ** 3

    def k_shape(self, n) -> float:
        """
        SFS EN 1995-1-1 §6.1.8(6.15)
        @return: cross section factor for rectangular cross sections (circular sections not included)
        """
        return min(1.3, 1 + 0.05 * self.section.get_H(self.loc[n]) / self.section.B)

    def sigma_alphad(self):
        pass

    def beta_c(self):
        if self.material.type == 'solid_timber':
            return 0.2
        return 0.1

    def radius_of_gyration(self):
        supp_y = self.lateral_support_y.copy()
        supp_z = self.lateral_support_z.copy()
        supp_y.append(0.0)
        supp_y.append(1.0)
        supp_z.append(0.0)
        supp_z.append(1.0)
        supp_y.sort()
        supp_z.sort()

        sig_md_list = [self.sigma_md(i) for i in range(self.nsect())]
        iy = []
        iz = []
        for i in range(len(supp_y) - 1):
            sig_md_in_seg = []
            locs = []
            for loc in self.loc:
                if supp_y[i] <= loc <= supp_y[i+1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iy.append(locs[sig_md_in_seg.index(min(sig_md_in_seg))])

        for i in range(len(supp_z) - 1):
            sig_md_in_seg = []
            locs = []
            for loc in self.loc:
                if supp_z[i] <= loc <= supp_z[i+1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iz.append(locs[sig_md_in_seg.index(min(sig_md_in_seg))])

        y = [sqrt(self.section.get_I(i)[0] / self.section.get_A(i)) for i in iy]
        z = [sqrt(self.section.get_I(i)[1] / self.section.get_A(i)) for i in iz]

        return np.asarray([np.asarray(y),
                           np.asarray(z)], dtype=object)

    def Ln(self):
        # Beta arvot saatu Liimapuukäsikirja Osa 3 sivu 3-35
        # TODO arvot tarkistettu vain liimapuulle.
        # TODO tarkista betan arvot palkille ja pilarille erikseen Lpkk3 (3-39)
        beta = [[], []]

        # TODO vaihtuvan poikkileikkauksen sauvassa nurjahduspituus tarkistettava. nyt vakio poikkileikkaus
        def form_beta(lc: list, index: int, beam=False):
            if len(lc) == 1:
                if beam and index == 1:
                    beta[index] = 1
                elif self.Sj1 >= 1e15 and self.Sj2 >= 1e15:
                    beta[index] = 0.7
                elif (self.Sj1 >= 1e15 and self.Sj2 <= 1e3) or (self.Sj2 >= 1e15 and self.Sj1 <= 1e3):
                    beta[index] = 0.85
                else:
                    beta[index] = 1
            else:
                for i in range(len(lc)):
                    if beam and index == 1:
                        beta[index].append(1)
                    elif i == 0:
                        if index == 1:
                            beta[index].append(0.85)
                        elif self.Sj1 >= 1e15:
                            beta[index].append(0.85)
                        else:
                            beta[index].append(1)
                    elif i == len(lc) - 1:
                        if index == 1:
                            beta[index].append(0.85)
                        elif self.Sj2 >= 1e15:
                            beta[index].append(0.85)
                        else:
                            beta[index].append(1)
                    else:
                        beta[index].append(1)

        if self.beta is not None:
            if isinstance(self.beta, float):
                return np.asarray([self.beta * self.L_c[0], self.beta * self.L_c[1]], dtype=object)
            elif isinstance(self.beta, list):
                return np.asarray([self.beta[0] * self.L_c[0], self.beta[1] * self.L_c[1]], dtype=object)
        if self.type == 'column':
            form_beta(self.L_c[0], 0)
            form_beta(self.L_c[1], 1)

            return np.asarray([beta[0] * self.L_c[0], beta[1] * self.L_c[1]], dtype=object)
        elif self.type == 'rigid-beam':
            form_beta(self.L_c[0], 0, True)
            form_beta(self.L_c[1], 1, True)

            return np.asarray([beta[0] * self.L_c[0], beta[1] * self.L_c[1]], dtype=object)

        if self.Sj1 >= 1e15 and self.Sj2 >= 1e15:
            beta = 0.7
        elif (self.Sj1 >= 1e15 and self.Sj2 <= 1e3) or (self.Sj2 >= 1e15 and self.Sj1 <= 1e3):
            beta = 0.85
        else:
            beta = 1
        return [beta * self.L_c[0], beta * self.L_c[1]]

    def lamda(self):
        return [self.Ln()[0] / self.radius_of_gyration()[0],
                self.Ln()[1] / self.radius_of_gyration()[1]]


    def lamda_rel(self):
        return [(self.lamda()[0] / pi) * sqrt(self.material.fc0k / self.material.E005),
                (self.lamda()[1] / pi) * sqrt(self.material.fc0k / self.material.E005)]


    def sigma(self):
        return [pi ** 2 * (self.material.E0mean / self.lamda()[0]),
                pi ** 2 * (self.material.E0mean / self.lamda()[1])]


    def k(self):
        return [0.5 * (1 + self.beta_c() * (self.lamda_rel()[0] - 0.3) + np.power(self.lamda_rel()[0], 2)),
                0.5 * (1 + self.beta_c() * (self.lamda_rel()[1] - 0.3) + np.power(self.lamda_rel()[1], 2))]


    def k_c(self):
        # SFS-EN 1995-1-1 § 6.3.2 (6.26)
        return [np.power(self.k()[0] + np.power(np.power(self.k()[0], 2) - np.power(self.lamda_rel()[0], 2), 0.5), -1),
                np.power(self.k()[1] + np.power(np.power(self.k()[1], 2) - np.power(self.lamda_rel()[1], 2), 0.5), -1)]


    def l_ef(self, i, yz, n):
        # SFS-EN 1995-1-1 § 6.3.3 (6.1)
        # TODO SFS 1995-1-1 s. 43 Taulukko 6.1
        # TODO Pistekuorma ja vakiomomentti tilanne puuttuu

        alpha = 0.9
        # if i == 0:
        #
        # elif i == len(self.L_c) - 1:
        #     pass
        # else:
        #     alpha = 0.9
        if self.edge_load == 'compression':
            return (self.L_c[yz][i] * alpha) + 2 * self.section.get_H(n)
        elif self.edge_load == 'tension':
            return (self.L_c[yz][i] * alpha) - 0.5 * self.section.get_H(n)

        return self.L_c[yz][i] * alpha


    def sigma_mcrit(self, i, yz, loc):
        '''
        For member undergoing bending relative to it's stronger axis: SFS-EN 1995-1-1 § 6.3.3 (6.31)
        @return: Single value of critical stress during lateral buckling check
        '''
        try:
            return pi * sqrt(self.material.E005 * self.section.get_I(loc)[1] * self.material.G005 * self.I_tor(loc)) / (self.l_ef(i, yz, loc) * self.section.get_Wel(loc)[0])
        except:
            return ((0.78 * self.section.B ** 2) / (self.section.H * self.l_ef(i, yz, loc))) * self.material.E005


    def lamda_relm(self, i, yz, n):
        '''
        SFS-EN 1995-1-1 § 6.3.3 (6.30)
        @return: Relative slenderness for bending problems      TODO taarvitaanko sauvan molempien akseleiden suhteen?
        '''
        return np.sqrt(self.material.fmk / self.sigma_mcrit(i, yz, n))


    def k_crit(self, i, yz, n):
        # SFS-EN 1995-1-1 § 6.3.3 (6.34)
        if self.lamda_relm(i, yz, n) <= 0.75:
            return 1
        elif 0.75 < self.lamda_relm(i, yz, n) <= 1.4:
            return 1.56 - 0.75 * self.lamda_relm(i, yz, n)
        else:
            return 1 / self.lamda_relm(i, yz, n) ** 2

    def check_single_tapered_beam_bending_tension(self):
        return round(self.sigma_m_alpha_d() / (self.k_m_alpha() * self.fmd), PREC)

    def check_tapered_beam_bending_tension(self):
        try:
            sig = self.section.k_l * self.sigma_md(self.loc.index(self.section.xap))[0]
        except:
            minusloc = None
            for loc in self.loc:
                if minusloc is None:
                    minusloc = loc
                    continue
                diff1 = self.section.xap - loc
                diff2 = self.section.xap - minusloc
                if diff2 > diff1 >= 0:
                    minusloc = loc
            plusloc = None
            for loc in self.loc[::-1]:
                if plusloc is None:
                    plusloc = loc
                    continue
                diff1 = loc - self.section.xap
                diff2 = plusloc - self.section.xap
                if diff2 > diff1 >= 0:
                    plusloc = loc
            sig = self.section.k_l * (
                        self.sigma_md(self.loc.index(minusloc))[0] + self.sigma_md(self.loc.index(plusloc))[
                    0]) / 2

        return round(sig / (self.k_r() * self.fmd), PREC)

    def check_perpendicular_to_grain_tension(self):
        print(f'sig t90d {self.sigma_t90d()}')
        print(f'k dis {self.k_dis()}')
        print(f'k vol {self.k_vol()}')
        print(f'f t90d {self.ft90d}')
        return self.sigma_t90d() / (self.k_dis() * self.k_vol() * self.ft90d)

    def check_normal_force(self, n):
        # Normal force utilization ratio along the grain:           SFS-EN 1995-1-1 § 6.1.2 and 6.1.4
        if self.ned[n] > 0:
            UN = round(abs(self.sigma_t0d(n)) / self.ft0d, 5)
        else:
            UN = round(abs(self.sigma_t0d(n)) / self.fc0d, 5)
        r = UN
        return r
#        return print('Normal force resistance capacity utilized until', UN * 100, '%')


    def check_perpendicular_tension(self):
        # Member tension perpendicular to the grain:                SFS-EN 1995-1-1 § 6.1.3
        # TODO syyn suuntaan vastaset voimat
        pass


    def check_perpendicular_compression(self, location=False):
        # Member compression perpendicular to the grain:            SFS-EN 1995-1-1 § 6.1.5

        #TODO ei laske oikein yhdessä pisteessä kahdelta suunnalta puristettua palkkia
        max_sigma = 0
        node = None
        l1 = None
        for n in self.nodes.values():
            n_sig = self.sigma_c90d(n)
            if isinstance(n_sig, list):
                n_sig = min(n_sig)
            if abs(n_sig) > abs(max_sigma):
                max_sigma = n_sig
                node = n
        if node:
            for n in self.nodes.values():
                if len(n.parents) == 1 or n == node:
                    continue

                other_members = list(filter(lambda v: v.member != self if True else False, n.parents.copy()))

                if len(other_members) == 0:
                    continue

                start_node = n.coord
                end_node = node.coord
                l11 = np.linalg.norm(end_node - start_node)
                if l1 is None:
                    l1 = l11
                    self.update_kc90(l1, other_members[0].cross_section.get_H(0))
                elif l11 < l1:
                    l1 = l11
                    # TODO check get_H
                    self.update_kc90(l1, other_members[0].cross_section.get_H(0))

        if not hasattr(self, 'k_c90'):
            return 0

        r = abs(round(max_sigma / (self.k_c90 * self.fc90d), 3))
        if location:
            return [r, node]
        return r


    def check_shear_force(self, n):
        # Shear force utilization ratio:                            SFS-EN 1995-1-1 § 6.1.7

        if self.type in ("beam", "rigid-beam"):
            start = list(filter(lambda v: v.member != self if True else False, self.parent.n1.parents.copy()))
            other_h_start = start[0].cross_section.H
            sumH_start = (other_h_start + self.section.get_H(0)) / self.length
            loc = self.loc[n]
            if loc < sumH_start:
                for i in range(len(self.loc)):
                    try:
                        if self.loc[n+i+1] > sumH_start:
                            n = n+i
                            break
                    except:
                        pass

            end = list(filter(lambda v: v.member != self if True else False, self.parent.n2.parents.copy()))
            other_h_end = end[0].cross_section.H
            sumH_end = (other_h_end + self.section.get_H(0)) / self.length
            if loc > 1 - sumH_end:
                for i in range(len(self.loc)):
                    try:
                        if self.loc[n-i-1] < 1 - sumH_end:
                            n = n-i
                            break
                    except:
                        pass

        UVy = round(abs(self.tau_d(n)[0]) / self.fvd, 5)
        UVz = round(abs(self.tau_d(n)[1]) / self.fvd, 5)
        r = max(UVy, UVz)
        return r
#        return print('Shear force resistance capacity utilized until', UV * 100, '%')

    def check_bending_moment(self, n):
        # Bending moment utilization ratio: TODO en:ssä f_md:lle arvot 2 eri suuntaan    SFS-EN 1995-1-1 § 6.1.6
        UM1 = round(abs(self.sigma_md(n)[0] / self.fmd + k_m * self.sigma_md(n)[1] / self.fmd), 5)
        UM2 = round(abs(self.sigma_md(n)[1] / self.fmd + k_m * self.sigma_md(n)[0] / self.fmd), 5)
        r = max(UM1, UM2)
        return r

    def check_torsion(self, n):
        # Torsion utilization ratio:                                SFS-EN 1995-1-1 § 6.1.8
        UTx = round(self.tau_tord(n) / (self.k_shape(n) * self.fvd), 5)
        r = UTx
        return r

    def check_inclined_compression(self):
        # TODO vinossa vaikuttavat puristusjännitykset              SFS-EN 1995-1-1 § 6.2.2
        pass
    # def alpha(self):
    #     return (sauvan.loppupisteen.korkeus - sauvan.alkupisteen.korkeus) / sauvan.pituus

    # def km_alpha(self):
    #     if self.ned >= 0:
    #         return np.sqrt(1 + (self.fmd * tan(alpha) / 0.75 / self.fvd) ** 2 + (self.fmd * (tan(alpha) ** 2) / self.ft90d) ** 2) ** -1
    #     else:
    #         return np.sqrt(1 + (self.fmd * tan(alpha) / 1.5 / self.fvd) ** 2 + (self.fmd * (tan(alpha) ** 2) / self.ft90d) ** 2) ** -1

    def check_bending_tension(self, n):
        #                                                           SFS-EN 1995-1-1 § 6.2.3
        UM1 = round(self.sigma_t0d(n) / self.ft0d + self.sigma_md(n)[0] / self.fmd + k_m * self.sigma_md(n)[1] / self.fmd, 5)
        UM2 = round(self.sigma_t0d(n) / self.ft0d + k_m * self.sigma_md(n)[0] / self.fmd + self.sigma_md(n)[1] / self.fmd, 5)
        # TODO pitääkö ei oletetut normaalivoiman suunnat rajata pois?
        r = max(UM1, UM2)
        return r

    def check_bending_compression(self, n):
        #                                                           SFS-EN 1995-1-1 § 6.2.4
        UM1 = round((self.sigma_c0d(n) / self.fc0d) ** 2 + self.sigma_md(n)[0] / self.fmd + k_m * self.sigma_md(n)[1] / self.fmd, 5)
        UM2 = round((self.sigma_c0d(n) / self.fc0d) ** 2 + k_m * self.sigma_md(n)[0] / self.fmd + self.sigma_md(n)[1] / self.fmd, 5)
        # TODO pitääkö ei oletetut normaalivoiman suunnat rajata pois?
        r = max(UM1, UM2)
        return r

    def check_buckling(self):
        # Utilization ratio against member buckling:                SFS-EN 1995-1-1 § 6.3.2
        sigma_c0d_list = [[], []]
        sigma_myd_list = [[], []]
        sigma_mzd_list = [[], []]
        def form_max_sigma_list(supports, index):
            supp = supports.copy()
            supp.append(0.0)
            supp.append(1.0)
            supp.sort()
            for i in range(len(supp) - 1):
                max_sigma_c0d_list = []
                max_sigma_myd_list = []
                max_sigma_mzd_list = []
                for loc in self.loc:
                    if supp[i] <= loc <= supp[i+1]:
                        max_sigma_c0d_list.append(self.sigma_c0d(self.loc.index(loc)))
                        max_sigma_myd_list.append(self.sigma_md(self.loc.index(loc))[0])
                        max_sigma_mzd_list.append(self.sigma_md(self.loc.index(loc))[1])
                max_sigma_c0d = min(max_sigma_c0d_list)
                max_sigma_myd = max(max_sigma_myd_list, key=abs)
                max_sigma_mzd = max(max_sigma_mzd_list, key=abs)

                sigma_c0d_list[index].append(max_sigma_c0d)
                sigma_myd_list[index].append(max_sigma_myd)
                sigma_mzd_list[index].append(max_sigma_mzd)

        form_max_sigma_list(self.lateral_support_y, 0)
        form_max_sigma_list(self.lateral_support_z, 1)

        sc0 = np.array(sigma_c0d_list[0])
        sc1 = np.array(sigma_c0d_list[1])

        my0 = np.array(sigma_myd_list[0])
        my1 = np.array(sigma_myd_list[1])
        mz0 = np.array(sigma_mzd_list[0])
        mz1 = np.array(sigma_mzd_list[1])

        UB = [abs(sc0 / (self.k_c()[0] * self.fc0d) + my0 / self.fmd + k_m * mz0 / self.fmd),
              abs(sc1 / (self.k_c()[1] * self.fc0d) + k_m * my1 / self.fmd + mz1 / self.fmd)]
        UB_fixed = [UB[0] * (sc0 < 0), UB[1] * (sc1 < 0)]

        return [max(UB_fixed[0]), max(UB_fixed[1])]

    def check_LT_buckling(self):
        # Utilization ratio against member lateral buckling:        SFS-EN 1995-1-1 § 6.3.3 (3) and (6)
        # TODO vaihtuvan poikkileikkausken tarkistelu
        # mitä k_criittisen arvoa käytetään vaihtuvalle poikkileikkaukselle
        ULT_list = []
        supp = self.lateral_support_z.copy()
        supp.append(0.0)
        supp.append(1.0)
        supp.sort()
        for i in range(len(supp) - 1):
            # max_sigma_myd_list = []
            # loc_list = []
            for loc in self.loc:
                if supp[i] <= loc <= supp[i+1]:
                    # max_sigma_myd_list.append(self.sigma_md(self.loc.index(loc))[0])
                    # loc_list.append(loc)
                    if max(self.myed, key=abs) != 0 and abs(max(self.ned, key=abs)) <= 0.1:
                        ULT_list.append(self.sigma_md(self.loc.index(loc))[0] / (self.k_crit(i, 1, loc) * self.fmd))
                    else:
                        ULT_list.append(
                            ((self.sigma_md(self.loc.index(loc))[0] / (self.k_crit(i, 1, loc) * self.fmd)) ** 2) +
                            self.sigma_c0d(self.loc.index(loc)) / (self.k_c()[1][i] * self.fc0d))

            # max_sigma_myd = max(max_sigma_myd_list, key=abs)
            # print(loc_list)
            # max_loc = loc_list[max_sigma_myd_list.index(max_sigma_myd)]
            # if max(self.myed, key=abs) != 0 and abs(max(self.ned, key=abs)) <= 0.1:
            #     ULT_list.append(max_sigma_myd / (self.k_crit(i, 1, max_loc) * self.fmd))
            # else:
            #     for loc in loc_list:
            #         ULT_list.append((self.sigma_md(self.loc.index(loc))[0] / (self.k_crit(i, 1, loc) * self.fmd) ** 2) +
            #                         self.sigma_c0d(self.loc.index(loc)) / (self.k_c()[1][i] * self.fc0d))
                # TODO sigma_cod positiivinen tai negatiivinen
        i = ULT_list.index(max(ULT_list))

        return max(ULT_list)

    def check_apex_bending_stress(self):
        # EN 1995-1-1 (6.41)
        sigma = self.k_l() * self.M_apd() / self.W_ap()[0]
        return sigma / (self.k_r() * self.fmd)

    def check_apex_perpendicular_tension(self):
        # EN 1995-1-1 (6.50)
        # print(f'sigma t90d {self.sigma_t90d()}')
        # print(f'k_dis {self.k_dis()}')
        # print(f'k_vol {self.k_vol()}')
        return self.sigma_t90d() / (self.k_dis() * self.k_vol() * self.ft90d)

    def check_apex_shear_perpendicular_tension_combined(self):
        # EN 1995-1-1 (6.53)
        # taking maximum shear force instead of reduced shear force
        # TODO palataan myöhemmin

        if isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki)):
            closesttodownx1 = None
            closesttodownx2 = None
            for loc in self.loc:
                if closesttodownx1 is None:
                    closesttodownx1 = loc
                    continue
                if closesttodownx2 is None:
                    closesttodownx2 = loc
                    continue
                if abs((self.section.dwn_x1 / self.length) - loc) < abs(closesttodownx1 - loc):
                    closesttodownx1 = loc
                if abs(self.section.dwn_x2 / self.length - loc) < abs(closesttodownx2 - loc):
                    closesttodownx2 = loc

            tau_d_dw_x1 = self.tau_d(self.loc.index(closesttodownx1))
            tau_d_dw_x2 = self.tau_d(self.loc.index(closesttodownx2))
            # TODO miksi metku laskee leikkausjännityksen z-suunnassa?
            tau = abs(min(tau_d_dw_x1[1], tau_d_dw_x2[1], key=abs))

        elif isinstance(self.section, HarjaPalkki):
            closesttoleft = None
            closesttoright = None
            for loc in self.loc:
                if closesttoleft is None:
                    closesttoleft = loc
                    continue
                if closesttoright is None:
                    closesttoright = loc
                    continue
                if abs(self.section.xap - (self.section.hap / 2 / self.length) - loc) < abs(closesttoleft - loc):
                    closesttoleft = loc
                if abs(self.section.xap + (self.section.hap / 2 / self.length) - loc) < abs(closesttoright - loc):
                    closesttoright = loc

            tau_d_dw_x1 = self.tau_d(self.loc.index(closesttoleft))
            tau_d_dw_x2 = self.tau_d(self.loc.index(closesttoright))
            # TODO miksi metku laskee leikkausjännityksen z-suunnassa?
            tau = abs(min(tau_d_dw_x1[1], tau_d_dw_x2[1], key=abs))

        return (tau / self.fvd) + self.sigma_t90d() / (self.k_dis() * self.k_vol() * self.ft90d)


    def add_section(self, ned=0.0, myed=0.0, mzed=0.0, vyed=0.0, vzed=0.0,
                    ted=0.0, loc=0.0):
        """ Adds new section with internal forces and their location """

        self.ned.append(ned)
        self.myed.append(myed)
        self.mzed.append(mzed)
        self.vyed.append(vyed)
        self.vzed.append(vzed)
        self.ted.append(ted)
        self.loc.append(loc)

    def clear_sections(self):
        self.ned.clear()
        self.myed.clear()
        self.mzed.clear()
        self.vyed.clear()
        self.vzed.clear()
        self.ted.clear()
        self.loc.clear()

    def check_section(self, n=0, verb=False):
        """ Verify resistance of section 'n'
            returns a list of utilization ratios
        """

        r = self.section_resistance(n=n, verb=verb)
        return r


    def section_resistance(self, n=0, axis='y', return_list=True, verb=False):
        """ Calculates resistance of cross-section
            Checks the following:
                Axial force
                Shear force
                Bending moment
                Torsion
                Interaction of tension and bending moment
                Interaction of compression and bending moment
            output: r .. maximum value of different utilization ratios
        """
        if axis == 'y':
            idx = 0
        else:
            idx = 1

        UN = self.check_normal_force(n)
        UNcp = self.check_perpendicular_compression()
        UV = self.check_shear_force(n)
        UM = self.check_bending_moment(n)
        UT = self.check_torsion(n)
        UBT = self.check_bending_tension(n)
        UBC = self.check_bending_compression(n)
        # UB = self.check_buckling()
        # ULT = self.check_LT_buckling()


        # return print('Normaalivoiman käyttöaste =', UN * 100, '%')

        if return_list:
            return [UN, UV, UM, UT, UBT, UBC, UNcp]
        # else:
        #     return max([UN, UV, UM, UT, UBT, UBC])

    #
    #
    #     if self.C < 3:
    #         MNRd = self.moment_axial_force_interact(UN, self.MRd[idx])
    #         if MNRd > 0.0:
    #             UMN = abs(self.Med) / MNRd
    #         else:
    #             UMN = UN + UM
    #             # UMN = INFEASIBLE
    #     else:
    #         UMN = UN + UM
    #
    #     if return_list:
    #         r = [UN, UV, UM, UMN]
    #     else:
    #         r = max([UN, UV, UM, UMN])
    #
    #     if verb:
    #         print("Cross-section design: " + self.__repr__() + " S" + str(self.fy))
    #         print("NEd = {0:4.2f} kN; NRd = {1:4.2f} kN => UN = {2:4.2f}".format(self.Ned * 1e-3, self.NRd * 1e-3, UN))
    #         print("VEd = {0:4.2f} kN; VRd = {1:4.2f} kN => UV = {2:4.2f}".format(self.Ved * 1e-3, self.VRd * 1e-3, UV))
    #         print(
    #             "MEd = {0:4.2f} kNm; MRd = {1:4.2f} kNm => UM = {2:4.2f}".format(self.Med * 1e-6, self.MRd[idx] * 1e-6,
    #                                                                              UM))
    #
    #     return r

    def check_cross_section(self):
        """ Checks if cross-section is strong enough.
            Gives boolean value for self.is_strong_enough
            Saves list of stress ratios to self.r
        """

        self.r = np.zeros_like(self.r) # [0, 1, 2, .. 8]
        # Cross-sectional stress ratios in member's nodes' locations
        self.r[:7] = self.check_sections()
        # Buckling about y - and z - axis
        buckling_r = self.check_buckling()
        self.r[7] = buckling_r[0]
        self.r[8] = buckling_r[1]
        # Lateral-Torsional buckling
        self.r[9] = self.check_LT_buckling()

        if max(self.r) > 1.0:
            self.is_strong_enough = False
        else:
            self.is_strong_enough = True

    def check_sections(self, class1or2=True, verb=False):
        """ Verify resistance of all sections """
        r = np.zeros_like(self.check_section())
        for n in range(self.nsect()):
            r = np.vstack((r, self.check_section(n)))
            # r.append(self.check_section(n))
        return np.max(r, axis=0)


    def nsect(self):
        """ Number of sections """
        return len(self.ned)

    def add_neighbours_support(self):
        for node in self.nodes.values():
            if node is self.parent.n1 or node is self.parent.n2:
                continue
            other_members = list(filter(lambda v: v.member != self if True else False, node.parents.copy()))
            if len(other_members) == 0:
                continue

            dist_to_start = np.linalg.norm(self.parent.coordinates[0] - node.coord)
            local_coord = dist_to_start / self.length

            self.lateral_support_y.append(local_coord)
            self.lateral_support_z.append(local_coord)

        self.calc_l_c()

    def calc_l_c(self):
        self.lateral_support_y.sort()
        self.lateral_support_z.sort()

        Lcy = []
        Lcz = []
        for i in range(len(self.lateral_support_y) + 1):
            if len(self.lateral_support_y) == 0:
                Lcy.append(self.length)
            elif i == 0:
                Lcy.append(self.length * self.lateral_support_y[i])
            elif i == len(self.lateral_support_y):
                Lcy.append(self.length * (1 - self.lateral_support_y[i - 1]))
            else:
                Lcy.append(self.length * (self.lateral_support_y[i] - self.lateral_support_y[i - 1]))
        for i in range(len(self.lateral_support_z) + 1):
            if len(self.lateral_support_z) == 0:
                Lcz.append(self.length)
            elif i == 0:
                Lcz.append(self.length * self.lateral_support_z[i])
            elif i == len(self.lateral_support_z):
                Lcz.append(self.length * (1 - self.lateral_support_z[i - 1]))
            else:
                Lcz.append(self.length * (self.lateral_support_z[i] - self.lateral_support_z[i - 1]))
        self.L_c = np.asarray([np.asarray(Lcy),
                               np.asarray(Lcz)], dtype=object)

    def __repr__(self):
        return self.parent.__repr__()

    # TODO lisää femiin sectio
    # def add_section(self, fem_model):
    #     """ Units: mm A and I_y are in mm**2 and mm**4, respectively"""
    #     s = fem.BeamSection(self.cross_section.A, self.cross_section.I[0])
    #     fem_model.add_section(s)

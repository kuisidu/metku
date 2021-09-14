"""
Timber member

@author Viktor Haimi
"""
from math import log, tan, pi, sqrt, isclose, sin, cos, radians
import numpy as np

try:
    import frame2d
    from materials.timber_data import Timber, T
    from eurocodes.en1995.en1995_1_1 import kmod, kdef, k_to_d, get_gammaM
    from eurocodes.en1995.constants import k_m
    from sections.timber.timber_section import *
    from framefem.framefem import FEMNode
    from framefem.elements import EBBeam, EBSemiRigidBeam, Rod
    from sections.steel.ISection import HEA
except:
    import metku.frame2d as frame2d
    from metku.materials.timber_data import Timber, T
    from metku.eurocodes.en1995.en1995_1_1 import kmod, kdef, k_to_d, get_gammaM
    from metku.eurocodes.en1995.constants import k_m
    from metku.sections.timber.timber_section import *
    from metku.framefem.framefem import FEMNode
    from metku.framefem.elements import EBBeam, EBSemiRigidBeam, Rod
    from metku.sections.steel.ISection import HEA

PREC = 5

import sys

if sys.version_info[0] == 3 and sys.version_info[1] < 9:
    from typing import List, Tuple, Dict
    list_array = List[np.ndarray]
    #list_float_FEMNode = List[[float,FEMNode]]
    list_float_FEMNode = Tuple[float,FEMNode]
    list_out = List[float]
    llist_out = List[List[float]]
    lllist_out = List[List[List[float]]]
    
    dict_str_float = Dict[str, float]
    
    tuple_out = Tuple
else:
    list_array = list[np.ndarray]
    list_float_FEMNode = list[float,FEMNode]
    list_out = list[float]
    llist_out = list[list[float]]
    lllist_out = list[list[list[float]]]
    
    dict_str_float = dict[str: float]
    
    tuple_out = tuple

class TimberMember:
    def __init__(self, profile: (TimberSection, TaperedSection), material: Timber, length: float, ldc: str, sc: int,
                 mtype: str = 'beam', varnished: bool = False, Sj1: int = np.inf, Sj2: int = np.inf, nodes: dict = None,
                 parent=None, lateral_support_y: list = None, lateral_support_z: list = None,
                 edge_load: str = 'compression', beta: (list, float) = None, k_vol: float=None):
        """
        Initialization method for timber member

        @param profile: poikkileikkaus luokkaa TimberSection
        @param material: sauvan materiaali muodossa: C14..D70 / GL24c..GL32h
        @param length: sauvan pituus [mm]
        @param ldc: Load duration class. options: 'perm' = permanent, 'lt' = long-term, 'mt' = medium-term,
                    'st' = short-term, 'inst' = instantaneous
        @param sc: Service class.
        @param mtype: member type: 'column', 'beam'
        @param varnished: lakkaus
        @param Sj1: Sauvan alkupään jäykkyys
        @param Sj2: Sauvan loppupään jäykkyys
        @param nodes:
        @param parent: TimberFrameMember
        @param lateral_support_y: locations of the supports in global y direction (This has nothing to do with FEM)
        @param lateral_support_z: locations of the supports in z direction (This has nothing to do with FEM)
        @param edge_load: 'compression' or 'tension'
        @param beta: Eulerin nurjahduskerroin
        @param k_vol: pakotettu k_vol
        """

        self.material = material
        self.ldc = ldc
        self.sc = sc
        self.varnished = varnished
        self.length = length

        self.type = mtype
        self.Sj1 = Sj1
        self.Sj2 = Sj2
        self.beta = beta
        self.perm_k_vol = k_vol

        self.ned = []
        self.myed = []
        self.mzed = []
        self.vyed = []
        self.vzed = []
        self.ted = []
        self.loc = []
        self.utilization = []

        self.nodes = nodes
        self.parent = parent

        self.lateral_support_y = lateral_support_y
        self.lateral_support_z = lateral_support_z
        self.L_c = np.asarray([np.asarray([]),
                               np.asarray([])], dtype=object)

        self.edge_load = edge_load

        self.kmod = kmod(self.material.type, self.sc, self.ldc)

        self.section = profile

        self.R = 0

        self.gammaM = get_gammaM(self.material.type)

    @property
    def k_L(self):
        # SFS-EN 1995-1-1 (3.4)
        k_L = 1.0
        if self.material.type == 'lvl':
            if self.NEd > 0 and self.length != 3000:
                k_L = min(1.1, (3000 / self.length) ** (self.material.s / 2))
        return k_L

    def get_kdef(self):
        return kdef(self.material.type, self.sc)

    @property
    def R(self) -> (int, float):
        return self.section.R

    @R.setter
    def R(self, val: (int, float)):
        self.section.R = val

    @property
    def fmd(self) -> float:
        if self.material.type == 'lvl':
            if self.material.direction == 'flat':
                return self.section.k_h * k_to_d(self.material.fm0flatk, self.kmod, self.gammaM)
        return self.section.k_h * k_to_d(self.material.fmk, self.kmod, self.gammaM)

    # ilman k_h kerrointa
    @property
    def unmodified_fmd(self) -> float:
        return k_to_d(self.material.fmk, self.kmod, self.gammaM)

    @property
    def ft0d(self) -> float:
        if self.material.type == 'solid_timber' or self.material.type == 'glt':
            return self.section.k_h * k_to_d(self.material.ft0k, self.kmod, self.gammaM)
        elif self.material.type == 'lvl':
            return self.k_L * k_to_d(self.material.ft0k, self.kmod, self.gammaM)

    # ilman k_h kerrointa
    @property
    def unmodified_ft0d(self) -> float:
        return k_to_d(self.material.ft0k, self.kmod, self.gammaM)

    @property
    def ft90d(self) -> float:
        if self.material.type == 'lvl':
            if self.material.direction == 'edge':
                return k_to_d(self.material.ft90edgek, self.kmod, self.gammaM)
            raise Exception('lvl does not have ft90d, try ft90edged')
        return k_to_d(self.material.ft90k, self.kmod, self.gammaM)

    @property
    def fc0d(self) -> float:
        return k_to_d(self.material.fc0k, self.kmod, self.gammaM)

    @property
    def fc90d(self) -> float:
        if self.material.type == 'lvl':
            if self.material.direction == 'flat':
                return k_to_d(self.material.fc90flatk, self.kmod, self.gammaM)
            return k_to_d(self.material.fc90edgek, self.kmod, self.gammaM)
        return k_to_d(self.material.fc90k, self.kmod, self.gammaM)

    @property
    def fvd(self) -> float:
        return k_to_d(self.material.fvk, self.kmod, self.gammaM)

    @property
    def ft90edged(self) -> float:
        if self.material.type == 'glt' or self.material.type == 'solid_timber':
            raise Exception(f'{self.material.type} does not have ft90edged, try ft90d')
        return k_to_d(self.material.ft90edgek, self.kmod, self.gammaM)

    @property
    def fr0d(self) -> float:
        if self.material.type == 'glt' or self.material.type == 'solid_timber':
            raise Exception(f'{self.material.type} does not have fr0d')
        return k_to_d(self.material.fr0k, self.kmod, self.gammaM)

    @property
    def MEdY(self) -> float:
        return max(abs(np.array(self.myed)))

    @property
    def MEdZ(self) -> float:
        return max(abs(np.array(self.mzed)))

    @property
    def NEd(self) -> float:
        """ Maximum axial force (absolute value) """
        return max(abs(np.array(self.ned)))

    @property
    def A(self) -> (int, float):
        return self.section.get_A(0)

    def get_A(self, at: float) -> (int, float):
        return self.section.get_A(at)

    @property
    def Iy(self) -> float:
        return self.section.I[0]

    @Iy.setter
    def Iy(self, val: float):
        self.section.I[0] = val

    @property
    def Iz(self) -> float:
        return self.section.I[1]

    @Iz.setter
    def Iz(self, val: float):
        self.section.I[1] = val

    @property
    def E0mean(self) -> (int, float):
        return self.material.E0mean

    @E0mean.setter
    def E0mean(self, val: (int, float)):
        self.material.E0mean = val

    @property
    def h(self) -> float:
        """
        Property, returns member's cross-section's height
        """
        return self.section.get_H(0)

    @h.setter
    def h(self, val: (int, float)):
        """
        Sets cross-section's height to given value.
        """
        self.section.H = val

    @property
    def b(self) -> float:
        """
        Property, returns member's cross-section's breadth
        """
        return self.section.B

    @b.setter
    def b(self, val: (int, float)):
        """
        Sets cross-section's height to given value.
        """
        self.section.B = val

    def sigma_t0d(self, n: int) -> float:
        """
        EN 1995-1-1 6.1.2
        @param n: indeksi
        @return:
        """
        return self.ned[n] / self.section.get_A(self.loc[n])

    def sigma_c0d(self, n: int) -> float:
        """
        EN 1995-1-1 6.1.4
        @param n: indeksi
        @return:
        """
        return self.ned[n] / self.section.get_A(self.loc[n])

    def sigma_m_alpha_d(self) -> float:
        # SFS 1995-1-1 (6.37)
        list = [self.myed[self.loc.index(i)] / self.section.get_Wel(i)[0] for i in self.loc]
        min_val = min(list)
        max_val = max(list)
        abs_max = max(abs(min_val), abs(max_val))
        if abs_max == abs(min_val):
            return min_val
        else:
            return max_val

    def sigma_t90d(self, kaava6_55:bool=False) -> float:
        if kaava6_55:
            p = 0
            for load in self.parent.loads:
                if load.direction == 'y' and load.values[0] == load.values[1] and self.edge_load == 'compression':
                    p = load.values[0]
            # SFS 1995-1-1 (6.55)
            return self.k_p() * self.M_apd() / self.W_ap()[0] - 0.6 * (p / self.section.B)
        # SFS 1995-1-1 (6.54)
        return self.k_p() * self.M_apd() / self.W_ap()[0]

    def M_apd(self) -> float:
        """
        harjavyöhykkeen momentti
        jos solmu löytyy keskeltä niin käytetään sitä, muuten lasketaan kahden viereisen solmun keskiarvo
        @return:
        """
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

    def W_ap(self) -> np.ndarray:
        """
        harjavyöhykkeen kimmoinen taivutusvastus
        @return:
        """
        return self.section.get_Wel(self.section.xap)

    def k_r(self) -> float:
        # SFS 1995-1-1 (6.49)
        # TODO Mahapalkki kaarevien ryhmään, vaatii mahapalkin uudelleen määrittämisen
        if isinstance(self.section, (HarjaPalkki, PulpettiPalkki, MahaPalkki)):
            return 1
        elif isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki)):
            if self.section.r_in * 1000 / self.section.t >= 240:
                return 1
            else:
                return 0.76 + 0.001 * (self.section.r_in * 1000 / self.section.t)
        else:
            raise NotImplementedError

    def V(self) -> float:
        # SFS 1995-1-1 (Kuva 6.9)
        # V <= 2/3 * koko palkin tilavuus
        if isinstance(self.section, (HarjaPalkki, MahaPalkki)):
            return self.section.B * np.power(self.section.hap, 2)
        elif isinstance(self.section, KaarevaHarjaPalkki):

            # l1 = self.section.h0 * cos(self.alfa())
            # l2 = (self.section.up_y1 - self.section.dwn_y1) / cos(self.alfa())
            # l3 = self.section.up_x1 / cos(self.alfa())
            # V_leg = self.section.B * l3 * (l1 + l2) * 0.5
            # r = self.V_tot() - 2 * V_leg

            from scipy import integrate
            up = integrate.quad(self.section.f1, 0, self.section.up_x1)
            dwn = integrate.quad(self.section.f2, 0, self.section.up_x1)
            intlegs = (up[0] - dwn[0]) * 2 * self.section.B

            r = self.V_tot() - intlegs

            if r < 0:
                print(f'harja-alueen tilavuudessa laskuvirhe, epätarkkojen oletusten vuoksi')
                # if intlegs < 0:
                #     print(self.section.intf1)
                #     print(self.section.intf2)
                #     print(self.section.intf3)
                #     print(self.section.intf4)
                #     print(self.section.intf6)
            return r
        elif isinstance(self.section, KaarevaPalkki):
            return self.section.B * self.section.hap * (self.section.r_in + self.section.hap * 0.5) * ((self.alfa() + self.section.beta) * 0.5)
        else:
            raise NotImplementedError

    def V_tot(self) -> float:
        if isinstance(self.section, (HarjaPalkki, MahaPalkki)):
            return self.section.B * self.length * (self.section.hap - self.length * np.tan(self.alfa()) / 4)
        elif isinstance(self.section, (KaarevaPalkki, KaarevaHarjaPalkki)):
            return self.section.A_face() * self.section.B
        elif isinstance(self.section, PulpettiPalkki):
            return (self.section.H0 + self.section.Hap) / 2 * self.length * self.section.B
        else:
            return self.section.A * self.length

    def k_vol(self) -> float:
        # SFS 1995-1-1 (6.51)
        if self.perm_k_vol is not None:
            return self.perm_k_vol
        if self.material.type == 'solid timber':
            return 1.0
        v = self.V()
        if v > (2 / 3) * self.V_tot():
            v = (2 / 3) * self.V_tot()
        k_vol = np.power(1e7 / v, 0.2)
        return k_vol

    def k_dis(self) -> float:
        # SFS 1995-1-1 (6.52)
        if isinstance(self.section, (HarjaPalkki, PulpettiPalkki, MahaPalkki)):
            return 1.4
        elif isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki)):
            return 1.7

    def k_p(self) -> float:
        # SFS 1995-1-1 (6.56)
        if isinstance(self.section, (HarjaPalkki, MahaPalkki)):
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

    def k_l(self) -> float:
        # SFS 1995-1-1 (6.43)
        k1 = 1 + 1.4 * np.tan(self.alfa()) + 5.4 * np.tan(self.alfa()) ** 2
        if isinstance(self.section, (HarjaPalkki, MahaPalkki)):
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

    def k_m_alpha(self) -> float:
        # SFS 1995-1-1 (6.39) ja (6.40)
        sigma_malphad = self.sigma_m_alpha_d()
        if not self.section.flipped:
            if sigma_malphad >= 0:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (4 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1 / 2)
            else:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (2 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1 / 2)
        else:
            if sigma_malphad <= 0:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (4 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1 / 2)
            else:
                kma = np.power(1 + np.power(np.tan(self.alfa()) * (2 / 3) * self.fmd / self.fvd, 2) +
                               np.power(np.power(np.tan(self.alfa()), 2) * self.fmd / self.fc90d, 2), -1 / 2)
        return kma

    def alfa(self) -> float:
        """
        Palkin yläpinnan kaltevuuskulma
        @return: radiaaneina
        """
        if isinstance(self.section, PulpettiPalkki):
            a = np.arctan((self.section.hap - self.section.h0) / self.length)
        elif isinstance(self.section, HarjaPalkki):
            a = np.arctan((self.section.hap - self.section.h0) / (self.length / 2))
        elif isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki, MahaPalkki)):
            a = self.section.alpha
        else:
            a = 0
        return a

    def sigma_c90d(self, node: FEMNode) -> (list_out, float):
        """
        SFS 1995-1-1 6.1.5
        Rajoitusehdoissa pitää olla määritettynä H > B, muuten ohjelma asettaa pilarin "väärinpäin"

        @param node: node where to check sigma c90d
        @return:     0 if there is no other member on this node,
                     float if there is only one other member on this node
                     list of floats if there is more than one members on this node
        """

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
                # Rajoitusehdoissa pitää olla määritettynä H > B, muuten ohjelma asettaa pilarin "väärinpäin"
                lc90ef = om.cross_section.H + min(30, 0, om.cross_section.H) + min(30, om.cross_section.H)
            else:
                lc90ef = om.cross_section.H + min(30, om.cross_section.H) + min(30, om.cross_section.H)

            # SFS-EN 1995-1-1 § 6.1.5
            Aef = lc90ef * min(om.cross_section.B, self.section.B)
            sigma = om.member.ned[other_idx] / Aef
            return sigma

        if len(other_members) == 1:
            om = other_members[0]
            return get_sigma(om)
        else:
            return [get_sigma(om) for om in other_members]

    # def kc_perpendicular(self, l: (int, float), a: (int, float), l1: (int, float)) -> float:
    #     """
    #     Tukipainekerroin: RIL 205-1-2009 B.5.2a,
    #     @param l:   sen sauvan kosketuspinnan pituus, joka rasitetaan syysuuntaa vastaan esim. pilarin päällä oleva palkki
    #     @param a:   kosketuspinnan ulkopuolinen, puristetun palkin päätyyn asti ulottuva mitta, tai seuraavaan puristuspintaan
    #     @param l1:  puristuspintojen lähin välimatka toisiinsa nähden
    #     @return:    tukipainekerroin
    #     """
    #     self.update_kc90(l1, l)
    #     lc90ef = l + min(30, a, l) + min(30, l, l1 / 2)
    #     return lc90ef / l * self.k_c90

    def update_kc90(self, l1: (int, float), l: (int, float)):
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

    def sigma_md(self, n: int) -> list_out:
        """
        EN-1995-1-1 6.1.6
        @param n: indeksi
        @return:
        """
        Wel = self.section.get_Wel(self.loc[n])
        r = [self.myed[n] / Wel[0],
             self.mzed[n] / Wel[1]]
        return r

    def get_max_sigma_per_segment(self) -> lllist_out:
        """
        hakee suurimman sigma myd ja mzd per sivuttaistuetun palkin jokaiselle tukemattomalle palkin välille.
        @return:
        """
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
                if supp_y[i] <= loc <= supp_y[i + 1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iy.append(max(sig_md_in_seg, key=lambda x: abs(x[0])))

        for i in range(len(supp_z) - 1):
            sig_md_in_seg = []
            locs = []
            for loc in self.loc:
                if supp_z[i] <= loc <= supp_z[i + 1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iz.append(max(sig_md_in_seg, key=lambda x: abs(x[0])))
        return [iy, iz]

    def tau_d(self, n: int) -> list_out:
        """
        Poikkileikkauksessa vaikuttava työntövoima muutetaan leikkausjännitykseski Jourawskin kaavaa hyväksikäyttäen
        @return: palauttaa suorakaidepoikkileikkauksen leikkausjännityksen arvon, ottaa huomioon k_cr kertoimen
        """
        k_cr = 1.0
        if self.material.type == 'solid_timber' or self.material.type == 'glt' or self.sc == 1:
            if not self.varnished:
                k_cr = 0.67

        return [(3 / 2 * self.vyed[n]) / (self.section.get_H(self.loc[n]) * self.section.B * k_cr),
                (3 / 2 * self.vzed[n]) / (self.section.get_H(self.loc[n]) * self.section.B * k_cr)]

    def tau_tord(self, n: int) -> float:
        """
        EN-1995-1-1 6-1-8
        Alpha and beta are based on the theoretical tables given in: T. Salmi 2000, Lujuusopin perusteet, p. 258
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
            beta = -0.0007 * section_height_to_width_ratio ** 2 + 0.0193 * section_height_to_width_ratio + 0.1923
        Wv = beta * self.section.get_H(self.loc[n]) * self.section.B ** 2
        return self.ted[n] / Wv

    def I_tor(self, loc: float) -> float:
        """
        @return: torsional constant for rectangular cross sections
        """
        ratio = self.section.get_H(loc) / self.section.B
        if ratio < 1:
            ratio = 1 / ratio
        if 1 < ratio < 1.05:
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
            alpha = -0.0022 * ratio ** 2 + 0.425 * ratio + 0.1116

        return alpha * self.section.get_H(loc) * self.section.B ** 3

    def k_shape(self, n: int) -> float:
        """
        SFS EN 1995-1-1 §6.1.8(6.15)
        @return: cross section factor for rectangular cross sections (circular sections not included)
        """
        return min(1.3, 1 + 0.05 * self.section.get_H(self.loc[n]) / self.section.B)

    def beta_c(self) -> float:
        """
        EN 1995-1-1 6.3.2 (6.29)
        @return:
        """
        if self.material.type == 'solid_timber':
            return 0.2
        return 0.1

    def radius_of_gyration(self, get_locs: bool = False) -> (np.ndarray, list):
        """
        jäyhyyssäteen arvo etsitään suurimman taivutusjännityksen alueelta
        @param get_locs:
        @return:
        """
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
                if supp_y[i] <= loc <= supp_y[i + 1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iy.append(locs[sig_md_in_seg.index(max(sig_md_in_seg))])

        for i in range(len(supp_z) - 1):
            sig_md_in_seg = []
            locs = []
            for loc in self.loc:
                if supp_z[i] <= loc <= supp_z[i + 1]:
                    sig_md_in_seg.append(sig_md_list[self.loc.index(loc)])
                    locs.append(loc)
            iz.append(locs[sig_md_in_seg.index(max(sig_md_in_seg))])

        y = [sqrt(self.section.get_I(i)[0] / self.section.get_A(i)) for i in iy]
        z = [sqrt(self.section.get_I(i)[1] / self.section.get_A(i)) for i in iz]
        if get_locs:
            return [np.asarray([np.asarray(y),
                                np.asarray(z)], dtype=object), [iy, iz]]
        else:
            return np.asarray([np.asarray(y),
                               np.asarray(z)], dtype=object)

    def Ln(self) -> (llist_out, np.ndarray):
        """
        Eulerin teoriaan perustuva nurjahduspituus
        @return:
        """
        # Beta arvot saatu Liimapuukäsikirja Osa 3 sivu 3-35
        # TODO tarkista betan arvot palkille ja pilarille erikseen Lpkk3 (3-39)
        beta = [[], []]

        def form_beta(lc: list, index: int, beam=False):
            startsupp = None
            endsupp = None
            for supp in self.parent.frame2d.supports.values():
                if list(self.parent.n1.coord) == supp.coordinate:
                    startsupp = supp
                elif list(self.parent.n2.coord) == supp.coordinate:
                    endsupp = supp

            if len(lc) == 1:
                if beam:
                    if self.Sj1 < 1e3 and self.Sj2 < 1e3:
                        beta[index] = 1
                    elif (self.Sj1 >= 1e15 and self.Sj2 <= 1e3) or (self.Sj2 >= 1e15 and self.Sj1 <= 1e3):
                        beta[index] = 0.9
                    elif self.Sj1 >= 1e15 and self.Sj2 >= 1e15:
                        beta[index] = 0.8
                    else:
                        beta[index] = 0.85

                elif startsupp is not None and endsupp is not None:
                    if isinstance(startsupp, frame2d.FixedSupport) and isinstance(endsupp, frame2d.FixedSupport):
                        beta[index] = 0.7
                    elif isinstance(startsupp, frame2d.XYHingedSupport) and isinstance(endsupp, frame2d.XYHingedSupport):
                        beta[index] = 1
                    elif (isinstance(startsupp, frame2d.FixedSupport) or isinstance(endsupp, frame2d.FixedSupport) and
                          isinstance(startsupp, frame2d.XYHingedSupport) or isinstance(endsupp, frame2d.XYHingedSupport)):
                        beta[index] = 0.85
                    elif isinstance(startsupp, frame2d.YHingedSupport) or isinstance(endsupp, frame2d.YHingedSupport):
                        if isinstance(startsupp, frame2d.FixedSupport) or isinstance(endsupp, frame2d.FixedSupport):
                            beta[index] = 2.2
                        elif isinstance(startsupp, frame2d.XYHingedSupport) or isinstance(endsupp, frame2d.XYHingedSupport):
                            beta[index] = 2
                    elif isinstance(startsupp, frame2d.XHingedSupport) or isinstance(endsupp, frame2d.XHingedSupport):
                        if isinstance(startsupp, frame2d.FixedSupport) or isinstance(endsupp, frame2d.FixedSupport):
                            beta[index] = 0.85
                        elif isinstance(startsupp, frame2d.XYHingedSupport) or isinstance(endsupp, frame2d.XYHingedSupport):
                            beta[index] = 1
                elif isinstance(startsupp, frame2d.FixedSupport) or isinstance(endsupp, frame2d.FixedSupport):
                    beta[index] = 2.2
                elif isinstance(startsupp, frame2d.XYHingedSupport) or isinstance(endsupp, frame2d.XYHingedSupport):
                    beta[index] = 2

            else:
                # jos sauvat on jaettu useampaan osaan/segmenttiin
                for i in range(len(lc)):
                    # palkin segmentin z suunnassa on kiinteä beta=1 nurjahduskerroin
                    if beam and index == 1:
                        beta[index].append(1)
                    # sauvan alkupäästä ensimmäinen segmentti beta=0.8, sitä seuraavat segmentit beta=1, jne.
                    # beta arvot annetaan tukien vapausasteiden ja sauvan päiden kiertymisjäykkyyksien mukaan y- ja z-suunnissa
                    elif i == 0:
                        if isinstance(startsupp, frame2d.FixedSupport) or self.Sj1 >= 1e15:
                            beta[index].append(0.8)
                        elif isinstance(startsupp, frame2d.XYHingedSupport) or self.Sj1 <= 1e3:
                            beta[index].append(1)
                        elif 1e3 < self.Sj1 < 1e15:
                            beta[index].append(0.85)
                        else:
                            beta[index].append(1)
                    elif i == len(lc) - 1:
                        # sauvan alkupäästä viimeisessä segmentissä beta=0.8, sitä seuraavat segmentit beta=1, jne.
                        if isinstance(endsupp, frame2d.FixedSupport) or self.Sj2 >= 1e15:
                            beta[index].append(0.8)
                        elif isinstance(endsupp, frame2d.XYHingedSupport) or self.Sj2 <= 1e3:
                            beta[index].append(1)
                        elif 1e3 < self.Sj2 < 1e15:
                            beta[index].append(0.85)
                        else:
                            beta[index].append(1)
                    else:
                        # sauvojen välisegmenteille lokaalissa z-suunnassa beta=1, mutta y-suunnassa beta=0.9
                        if index == 0:
                            beta[index].append(0.9)
                        else:
                            beta[index].append(1)

        if self.beta is not None:
            if isinstance(self.beta, (int, float)):
                return np.asarray([self.beta * self.L_c[0], self.beta * self.L_c[1]], dtype=object)
            elif isinstance(self.beta, list):
                return np.asarray([self.beta[0] * self.L_c[0], self.beta[1] * self.L_c[1]], dtype=object)
        if self.type == 'column':
            form_beta(self.L_c[0], 0)
            form_beta(self.L_c[1], 1)
            self.beta = beta
            return np.asarray([beta[0] * self.L_c[0], beta[1] * self.L_c[1]], dtype=object)
        elif self.type in ('rigid-beam', 'beam'):
            form_beta(self.L_c[0], 0, True)
            form_beta(self.L_c[1], 1, True)
            self.beta = beta
            return np.asarray([beta[0] * self.L_c[0], beta[1] * self.L_c[1]], dtype=object)

        if self.Sj1 >= 1e15 and self.Sj2 >= 1e15:
            beta = 0.7
        elif (self.Sj1 >= 1e15 and self.Sj2 <= 1e3) or (self.Sj2 >= 1e15 and self.Sj1 <= 1e3):
            beta = 0.85
        else:
            beta = 1
        return [beta * self.L_c[0], beta * self.L_c[1]]

    def lamda(self) -> list_array:
        """
        EN 1995-1-1 6.3.2
        @return:
        """
        i = self.radius_of_gyration()
        ln = self.Ln()
        return [ln[0] / i[0],
                ln[1] / i[1]]

    def lamda_rel(self) -> list_array:
        """
        EN 1995-1-1 6.3.2 (6.21, 6.22)
        @return:
        """
        lam = self.lamda()
        return [(lam[0] / pi) * sqrt(self.material.fc0k / self.material.E005),
                (lam[1] / pi) * sqrt(self.material.fc0k / self.material.E005)]

    def sigma(self) -> list_array:
        """
        sigma kriittinen
        @return:
        """
        lam = self.lamda()
        return [pi ** 2 * (self.material.E0mean / (lam[0] ** 2)),
                pi ** 2 * (self.material.E0mean / (lam[1] ** 2))]

    def k(self) -> list_array:
        """
        EN 1995-1-1 6.3.2 (6.27, 6.28)
        @return:
        """
        l_rel = self.lamda_rel()
        return [0.5 * (1 + self.beta_c() * (l_rel[0] - 0.3) + np.power(l_rel[0], 2)),
                0.5 * (1 + self.beta_c() * (l_rel[1] - 0.3) + np.power(l_rel[1], 2))]

    def k_c(self) -> list_array:
        """
        EN 1995-1-1 6.3.2 (6.25, 6.26)
        @return:
        """
        # SFS-EN 1995-1-1 § 6.3.2 (6.26)
        k = self.k()
        l_rel = self.lamda_rel()
        r = [np.power(k[0] + np.power(np.power(k[0], 2) - np.power(l_rel[0], 2), 0.5), -1),
             np.power(k[1] + np.power(np.power(k[1], 2) - np.power(l_rel[1], 2), 0.5), -1)]
        return r

    def l_ef(self, i: int, yz: int, loc: float) -> float:
        # SFS-EN 1995-1-1 § 6.3.3 (6.1)
        # TODO SFS 1995-1-1 s. 43 Taulukko 6.1
        # TODO Pistekuorma ja vakiomomentti tilanne puuttuu

        alpha = 0.9

        if self.edge_load == 'compression':
            return (self.L_c[yz][i] * alpha) + 2 * self.section.get_H(loc)
        elif self.edge_load == 'tension':
            return (self.L_c[yz][i] * alpha) - 0.5 * self.section.get_H(loc)

        return self.L_c[yz][i] * alpha

    def sigma_mcrit(self, i: int, yz: int, loc: float) -> float:
        """
        For member undergoing bending relative to it's stronger axis: SFS-EN 1995-1-1 § 6.3.3 (6.31)
        @return: Single value of critical stress during lateral buckling check
        """
        try:
            return pi * sqrt(self.material.E005 * self.section.get_I(loc)[1] * self.material.G005 * self.I_tor(loc)) / (
                        self.l_ef(i, yz, loc) * self.section.get_Wel(loc)[0])
        except:
            return ((0.78 * self.section.B ** 2) / (
                        self.section.get_H(loc) * self.l_ef(i, yz, loc))) * self.material.E005

    def lamda_relm(self, i: int, yz: int, loc: float) -> float:
        """
        SFS-EN 1995-1-1 § 6.3.3 (6.30)
        @return: Relative slenderness for bending problems
        """
        return np.sqrt(self.material.fmk / self.sigma_mcrit(i, yz, loc))

    def k_crit(self, i: int, yz: int, loc) -> float:
        # SFS-EN 1995-1-1 § 6.3.3 (6.34)
        if self.lamda_relm(i, yz, loc) <= 0.75:
            return 1
        elif 0.75 < self.lamda_relm(i, yz, loc) <= 1.4:
            return 1.56 - 0.75 * self.lamda_relm(i, yz, loc)
        else:
            return 1 / self.lamda_relm(i, yz, loc) ** 2

    def check_single_tapered_beam_bending_tension(self) -> float:
        # SFS 1995-1-1 (6.38)
        return self.sigma_m_alpha_d() / (self.k_m_alpha() * self.fmd)

    def check_tapered_beam_bending_tension(self) -> float:
        """
        EN 1995-1-1 6.4.3 (6.41) kuva 6.9
        @return: APXBS käyttöasteen
        """
        try:
            sig = self.k_l() * self.sigma_md(self.loc.index(self.section.xap))[0]
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
            sig = self.k_l() * (
                    self.sigma_md(self.loc.index(minusloc))[0] + self.sigma_md(self.loc.index(plusloc))[
                0]) / 2

        return sig / (self.k_r() * self.fmd)

    def check_perpendicular_to_grain_tension(self) -> float:
        """
        EN 1995-1-1 6.4.3 (6.50)
        @return: APXPT käyttöasteen
        """
        return self.sigma_t90d() / (self.k_dis() * self.k_vol() * self.ft90d)

    def check_normal_force(self, n: int) -> float:
        # Normal force utilization ratio along the grain:           SFS-EN 1995-1-1 § 6.1.2 and 6.1.4
        if self.ned[n] > 0:
            UN = abs(self.sigma_t0d(n)) / self.ft0d
        else:
            UN = abs(self.sigma_t0d(n)) / self.fc0d
        r = UN
        return r

    # def check_perpendicular_tension(self):
    #     # Member tension perpendicular to the grain:                SFS-EN 1995-1-1 § 6.1.3
    #     # TODO syyn suuntaan vastaset voimat
    #     pass

    def check_perpendicular_compression(self, location: bool = False) -> (list_float_FEMNode, float):
        # Member compression perpendicular to the grain:            SFS-EN 1995-1-1 § 6.1.5

        # TODO ei laske oikein yhdessä pisteessä kahdelta suunnalta puristettua palkkia
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
                elif l11 < l1:
                    l1 = l11
                self.update_kc90(l1, other_members[0].cross_section.get_H(0))

        if not hasattr(self, 'k_c90'):
            if location:
                return [0, node]
            return 0

        r = abs(max_sigma / (self.k_c90 * self.fc90d))
        if location:
            return [r, node]
        return r

    def check_shear_force(self, n: int, reduce: bool = True) -> float:
        """
        Laskee leikkausjännityksen käyttöasteen annetun indeksin pisteessä

        Jos pienennys on True niin mitoittavan leikkausvoiman pienennys simuloidaan RIL 205-1-2009
        kohdan 6.1.7 mukaisesti palkin päissä.
        Jos annettu indeksi on pienennysalueen sisällä niin etsitään indeksi joka on lähimpänä pienennysalueen reunaa,
        mutta silti pienennyksalueen sisällä

        @param n: indeksi
        @param reduce: pienennys
        @return: käyttöaste
        """
        # Shear force utilization ratio:                            SFS-EN 1995-1-1 § 6.1.7

        if reduce:
            if self.type in ("beam", "rigid-beam"):
                start = list(filter(lambda v: v.member != self if True else False, self.parent.n1.parents.copy()))
                if len(start) == 0:
                    other_h_start = 0
                else:
                    other_h_start = start[0].cross_section.H / 2
                sumH_start = (other_h_start + self.section.get_H(0)) / self.length
                loc = self.loc[n]
                if loc < sumH_start:
                    for i in range(len(self.loc)):
                        try:
                            if self.loc[n + i + 1] > sumH_start:
                                n = n + i
                                break
                        except:
                            pass

                end = list(filter(lambda v: v.member != self if True else False, self.parent.n2.parents.copy()))
                if len(end) == 0:
                    other_h_end = 0
                else:
                    other_h_end = end[0].cross_section.H / 2
                sumH_end = (other_h_end + self.section.get_H(0)) / self.length
                if loc > 1 - sumH_end:
                    for i in range(len(self.loc)):
                        try:
                            if self.loc[n - i - 1] < 1 - sumH_end:
                                n = n - i
                                break
                        except:
                            pass

        UVy = abs(self.tau_d(n)[0]) / self.fvd
        UVz = abs(self.tau_d(n)[1]) / self.fvd
        r = max(UVy, UVz)
        return r

    def check_bending_moment(self, n: int) -> float:
        # Bending moment utilization ratio:                         SFS-EN 1995-1-1 § 6.1.6
        UM1 = abs(self.sigma_md(n)[0] / self.fmd + k_m * self.sigma_md(n)[1] / self.fmd)
        UM2 = abs(self.sigma_md(n)[1] / self.fmd + k_m * self.sigma_md(n)[0] / self.fmd)
        r = max(UM1, UM2)
        return r

    def check_torsion(self, n: int) -> float:
        # Torsion utilization ratio:                                SFS-EN 1995-1-1 § 6.1.8
        UTx = self.tau_tord(n) / (self.k_shape(n) * self.fvd)
        r = UTx
        return r

    def check_bending_tension(self, n: int) -> float:
        #                                                           SFS-EN 1995-1-1 § 6.2.3
        # TODO z-suuntainen taivutusjännitys tarkistettava merkkien osalta (deplanaatiota ei ole otettu huomioon)
        UM1 = self.sigma_t0d(n) / self.ft0d + abs(self.sigma_md(n)[0]) / self.fmd + k_m * self.sigma_md(n)[1] / self.fmd
        UM2 = self.sigma_t0d(n) / self.ft0d + k_m * abs(self.sigma_md(n)[0]) / self.fmd + self.sigma_md(n)[1] / self.fmd
        r = max(UM1, UM2)
        return r

    def check_bending_compression(self, n: int) -> float:
        #                                                           SFS-EN 1995-1-1 § 6.2.4
        # TODO katso ylemää
        UM1 = (self.sigma_c0d(n) / self.fc0d) ** 2 + abs(self.sigma_md(n)[0]) / self.fmd + k_m * self.sigma_md(n)[1] / self.fmd
        UM2 = (self.sigma_c0d(n) / self.fc0d) ** 2 + k_m * abs(self.sigma_md(n)[0]) / self.fmd + self.sigma_md(n)[1] / self.fmd
        r = max(UM1, UM2)
        return r

    def check_buckling(self) -> list_out:
        # Utilization ratio against member buckling:                SFS-EN 1995-1-1 § 6.3.2 (6.23, 6.24))
        sigma_c0d_list = [[], []]
        sigma_myd_list = [[], []]
        sigma_mzd_list = [[], []]

        def form_max_sigma_list(supports: list_out, index: int):
            supp = supports.copy()
            supp.append(0.0)
            supp.append(1.0)
            supp.sort()
            for i in range(len(supp) - 1):
                max_sigma_c0d_list = []
                max_sigma_myd_list = []
                max_sigma_mzd_list = []
                for loc in self.loc:
                    if supp[i] <= loc <= supp[i + 1]:
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

        UB = [(abs(sc0) / (self.k_c()[0] * self.fc0d) + abs(my0) / self.fmd + k_m * mz0 / self.fmd),
              (abs(sc1) / (self.k_c()[1] * self.fc0d) + k_m * abs(my1) / self.fmd + mz1 / self.fmd)]
        UB_fixed = [UB[0] * (sc0 < 0), UB[1] * (sc1 < 0)]

        return [max(UB_fixed[0]), max(UB_fixed[1])]

    def check_LT_buckling(self, get_full: bool = False) -> (float, dict_str_float):
        # Utilization ratio against member lateral buckling:        SFS-EN 1995-1-1 § 6.3.3 (3) and (6)
        ULT_list = []
        supp = self.lateral_support_z.copy()
        supp.append(0.0)
        supp.append(1.0)
        supp.sort()
        for i in range(len(supp) - 1):
            for loc in self.loc:
                if supp[i] <= loc <= supp[i + 1]:
                    if max(self.myed, key=abs) != 0 and abs(max(self.ned, key=abs)) <= 0.1:
                        ULT_list.append({'r': self.sigma_md(self.loc.index(loc))[0] / (self.k_crit(i, 1, loc) * self.fmd), 'kaava': 33, 'loc': loc, 'i': i})
                    else:
                        ULT_list.append({
                            'r': ((self.sigma_md(self.loc.index(loc))[0] / (self.k_crit(i, 1, loc) * self.fmd)) ** 2) +
                            abs(self.sigma_c0d(self.loc.index(loc))) / (self.k_c()[1][i] * self.fc0d), 'kaava': 35, 'loc': loc, 'i': i})

                # TODO sigma_cod positiivinen tai negatiivinen
        if get_full:
            return max(ULT_list, key=lambda x: x['r'])
        return max(ULT_list, key=lambda x: x['r'])['r']

    def bracing(self, n: int=None, ks=None) -> (tuple_out[float, float, float], tuple_out[float, float]):
        """
        # SFS EN 1995-1-1 9.2.5
        kaavat 9.34, 9.35, 9.36, 9.37

        @param n: kehien määrä
        @param ks: ks
        @return: if n is given returns C, Fd and qd, and if its not given returns C and Fd
        """

        # Z-suunnan stabiloivan tuen voima ja jousijäykkyys 1.muoto
        m = len(self.lateral_support_z) + 1
        a = self.parent.length / m
        if ks is None:
            ks = 2 + 2 * cos(radians(180 / m))
            #print(f'ks {ks}')
        Nd_list = []
        supp = self.lateral_support_z.copy()
        supp.append(0.0)
        supp.append(1.0)
        supp.sort()
        for i in range(len(supp) - 1):
            for loc in self.loc:
                if supp[i] <= loc <= supp[i + 1]:
                    Md = self.myed[self.loc.index(loc)]
                    Nd_temp = (1 - self.k_crit(i, 1, loc)) * (Md / self.section.get_H(loc))
                    Nd_list.append(Nd_temp)

        Nd = max(Nd_list)
        C = Nd / a
        if self.material.type == 'solid_timber':
            Fd = Nd / 50
        else:
            Fd = Nd / 80

        if n is not None:
            qd = min([1, sqrt(15/self.parent.length)]) * ((n * Nd) / (50 * self.parent.length))
            return C, Fd, qd
        return C, Fd

    def check_apex_bending_stress(self) -> float:
        # EN 1995-1-1 (6.41)
        sigma = self.k_l() * self.M_apd() / self.W_ap()[0]
        return sigma / (self.k_r() * self.fmd)

    def check_apex_perpendicular_tension(self) -> float:
        # EN 1995-1-1 (6.50)
        return self.sigma_t90d() / (self.k_dis() * self.k_vol() * self.ft90d)

    def check_apex_shear_perpendicular_tension_combined(self, kaava6_55:bool=False) -> float:
        # EN 1995-1-1 (6.53)
        #
        n = 0
        start = list(filter(lambda v: v.member != self if True else False, self.parent.n1.parents.copy()))
        if len(start) == 0:
            other_h_start = 0
        else:
            other_h_start = start[0].cross_section.H / 2
        sumH_start = (other_h_start + self.section.get_H(0)) / self.length
        loc = self.loc[n]
        if loc < sumH_start:
            for i in range(len(self.loc)):
                try:
                    if self.loc[n + i + 1] > sumH_start:
                        n = n + i
                        break
                except:
                    pass

        tau1 = self.tau_d(n)[1]
        n = self.nsect() - 1
        end = list(filter(lambda v: v.member != self if True else False, self.parent.n2.parents.copy()))
        if len(end) == 0:
            other_h_end = 0
        else:
            other_h_end = end[0].cross_section.H / 2
        sumH_end = (other_h_end + self.section.get_H(0)) / self.length
        if loc > 1 - sumH_end:
            for i in range(len(self.loc)):
                try:
                    if self.loc[n - i - 1] < 1 - sumH_end:
                        n = n - i
                        break
                except:
                    pass
        tau2 = self.tau_d(n)[1]
        tau = max(tau1, tau2)
        # if isinstance(self.section, (KaarevaHarjaPalkki, KaarevaPalkki)):
        #     closesttodownx1 = None
        #     closesttodownx2 = None
        #     for loc in self.loc:
        #         if closesttodownx1 is None:
        #             closesttodownx1 = loc
        #             continue
        #         if closesttodownx2 is None:
        #             closesttodownx2 = loc
        #             continue
        #         if abs((self.section.dwn_x1 / self.length) - loc) < abs(closesttodownx1 - loc):
        #             closesttodownx1 = loc
        #         if abs(self.section.dwn_x2 / self.length - loc) < abs(closesttodownx2 - loc):
        #             closesttodownx2 = loc
        #
        #     tau_d_dw_x1 = self.tau_d(self.loc.index(closesttodownx1))
        #     tau_d_dw_x2 = self.tau_d(self.loc.index(closesttodownx2))
        #
        #     tau = abs(min(tau_d_dw_x1[1], tau_d_dw_x2[1], key=abs))
        #
        # elif isinstance(self.section, (HarjaPalkki, MahaPalkki)):
        #     closesttoleft = None
        #     closesttoright = None
        #     for loc in self.loc:
        #         if closesttoleft is None:
        #             closesttoleft = loc
        #             continue
        #         if closesttoright is None:
        #             closesttoright = loc
        #             continue
        #         if abs(self.section.xap - (self.section.hap / 2 / self.length) - loc) < abs(closesttoleft - loc):
        #             closesttoleft = loc
        #         if abs(self.section.xap + (self.section.hap / 2 / self.length) - loc) < abs(closesttoright - loc):
        #             closesttoright = loc
        #
        #     tau_d_dw_x1 = self.tau_d(self.loc.index(closesttoleft))
        #     tau_d_dw_x2 = self.tau_d(self.loc.index(closesttoright))
        #
        #     tau = abs(min(tau_d_dw_x1[1], tau_d_dw_x2[1], key=abs))
        # else:
        #     raise NotImplementedError
        return (tau / self.fvd) + self.sigma_t90d(kaava6_55) / (self.k_dis() * self.k_vol() * self.ft90d)

    def add_section(self, ned: float = 0.0, myed: float = 0.0, mzed: float = 0.0, vyed: float = 0.0, vzed: float = 0.0,
                    ted: float = 0.0, loc: float = 0.0):
        """ Adds new section with internal forces and their location """

        self.ned.append(ned)
        self.myed.append(myed)
        self.mzed.append(mzed)
        self.vyed.append(vyed)
        self.vzed.append(vzed)
        self.ted.append(ted)
        self.loc.append(loc)
        self.utilization.append({})

    def clear_sections(self):
        self.ned.clear()
        self.myed.clear()
        self.mzed.clear()
        self.vyed.clear()
        self.vzed.clear()
        self.ted.clear()
        self.loc.clear()
        self.utilization.clear()

    def check_section(self, n: int = 0, reduce_uv: bool = True, verb: bool = False) -> list_out:
        """ Verify resistance of section 'n'
            returns a list of utilization ratios
        """

        r = self.section_resistance(n=n, reduce_uv=reduce_uv)
        self.utilization[n] = {'un': r[0], 'uv': r[1], 'um': r[2], 'ut': r[3], 'ubt': r[4], 'ubc': r[5]}
        return r

    def section_resistance(self, n: int = 0, return_list: bool = True, reduce_uv: bool = True) -> (list_out, float):
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

        UN = self.check_normal_force(n)
        UV = self.check_shear_force(n, reduce=reduce_uv)
        UM = self.check_bending_moment(n)
        UT = self.check_torsion(n)
        UBT = self.check_bending_tension(n)
        UBC = self.check_bending_compression(n)

        if return_list:
            return [UN, UV, UM, UT, UBT, UBC]
        else:
            return max([UN, UV, UM, UT, UBT, UBC])

    def check_cross_section(self):
        """ Checks if cross-section is strong enough.
            Gives boolean value for self.is_strong_enough
            Saves list of stress ratios to self.r
        """

        self.r = []
        # Cross-sectional stress ratios in member's nodes' locations
        self.r.extend(self.check_sections())
        # Buckling about y - and z - axis
        buckling_r = self.check_buckling()
        self.r.append(buckling_r[0])
        self.r.append(buckling_r[1])
        # Lateral-Torsional buckling
        self.r.append(self.check_LT_buckling())

        if max(self.r) > 1.0:
            self.is_strong_enough = False
        else:
            self.is_strong_enough = True

    def check_sections(self, reduce_uv: bool = True) -> np.ndarray:
        """ Verify resistance of all sections """
        r = np.zeros_like(self.check_section())
        for n in range(self.nsect()):
            r = np.vstack((r, self.check_section(n, reduce_uv=reduce_uv)))
            # r.append(self.check_section(n))
        return np.max(r, axis=0)

    def nsect(self) -> int:
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

    def design(self) -> float:
        r = max(self.check_sections())
        r = max(r, np.max(self.check_buckling()))
        if self.type in ('beam', 'rigid-beam'):
            r = max(r, self.check_LT_buckling())
            r = max(r, self.check_perpendicular_compression())

        return r

    def __repr__(self):
        return self.parent.__repr__()

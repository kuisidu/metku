"""
Timber section parameters, including timber fire design calculations (simple method)

@author: Viktor Haimi
"""

from abc import ABC, abstractmethod
import numpy as np
from scipy import integrate
try:
    from metku.eurocodes.en1995 import en1995_1_1, en1995_1_2
    from metku.eurocodes.en1995.fire_protection import *
    from metku.materials.timber_data import Timber, T
except:
    from eurocodes.en1995 import en1995_1_1, en1995_1_2
    from eurocodes.en1995.fire_protection import *
    from materials.timber_data import Timber, T
import matplotlib.pyplot as plt
import math

import sys

if sys.version_info[0] == 3 and sys.version_info[1] < 9:
    from typing import List
    list_out = List[float]
else:
    list_out = list[float]


class TimberSection:
    def __init__(self, B: (int, float), H: (int, float), material: T=T.GL30c,
                 fire_protection_right: FireProtection = None, fire_protection_left: FireProtection = None,
                 fire_protection_top: FireProtection = None, fire_protection_bottom: FireProtection = None,
                 fire_protection_generic: FireProtection = None, fire_protection_sides: int = 3,
                 sides_on_fire: int = 3, lvldir: str = 'edge'):
        """
        Tasakorkea ja -leveä suorakaiteinen sauva

        @param B: leveys [mm]
        @param H: korkeus [mm]
        @param material: materiaali
        @param fire_protection_right: palosuoja joka laitetaan oikealle, yliajaa fire_protection_genericin
        @param fire_protection_left: palosuoja joka laitetaan vasemmalle, yliajaa fire_protection_genericin
        @param fire_protection_top: palosuoja joka laitetaan ylös, yliajaa fire_protection_genericin
        @param fire_protection_bottom: palosuoja joka laitetaan alas, yliajaa fire_protection_genericin
        @param fire_protection_generic: geneerinen palosuoja joka annetaan monelle sivulle samaan aikaan
        @param fire_protection_sides: jos antaa 3 niin fire_protection_top jää tyhjäksi.
                                      jos antaa jonkun muun kuin 3 niin suojaa kaikki neljä sivua fire_protection_genericillä
        @param sides_on_fire: jos antaa 4 niin poltetaan kaikilta neljältä sivulta,
                              jos antaa jonkun muun kuin 4 niin poltetaan kolmelta sivulta. (ylhäältä ei polteta)
        @param lvldir: lvl:n suuntaus. lappeellaan tai syrjällään. 'edge' tai 'flat'
        """
        self._H = H
        self._B = B

        self.timber_member = None

        self.beta_0 = None
        self.beta_n = None
        self.R = 0

        self.__fire_protection_generic = fire_protection_generic
        self.fire_protection_sides = fire_protection_sides
        self.sides_on_fire = sides_on_fire
        if fire_protection_sides == 3:
            self.fire_protection_right = fire_protection_right if fire_protection_right else fire_protection_generic
            self.fire_protection_left = fire_protection_left if fire_protection_left else fire_protection_generic
            self.fire_protection_bottom = fire_protection_bottom if fire_protection_bottom else fire_protection_generic
            self.fire_protection_top = fire_protection_top
        else:
            self.fire_protection_right = fire_protection_right if fire_protection_right else fire_protection_generic
            self.fire_protection_left = fire_protection_left if fire_protection_left else fire_protection_generic
            self.fire_protection_bottom = fire_protection_bottom if fire_protection_bottom else fire_protection_generic
            self.fire_protection_top = fire_protection_top if fire_protection_top else fire_protection_generic

        self.lvldir = lvldir
        self.material = material

    @property
    def material(self) -> Timber:
        return self.__material

    def set_beta_0_and_n(self):
        """
        laskee ja asettaa hiiltymisnopeudet beta_0:n ja beta_n:n
        SFS-EN1995-1-2 Taulukko 3.1 Tai RIL-205-2-2009 Taulukon 3.2
        @return: None
        """
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self._B, self._H))

    @material.setter
    def material(self, val: T):
        self.__material = Timber(val)
        if self.__material.type == 'lvl':
            self.__material.direction = self.lvldir
        self.set_beta_0_and_n()

        # RIL-205-2-2009 3.4.3
        # Asetetaan palosuojille millainen materiaali suojataan jotta voidaan laskea ta ja tf
        if self.fire_protection_right:
            self.fire_protection_right.protected_material(self.material)
            if isinstance(self.fire_protection_right, (GypsumPlasterboardF, GypsumPlasterboardFF, GypsumPlasterboardAF)):
                self.fire_protection_right.set_tf_for_gypsumF()
            self.fire_protection_right.set_ta()
        if self.fire_protection_left:
            self.fire_protection_left.protected_material(self.material)
            if isinstance(self.fire_protection_left, (GypsumPlasterboardF, GypsumPlasterboardFF, GypsumPlasterboardAF)):
                self.fire_protection_left.set_tf_for_gypsumF()
            self.fire_protection_left.set_ta()
        if self.fire_protection_bottom:
            self.fire_protection_bottom.protected_material(self.material)
            if isinstance(self.fire_protection_bottom, (GypsumPlasterboardF, GypsumPlasterboardFF, GypsumPlasterboardAF)):
                self.fire_protection_bottom.set_tf_for_gypsumF()
            self.fire_protection_bottom.set_ta()
        if self.fire_protection_top:
            self.fire_protection_top.protected_material(self.material)
            if isinstance(self.fire_protection_top, (GypsumPlasterboardF, GypsumPlasterboardFF, GypsumPlasterboardAF)):
                self.fire_protection_top.set_tf_for_gypsumF()
            self.fire_protection_top.set_ta()

    @property
    def H(self) -> float:
        h = self._H
        h -= en1995_1_2.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= en1995_1_2.d_char_n(self.fire_protection_top, self.R, self.beta_n)
        if h <= 0:
            h = 1e-5
        return h

    @H.setter
    def H(self, val: (int, float)):
        self._H = val

    def get_H(self, x: (int, float)) -> (int, float):
        return self.H

    def get_A(self, x: (int, float)) -> (int, float):
        return self.A

    @property
    def B(self) -> float:
        b = self._B
        b -= en1995_1_2.d_char_n(self.fire_protection_right, self.R, self.beta_n)
        b -= en1995_1_2.d_char_n(self.fire_protection_left, self.R, self.beta_n)
        if b <= 0:
            b = 1e-5
        return b

    @B.setter
    def B(self, val: (int, float)):
        self._B = val

    @property
    def I(self) -> np.ndarray:
        return np.asarray([(self.B * self.H ** 3) / 12, (self.H * self.B ** 3) / 12])

    def get_I(self, x: (int, float)) -> np.ndarray:
        return self.I

    @property
    def Wel(self) -> np.ndarray:
        return np.asarray([(self.B * self.H ** 2) / 6, (self.H * self.B ** 2) / 6])

    def get_Wel(self, x: (int, float)) -> np.ndarray:
        return self.Wel

    @property
    def k_h(self) -> float:
        # SFS-EN1995-1-1 3.2, 3.3, 3.4
        # TODO pitää kysyä samilta
        # toimii kaikille rajatiloille, jos halutaan poistaa palolta niin pitää muokata koodia
        Med = max(self.timber_member.MEdY, self.timber_member.MEdZ)
        if self.material.rhok <= 700 and self.material.type == 'solid_timber':
            if self.timber_member.NEd > 0 and self.B < 150 and Med != 0 and self.H < 150:
                return min(min(1.3, (150 / self.B) ** 0.2), min(1.3, (150 / self.H) ** 0.2))
            elif self.timber_member.NEd > 0 and self.B < 150:
                return min(1.3, (150 / self.B) ** 0.2)
            elif Med != 0 and self.H < 150:
                return min(1.3, (150 / self.H) ** 0.2)
            else:
                return 1.0
        elif self.material.type == 'glt':
            if self.timber_member.NEd > 0 and self.B < 600 and Med != 0 and self.H < 600:
                return min(min(1.1, (600 / self.B) ** 0.1), min(1.1, (600 / self.H) ** 0.1))
            elif self.timber_member.NEd > 0 and self.B < 600:
                return min(1.1, (600 / self.B) ** 0.1)
            elif Med != 0 and self.H < 600:
                return min(1.1, (600 / self.H) ** 0.1)
            else:
                return 1.0
        elif self.material.type == 'lvl':
            if Med != 0 and self.H != 300:
                return min(1.2, (300 / self.H) ** self.material.s)
            else:
                return 1.0
        else:
            return 1.0

    @property
    def A(self) -> (int, float):
        '''
        RIL-205-2-2009 4.1.1S 4.1.2S
        Laskee palolle altistuneen rakenteen poikkileikkauksen tehollisen pinta-alan palonkestoajan funktiona
        @return: tehollinen pinta-ala A_ef [mm^2]
        '''

        b = self.B
        h = self.H

        if b <= 0 or h <= 0:
            return 0
        return b*h

    @property
    def fire_protection_generic(self) -> FireProtection:
        return self.__fire_protection_generic

    @fire_protection_generic.setter
    def fire_protection_generic(self, val: FireProtection):
        self.__fire_protection_generic = val
        if self.fire_protection_sides == 3:
            self.fire_protection_right = val
            self.fire_protection_left = val
            self.fire_protection_bottom = val
        else:
            self.fire_protection_right = val
            self.fire_protection_left = val
            self.fire_protection_bottom = val
            self.fire_protection_top = val

    def get_R(self) -> (int, float):
        return self.R

    def __str__(self):
        return f'TimberSection({round(self._B, 0)}X{round(self._H, 0)})'

    def __repr__(self):
        return f"TimberSection {self.material.timber_type} {round(self.H, 0)}X{round(self.B, 0)}"


class TaperedSection(TimberSection, ABC):
    def __init__(self, B: (int, float), H: (int, float), material: T=T.C14,
                 fire_protection_right: FireProtection = None, fire_protection_left: FireProtection = None,
                 fire_protection_top: FireProtection = None, fire_protection_bottom: FireProtection = None,
                 fire_protection_generic: FireProtection = None, fire_protection_sides=3, sides_on_fire=3, lvldir='edge'):
        super().__init__(B, H, material, fire_protection_right, fire_protection_left,
                         fire_protection_top, fire_protection_bottom, fire_protection_generic, fire_protection_sides,
                         sides_on_fire, lvldir)
        self.small_elements = []
        self.num_elements = -1

    def create_elements(self):
        # TODO ei laske oiken jos nodet ei ole tasavälein
        if len(self.small_elements) > 0:
            raise Exception('Use update_elements to update')
        if self.num_elements < 0:
            raise Exception('num_elements not updated')
        elem_l = 1 / self.num_elements
        for i in range(self.num_elements):
            a = self.get_A(i * elem_l + 0.5 * elem_l)
            iy = self.get_I(i * elem_l + 0.5 * elem_l)
            self.small_elements.append(SmallElementSection(a, iy))

    def update_elements(self):
        elem_l = 1 / self.num_elements
        for i in range(len(self.small_elements)):
            a = self.get_A(i * elem_l + 0.5 * elem_l)
            iy = self.get_I(i * elem_l + 0.5 * elem_l)
            self.small_elements[i].A = a
            self.small_elements[i].I = iy

    def get_H(self, x: (int, float)) -> float:
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """

        h = self.get_unmodified_h(x)

        h -= en1995_1_2.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= en1995_1_2.d_char_n(self.fire_protection_top, self.R, self.beta_n)
        if h <= 0:
            h = 1e-5
        return h

    def get_A(self, x: (int, float)) -> float:
        """
        Calculates Area of cross section

        @param x: local coordinate 0 ... 1
        @return: area at location
        """
        b = self.B
        h = self.get_H(x)

        if b <= 0 or h <= 0:
            return 1e-5
        return b * h

    def get_Wel(self, x: (int, float)) -> np.ndarray:
        return np.asarray([(self.B * self.get_H(x) ** 2) / 6, (self.get_H(x) * self.B ** 2) / 6])

    def get_I(self, x: (int, float)) -> np.ndarray:
        return np.asarray([(self.B * self.get_H(x) ** 3) / 12, (self.get_H(x) * self.B ** 3) / 12])

    @abstractmethod
    def get_unmodified_h(self, x: (int, float)) -> float:
        pass


class PulpettiPalkki(TaperedSection):
    def __init__(self, B, H0, Hap, material=T.C14, flipped=False, fire_protection_right: FireProtection = None,
                 fire_protection_left: FireProtection = None, fire_protection_top: FireProtection = None,
                 fire_protection_bottom: FireProtection = None, fire_protection_generic: FireProtection = None,
                 fire_protection_sides=3, sides_on_fire=3, lvldir='edge'):
        """
        Pulpettipalkki

        @param B: paksuus [mm]
        @param H0: palkin alkupään korkeus [mm]
        @param Hap: palkin loppupään korkeus[mm]
        @param material: materiaali
        @param flipped: ylösalaisin
        @param fire_protection_right: palosuoja joka laitetaan oikealle, yliajaa fire_protection_genericin
        @param fire_protection_left: palosuoja joka laitetaan vasemmalle, yliajaa fire_protection_genericin
        @param fire_protection_top: palosuoja joka laitetaan ylös, yliajaa fire_protection_genericin
        @param fire_protection_bottom: palosuoja joka laitetaan alas, yliajaa fire_protection_genericin
        @param fire_protection_generic: geneerinen palosuoja joka annetaan monelle sivulle samaan aikaan
        @param fire_protection_sides: jos antaa 3 niin fire_protection_top jää tyhjäksi.
                                      jos antaa jonkun muun kuin 3 niin suojaa kaikki neljä sivua fire_protection_genericillä
        @param sides_on_fire: jos antaa 4 niin poltetaan kaikilta neljältä sivulta,
                              jos antaa jonkun muun kuin 4 niin poltetaan kolmelta sivulta. (ylhäältä ei polteta)
        @param lvldir: lvl:n suuntaus. lappeellaan tai syrjällään. 'edge' tai 'flat'
        """
        self.__material = None
        self.h0 = H0
        self.hap = Hap
        self.flipped = flipped
        super().__init__(B, H0, material=material, fire_protection_right=fire_protection_right,
                         fire_protection_left=fire_protection_left,
                         fire_protection_top=fire_protection_top,
                         fire_protection_bottom=fire_protection_bottom,
                         fire_protection_generic=fire_protection_generic,
                         fire_protection_sides=fire_protection_sides,
                         sides_on_fire=sides_on_fire, lvldir=lvldir)

    @property
    def H0(self) -> float:
        return self.h0

    @property
    def Hap(self) -> float:
        return self.hap

    @H0.setter
    def H0(self, val):
        self.h0 = val
        self.update_elements()

    @Hap.setter
    def Hap(self, val):
        self.hap = val
        self.update_elements()

    def get_unmodified_h(self, x: (int, float)) -> float:
        k = self.hap - self.h0
        h = self.h0 + k * x
        return h

    def __str__(self):
        return f'Pulpettipalkki({round(self._B, 0)}X{round(self.H0, 0)}X{round(self.Hap, 0)})'


class HarjaPalkki(TaperedSection):
    def __init__(self, B, H0, Hap, xap=0.5, material=T.C14, fire_protection_right: FireProtection = None,
                 fire_protection_left: FireProtection = None, fire_protection_top: FireProtection = None,
                 fire_protection_bottom: FireProtection = None, fire_protection_generic: FireProtection = None,
                 fire_protection_sides=3, sides_on_fire=3, lvldir='edge'):
        """
        Harhjapalkki

        @param B: paksuus [mm]
        @param H0: palkin alkupään korkeus [mm]
        @param Hap:
        @param xap:
        @param material:
        @param fire_protection_right: palosuoja joka laitetaan oikealle, yliajaa fire_protection_genericin
        @param fire_protection_left: palosuoja joka laitetaan vasemmalle, yliajaa fire_protection_genericin
        @param fire_protection_top: palosuoja joka laitetaan ylös, yliajaa fire_protection_genericin
        @param fire_protection_bottom: palosuoja joka laitetaan alas, yliajaa fire_protection_genericin
        @param fire_protection_generic: geneerinen palosuoja joka annetaan monelle sivulle samaan aikaan
        @param fire_protection_sides: jos antaa 3 niin fire_protection_top jää tyhjäksi.
                                      jos antaa jonkun muun kuin 3 niin suojaa kaikki neljä sivua fire_protection_genericillä
        @param sides_on_fire: jos antaa 4 niin poltetaan kaikilta neljältä sivulta,
                              jos antaa jonkun muun kuin 4 niin poltetaan kolmelta sivulta. (ylhäältä ei polteta)
        @param lvldir: lvl:n suuntaus. lappeellaan tai syrjällään. 'edge' tai 'flat'
        """
        self.__material = None
        self.h0 = H0
        self.hap = Hap
        self.xap = xap
        super().__init__(B, H0, material=material, fire_protection_right=fire_protection_right,
                         fire_protection_left=fire_protection_left,
                         fire_protection_top=fire_protection_top,
                         fire_protection_bottom=fire_protection_bottom,
                         fire_protection_generic=fire_protection_generic,
                         fire_protection_sides=fire_protection_sides,
                         sides_on_fire=sides_on_fire, lvldir=lvldir)

    def get_unmodified_h(self, x: (int, float)) -> float:
        """
        Defines harjapalkki's height
        @param x: Harjapalkki's longitudinal coordinate
        @param h0: Section height at the beginning and at the end of the beam
        @param hap: Section height at the paramount point of the beam
        @param xap: X coordinate of the highest cross section of the beam
        @param L: Beam span length (support node to support node)
        @return: Height
        """

        h = (((self.hap - self.h0) / self.xap) * x + self.h0 - max(
            ((self.hap - self.h0) / self.xap) * x + self.h0 - ((self.h0 - self.hap) / (1 - self.xap) * x + self.hap +
                                                               ((self.hap - self.h0) / (1 - self.xap)) * self.xap), 0))
        return h

    @property
    def H0(self) -> float:
        return self.h0

    @property
    def Hap(self) -> float:
        return self.hap

    @property
    def Xap(self) -> float:
        return self.xap

    @H0.setter
    def H0(self, val: (int, float)):
        self.h0 = val
        self.update_elements()

    @Hap.setter
    def Hap(self, val: (int, float)):
        self.hap = val
        self.update_elements()

    @Xap.setter
    def Xap(self, val: (int, float)):
        self.xap = val
        self.update_elements()

    def __str__(self):
        return f'Harjapalkki({round(self._B, 0)}X{round(self.h0, 0)}X{round(self.hap, 0)})'


class MahaPalkki(TaperedSection):
    def __init__(self, B, H0, alpha, xap=0.5, material=T.C14, fire_protection_right: FireProtection = None,
                 fire_protection_left: FireProtection = None, fire_protection_top: FireProtection = None,
                 fire_protection_bottom: FireProtection = None, fire_protection_generic: FireProtection = None,
                 fire_protection_sides=3, sides_on_fire=3, lvldir='edge'):
        """
        Mahapalkki harjapalkkitapainen

        @param B: paksuus [mm]
        @param H0: palkin alkupään korkeus [mm]
        @param alpha:
        @param xap:
        @param material:
        @param fire_protection_right: palosuoja joka laitetaan oikealle, yliajaa fire_protection_genericin
        @param fire_protection_left: palosuoja joka laitetaan vasemmalle, yliajaa fire_protection_genericin
        @param fire_protection_top: palosuoja joka laitetaan ylös, yliajaa fire_protection_genericin
        @param fire_protection_bottom: palosuoja joka laitetaan alas, yliajaa fire_protection_genericin
        @param fire_protection_generic: geneerinen palosuoja joka annetaan monelle sivulle samaan aikaan
        @param fire_protection_sides: jos antaa 3 niin fire_protection_top jää tyhjäksi.
                                      jos antaa jonkun muun kuin 3 niin suojaa kaikki neljä sivua fire_protection_genericillä
        @param sides_on_fire: jos antaa 4 niin poltetaan kaikilta neljältä sivulta,
                              jos antaa jonkun muun kuin 4 niin poltetaan kolmelta sivulta. (ylhäältä ei polteta)
        @param lvldir: lvl:n suuntaus. lappeellaan tai syrjällään. 'edge' tai 'flat'
        """
        self.__material = None
        self.h0 = H0
        self.alpha = math.radians(alpha)
        self.beta = None
        self.parent = None
        self.xap = xap
        super().__init__(B, H0, material=material, fire_protection_right=fire_protection_right,
                         fire_protection_left=fire_protection_left,
                         fire_protection_top=fire_protection_top,
                         fire_protection_bottom=fire_protection_bottom,
                         fire_protection_generic=fire_protection_generic,
                         fire_protection_sides=fire_protection_sides,
                         sides_on_fire=sides_on_fire, lvldir=lvldir)

    def get_unmodified_h(self, x: (int, float)) -> float:
        x *= self.parent.length

        ya = self.parent.length * math.tan(self.alpha)
        self.beta = np.arctan(ya / (self.parent.length - self.xap * self.parent.length))

        h = ((math.tan(self.alpha) * x + self.h0) / math.cos(self.alpha) -
             max(0, (x * math.tan(self.beta) - self.xap * self.parent.length * math.tan(self.beta))) / math.cos(self.alpha))
        return h

    @property
    def H0(self) -> float:
        return self.h0

    @property
    def Hap(self) -> float:
        return self.get_H(self.xap)

    @property
    def hap(self) -> float:
        return self.Hap

    @property
    def Xap(self) -> float:
        return self.xap

    @H0.setter
    def H0(self, val: (int, float)):
        self.h0 = val
        self.update_elements()

    @Xap.setter
    def Xap(self, val: (int, float)):
        self.xap = val
        self.update_elements()

    def __str__(self):
        return f'MahaPalkki({round(self._B, 0)}X{round(self.h0, 0)}X{round(self.hap, 0)})'


class KaarevaPalkki(TaperedSection):
    def __init__(self, B, H0, alpha, beta, r_in, material=T.C14, t=45, fire_protection_right: FireProtection = None,
                 fire_protection_left: FireProtection = None, fire_protection_top: FireProtection=None,
                 fire_protection_bottom: FireProtection = None, fire_protection_generic: FireProtection=None,
                 fire_protection_sides=3, sides_on_fire=3, lvldir='edge'):
        """
        Kaarevapalkki

        @param B: paksuus [mm]
        @param H0: palkin alkupään korkeus [mm]
        @param alpha:
        @param beta:
        @param r_in:
        @param material:
        @param t:
        @param fire_protection_right: palosuoja joka laitetaan oikealle, yliajaa fire_protection_genericin
        @param fire_protection_left: palosuoja joka laitetaan vasemmalle, yliajaa fire_protection_genericin
        @param fire_protection_top: palosuoja joka laitetaan ylös, yliajaa fire_protection_genericin
        @param fire_protection_bottom: palosuoja joka laitetaan alas, yliajaa fire_protection_genericin
        @param fire_protection_generic: geneerinen palosuoja joka annetaan monelle sivulle samaan aikaan
        @param fire_protection_sides: jos antaa 3 niin fire_protection_top jää tyhjäksi.
                                      jos antaa jonkun muun kuin 3 niin suojaa kaikki neljä sivua fire_protection_genericillä
        @param sides_on_fire: jos antaa 4 niin poltetaan kaikilta neljältä sivulta,
                              jos antaa jonkun muun kuin 4 niin poltetaan kolmelta sivulta. (ylhäältä ei polteta)
        @param lvldir: lvl:n suuntaus. lappeellaan tai syrjällään. 'edge' tai 'flat'
        """
        self.h0 = H0
        self.parent = None
        self.alpha = np.pi * alpha / 180
        self.beta = np.pi * beta / 180
        self.r_in = r_in
        self.t = t

        super().__init__(B, H0, material=material, fire_protection_right=fire_protection_right,
                         fire_protection_left=fire_protection_left,
                         fire_protection_top=fire_protection_top,
                         fire_protection_bottom=fire_protection_bottom,
                         fire_protection_generic=fire_protection_generic,
                         fire_protection_sides=fire_protection_sides,
                         sides_on_fire=sides_on_fire, lvldir=lvldir)

        self.xap = 0.5
        self.hap = None

    def update_circles(self):
        self.L_y = 0.5 * self.parent.length * np.tan(self.beta)
        self.U_y = 0.5 * self.parent.length * (np.tan(self.alpha) + np.tan(self.beta)) - self.L_y + self.h0

        self.x_circle = self.parent.length / 2
        self.y_circle_up = -self.r_in * np.cos(self.alpha) + (0.5 * self.parent.length - self.r_in * np.sin(self.alpha)) * np.tan(
            self.alpha) + self.h0
        self.y_circle_dw = (0.5 * self.parent.length - self.r_in * np.sin(self.beta)) * np.tan(self.beta) - self.r_in * np.cos(
            self.beta)

        self.up_x1, self.up_y1 = [self.parent.length / 2 - self.r_in * np.sin(self.alpha),
                                  (self.parent.length / 2 - self.r_in * np.sin(self.alpha)) * np.tan(self.alpha) + self.h0]
        self.dwn_x1, self.dwn_y1 = [self.parent.length / 2 - self.r_in * np.sin(self.beta),
                                    (self.parent.length / 2 - self.r_in * np.sin(self.beta)) * np.tan(self.beta)]
        self.up_x2, self.up_y2 = [self.parent.length / 2 - self.r_in * np.sin(self.alpha) + 2 * (self.r_in * np.sin(self.alpha)),
                                  (self.parent.length / 2 - self.r_in * np.sin(self.alpha)) * np.tan(self.alpha) + self.h0]
        self.dwn_x2, self.dwn_y2 = [self.parent.length / 2 - self.r_in * np.sin(self.beta) + 2 * (self.r_in * np.sin(self.beta)),
                                    (self.parent.length / 2 - self.r_in * np.sin(self.beta)) * np.tan(self.beta)]
        self.hap = self.get_H(0.5)

    def f1(self, x: float) -> float:
        return np.tan(self.alpha) * x + self.h0

    def f2(self, x: float) -> float:
        return np.tan(self.beta) * x

    def f3(self, x: float) -> float:
        return (self.parent.length - x) * np.tan(self.alpha) + self.h0

    def f4(self, x: float) -> float:
        return 2 * self.L_y - np.tan(self.beta) * x

    def f5(self, x: float) -> float:
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_up

    def f6(self, x: float) -> float:
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_dw

    @property
    def H0(self) -> float:
        return self.h0

    @property
    def Hap(self) -> float:
        return self.hap

    @H0.setter
    def H0(self, val: (int, float)):
        self.h0 = val
        self.update_circles()
        self.update_elements()

    @property
    def Alpha(self) -> float:
        return self.alpha

    @Alpha.setter
    def Alpha(self, val: (int, float)):
        self.alpha = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def Beta(self) -> float:
        return self.beta

    @Beta.setter
    def Beta(self, val: (int, float)):
        self.beta = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def R_in(self) -> float:
        return self.r_in

    @R_in.setter
    def R_in(self, val: (int, float)):
        self.r_in = val
        self.update_circles()
        self.update_elements()

    def get_unmodified_h(self, x: (int, float)) -> float:
        x *= self.parent.length
        if 0 <= x < self.up_x1:
            h = self.f1(x) - self.f2(x)
        elif self.up_x1 <= x < self.dwn_x1:
            h = self.f5(x) - self.f2(x)
        elif self.dwn_x1 <= x < self.dwn_x2:
            h = self.f5(x) - self.f6(x)
        elif self.dwn_x2 <= x < self.up_x2:
            h = self.f5(x) - self.f4(x)
        else:
            h = self.f3(x) - self.f4(x)
        return h

    def A_face(self):
        if1 = integrate.quad(self.f1, 0, self.up_x1)
        if2 = integrate.quad(self.f2, 0, self.dwn_x1)
        if3 = integrate.quad(self.f3, self.up_x2, self.parent.length)
        if4 = integrate.quad(self.f4, self.dwn_x2, self.parent.length)
        if5 = integrate.quad(self.f5, self.up_x1, self.up_x2)
        if6 = integrate.quad(self.f6, self.dwn_x1, self.dwn_x2)
        a = if1[0] + if3[0] + if5[0] - if2[0] - if4[0] - if6[0]
        return a

    def __str__(self):
        return f'KaarevaPalkki({round(self._B, 0)}X{round(self.H0, 0)}X{round(self.Alpha, 0)}X{round(self.Beta, 0)}X{round(self.R_in, 0)})'


class KaarevaHarjaPalkki(TaperedSection):
    def __init__(self, B, H0, alpha, beta, r_in, material=T.C14, t=45, fire_protection_right: FireProtection=None, fire_protection_left: FireProtection=None,
                 fire_protection_top: FireProtection=None, fire_protection_bottom: FireProtection=None,
                 fire_protection_generic: FireProtection=None, fire_protection_sides=3, sides_on_fire=3, lvldir='edge'):
        """
        Kaarevaharjapalkki

        @param B: paksuus [mm]
        @param H0: palkin alkupään korkeus [mm]
        @param alpha:
        @param beta:
        @param r_in:
        @param material:
        @param t:
        @param fire_protection_right: palosuoja joka laitetaan oikealle, yliajaa fire_protection_genericin
        @param fire_protection_left: palosuoja joka laitetaan vasemmalle, yliajaa fire_protection_genericin
        @param fire_protection_top: palosuoja joka laitetaan ylös, yliajaa fire_protection_genericin
        @param fire_protection_bottom: palosuoja joka laitetaan alas, yliajaa fire_protection_genericin
        @param fire_protection_generic: geneerinen palosuoja joka annetaan monelle sivulle samaan aikaan
        @param fire_protection_sides: jos antaa 3 niin fire_protection_top jää tyhjäksi.
                                      jos antaa jonkun muun kuin 3 niin suojaa kaikki neljä sivua fire_protection_genericillä
        @param sides_on_fire: jos antaa 4 niin poltetaan kaikilta neljältä sivulta,
                              jos antaa jonkun muun kuin 4 niin poltetaan kolmelta sivulta. (ylhäältä ei polteta)
        @param lvldir: lvl:n suuntaus. lappeellaan tai syrjällään. 'edge' tai 'flat'
        """
        self.h0 = H0
        self.parent = None
        self.alpha = np.pi * alpha / 180
        self.beta = np.pi * beta / 180
        self.r_in = r_in
        self.t = t
        super().__init__(B, H0, material=material, fire_protection_right=fire_protection_right,
                         fire_protection_left=fire_protection_left,
                         fire_protection_top=fire_protection_top,
                         fire_protection_bottom=fire_protection_bottom,
                         fire_protection_generic=fire_protection_generic,
                         fire_protection_sides=fire_protection_sides,
                         sides_on_fire=sides_on_fire, lvldir=lvldir)
        self.hap = None
        self.xap = 0.5

    def update_circles(self):
        self.L_y = 0.5 * self.parent.length * np.tan(self.beta)
        self.U_y = 0.5 * self.parent.length * (np.tan(self.alpha) + np.tan(self.beta)) - self.L_y + self.h0

        self.x_circle = self.parent.length / 2
        self.y_circle_up = -self.r_in * np.cos(self.alpha) + (0.5 * self.parent.length - self.r_in * np.sin(self.alpha)) * np.tan(
            self.alpha) + self.h0
        self.y_circle_dw = (0.5 * self.parent.length - self.r_in * np.sin(self.beta)) * np.tan(self.beta) - self.r_in * np.cos(
            self.beta)

        self.up_x1, self.up_y1 = [self.parent.length / 2 - self.r_in * np.sin(self.alpha),
                                  (self.parent.length / 2 - self.r_in * np.sin(self.alpha)) * np.tan(self.alpha) + self.h0]
        self.dwn_x1, self.dwn_y1 = [self.parent.length / 2 - self.r_in * np.sin(self.beta),
                                    (self.parent.length / 2 - self.r_in * np.sin(self.beta)) * np.tan(self.beta)]
        self.up_x2, self.up_y2 = [self.parent.length / 2 + self.r_in * np.sin(self.alpha),
                                  (self.parent.length / 2 - self.r_in * np.sin(self.alpha)) * np.tan(self.alpha) + self.h0]
        self.dwn_x2, self.dwn_y2 = [self.parent.length / 2 + self.r_in * np.sin(self.beta),
                                    (self.parent.length / 2 - self.r_in * np.sin(self.beta)) * np.tan(self.beta)]
        self.hap = self.get_H(0.5)

    def f1(self, x: float) -> float:
        return np.tan(self.alpha) * x + self.h0

    def f2(self, x: float) -> float:
        return np.tan(self.beta) * x

    def f3(self, x: float) -> float:
        return (self.parent.length - x) * np.tan(self.alpha) + self.h0

    def f4(self, x: float) -> float:
        return 2 * self.L_y - np.tan(self.beta) * x

    def f5(self, x: float) -> float:
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_up

    def f6(self, x: float) -> float:
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_dw

    @property
    def H(self) -> float:
        return self.get_H(0)

    @property
    def H0(self) -> float:
        return self.h0

    @property
    def Hap(self) -> float:
        return self.hap

    @H0.setter
    def H0(self, val: (int, float)):
        self.h0 = val
        self.update_circles()
        self.update_elements()

    @property
    def Alpha(self) -> float:
        return self.alpha

    @Alpha.setter
    def Alpha(self, val: (int, float)):
        self.alpha = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def Beta(self) -> float:
        return self.beta

    @Beta.setter
    def Beta(self, val: (int, float)):
        self.beta = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def R_in(self) -> float:
        return self.r_in

    @R_in.setter
    def R_in(self, val: (int, float)):
        self.r_in = val
        self.update_circles()
        self.update_elements()

    def get_unmodified_h(self, x: (int, float)) -> float:
        x *= self.parent.length
        if 0 <= x < self.dwn_x1:
            h = self.f1(x) - self.f2(x)
        elif self.dwn_x1 <= x < self.parent.length / 2:
            h = self.f1(x) - self.f6(x)
        elif self.parent.length / 2 <= x < self.dwn_x2:
            h = self.f3(x) - self.f6(x)
        else:
            h = self.f3(x) - self.f4(x)
        return h

    def A_face(self):
        if1 = integrate.quad(self.f1, 0, self.parent.length/2)
        if2 = integrate.quad(self.f2, 0, self.dwn_x1)
        if3 = integrate.quad(self.f3, self.parent.length/2, self.parent.length)
        if4 = integrate.quad(self.f4, self.dwn_x2, self.parent.length)
        if6 = integrate.quad(self.f6, self.dwn_x1, self.dwn_x2)
        self.intf1 = if1[0]
        self.intf2 = if2[0]
        self.intf3 = if3[0]
        self.intf4 = if4[0]
        self.intf6 = if6[0]
        a = if1[0] + if3[0] - if2[0] - if4[0] - if6[0]

        return a

    def __str__(self):
        return f'KaarevaHarjaPalkki({round(self._B, 0)}X{round(self.H0, 0)}X{round(self.Alpha, 0)}X{round(self.Beta, 0)}X{round(self.R_in, 0)}X{round(self.Hap, 0)})'


class SmallElementSection:
    def __init__(self, A: (int, float), I: list_out):
        self.A = A
        self.I = I

    def __repr__(self):
        return f'SmallElementSection A:{self.A}, I:{self.I}'


if __name__ == '__main__':
    import frame2d.frame2d as f2d
    from time import time

    ts = TimberSection(30, 50)
    pp = PulpettiPalkki(200, 300, 500)
    hp = HarjaPalkki(20, 30, 100)
    mp = MahaPalkki(100, 200, 3)
    kp = KaarevaPalkki(100, 600, 20, 15, 15000)
    khp = KaarevaHarjaPalkki(100, 600, 20, 15, 15000)
    obj = type('', (), {})()
    obj.length = 10000
    khp.parent = obj
    kp.parent = obj
    khp.update_circles()
    kp.update_circles()

    h = []
    i_l = []
    start = time()
    for i in np.arange(1 + 0.001, step=0.001):
        # print(f'{i}   h {t.get_H(i)}')
        #print(pp.get_A(i))
        h.append(khp.f4(i))
        i_l.append(i)
        # print(f'H: {t.get_H(i)}')
        # print(f'A: {t.get_A(i)}')
        # print(f'I: {t.get_I(i)}')
        # print(f'Wel: {t.get_Wel(i)}')
    print(f'aikaa kulunut {round(time() - start, 20)} s')
    print(f'laskentakertoja {len(h)}')
    plt.plot(i_l, h)
    plt.show()

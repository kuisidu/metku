"""
Timber section parameters, including timber fire design calculations (simple method)

@author: Viktor Haimi
"""

from abc import ABC
import numpy as np
from eurocodes.en1995 import en1995_1_1, en1995_1_2
from eurocodes.en1995.fire_protection import *
from materials.timber_data import Timber, T
import matplotlib.pyplot as plt
from typing import Type, TypeVar

Num = TypeVar('Num', int, float)

class TimberSection():
    def __init__(self, B: Num, H: Num, Ned: Num=0.0, Med: Num=0.0, Ved: Num=0.0, design_code=en1995_1_1,
                 fire_protection_right=None, fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):
        """

        @param material: materiaali
        @param B: leveys
        @param H: korkeus
        @param Ned: ?
        @param Med: ?
        @param Ved: ?
        @param design_code: kantsis olla 1995_1_1
        @param fire_protection_right: palosuoja joka laitetaan oikealle, yliajaa fire_protection_genericin
        @param fire_protection_left: palosuoja joka laitetaan vasemmalle, yliajaa fire_protection_genericin
        @param fire_protection_top: palosuoja joka laitetaan ylös, yliajaa fire_protection_genericin
        @param fire_protection_bottom: palosuoja joka laitetaan alas, yliajaa fire_protection_genericin
        @param fire_protection_generic: geneerinen palosuoja joka annetaan monelle sivulle samaan aikaan
        @param fire_protection_sides: jos antaa 3 niin fire_protection_top jää tyhjäksi.
                                      jos antaa jonkun muun kuin 3 niin suojaa kaikki neljä sivua fire_protection_genericillä
        @param sides_on_fire: jos antaa 4 niin poltetaan kaikilta neljältä sivulta,
                              jos antaa jonkun muun kuin 4 niin poltetaan kolmelta sivulta. (ylhäältä ei polteta)
        """
        self.__material = None
        self.__H = H
        self.__B = B
        print(self.__B)
        self.Ned = Ned
        self.Med = Med
        self.Ved = Ved
        self.code = design_code

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


    @property
    def material(self):
        return self.__material

    def set_beta_0_and_n(self):
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self.__B, self.__H))

    @material.setter
    def material(self, val):
        self.__material = val
        self.set_beta_0_and_n()

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
    def H(self):
        h = self.__H
        h -= self.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= self.d_char_n(self.fire_protection_top, self.R, self.beta_n)
        return h

    @H.setter
    def H(self, val):
        self.__H = val

    def get_H(self, at):
        return self.H

    def get_A(self, at):
        return self.A

    @property
    def B(self):
        b = self.__B
        b -= self.d_char_n(self.fire_protection_right, self.R, self.beta_n)
        # TODO nurkkapilarit tarkistettava
        b -= self.d_char_n(self.fire_protection_left, self.R, self.beta_n)
        return b

    @B.setter
    def B(self, val):
        self.__B = val

    @property
    def I(self):
        return np.asarray([(self.B * self.H ** 3) / 12, (self.H * self.B ** 3) / 12])

    def get_I(self, at):
        return self.I

    @property
    def Wel(self):
        return np.asarray([(self.B * self.H ** 2) / 6, (self.H * self.B ** 2) / 6])

    def get_Wel(self, at):
        return self.Wel

    @property
    def k_h(self):
        # TODO pitää kysyä samilta
        if self.material.rhok <= 700 and self.material.type == 'solid_timber':
            if self.Ned > 0 and self.B < 150 and self.Med != 0 and self.H < 150:
                return min(min(1.3, (150 / self.B) ** 0.2), min(1.3, (150 / self.H) ** 0.2))
            elif self.Ned > 0 and self.B < 150:
                return min(1.3, (150 / self.B) ** 0.2)
            elif self.Med != 0 and self.H < 150:
                return min(1.3, (150 / self.H) ** 0.2)
            else:
                return 1.0
        elif self.material.type == 'glt':
            if self.Ned > 0 and self.B < 600 and self.Med != 0 and self.H < 600:
                return min(min(1.1, (600 / self.B) ** 0.1), min(1.1, (600 / self.H) ** 0.1))
            elif self.Ned > 0 and self.B < 600:
                return min(1.1, (600 / self.B) ** 0.1)
            elif self.Med != 0 and self.H < 600:
                return min(1.1, (600 / self.H) ** 0.1)
            else:
                return 1.0
        elif self.material.type == 'lvl':
            if self.Med != 0 and self.H != 300:
                return min(1.2, (300 / self.H) ** self.material.s)
            else:
                return 1.0
        else:
            return 1.0

    def __repr__(self):
        return f"{self.material.timber_type} {self.H:.0f}X{self.B:.0f}"


    #########################  SFS-EN 1995-1-2 FIRE DESIGN ##########################

    @property
    def A(self):
        '''
        Laskee palolle altistuneen rakenteen poikkileikkauksen tehollisen pinta-alan palonkestoajan funktiona
        @return: tehollinen pinta-ala A_ef [mm^2]
        '''

        b = self.B
        h = self.H

        if b <= 0 or h <= 0:
            return 0
        return b*h

    @property
    def fire_protection_generic(self):
        return self.__fire_protection_generic

    @fire_protection_generic.setter
    def fire_protection_generic(self, val):
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

    def d_char_n(self, fp: Type[FireProtection], t: Num, beta_n: float) -> float:
        """
        Returns the burned depth of single side
        @param fp: fire protection on the side
        @param t: how long the fire has burned
        @param beta_n: charring speed
        @return: depth of charring
        """

        k0 = 1.0
        d0 = 7
        if fp:
            if fp.t_ch <= 20:
                if t < 20:
                    k0 = t / 20
            else:
                if t <= fp.t_ch:
                    k0 = t/fp.t_ch
        else:
            if t < 20:
                k0 = t / 20

        if not fp:
            return t * beta_n + k0 * d0
        if t < fp.t_ch:
            return 0
        if fp.t_ch < fp.t_f:
            if t < fp.t_f:
                return (t - fp.t_ch) * 0.73 * beta_n + k0 * d0
            elif fp.t_f < t < fp.t_a:
                return (fp.t_f - fp.t_ch) * 0.73 * beta_n + (t - fp.t_f) * 2 * beta_n + k0 * d0
            else:
                return (fp.t_f - fp.t_ch) * 0.73 * beta_n + (fp.t_a - fp.t_f) * 2 * beta_n + (t - fp.t_a) * beta_n + k0 * d0
        if t < fp.t_a:
            dt = t - fp.t_f
            return 2 * beta_n * dt + k0 * d0
        else:
            t_fa = fp.t_a - fp.t_f
            dt = t - fp.t_a
            return 2 * beta_n * t_fa + beta_n * dt + k0 * d0

    def get_R(self):
        return self.R

    def __str__(self):
        return f'TimberSection({self.__B}X{self.__H})'


class TaperedSection(ABC):
    def get_A(self, at):
        pass

    def get_H(self, at):
        pass

    def get_I(self, at):
        pass

    def get_Wel(self, at):
        pass

class PulpettiPalkki(TimberSection, TaperedSection):
    def __init__(self, B, H, flipped=False, fire_protection_right=None, fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):
        self.__H = H
        self.__B = B
        self.h0 = H[0]
        self.hap = H[1]
        self.flipped = flipped
        super().__init__(B, H, fire_protection_right=fire_protection_right, fire_protection_left=fire_protection_left,
                                             fire_protection_top=fire_protection_top, fire_protection_bottom=fire_protection_bottom,
                                             fire_protection_generic=fire_protection_generic, fire_protection_sides=fire_protection_sides,
                                             sides_on_fire=sides_on_fire)

    @property
    def H(self):
        return self.get_H(0)
        raise Exception('Use get_H(at) method to get height at specific point')

    @H.setter
    def H(self, val):
        if val is not list:
            raise Exception('value must be a list of two numbers [h1, h2]')
        self.__H = val

    def get_H(self, x: Num):
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """
        if not isinstance(self.__H, list):
            return self.__H
        k = self.__H[1] - self.__H[0]
        h = self.__H[0] + k * x
        h -= self.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= self.d_char_n(self.fire_protection_top, self.R, self.beta_n)
        return h

    @property
    def A(self):
        return self.get_A(0)
        raise Exception('Use get_A(at) method to get area at specific point')

    def get_A(self, x):
        """
        Calculates A

        @param x: local coordinate 0 ... 1
        @return: area at location
        """
        b = self.B
        h = self.get_H(x)

        if b <= 0 or h <= 0:
            return 0
        return b * h

    @property
    def I(self):
        return self.get_I(0)
        raise Exception('Use get_I(at) method to get I at specific point')

    def get_I(self, x):
        return np.asarray([(self.B * self.get_H(x) ** 3) / 12, (self.get_H(x) * self.B ** 3) / 12])

    @property
    def Wel(self):
        return self.get_Wel(0)
        raise Exception('Use get_Wel(at) method to get Wel at specific point')

    def get_Wel(self, x):
        return np.asarray([(self.B * self.get_H(x) ** 2) / 6, (self.get_H(x) * self.B ** 2) / 6])

    def set_beta_0_and_n(self):
        # TODO betat x:n suhteen
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self.__B, self.__H[0]))

    def __str__(self):
        return f'Pulpettipalkki({self.__B}X{self.__H[0]}X{self.__H[1]})'

class HarjaPalkki(TimberSection, TaperedSection):
    def __init__(self, B, H, k_p=0.02, k_l=1.13, flipped=False, fire_protection_right=None, fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):
        self.__material = None
        self.__H = H
        self.__B = B
        if len(self.__H) == 3:
            self.h0 = H[0]
            self.hap = H[1]
            self.xap = H[2]
        elif len(self.__H) == 2:
            self.h0 = H[0]
            self.hap = H[1]
            self.xap = 0.5
        self.k_p = k_p
        self.k_l = k_l
        self.flipped = flipped
        super().__init__(B, H, fire_protection_right=fire_protection_right, fire_protection_left=fire_protection_left,
                        fire_protection_top=fire_protection_top, fire_protection_bottom=fire_protection_bottom,
                        fire_protection_generic=fire_protection_generic, fire_protection_sides=fire_protection_sides,
                        sides_on_fire=sides_on_fire)



    def hapa(self, x):
        """
        Defines harjapalkki
        @param x: Harjapalkki's longitudinal coordinate
        @param h0: Section height at the beginning and at the end of the beam
        @param hap: Section height at the paramount point of the beam
        @param xap: X coordinate of the highest cross section of the beam
        @param L: Beam span length (support node to support node)
        @return:
        """

        if not isinstance(self.__H, list):
            return self.__H
        h = (((self.hap - self.h0) / self.xap) * x + self.h0 - max(
            ((self.hap - self.h0) / self.xap) * x + self.h0 - ((self.h0 - self.hap) / (1 - self.xap) * x + self.hap +
                                                               ((self.hap - self.h0) / (1 - self.xap)) * self.xap), 0))
        return h

    @property
    def H(self):
        return self.get_H(0)
        raise Exception('Use get_H(at) method to get height at specific point')

    @H.setter
    def H(self, val):
        if val is not list:
            raise Exception('value must be a list of two numbers [h1, h2]')
        self.__H = val

    def get_H(self, x: Num):
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """
        if not isinstance(self.__H, list):
            return self.__H

        h = self.hapa(x)

        h -= self.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= self.d_char_n(self.fire_protection_top, self.R, self.beta_n)
        return h

    def get_A(self, x):
        """
        Calculates A

        @param x: local coordinate 0 ... 1
        @return: area at location
        """
        b = self.B
        h = self.get_H(x)

        if b <= 0 or h <= 0:
            return 0
        return b * h

    def get_Wel(self, x):
        return np.asarray([(self.B * self.get_H(x) ** 2) / 6, (self.get_H(x) * self.B ** 2) / 6])

    def get_I(self, x):
        return np.asarray([(self.B * self.get_H(x) ** 3) / 12, (self.get_H(x) * self.B ** 3) / 12])

    def set_beta_0_and_n(self):
        # TODO betat x:n suhteen
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self.__B, self.__H[0]))

    def __str__(self):
        return f'Harjapalkki({self.__B}X{self.__H[0]}X{self.__H[1]})'

class KaannettyHarjaPalkki:
    def __init__(self, B, H):
        pass


class KaarevaPalkki:
    def __init__(self, B, H):
        pass


class KaarevaHarjaPalkki(TimberSection, TaperedSection):
    def __init__(self, B, H, k_p=0.03, k_l=1.113, r_in=13, t=45, fire_protection_right=None, fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):
        self.__H = H
        self.__B = B
        if len(self.__H) == 4:
            self.h0 = H[0]
            self.hap = H[1]
            self.hp = H[2]
            self.xap = H[3]
        elif len(self.__H) == 3:
            self.h0 = H[0]
            self.hap = H[1]
            self.hp = H[2]
            self.xap = 0.5
        self.k_p = k_p
        self.k_l = k_l
        self.r_in = r_in
        self.t = t
        super().__init__(B, H, fire_protection_right=fire_protection_right, fire_protection_left=fire_protection_left,
                        fire_protection_top=fire_protection_top, fire_protection_bottom=fire_protection_bottom,
                        fire_protection_generic=fire_protection_generic, fire_protection_sides=fire_protection_sides,
                        sides_on_fire=sides_on_fire)


    def kahapa(self, x):
        """
        Defines harjapalkki
        @param x: Harjapalkki's longitudinal coordinate
        @param h0: Section height at the beginning and at the end of the beam
        @param hap: Section height at the paramount point of the beam
        @param xap: X coordinate of the highest cross section of the beam
        @param L: Beam span length (support node to support node)
        @return:
        """
        if len(self.__H) == 4:
            h0 = self.__H[0]
            hap = self.__H[1]
            hp = self.__H[2]
            xap = self.__H[3]
        elif len(self.__H) == 3:
            h0 = self.__H[0]
            hap = self.__H[1]
            hp = self.__H[2]
            xap = 0.5
        else:
            return self.__H
        # TODO laske korkeus oikein
        # h = (((hap - h0) / xap) * x + h0 - max(
        #     ((hap - h0) / xap) * x + h0 - ((h0 - hap) / (1 - xap) * x + hap + ((hap - h0) / (1 - xap)) * xap), 0)) -\
        #     (2 * hp * x ** 2 + 2 * hp * x)
        h = (((self.hap - self.h0) / self.xap) * x + self.h0 - max(
            ((self.hap - self.h0) / self.xap) * x + self.h0 - ((self.h0 - self.hap) / (1 - self.xap) * x + self.hap +
                                                               ((self.hap - self.h0) / (1 - self.xap)) * self.xap), 0))
        return h

    @property
    def H(self):
        return self.get_H(0)
        raise Exception('Use get_H(at) method to get height at specific point')

    @H.setter
    def H(self, val):
        if val is not list:
            raise Exception('value must be a list of two numbers [h1, h2]')
        self.__H = val

    def get_H(self, x: Num):
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """
        if not isinstance(self.__H, list):
            return self.__H

        h = self.kahapa(x)

        h -= self.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= self.d_char_n(self.fire_protection_top, self.R, self.beta_n)
        return h

    def get_A(self, x):
        """
        Calculates A

        @param x: local coordinate 0 ... 1
        @return: area at location
        """
        b = self.B
        h = self.get_H(x)

        if b <= 0 or h <= 0:
            return 0
        return b * h

    def get_Wel(self, x):
        return np.asarray([(self.B * self.get_H(x) ** 2) / 6, (self.get_H(x) * self.B ** 2) / 6])

    def get_I(self, x):
        return np.asarray([(self.B * self.get_H(x) ** 3) / 12, (self.get_H(x) * self.B ** 3) / 12])

    def set_beta_0_and_n(self):
        # TODO betat x:n suhteen
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self.__B, self.__H[0]))

    def __str__(self):
        return f'KaarevaHarjaPalkki({self.__B}X{self.__H[0]}X{self.__H[1]})'

class SmallElementSection:
    def __init__(self, A, I):
        self.A = A
        self.I = I

if __name__ == '__main__':
    '''
    TODO tehollinen pinta-ala pienenee millä tahansa palonketoajalla 
    Katso RIL 205-2-2009 lauseke 4.1.1 s. 32
    
    '''
    t = KaarevaHarjaPalkki(240, [540, 1100, 300])
    t.material = Timber()
    for i in np.arange(1 + 0.1, step=0.1):
        print(t.get_H(i))
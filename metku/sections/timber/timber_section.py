"""
Timber section parameters, including timber fire design calculations (simple method)

@author: Viktor Haimi
"""

from abc import ABC, abstractmethod
import numpy as np
from scipy import integrate
from eurocodes.en1995 import en1995_1_1, en1995_1_2
from eurocodes.en1995.fire_protection import *
from materials.timber_data import Timber, T
import matplotlib.pyplot as plt
from typing import Type, TypeVar

Num = TypeVar('Num', int, float)

class TimberSection():
    def __init__(self, B: Num, H: Num, material=T.C14, Ned: Num=0.0, Med: Num=0.0, Ved: Num=0.0, design_code=en1995_1_1,
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
        self.__H = H
        self.__B = B

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

        self.material = material

    @property
    def material(self):
        return self.__material

    def set_beta_0_and_n(self):
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self.__B, self.__H))

    @material.setter
    def material(self, val):
        self.__material = Timber(val)
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
    def __init__(self):
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

    @abstractmethod
    def get_A(self, at):
        pass

    @abstractmethod
    def get_H(self, at):
        pass

    @abstractmethod
    def get_I(self, at):
        pass

    @abstractmethod
    def get_Wel(self, at):
        pass


class PulpettiPalkki(TimberSection, TaperedSection):
    def __init__(self, B, H0, Hap, material=T.C14, flipped=False, fire_protection_right=None, fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):
        self.__H = H0
        self.__B = B
        self.h0 = H0
        self.hap = Hap
        self.flipped = flipped
        TimberSection.__init__(self, B, H0, material=material, fire_protection_right=fire_protection_right, fire_protection_left=fire_protection_left,
                         fire_protection_top=fire_protection_top, fire_protection_bottom=fire_protection_bottom,
                         fire_protection_generic=fire_protection_generic, fire_protection_sides=fire_protection_sides,
                         sides_on_fire=sides_on_fire)
        TaperedSection.__init__(self)

    @property
    def H(self):
        return self.get_H(0)
        raise Exception('Use get_H(at) method to get height at specific point')

    @H.setter
    def H(self, val):
        if val is not list:
            raise Exception('value must be a list of two numbers [h1, h2]')
        self.__H = val

    @property
    def H0(self):
        return self.h0

    @property
    def Hap(self):
        return self.hap

    @H0.setter
    def H0(self, val):
        self.h0 = val
        self.update_elements()

    @Hap.setter
    def Hap(self, val):
        self.hap = val
        self.update_elements()

    def get_H(self, x: Num):
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """

        k = self.hap - self.h0
        h = self.h0 + k * x
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
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self.__B, self.h0))

    def __str__(self):
        return f'Pulpettipalkki({self.__B}X{self.h0}X{self.hap})'

class HarjaPalkki(TimberSection, TaperedSection):
    def __init__(self, B, H0, Hap, xap=0.5, material=T.C14, k_p=0.02, k_l=1.13, flipped=False, fire_protection_right=None, fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):
        self.__material = None
        self.__H = H0
        self.__B = B
        self.h0 = H0
        self.hap = Hap
        self.xap = xap
        self.k_p = k_p
        self.k_l = k_l
        self.flipped = flipped
        TimberSection.__init__(self, B, H0, material=material, fire_protection_right=fire_protection_right, fire_protection_left=fire_protection_left,
                        fire_protection_top=fire_protection_top, fire_protection_bottom=fire_protection_bottom,
                        fire_protection_generic=fire_protection_generic, fire_protection_sides=fire_protection_sides,
                        sides_on_fire=sides_on_fire)
        TaperedSection.__init__(self)


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

    @property
    def H0(self):
        return self.h0

    @property
    def Hap(self):
        return self.hap

    @property
    def Xap(self):
        return self.xap

    @H0.setter
    def H0(self, val):
        self.h0 = val
        self.update_elements()

    @Hap.setter
    def Hap(self, val):
        self.hap = val
        self.update_elements()

    @Xap.setter
    def Xap(self, val):
        self.xap = val
        self.update_elements()

    def get_H(self, x: Num):
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """

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
                                                             self.material.hardness, min(self.__B, self.h0))

    def __str__(self):
        return f'Harjapalkki({self.__B}X{self.h0}X{self.hap})'


class KaarevaPalkki(TimberSection, TaperedSection):
    def __init__(self, B, H0, alpha, beta, r_in, material=T.C14, t=45, fire_protection_right=None,
                 fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):

        self.__H = H0
        self.__B = B
        self.h0 = H0
        self.parent = None
        self.alpha = np.pi * alpha / 180
        self.beta = np.pi * beta / 180
        self.r_in = r_in
        self.t = t

        TimberSection.__init__(self, B, H0, material=material, fire_protection_right=fire_protection_right,
                         fire_protection_left=fire_protection_left,
                         fire_protection_top=fire_protection_top, fire_protection_bottom=fire_protection_bottom,
                         fire_protection_generic=fire_protection_generic, fire_protection_sides=fire_protection_sides,
                         sides_on_fire=sides_on_fire)
        TaperedSection.__init__(self)

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

    def f1(self, x):
        return np.tan(self.alpha) * x + self.h0

    def f2(self, x):
        return np.tan(self.beta) * x

    def f3(self, x):
        return (self.parent.length - x) * np.tan(self.alpha) + self.h0

    def f4(self, x):
        return 2 * self.L_y - np.tan(self.beta) * x

    def f5(self, x):
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_up

    def f6(self, x):
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_dw

    @property
    def H(self):
        return self.get_H(0)

    @H.setter
    def H(self, val):
        if val is not list:
            raise Exception('value must be a list of five or four numbers [h1, h2]')
        self.__H = val

    @property
    def H0(self):
        return self.h0

    @property
    def Hap(self):
        return self.hap

    @H0.setter
    def H0(self, val):
        self.h0 = val
        self.update_circles()
        self.update_elements()

    @property
    def Alpha(self):
        return self.alpha

    @Alpha.setter
    def Alpha(self, val):
        self.alpha = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def Beta(self):
        return self.beta

    @Beta.setter
    def Beta(self, val):
        self.beta = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def R_in(self):
        return self.r_in

    @R_in.setter
    def R_in(self, val):
        self.r_in = val
        self.update_circles()
        self.update_elements()

    def get_H(self, x: Num):
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """

        x *= self.parent.length
        if 0 <= x < self.up_x1:
            h = self.f1(x) - self.f2(x)
        elif self.up_x1 <= x < self.dwn_x1:
            h = self.f5(x) - self.f2(x)
        elif self.dwn_x1 <= x < self.dwn_x2:
            h = self.f5(x) - self.f6(x)
        elif self.dwn_x2 <= x < self.up_x2:
            h = self.f5(x) - self.f4(x)
        elif self.up_x2 <= x:
            h = self.f3(x) - self.f4(x)

        h -= self.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= self.d_char_n(self.fire_protection_top, self.R, self.beta_n)
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

    def plot(self):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Arc

        x1 = [0, self.up_x1]
        x2 = [self.up_x2, self.parent.length]
        x3 = [0, self.dwn_x1]
        x4 = [self.dwn_x2, self.parent.length]
        y1 = [self.h0,  self.up_y1]
        y2 = [self.up_y2, self.h0]
        y3 = [0, self.dwn_y1]
        y4 = [self.dwn_y2, 0]

        cx1, cy1 = self.x_circle, self.y_circle_up
        cx2, cy2 = self.x_circle, self.y_circle_dw
        a, b = 2 * self.r_in, 2 * self.r_in
        # Initialise the figure and axes.
        fig, ax = plt.subplots(1, figsize=(12, 4))

        # Set the title for the figure
        fig.suptitle('Kaareva palkki', fontsize=15)

        # Draw all the lines in the same plot, assigning a label for each one to be
        # shown in the legend.
        ax.plot(x1, y1, color="b")
        ax.plot(x2, y2, color="b")
        ax.plot(x3, y3, color="b")
        ax.plot(x4, y4, color="b")

        ax.add_patch(Arc((cx1, cy1), a, b, theta1=90-180*self.alpha/np.pi, theta2=90+180*self.alpha/np.pi, edgecolor='b', lw=1.5))
        ax.add_patch(Arc((cx2, cy2), a, b, theta1=90-180*self.beta/np.pi, theta2=90+180*self.beta/np.pi, edgecolor='b', lw=1.5))

        plt.xlim(-1000, self.parent.length + 1000)
        plt.ylim(-1000, 5000)
        ax.set_aspect(0.7)
        plt.grid(linestyle='--')

        plt.show()

    def set_beta_0_and_n(self):
        # TODO betat x:n suhteen
        self.beta_0, self.beta_n = en1995_1_2.charring_speed(self.material.rhok, self.material.type,
                                                             self.material.hardness, min(self.__B, self.__H))


class KaarevaHarjaPalkki(TimberSection, TaperedSection):
    def __init__(self, B, H0, alpha, beta, r_in, material=T.C14, t=45, fire_protection_right=None, fire_protection_left=None,
                 fire_protection_top=None, fire_protection_bottom=None,
                 fire_protection_generic=None, fire_protection_sides=3, sides_on_fire=3):

        self.__H = H0
        self.__B = B
        self.h0 = H0
        self.parent = None
        self.alpha = np.pi * alpha / 180
        self.beta = np.pi * beta / 180
        self.r_in = r_in
        self.t = t
        TimberSection.__init__(self, B, H0, material=material, fire_protection_right=fire_protection_right, fire_protection_left=fire_protection_left,
                        fire_protection_top=fire_protection_top, fire_protection_bottom=fire_protection_bottom,
                        fire_protection_generic=fire_protection_generic, fire_protection_sides=fire_protection_sides,
                        sides_on_fire=sides_on_fire)
        TaperedSection.__init__(self)
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
        self.up_x2, self.up_y2 = [self.parent.length / 2 - self.r_in * np.sin(self.alpha) + 2 * (self.r_in * np.sin(self.alpha)),
                                  (self.parent.length / 2 - self.r_in * np.sin(self.alpha)) * np.tan(self.alpha) + self.h0]
        self.dwn_x2, self.dwn_y2 = [self.parent.length / 2 - self.r_in * np.sin(self.beta) + 2 * (self.r_in * np.sin(self.beta)),
                                    (self.parent.length / 2 - self.r_in * np.sin(self.beta)) * np.tan(self.beta)]
        self.hap = self.get_H(0.5)

    def f1(self, x):
        return np.tan(self.alpha) * x + self.h0

    def f2(self, x):
        return np.tan(self.beta) * x

    def f3(self, x):
        return (self.parent.length - x) * np.tan(self.alpha) + self.h0

    def f4(self, x):
        return 2 * self.L_y - np.tan(self.beta) * x

    def f5(self, x):
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_up

    def f6(self, x):
        return np.sqrt(self.r_in ** 2 - (self.parent.length / 2 - x) ** 2) + self.y_circle_dw

    @property
    def H(self):
        return self.get_H(0)
        raise Exception('Use get_H(at) method to get height at specific point')

    @H.setter
    def H(self, val):
        if val is not list:
            raise Exception('value must be a list of five or four numbers [h1, h2]')
        self.__H = val

    @property
    def H0(self):
        return self.h0

    @property
    def Hap(self):
        return self.hap

    @H0.setter
    def H0(self, val):
        self.h0 = val
        self.update_circles()
        self.update_elements()

    @property
    def Alpha(self):
        return self.alpha

    @Alpha.setter
    def Alpha(self, val):
        self.alpha = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def Beta(self):
        return self.beta

    @Beta.setter
    def Beta(self, val):
        self.beta = np.pi * val / 180
        self.update_circles()
        self.update_elements()

    @property
    def R_in(self):
        return self.r_in

    @R_in.setter
    def R_in(self, val):
        self.r_in = val
        self.update_circles()
        self.update_elements()

    def get_H(self, x: Num):
        """
        Getter for H value

        @param x: local coordinate 0 ... 1
        @return: height at location
        """

        x *= self.parent.length
        if 0 <= x < self.dwn_x1:
            h = self.f1(x) - self.f2(x)
        elif self.dwn_x1 <= x < self.parent.length / 2:
            h = self.f1(x) - self.f6(x)
        elif self.parent.length / 2 <= x < self.dwn_x2:
            h = self.f3(x) - self.f6(x)
        elif self.dwn_x2 <= x:
            h = self.f3(x) - self.f4(x)

        h -= self.d_char_n(self.fire_protection_bottom, self.R, self.beta_n)
        if self.sides_on_fire == 4:
            h -= self.d_char_n(self.fire_protection_top, self.R, self.beta_n)
        return h

    def A_face(self):
        if1 = integrate.quad(self.f1, 0, self.parent.length/2)
        if2 = integrate.quad(self.f2, 0, self.dwn_x1)
        if3 = integrate.quad(self.f3, self.parent.length/2, self.parent.length)
        if4 = integrate.quad(self.f4, self.dwn_x2, self.parent.length)
        if6 = integrate.quad(self.f6, self.dwn_x1, self.dwn_x2)
        a = if1[0] + if3[0] - if2[0] - if4[0] - if6[0]

        return a

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
                                                             self.material.hardness, min(self.__B, self.__H))

    def plot(self):

        import matplotlib.pyplot as plt
        from matplotlib.patches import Arc

        x1 = [0, self.parent.length/2, self.parent.length]
        x2 = [0, self.dwn_x1]
        x3 = [self.dwn_x2, self.parent.length]
        y1 = [self.h0,  self.U_y, self.h0]
        y2 = [0, self.dwn_y1]
        y3 = [self.dwn_y2, 0]

        cx2, cy2 = self.x_circle, self.y_circle_dw
        a, b = 2 * self.r_in, 2 * self.r_in
        # Initialise the figure and axes.
        fig, ax = plt.subplots(1, figsize=(12, 4))

        # Set the title for the figure
        fig.suptitle('Kaareva harjapalkki', fontsize=15)

        # Draw all the lines in the same plot, assigning a label for each one to be
        # shown in the legend.
        ax.plot(x1, y1, color="b")
        ax.plot(x2, y2, color="b")
        ax.plot(x3, y3, color="b")

        ax.add_patch(Arc((cx2, cy2), a, b, theta1=90-180*self.beta/np.pi, theta2=90+180*self.beta/np.pi, edgecolor='b', lw=1.5))

        plt.xlim(-100, self.parent.length + 100)
        plt.ylim(-100, 5000)
        ax.set_aspect(1)

        plt.grid(linestyle='--')

        plt.show()

    def __str__(self):
        return f'KaarevaHarjaPalkki({self.__B}X{self.H0}X{self.Hap})'


class SmallElementSection:
    def __init__(self, A, I):
        self.A = A
        self.I = I

    def __repr__(self):
        return f'SmallElementSection A:{self.A}, I:{self.I}'


if __name__ == '__main__':
    '''
    TODO tehollinen pinta-ala pienenee millä tahansa palonketoajalla 
    Katso RIL 205-2-2009 lauseke 4.1.1 s. 32
    '''
      ############### KAAREVAN PALKKIN SEKÄ HARJAPALKIN MÄÄRITYS ################

    import frame2d.frame2d as f2d
    t = KaarevaHarjaPalkki(240, 600, 20, 15, 10000)
    obj = type('', (), {})()
    obj.length = 16000
    t.parent = obj
    t.update_circles()
    t.plot()

    for i in np.arange(1 + 0.01, step=0.01):
        print(f'i {i} h {t.get_H(i)}')

    # i = integrate.quad(t.get_H, 0, 1)
    # #integrate.quad_explain()
    # print('integrate', i)
    #print(t.V_tot())

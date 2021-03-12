"""

"""
import numpy as np
from eurocodes.en1995 import en1995_1_1, en1995_1_2
from typing import Type, TypeVar

Num = TypeVar('Num', int, float)


class FireProtection():
    def __init__(self, t_ch: Num, b: Num, t_f: Num=None, t_a: Num=None):
        self.h_p = b
        self.t_ch = t_ch
        self.t_f = t_f
        self.t_a = t_a
        self.prot_mat = None
        self.prot_mat_beta_0 = None
        self.prot_mat_beta_n = None

    def protected_material(self, mat):
        beta_0, beta_n = en1995_1_2.charring_speed(mat.rhok, mat.type, mat.hardness, self.h_p)
        self.prot_mat = mat
        self.prot_mat_beta_0 = beta_0
        self.prot_mat_beta_n = beta_n


    def set_ta(self):
        t_a = min(2 * self.t_f, 25 / 2 * self.prot_mat_beta_n + self.t_f)
        self.t_a = t_a

    def set_tf_for_gypsumF(self):
        self.t_f = 12.5 / (0.73 * self.prot_mat_beta_n) + self.t_ch


class GypsumPlasterboardA(FireProtection):
    """
    Parameter
    """
    def __init__(self, sauma):
        if sauma < 2:
            super().__init__(22, 13, 22)
        else:
            super().__init__(13, 13, 13)


class GypsumPlasterboardF(FireProtection):
    def __init__(self, sauma):
        if sauma < 2:
            super().__init__(28, 15)
        else:
            super().__init__(19, 15)


class GypsumPlasterboardH(FireProtection):
    def __init__(self, sauma):
        if sauma < 2:
            super().__init__(11, 9, 11)
        else:
            super().__init__(2, 9, 2)


class GypsumPlasterboardAA(FireProtection):
    def __init__(self, sauma):
        if sauma < 2:
            super().__init__(40, 26, 40)
        else:
            super().__init__(31, 26, 31)


class GypsumPlasterboardFF(FireProtection):
    def __init__(self, sauma):
        if sauma < 2:
            super().__init__(61, 30)
        else:
            super().__init__(52, 30)


class GypsumPlasterboardHH(FireProtection):
    def __init__(self, sauma):
        if sauma < 2:
            super().__init__(23, 18, 23)
        else:
            super().__init__(14, 18, 14)


class GypsumPlasterboardAF(FireProtection):
    def __init__(self, sauma):
        if sauma < 2:
            super().__init__(46, 28)
        else:
            super().__init__(37, 28)


class WoodenFireProtection(FireProtection):
    def __init__(self, material, h_p):
        beta_0, beta_n = en1995_1_2.charring_speed(material.rhok, material.type, material.hardness, h_p)
        t_ch = h_p / beta_0
        super().__init__(t_ch, h_p, t_ch)

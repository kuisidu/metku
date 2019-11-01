# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 23:12:31 2018

Rectangular hollow sections

@author: kmela
"""

import math
import numpy as np

try:
    from src.eurocodes.en1993 import en1993_1_1
    from src.sections.steel.steel_section import SteelSection
except:
    from eurocodes.en1993 import en1993_1_1
    from sections.steel.steel_section import SteelSection


class RHS(SteelSection):
    """ Rectangular hollow sections """

    def __init__(self, H, B, T, fy=355):
        """ Rectangular hollow sections
            
            H -- height
            B -- width
            T -- wall thickness
            fy -- yield strength
        """
        R = RHS_outer_radius(T)
        A = RHS_area(H, B, T, R)
        I = [0.0, 0.0]
        I[0] = RHS_second_moment(H, B, T, R)
        I[1] = RHS_second_moment(B, H, T, R)
        Ashear = RHS_shear_area(A, H, B)
        Au = RHS_paint_area(H, B, R)
        Wel = [0.0, 0.0]
        Wel[0] = RHS_elastic_modulus(H, I[0])
        Wel[1] = RHS_elastic_modulus(B, I[1])
        Wpl = [0.0, 0.0]
        Wpl[0] = RHS_plastic_modulus(H, B, R, T)
        Wpl[1] = RHS_plastic_modulus(B, H, R, T)

        super().__init__(fy,
                         A,
                         I,
                         Au,
                         Wpl,
                         Wel,
                         Ashear)

        self.H = H
        self.B = B
        self.T = T
        self.Iw = 0
        self.It = self.torsional_constant()

        self.imp_factor = [en1993_1_1.buckling_curve["c"],
                           en1993_1_1.buckling_curve["c"]]
        self.imp_factor_LT_gen = en1993_1_1.buckling_curve["c"]  # TODO
        self.imp_factor_LT = en1993_1_1.buckling_curve["c"]  # TODO

    def __repr__(self):
        return f"{type(self).__name__} {self.H:.0f}X{self.B:.0f}X{self.T:.0f}"

    def update_properties(self):
        """
        Updates cross-sectional properties
        """

        try:
            self.R = RHS_outer_radius(self.T)
            self.A = RHS_area(self.H, self.B, self.T, self.R)
            I = [0.0, 0.0]
            I[0] = RHS_second_moment(self.H, self.B, self.T, self.R)
            I[1] = RHS_second_moment(self.B, self.H, self.T, self.R)
            self.I = np.asarray(I)
            self.Ashear = RHS_shear_area(self.A, self.H, self.B)
            self.Au = RHS_paint_area(self.H, self.B, self.R)
            Wel = [0.0, 0.0]
            Wel[0] = RHS_elastic_modulus(self.H, I[0])
            Wel[1] = RHS_elastic_modulus(self.B, I[1])
            self.Wel = np.asarray(Wel)
            Wpl = [0.0, 0.0]
            Wpl[0] = RHS_plastic_modulus(self.H, self.B, self.R, self.T)
            Wpl[1] = RHS_plastic_modulus(self.B, self.H, self.R, self.T)
            self.Wpl = np.asarray(Wpl)
            self.It = self.torsional_constant()

        except AttributeError:
            pass

    def __setattr__(self, key, value):

        dim_vars = ['H', 'B', 'T']
        super().__setattr__(key, value)
        # Whenever dimensions change, update cross-sectional properties
        if key in dim_vars:
            self.update_properties()


    def torsional_constant(self):
        """
        Calculated according to EN 10219-2
        
        """

        # Rc .. corner radius at center line
        Rc = self.R - 0.5 * self.T
        h = 2 * ((self.B - self.T) + (self.H - self.T)) - 2 * Rc * (
                    4 - math.pi)
        Ah = (self.B - self.T) * (self.H - self.T) - Rc ** 2 * (4 - math.pi)
        K = 2 * Ah * self.T / h
        It = (self.T ** 3 * h / 3 + 2 * K * Ah) * 1e-4

        return It


    def flange_class(self, verb=False):
        """ Determine class of compressed flange """
        cf = self.B - 2 * self.R
        rf = cf / self.T
        cFlange = en1993_1_1.internal_part_in_compression(rf, self.eps)

        return cFlange

    def web_class_comp(self, verb=False):
        """ Determine class of compressed web """
        cw = self.H - 2 * self.R
        rw = cw / self.T
        cWeb = en1993_1_1.internal_part_in_compression(rw, self.eps)

        return cWeb

    def web_class_bend(self, verb=False):
        """ Determine class of web in bending """
        cw = self.H - 2 * self.R
        rw = cw / self.T
        cWeb = en1993_1_1.internal_part_in_bending(rw, self.eps)

        return cWeb

    def web_class_comp_bend(self, Ned, verb=False):

        return 2

    def moment_axial_force_interact(self, UN, MRd):
        """ Interaction rule for section resistance for combined
            axial force and bending
            
            input: UN .. NEd/NRd
                  MRd .. moment resistance
        """
        aw = min((self.A - 2 * self.B * self.T) / self.A, 0.5)

        if UN > 1.0:
            MNRd = 0.0
        else:
            MNRd = MRd * min((1 - UN) / (1 - 0.5 * aw), 1)

        return MNRd


class SHS(RHS):
    """ Square hollow sections
        Subclass of RHS
    """

    def __init__(self, H, T, fy=355):
        RHS.__init__(self, H, H, T, fy)


# Functions for computing RHS section properties
def RHS_area(H, B, T, R):
    """ cross-sectional area """
    A = 2 * (H - 2 * R) * T + 2 * (B - 2 * R) * T + math.pi * (
                R ** 2 - (R - T) ** 2)
    return A


def RHS_paint_area(H, B, R):
    """ circumference """
    Au = 2 * (H + B) - 8 * R + 2 * math.pi * R
    return Au


def RHS_shear_area(A, H, B):
    Ashear = A * B / (H + B)
    return Ashear


def RHS_second_moment(H, B, T, R):
    """ Second moment of area """
    I1 = 2 * (1 / 12 * T * (H - 2 * R) ** 3 + 1 / 12 * (B - 2 * R) * T ** 3 + (
                0.5 * (H - T)) ** 2 * (B - 2 * R) * T)
    IR = 4 * ((math.pi / 16 - 4 / 9 / math.pi) * R ** 4 + (
                0.5 * H - R + 4 * R / 3 / math.pi) ** 2 * math.pi * R ** 2 / 4)
    IRneg = 4 * ((math.pi / 16 - 4 / 9 / math.pi) * (R - T) ** 4 + (
                0.5 * H - R + 4 * (R - T) / 3 / math.pi) ** 2 * math.pi * (
                             R - T) ** 2 / 4)
    I1 = I1 + IR - IRneg

    return I1


def RHS_elastic_modulus(H, I1):
    Wel = I1 / (0.5 * H)

    return Wel


def RHS_plastic_modulus(H, B, R, T):
    Wp = 2 * ((0.5 * H - R) * T * (0.5 * H - R)) + (B - 2 * R) * T * (H - T)
    WR = 2 * math.pi / 4 * R ** 2 * 2 * (0.5 * H - R + 4 * R / 3 / math.pi)
    Wneg = 2 * math.pi / 4 * (R - T) ** 2 * 2 * (
                0.5 * H - R + 4 * (R - T) / 3 / math.pi)
    Wp = Wp + WR - Wneg

    return Wp


def RHS_outer_radius(T):
    if T <= 6.0:
        R = 2 * T
    elif T <= 10.0:
        R = 2.5 * T
    else:
        R = 3 * T
    return R

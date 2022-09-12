#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

"""
Formulas from
https://sections.arcelormittal.com/repository2/Sections/Sections_MB_ArcelorMittal_FR_EN_DE.pdf
"""

import numpy as np


def area(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """ Area of section [mm^2]"""
    A = 2 * tf * b + (h - 2 * tf) * tw + (4 - np.pi) * r ** 2
    return A


def shear_area_z(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """ Shear area of section [mm^2]"""
    A = area(h, b, tf, tw, r)
    avz = A - 2.0 * b * tf + (tw + 2 * r) * tf
    return avz


def shear_area_y(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """ Shear area of section [mm^2]"""
    avy = 2 * b * tf
    return avy


def perimeter(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """ Perimeter of cross-section """
    au = 2 * b + 2 * (b - 2 * r - tw) + 2 * (h - 2 * r) + 2 * np.pi * r
    return au


def second_moment_of_area_y(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """ Second moment of area, I_y is stronger direction [mm^4]"""
    iy = 1.0 / 12.0 * (b * h ** 3.0 - (b - tw) * (
            h - 2.0 * tf) ** 3.0) + 0.03 * r ** 4 + 0.2146 * r ** 2.0 * (
                 h - 2.0 * tf - 0.4468 * r) ** 2.0
    return iy


def second_moment_of_area_z(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """Second moment of area, I_y is stronger direction [mm^4]"""
    iz = 1.0 / 12.0 * (2.0 * tf * b ** 3.0 + (
            h - 2.0 * tf) * tw ** 3.0) + 0.03 * r ** 4.0 + 0.2146 * r ** 2.0 * (
                 tw + 0.4468 * r) ** 2.0
    return iz


def torsion_constant(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """
    Torsional constant [mm^4]
    https://www.steelconstruction.info/images/6/6f/Sci_p385.pdf
    page 76
    """
    alpha_1 = -0.042 + 0.2204 * tw / tf + 0.1355 * r / tf - 0.0865 * r * tw / tf ** 2-0.0725 * tw ** 2 / tf ** 2
    D_1 = ((tf + r) ** 2 + tw * (r + tw / 4)) / (2 * r + tf)
    It = (2 * b * tf ** 3 + (h - 2 * tf) * tw ** 3) / 3 + 2 * alpha_1 * D_1 ** 4 - 4 * 0.105 * tf ** 4

    return It


def warping_constant(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """Warping constant [mm^6]"""
    iw = (tf * b ** 3.0) / 24.0 * (h - tf) ** 2.0
    return iw


def plastic_section_modulus_y(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """Plastic section modulus [mm^3]"""
    wply = (tw * h ** 2.0) / 4.0 + (b - tw) * (h - tf) * tf + (
            4.0 - np.pi) / 2.0 * r ** 2.0 * (h - 2.0 * tf) + (
                   3.0 * np.pi - 10.0) / 3.0 * r ** 3.0
    return wply


def plastic_section_modulus_z(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """Plastic section modulus [mm^3]"""
    wplz = (b ** 2.0 * tf) / 2.0 + (
            h - 2.0 * tf) / 4.0 * tw ** 2.0 + r ** 3.0 * (
                   10.0 / 3.0 - np.pi) + (
                   2.0 - np.pi / 2.0) * tw * r ** 2.0
    return wplz


def elastic_section_modulus_y(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """Elastic section modulus [mm^3]"""
    wely = second_moment_of_area_y(h, b, tf, tw, r) / (0.5 * h)
    return wely


def elastic_section_modulus_z(h: float, b: float, tf: float, tw: float, r: float) -> float:
    """Elastic section modulus [mm^3]"""
    welz = second_moment_of_area_z(h, b, tf, tw, r) / (0.5 * b)
    return welz

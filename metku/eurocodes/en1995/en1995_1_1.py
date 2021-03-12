# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 12:48:16 2021

Design rules of EN 1995-1-1, Design of Timber Structures

@author: Viktor Haimi
"""

from eurocodes.en1995.constants import *

kmod_data = {'solid_timber': {1: {'perm': 0.6, 'lt': 0.7, "mt": 0.8, "st": 0.9, "inst": 1.1},
                              2: {'perm': 0.6, 'lt': 0.7, "mt": 0.8, "st": 0.9, "inst": 1.1},
                              3: {'perm': 0.5, 'lt': 0.55, "mt": 0.65, "st": 0.7, "inst": 0.9}
                              },
             'glt': {1: {'perm': 0.6, 'lt': 0.7, "mt": 0.8, "st": 0.9, "inst": 1.1},
                     2: {'perm': 0.6, 'lt': 0.7, "mt": 0.8, "st": 0.9, "inst": 1.1},
                     3: {'perm': 0.5, 'lt': 0.55, "mt": 0.65, "st": 0.7, "inst": 0.9}
                     },
             'lvl': {1: {'perm': 0.6, 'lt': 0.7, "mt": 0.8, "st": 0.9, "inst": 1.1},
                     2: {'perm': 0.6, 'lt': 0.7, "mt": 0.8, "st": 0.9, "inst": 1.1},
                     3: {'perm': 0.5, 'lt': 0.55, "mt": 0.65, "st": 0.7, "inst": 0.9}
                    }
             }

kdef_data = {'solid_timber': {1: 0.6, 2: 0.8, 3: 2.0},
             'glt': {1: 0.6, 2: 0.8, 3: 2.0},
             'lvl': {1: 0.6, 2: 0.8, 3: 2.0}
             }


def kmod(material: str, service_class: int, load_duration_class: str) -> float:
    """
    SFS-EN 1995-1-1 § 3, taulukko 3.1
    @param material:                valittu materiaali
    @param service_class:           käyttöluokka
    @param load_duration_class:     aikaluokka
    @return:                        palauttaa k_mod arvon
    """
    return kmod_data[material][service_class][load_duration_class]


def kdef(material: str, service_class: int) -> float:
    """
    SFS-EN 1995-1-1 § 3, taulukko 3.2
    @param material:        valittu materiaali
    @param service_class:   käyttöluokka
    @return:                palauttaa k_def arvon
    """
    return kdef_data[material][service_class]

# SFS-EN 1995-1-1 § 6.1.6
k_m = 0.7

def k_to_d(f_k: float, kmod: float, safety_factor: float) -> float:
    return (kmod * f_k) / safety_factor

def get_gammaM(timber_type):
    if timber_type == 'solid_timber':
        return gammaM_st
    elif timber_type == 'glt':
        return gammaM_glt
    elif timber_type == 'lvl':
        return gammaM_lvl
    else:
        return gammaM_other

def E_mean_fin():
    pass

def G_mean_fin():
    pass

def K_ser_fin():
    pass


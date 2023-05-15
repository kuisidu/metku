#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2023. Metku team.
#  All rights reserved.

from dataclasses import dataclass
import numpy as np


"""
@dataclass
class SIBaseUnit:
    name: str = None
    symbol: str = None
    value: float = 1
    factor: float = 1
    s: float = 0
    m: float = 0
    kg: float = 0
    A: float = 0
    K: float = 0
    mol: float = 0
    cd: float = 0

    @property
    def units(self):
        return self.s, self.m, self.kg, self.A, self.K, self.mol, self.cd

    def as_array(self):
        return np.array(self.units)

    def __add__(self, other):
        if isinstance(other, SIBaseUnit):
            return self.value + other.value

    def __mul__(self, other):
        if isinstance(other, SIBaseUnit):
            self.value *= other.value
            self.as_array()
        return self.value * self.as_array()

    def __rmul__(self, other):
        self.value *= other
        return self.value * self.as_array()
"""

@dataclass
class Unit:
    name: str
    quantity: str
    symbol: str
    base: float

    def __mul__(self, other):
        return other * self.base

    def __rmul__(self, other):
        return other * self.base

    def __str__(self):
        return self.symbol

    def __call__(self, unit):
        return f"{unit / self.base :.2f} {self.symbol}"

"""
@dataclass
class Length:
    mm: float = SIBaseUnit(name="millimeter", symbol="mm", m=1e-3)
    m: float = SIBaseUnit(name="millimeter", symbol="mm", m=1)
"""

@dataclass
class Force:
    N: float = Unit("Newton", "force", "N", 1)
    kN: float = Unit("kiloNewton", "force", "kN", 1e3)
    MN: float = Unit("megaNewton", "force", "MN", 1e6)

"""
@dataclass
class Mass:
    g: float = Unit("gram", "mass", "g", 1e-3)
    kg: float = Unit("kilogram", "mass", "kg", 1)


@dataclass
class Parameters:
    length: Length = Length.mm
    force: Force = Force.N
    mass: Mass = Mass.kg
"""

Settings = Parameters()

if __name__ == '__main__':
    mass = 100 * Length.m
    mass2 = 1000 * Length.mm
    print(mass)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

# @Filename : cross_section
# @Date : 2022
# @Project: metku
# @AUTHOR : Jaakko Huusko

from dataclasses import dataclass


@dataclass
class CrossSection:

    material: object = None
    profile: object = None

    # Forces
    N: float = 0
    Vy: float = 0
    Vz: float = 0
    MT: float = 0
    My: float = 0
    Mz: float = 0

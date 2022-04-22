# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
@author Viktor Haimi
"""

from enum import Enum


class LoadIDs(Enum):
    SLS_Characteristic = 1
    SLS_Frequent = 2
    SLS_Quasi_permanent = 3
    ULS = 4
    ACC = 5


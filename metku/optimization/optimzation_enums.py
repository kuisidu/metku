#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

from enum import Enum, auto, StrEnum

class ObjectiveTypeEnum(StrEnum):
    MIN = auto()
    MAX = auto()

class ObjectiveEnum(Enum):
    WEIGHT = 1
    COST = 2

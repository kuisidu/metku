#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

# @Filename : test_steemember
# @Date : 2022
# @Project: metku
# @AUTHOR : Jaakko Huusko
import unittest

from metku.structures.steel.steel_member import SteelMember
from metku.sections.steel.ISection import HEA
from metku.sections.steel.wi import WISection


class MyTestCase(unittest.TestCase):


    def test_Mcr(self):
        section = WISection(h=1000, tw=6, b=300, tf=16, fy=355)
        smem = SteelMember(section, length=10_000)
        mcr = smem.mcrit(C=[1.132, 0.459, 0.525], k=[1,1], za=500)
        self.assertAlmostEqual(mcr / 1e0, 601.8)

    def test_MbRd(self):
        section = WISection(h=1000, tw=6, b=300, tf=16, fy=355)
        smem = SteelMember(section, length=10_000)
        smem.add_section(myed=700e6)
        mcr = smem.mcrit(C=[1.132, 0.459, 0.525], k=[1,1], za=500)
        MbRd = smem.LT_buckling_strength(mcr,
                                         method="general")
        self.assertAlmostEqual(MbRd / 1e6, 402.8, delta=1)



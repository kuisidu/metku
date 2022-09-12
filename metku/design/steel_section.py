#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

# @Filename : steel_section
# @Date : 2022
# @Project: metku
# @AUTHOR : Jaakko Huusko

from dataclasses import dataclass

import metku.eurocodes.en1993.en1993_1_1 as EC3_1
from metku.design.cross_section import CrossSection


@dataclass
class SteelSection(CrossSection):

    @property
    def NcRd(self) -> float:
        """
        Compression strength of cross-section
        :return: float [N]
        """
        return EC3_1.compression_resistance(self.profile.A, self.material.fy)

    @property
    def NtRd(self) -> float:
        """
        Tensile strength of cross-section
        :return: float [N]
        """
        return EC3_1.tension_resistance(self.profile.A, self.material.fy)

    @property
    def VyRd(self) -> float:
        """
        Shear strength about y-axis
        :return:
        """
        return EC3_1.shear_resistance(self.profile.Avy, self.material.fy)

    @property
    def VzRd(self) -> float:
        """
        Shear strength about z-axis
        :return:
        """
        return EC3_1.shear_resistance(self.profile.Avz, self.material.fy)

    @property
    def MplyRd(self) -> float:
        """
        Plastinc bending moment strength about y-axis
        :return:
        """
        return EC3_1.bending_resistance(self.profile.Wply, self.material.fy)

    @property
    def MelyRd(self) -> float:
        """
        Elastic bending moment strength about y-axis
        :return:
        """
        return EC3_1.bending_resistance(self.profile.Wely, self.material.fy)

    @property
    def MplzRd(self) -> float:
        """
        Plastinc bending moment strength about y-axis
        :return:
        """
        return EC3_1.bending_resistance(self.profile.Wplz, self.material.fy)

    @property
    def MelzRd(self) -> float:
        """
        Elastic bending moment strength about y-axis
        :return:
        """
        return EC3_1.bending_resistance(self.profile.Welz, self.material.fy)

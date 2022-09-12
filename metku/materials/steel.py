# -*- coding: utf-8 -*-
#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright 2022. Metku team
#  All rights reserved.
#  Author(s): Jaakko Huusko

from dataclasses import dataclass, field

@dataclass
class SteelMaterial:
    """
    Steel material class

    Default values:
    ---------------
        fy: None, [MPa]
        fu: None, [MPa]
        fy_40: None, [MPa]
        nu: 0.3, []
        E: 210, [GPa]
        rho: 7850 [kg/m³]
    """
    fy: float # MPa
    fu: float # MPa
    fy_40: float = None # Yield strength when thickness over 40 mm
    nu: float = 0.3
    E: float = 210  # GPa
    rho: float = 7850  # kg / m³

    def __post_init__(self):
        if self.fy_40 is None:
            self.fy_40 = self.fy

class Steel:
    """
    Default steel materials
    Custom steel material can be created using SteelMaterial -class
    """
    # ---------- S235 -----------
    S235 = SteelMaterial(fy=235, fu=360)

    # ---------- S275 -----------
    S275 = SteelMaterial(fy=275, fu=430)

    # ---------- S355 -----------
    S355 = SteelMaterial(fy=355, fu=510)
    S355J2H = SteelMaterial(fy=355, fu=510)

    # ---------- S420 -----------
    S420 = SteelMaterial(fy=420, fu=500)


class StainlessSteel:
    """
    Default stainless steel materials
    Custom steel material can be created using SteelMaterial -class
    """
    EN14301_coldstrip = SteelMaterial(fy=230, fu=540, E=200)
    EN14301_hotsrip = SteelMaterial(fy=210, fu=520, E=200)
    EN14301_hotplate = SteelMaterial(fy=210, fu=520, E=200)
    EN14301_sections = SteelMaterial(fy=190, fu=500, E=200)


if __name__ == '__main__':
    s235 = Steel.S235
    print(help(Steel))

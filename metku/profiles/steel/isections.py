#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

from dataclasses import dataclass

import numpy as np

import metku.profiles.steel.isection_formulas as iformulas


@dataclass
class ISection:
    name: str = None
    h: float = 0
    b: float = 0
    tw: float = 0
    tf: float = 0
    r: float = 0

    def __str__(self):
        return f"ISections.{self.name}"

    def __repr__(self):
        return f"ISection({self.name}, {self.h}, {self.b}, {self.tw}, {self.tf}, {self.r})"

    def __copy__(self):
        return ISection(self.name, self.h, self.b, self.tw, self.tf, self.r)

    @property
    def dimensions(self):
        return self.h, self.b, self.tf, self.tw, self.r

    @property
    def A(self):
        return iformulas.area(*self.dimensions)

    @property
    def Avy(self):
        return iformulas.shear_area_y(*self.dimensions)

    @property
    def Avz(self):
        return iformulas.shear_area_z(*self.dimensions)

    @property
    def P(self):
        return iformulas.perimeter(*self.dimensions)

    @property
    def Iy(self):
        return iformulas.second_moment_of_area_y(*self.dimensions)

    @property
    def Iz(self):
        return iformulas.second_moment_of_area_z(*self.dimensions)

    @property
    def I(self):
        return np.array([self.Iy, self.Iz])

    @property
    def It(self):
        return iformulas.torsion_constant(*self.dimensions)

    @property
    def Iw(self):
        return iformulas.warping_constant(*self.dimensions)

    @property
    def Wply(self):
        return iformulas.plastic_section_modulus_y(*self.dimensions)

    @property
    def Wely(self):
        return iformulas.elastic_section_modulus_y(*self.dimensions)

    @property
    def Wplz(self):
        return iformulas.plastic_section_modulus_z(*self.dimensions)

    @property
    def Welz(self):
        return iformulas.elastic_section_modulus_z(*self.dimensions)

    @property
    def iy(self):
        return np.sqrt(self.Iy / self.A)

    @property
    def iz(self):
        return np.sqrt(self.Iz / self.A)


class ISections(ISection):
    # section properties according to:
    # https://eurocodeapplied.com/design/en1993/ipe-hea-heb-hem-design-properties

    # ----------- IPE -------------
    IPE80: ISection = ISection(name="IPE80", h=80.0, b=46.0, tw=3.8, tf=5.2, r=5.0)
    IPE100: ISection = ISection(name="IPE100", h=100.0, b=55.0, tw=4.1, tf=5.7, r=7.0)
    IPE120: ISection = ISection(name="IPE120", h=120.0, b=64.0, tw=4.4, tf=6.3, r=7.0)
    IPE140: ISection = ISection(name="IPE140", h=140.0, b=73.0, tw=4.7, tf=6.9, r=7.0)
    IPE160: ISection = ISection(name="IPE160", h=160.0, b=82.0, tw=5.0, tf=7.4, r=9.0)
    IPE180: ISection = ISection(name="IPE180", h=180.0, b=91.0, tw=5.3, tf=8.0, r=9.0)
    IPE200: ISection = ISection(name="IPE200", h=200.0, b=100.0, tw=5.6, tf=8.5, r=12.0)
    IPE220: ISection = ISection(name="IPE220", h=220.0, b=110.0, tw=5.9, tf=9.2, r=12.0)
    IPE240: ISection = ISection(name="IPE240", h=240.0, b=120.0, tw=6.2, tf=9.8, r=15.0)
    IPE270: ISection = ISection(name="IPE270", h=270.0, b=135.0, tw=6.6, tf=10.2, r=15.0)
    IPE300: ISection = ISection(name="IPE300", h=300.0, b=150.0, tw=7.1, tf=10.7, r=15.0)
    IPE330: ISection = ISection(name="IPE330", h=330.0, b=160.0, tw=7.5, tf=11.5, r=18.0)
    IPE360: ISection = ISection(name="IPE360", h=360.0, b=170.0, tw=8.0, tf=12.7, r=18.0)
    IPE400: ISection = ISection(name="IPE400", h=400.0, b=180.0, tw=8.6, tf=13.5, r=21.0)
    IPE450: ISection = ISection(name="IPE450", h=450.0, b=190.0, tw=9.4, tf=14.6, r=21.0)
    IPE500: ISection = ISection(name="IPE500", h=500.0, b=200.0, tw=10.2, tf=16.0, r=21.0)
    IPE550: ISection = ISection(name="IPE550", h=550.0, b=210.0, tw=11.1, tf=17.2, r=24.0)
    IPE600: ISection = ISection(name="IPE600", h=600.0, b=220.0, tw=12.0, tf=19.0, r=24.0)
    # ----------- HEA -------------
    HEA100: ISection = ISection(name="HEA100", h=96.0, b=100.0, tw=5.0, tf=8.0, r=12.0)
    HEA120: ISection = ISection(name="HEA120", h=114.0, b=120.0, tw=5.0, tf=8.0, r=12.0)
    HEA140: ISection = ISection(name="HEA140", h=133.0, b=140.0, tw=5.5, tf=8.5, r=12.0)
    HEA160: ISection = ISection(name="HEA160", h=152.0, b=160.0, tw=6.0, tf=9.0, r=15.0)
    HEA180: ISection = ISection(name="HEA180", h=171.0, b=180.0, tw=6.0, tf=9.5, r=15.0)
    HEA200: ISection = ISection(name="HEA200", h=190.0, b=200.0, tw=6.5, tf=10.0, r=18.0)
    HEA220: ISection = ISection(name="HEA220", h=210.0, b=220.0, tw=7.0, tf=11.0, r=18.0)
    HEA240: ISection = ISection(name="HEA240", h=230.0, b=240.0, tw=7.5, tf=12.0, r=21.0)
    HEA260: ISection = ISection(name="HEA260", h=250.0, b=260.0, tw=7.5, tf=12.5, r=24.0)
    HEA280: ISection = ISection(name="HEA280", h=270.0, b=280.0, tw=8.0, tf=13.0, r=24.0)
    HEA300: ISection = ISection(name="HEA300", h=290.0, b=300.0, tw=8.5, tf=14.0, r=27.0)
    HEA320: ISection = ISection(name="HEA320", h=310.0, b=300.0, tw=9.0, tf=15.5, r=27.0)
    HEA340: ISection = ISection(name="HEA340", h=330.0, b=300.0, tw=9.5, tf=16.5, r=27.0)
    HEA360: ISection = ISection(name="HEA360", h=350.0, b=300.0, tw=10.0, tf=17.5, r=27.0)
    HEA400: ISection = ISection(name="HEA400", h=390.0, b=300.0, tw=11.0, tf=19.0, r=27.0)
    HEA450: ISection = ISection(name="HEA450", h=440.0, b=300.0, tw=11.5, tf=21.0, r=27.0)
    HEA500: ISection = ISection(name="HEA500", h=490.0, b=300.0, tw=12.0, tf=23.0, r=27.0)
    HEA550: ISection = ISection(name="HEA550", h=540.0, b=300.0, tw=12.5, tf=24.0, r=27.0)
    HEA600: ISection = ISection(name="HEA600", h=590.0, b=300.0, tw=13.0, tf=25.0, r=27.0)
    HEA650: ISection = ISection(name="HEA650", h=640.0, b=300.0, tw=13.5, tf=26.0, r=27.0)
    HEA700: ISection = ISection(name="HEA700", h=690.0, b=300.0, tw=14.5, tf=27.0, r=27.0)
    HEA800: ISection = ISection(name="HEA800", h=790.0, b=300.0, tw=15.0, tf=28.0, r=30.0)
    HEA900: ISection = ISection(name="HEA900", h=890.0, b=300.0, tw=16.0, tf=30.0, r=30.0)
    HEA1000: ISection = ISection(name="HEA1000", h=990.0, b=300.0, tw=16.5, tf=31.0, r=30.0)
    # ----------- HEB -------------
    HEB100: ISection = ISection(name="HEB100", h=100.0, b=100.0, tw=6.0, tf=10.0, r=12.0)
    HEB120: ISection = ISection(name="HEB120", h=120.0, b=120.0, tw=6.5, tf=11.0, r=12.0)
    HEB140: ISection = ISection(name="HEB140", h=140.0, b=140.0, tw=7.0, tf=12.0, r=12.0)
    HEB160: ISection = ISection(name="HEB160", h=160.0, b=160.0, tw=8.0, tf=13.0, r=15.0)
    HEB180: ISection = ISection(name="HEB180", h=180.0, b=180.0, tw=8.5, tf=14.0, r=15.0)
    HEB200: ISection = ISection(name="HEB200", h=200.0, b=200.0, tw=9.0, tf=15.0, r=18.0)
    HEB220: ISection = ISection(name="HEB220", h=220.0, b=220.0, tw=9.5, tf=16.0, r=18.0)
    HEB240: ISection = ISection(name="HEB240", h=240.0, b=240.0, tw=10.0, tf=17.0, r=21.0)
    HEB260: ISection = ISection(name="HEB260", h=260.0, b=260.0, tw=10.0, tf=17.5, r=24.0)
    HEB280: ISection = ISection(name="HEB280", h=280.0, b=280.0, tw=10.5, tf=18.0, r=24.0)
    HEB300: ISection = ISection(name="HEB300", h=300.0, b=300.0, tw=11.0, tf=19.0, r=27.0)
    HEB320: ISection = ISection(name="HEB320", h=320.0, b=300.0, tw=11.5, tf=20.5, r=27.0)
    HEB340: ISection = ISection(name="HEB340", h=340.0, b=300.0, tw=12.0, tf=21.5, r=27.0)
    HEB360: ISection = ISection(name="HEB360", h=360.0, b=300.0, tw=12.5, tf=22.5, r=27.0)
    HEB400: ISection = ISection(name="HEB400", h=400.0, b=300.0, tw=13.5, tf=24.0, r=27.0)
    HEB450: ISection = ISection(name="HEB450", h=450.0, b=300.0, tw=14.0, tf=26.0, r=27.0)
    HEB500: ISection = ISection(name="HEB500", h=500.0, b=300.0, tw=14.5, tf=28.0, r=27.0)
    HEB550: ISection = ISection(name="HEB550", h=550.0, b=300.0, tw=15.0, tf=29.0, r=27.0)
    HEB600: ISection = ISection(name="HEB600", h=600.0, b=300.0, tw=15.5, tf=30.0, r=27.0)
    HEB650: ISection = ISection(name="HEB650", h=650.0, b=300.0, tw=16.0, tf=31.0, r=27.0)
    HEB700: ISection = ISection(name="HEB700", h=700.0, b=300.0, tw=17.0, tf=32.0, r=27.0)
    HEB800: ISection = ISection(name="HEB800", h=800.0, b=300.0, tw=17.5, tf=33.0, r=30.0)
    HEB900: ISection = ISection(name="HEB900", h=900.0, b=300.0, tw=18.5, tf=35.0, r=30.0)
    HEB1000: ISection = ISection(name="HEB1000", h=1000.0, b=300.0, tw=19.0, tf=36.0, r=30.0)
    # ----------- HEM -------------
    HEM100: ISection = ISection(name="HEM100", h=120.0, b=106.0, tw=12.0, tf=20.0, r=12.0)
    HEM120: ISection = ISection(name="HEM120", h=140.0, b=126.0, tw=12.5, tf=21.0, r=12.0)
    HEM140: ISection = ISection(name="HEM140", h=160.0, b=146.0, tw=13.0, tf=22.0, r=12.0)
    HEM160: ISection = ISection(name="HEM160", h=180.0, b=166.0, tw=14.0, tf=23.0, r=15.0)
    HEM180: ISection = ISection(name="HEM180", h=200.0, b=186.0, tw=14.5, tf=24.0, r=15.0)
    HEM200: ISection = ISection(name="HEM200", h=220.0, b=206.0, tw=15.0, tf=25.0, r=18.0)
    HEM220: ISection = ISection(name="HEM220", h=240.0, b=226.0, tw=15.5, tf=26.0, r=18.0)
    HEM240: ISection = ISection(name="HEM240", h=270.0, b=248.0, tw=18.0, tf=32.0, r=21.0)
    HEM260: ISection = ISection(name="HEM260", h=290.0, b=268.0, tw=18.0, tf=32.5, r=24.0)
    HEM280: ISection = ISection(name="HEM280", h=310.0, b=288.0, tw=18.5, tf=33.0, r=24.0)
    HEM300: ISection = ISection(name="HEM300", h=340.0, b=310.0, tw=21.0, tf=39.0, r=27.0)
    HEM320: ISection = ISection(name="HEM320", h=359.0, b=309.0, tw=21.0, tf=40.0, r=27.0)
    HEM340: ISection = ISection(name="HEM340", h=377.0, b=309.0, tw=21.0, tf=40.0, r=27.0)
    HEM360: ISection = ISection(name="HEM360", h=395.0, b=308.0, tw=21.0, tf=40.0, r=27.0)
    HEM400: ISection = ISection(name="HEM400", h=432.0, b=307.0, tw=21.0, tf=40.0, r=27.0)
    HEM450: ISection = ISection(name="HEM450", h=478.0, b=307.0, tw=21.0, tf=40.0, r=27.0)
    HEM500: ISection = ISection(name="HEM500", h=524.0, b=306.0, tw=21.0, tf=40.0, r=27.0)
    HEM550: ISection = ISection(name="HEM550", h=572.0, b=306.0, tw=21.0, tf=40.0, r=27.0)
    HEM600: ISection = ISection(name="HEM600", h=620.0, b=305.0, tw=21.0, tf=40.0, r=27.0)
    HEM650: ISection = ISection(name="HEM650", h=668.0, b=305.0, tw=21.0, tf=40.0, r=27.0)
    HEM700: ISection = ISection(name="HEM700", h=716.0, b=304.0, tw=21.0, tf=40.0, r=27.0)
    HEM800: ISection = ISection(name="HEM800", h=814.0, b=303.0, tw=21.0, tf=40.0, r=30.0)
    HEM900: ISection = ISection(name="HEM900", h=910.0, b=302.0, tw=21.0, tf=40.0, r=30.0)
    HEM1000: ISection = ISection(name="HEM1000", h=1008.0, b=302.0, tw=21.0, tf=40.0, r=30.0)


if __name__ == '__main__':
    ipe = ISections.IPE270

    print(ipe)

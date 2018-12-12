# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:43:48 2018

@author: huuskoj
"""

import unittest
from src.sections.steel.ISection import IPE, HEA


# Arcelormittal values [mm]
# A, Iy, Wely, Wply, iy, Avz, Iz, Welz, Wplz, iz, ss, It, Iw
IPE200 = [28.5e2, 1591.0e4, 162.0e3, 182.0e3, 8.23, 11.5e3,
          117.0e4, 23.4e3, 36.5e3, 2.23, 32.6, 4.11e4, 10.5e6]

class TestISection(unittest.TestCase):

    def test_area(self):
        ipe = IPE(200)
        print(ipe)
        self.assertAlmostEqual(ipe.A, IPE200[0], delta=IPE200[0]*0.01)

    def test_Iy(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.I[0], IPE200[1], delta=IPE200[1]*0.01)

    def test_Wely(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.Wel[0], IPE200[2], delta=IPE200[2]*0.01)

    def test_Wply(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.Wpl[0], IPE200[3], delta=IPE200[3]*0.01)

    """
    def test_iy(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe., IPE200[4])
    """

    def test_Avz(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.Ashear, IPE200[5], delta=IPE200[5]*0.01)

    def test_Iz(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.I[1], IPE200[6], delta=IPE200[6]*0.01)

    def test_Welz(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.Wel[1], IPE200[7], delta=IPE200[7]*0.01)

    def test_Wplz(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.Wpl[1], IPE200[8], delta=IPE200[8]*0.01)

    def test_It(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.It, IPE200[-2], delta=IPE200[-2]*0.01)

    def test_Iw(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe.Iw, IPE200[-1], delta=IPE200[-1]*0.01)








if __name__ == '__main__':
    unittest.main()

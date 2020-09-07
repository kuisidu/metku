# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:43:48 2018

@author: huuskoj
"""

import unittest
from ISection import IPE, HEA


# Values from Arcelormittal catalogue, units [mm]
# A, Iy, Wely, Wply, iy, Avz, Iz, Welz, Wplz, iz, ss, It, Iw
IPE200 = [28.5e2, 1943.0e4, 194.0e3, 221.0e3,  8.26, 14.0e2,
          142.0e4, 28.5e3, 44.6e3, 2.24, 36.7, 6.98e4, 13.0e9]
        
ipe = IPE(200)

class TestISection(unittest.TestCase):

    def test_area(self):
        self.assertAlmostEqual(ipe.A, IPE200[0], delta=IPE200[0]*0.01)

    def test_Iy(self):
        self.assertAlmostEqual(ipe.I[0], IPE200[1], delta=IPE200[1]*0.01)

    def test_Wely(self):
        self.assertAlmostEqual(ipe.Wel[0], IPE200[2], delta=IPE200[2]*0.01)

    def test_Wply(self):
        self.assertAlmostEqual(ipe.Wpl[0], IPE200[3], delta=IPE200[3]*0.01)

    """
    def test_iy(self):
        ipe = IPE(200)
        self.assertAlmostEqual(ipe., IPE200[4])
    """

    def test_Avz(self):
        self.assertAlmostEqual(ipe.Ashear, IPE200[5], delta=IPE200[5]*0.01)

    def test_Iz(self):
        self.assertAlmostEqual(ipe.I[1], IPE200[6], delta=IPE200[6]*0.01)

    def test_Welz(self):
        self.assertAlmostEqual(ipe.Wel[1], IPE200[7], delta=IPE200[7]*0.01)

    def test_Wplz(self):
        self.assertAlmostEqual(ipe.Wpl[1], IPE200[8], delta=IPE200[8]*0.01)

    def test_It(self):
        self.assertAlmostEqual(ipe.It, IPE200[-2], delta=IPE200[-2]*0.01)

    def test_Iw(self):
        self.assertAlmostEqual(ipe.Iw, IPE200[-1], delta=IPE200[-1]*0.01)


if __name__ == '__main__':
    unittest.main()

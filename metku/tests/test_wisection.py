#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

import unittest

from metku.sections.steel.wi import WISection


class MyTestCase(unittest.TestCase):
    WI = WISection(h=1000,
                   b=300,
                   tf=16,
                   tw=6,
                   weld_throat=4,
                   fy=355)

    def test_flange_class(self):
        self.assertEqual(self.WI.flange_class(), 3)

    def test_web_class(self):
        self.assertEqual(self.WI.web_class_bend(), 4)

    def test_hw(self):
        self.assertEqual(self.WI.hw, 968)

    def test_c(self):
        self.assertAlmostEqual(self.WI.cf_top, 141.3, delta=0.1)

    def test_bw(self):
        self.assertAlmostEqual(self.WI.cw, 956.7, delta=0.1)

    def test_effective_web(self):
        self.assertAlmostEqual(self.WI.effective_web(psi=-1), 312.4, delta=0.1)

    def test_A(self):
        self.assertAlmostEqual(self.WI.A, 15408, delta=0.1)

    def test_Aeff(self):
        self.assertAlmostEqual(self.WI.effective_area(psi=-1), 14412.6, delta=1)

    def test_web_neg(self):
        self.WI.effective_web(psi=-1)
        self.assertAlmostEqual(self.WI.hneg, 166, delta=0.2)

    def test_Ieffy(self):
        self.assertAlmostEqual(self.WI.effective_second_moment_area_y(psi_init=-1) / 1e6, 2692, delta=1)

    def test_be1(self):
        self.WI.effective_web(psi=-1)
        self.assertAlmostEqual(self.WI.eff_web['be'][0], 125, delta=0.1)

    def test_be2(self):
        self.WI.effective_web(psi=-1)
        self.assertAlmostEqual(self.WI.eff_web['be'][1], 187.4, delta=0.1)

    def test_Ashear(self):
        self.assertAlmostEqual(self.WI.Ashear[0], 6969.6, delta=0.1)

    def test_shear_buckling(self):
        self.assertAlmostEqual(self.WI.shear_buckling_resistance(stiffener_spacing=1000, end_post="rigid") / 1e3, 763.5,
                               delta=1)

    def test_Weffy(self):
        self.assertAlmostEqual(self.WI.W_eff_y / 1e3, 5173, delta=1)

    def test_Iz(self):
        self.assertAlmostEqual(self.WI.Iz / 1e4, 7202, delta=1)

    def test_Iw(self):
        self.assertAlmostEqual(self.WI.Iw / 1e9, 17433, delta=1)

    def test_It(self):
        self.assertAlmostEqual(self.WI.It / 1e3, 888.9, delta=1)

if __name__ == '__main__':
    unittest.main()

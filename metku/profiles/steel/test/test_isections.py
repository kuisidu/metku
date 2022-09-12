#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

import unittest

from metku.profiles.steel import ISections

profile = ISections.IPE80


class MyTestCase(unittest.TestCase):
    """
    Test values from
    https://eurocodeapplied.com/design/en1993/ipe-hea-heb-hem-design-properties
    """

    def test_h(self):
        self.assertEqual(profile.h, 80)

    def test_b(self):
        self.assertEqual(profile.b, 46)

    def test_tw(self):
        self.assertEqual(profile.tw, 3.8)

    def test_tf(self):
        self.assertEqual(profile.tf, 5.2)

    def test_r(self):
        self.assertEqual(profile.r, 5)

    def test_perimeter(self):
        self.assertAlmostEqual(profile.P / 1e3, 0.328, places=3)

    def test_area(self):
        self.assertAlmostEqual(profile.A, 764, places=0)

    def test_shear_area_z(self):
        self.assertAlmostEqual(profile.Avz, 358, places=0)

    def test_shear_area_y(self):
        self.assertAlmostEqual(profile.Avy, 478, places=0)

    def test_second_moment_y(self):
        self.assertAlmostEqual(profile.Iy / 1e6, 0.8014, places=4)

    def test_radius_of_gyration_y(self):
        self.assertAlmostEqual(profile.iy, 32.4, places=1)

    def test_elastic_modulus_y(self):
        self.assertAlmostEqual(profile.Wely / 1e3, 20.03, places=2)

    def test_plastic_modulus_y(self):
        self.assertAlmostEqual(profile.Wply / 1e3, 23.22, places=2)

    def test_second_moment_z(self):
        self.assertAlmostEqual(profile.Iz / 1e6, 0.08489, places=5)

    def test_radius_of_gyration_z(self):
        self.assertAlmostEqual(profile.iz, 10.5, places=1)

    def test_elastic_modulus_z(self):
        self.assertAlmostEqual(profile.Welz / 1e3, 3.691, places=3)

    def test_plastic_modulus_z(self):
        self.assertAlmostEqual(profile.Wplz / 1e3, 5.818, places=3)

    def test_torsion_constant(self):
        self.assertAlmostEqual(profile.It / 1e3, 6.727, places=2)

    def test_torsion_modulus(self):
        # self.assertAlmostEqual(profile.Wt, 1.770e3)
        raise NotImplementedError

    def test_warping_constant(self):
        self.assertAlmostEqual(profile.Iw / 1e6, 115.1, places=1)

    def test_warping_modulus(self):
        # self.assertAlmostEqual(profile.Ww, 1.770e3)
        raise NotImplementedError


if __name__ == '__main__':
    unittest.main()

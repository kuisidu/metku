#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

import unittest
from metku.design import SteelSection
from metku.materials.steel import Steel
from metku.profiles.steel import ISections

profile = ISections.IPE80
material = Steel.S355
cross_section = SteelSection(material=material, profile=profile)

class MyTestCase(unittest.TestCase):


    def test_NcRd(self):
        """
        https://eurocodeapplied.com/design/en1993/ipe-hea-heb-hem-design-properties
        :return:
        """
        self.assertAlmostEqual(cross_section.NcRd, 271.34e3, places=-1)  # add assertion here

    def test_NtRd(self):
        """
        https://eurocodeapplied.com/design/en1993/ipe-hea-heb-hem-design-properties
        :return:
        """
        self.assertAlmostEqual(cross_section.NtRd, 271.34e3, places=-1)  # add assertion here

    def test_VzRd(self):
        """
        https://eurocodeapplied.com/design/en1993/ipe-hea-heb-hem-design-properties
        :return:
        """
        self.assertAlmostEqual(cross_section.VzRd, 73.31e3, places=-1)  # add assertion here


    def test_VyRd(self):
        """
        https://eurocodeapplied.com/design/en1993/ipe-hea-heb-hem-design-properties
        :return:
        """
        self.assertAlmostEqual(cross_section.VyRd, 98.05e3, places=-1)  # add assertion here

    def test_MelRdy(self):
        self.assertAlmostEqual(cross_section.MelyRd*1e-6, 7.11, places=2)

    def test_MplRdy(self):
        self.assertAlmostEqual(cross_section.MplyRd*1e-6, 8.24, places=2)

    def test_MelRdz(self):
        self.assertAlmostEqual(cross_section.MelzRd*1e-6, 1.31, places=2)

    def test_MplRdz(self):
        self.assertAlmostEqual(cross_section.MplzRd*1e-6, 2.07, places=2)

if __name__ == '__main__':
    unittest.main()

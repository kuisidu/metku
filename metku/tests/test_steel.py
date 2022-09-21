#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

import unittest
from metku.materials.steel import Steel

class TestSteel(unittest.TestCase):

    def test_fy(self):
        s235 = Steel.S235
        self.assertEqual(s235.fy, 235)

    def test_E(self):
        s235 = Steel.S235
        self.assertEqual(s235.E, 210e3)


if __name__ == "__main__":
    unittest.main()

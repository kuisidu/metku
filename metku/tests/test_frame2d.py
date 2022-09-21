#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

import unittest
import metku.frame2d as f2d


"""
TESTS TODO:
- calculation result validation

"""

class MyTestCase(unittest.TestCase):

    def create_frame(self) -> f2d.Frame2D:
        frame = f2d.Frame2D()
        return frame


    def test_something(self):
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()

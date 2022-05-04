import unittest
from metku.materials.steel import Steel
from pint import UnitRegistry

ureg = UnitRegistry()


class TestSteel(unittest.TestCase):

    def test_fy(self):
        s235 = Steel.S235
        self.assertEqual(s235.fy, 235 * ureg.MPa)

    def test_E(self):
        s235 = Steel.S235
        self.assertEqual(s235.E, 210e3 * ureg.MPa)


if __name__ == "__main__":
    unittest.main()


import unittest
from frame2d import Frame2D, SteelBeam, SteelColumn, PointLoad, LineLoad, FixedSupport

# Coordinates of members, loads and supports
coord1 = [[0, 0], [0, 3.65]]
coord2 = [[0, 3.65], [0, 7.3]]
coord3 = [[7.3, 0], [7.3, 3.65]]
coord4 = [[7.3, 3.65], [7.3, 7.3]]
coord5 = [[0, 3.65], [7.3, 3.65]]
coord6 = [[0, 7.3], [7.3, 7.3]]
supp_coord1 = [0.0, 0]
supp_coord2 = [7.3, 0]
# Members
col1 = SteelColumn(coord1, profile='ipe 450')
col2 = SteelColumn(coord2, profile='ipe 400')
col3 = SteelColumn(coord3, profile='ipe 450')
col4 = SteelColumn(coord4, profile='ipe 400')
beam1 = SteelBeam(coord5, profile='ipe 360')
beam2 = SteelBeam(coord6, profile='ipe 300')
# Loads
load1 = PointLoad(coord1[1], [36, 0, 0])
load2 = PointLoad(coord2[1], [18, 0, 0])
load3 = LineLoad(beam1, [-62.5, -62.5], 'y')
load4 = LineLoad(beam2, [-38, -38], 'y')
# Create empty frame 'envelope'
frame = Frame2D(num_elements=10)
# Add members
frame.add(col1)
frame.add(col2)
frame.add(col3)
frame.add(col4)
frame.add(beam1)
frame.add(beam2)
# Add loads
frame.add(load1)
frame.add(load2)
frame.add(load3)
frame.add(load4)
# Add supports
frame.add(FixedSupport(supp_coord1))
frame.add(FixedSupport(supp_coord2))
# Generate frame
frame.generate()
frame.calculate()

# Values from Robot
ALPHA_CR = [4.87840e1, 1.13178e2, 2.53259e2, 4.21790e2]

# Difference allowed between Robot value
DELTA = 0.5

class TestFrame2D(unittest.TestCase):
    
    
    def test_Med(self):
         self.assertAlmostEqual(col1.med, -101.14, delta=DELTA)
         self.assertAlmostEqual(col2.med, -139.82, delta=DELTA)
         self.assertAlmostEqual(col3.med, -175.56, delta=DELTA)
         self.assertAlmostEqual(col4.med, -300.62, delta=DELTA)
         self.assertAlmostEqual(beam1.med, -300.62, delta=DELTA)
         self.assertAlmostEqual(beam2.med, 175.56, delta=DELTA)
         
    def test_Ved(self):
         self.assertAlmostEqual(col1.ved, -24.11, delta=DELTA)
         self.assertAlmostEqual(col2.ved, -72.19, delta=DELTA)
         self.assertAlmostEqual(col3.ved, -143.60, delta=DELTA)
         self.assertAlmostEqual(col4.ved, -238.51, delta=DELTA)
         self.assertAlmostEqual(beam1.ved, 78.11, delta=DELTA)
         self.assertAlmostEqual(beam2.ved, 90.19, delta=DELTA)
    
    def test_Ned(self):
         self.assertAlmostEqual(col1.ned, 351.55, delta=DELTA)
         self.assertAlmostEqual(col2.ned, 133.80, delta=DELTA)
         self.assertAlmostEqual(col3.ned, 90.19, delta=DELTA)
         self.assertAlmostEqual(col4.ned, -12.08, delta=DELTA)
         self.assertAlmostEqual(beam1.ned, 382.10, delta=DELTA)
         self.assertAlmostEqual(beam2.ned, 143.60, delta=DELTA)
    
    def test_false_profile(self):
        with self.assertRaises(TypeError):
            beam1.profile = 'iipe 100'
    
    #def test_alphacr(self):
    #    self.assertAlmostEqual(frame.alpha_cr[0]*1e3, ALPHA_CR[0], delta=DELTA)
    
    
if __name__ == '__main__':
    unittest.main()
    frame.bmd(5)
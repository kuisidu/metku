import sys
sys.path.append(".\src")

from src.frame2d.frame2d import Frame2D, SteelBeam, SteelColumn, PointLoad,\
     FixedSupport, XHingedSupport, LineLoad, YHingedSupport, Hinge

import random
import time
import matplotlib.pyplot as plt

# Coordinates of members, loads and supports
coord1 = [[0,1], [1,1]]
coord2 = [[1,0], [1,1]]
coord3 = [[0.0,0], [0.0, 2]]
supp_coord1 = [0.0,0]
supp_coord2 = [1,0]
#hinge_coord = [0.5, 1]

# Loads and members
col1 = SteelColumn(coord2)
col2 = SteelColumn(coord3)
beam1 = SteelBeam(coord1)

load1 = PointLoad(coord3[1], [50, 0,0])
load2 = LineLoad(beam1, [-50, -50], 'y')
# Create empty frame 'envelope'
frame = Frame2D(num_elements=5)

# Add members
frame.add(col1)
frame.add(col2)
frame.add(beam1)

# Add loads
frame.add(load1)
frame.add(load2)

# Add supports
frame.add(FixedSupport(supp_coord1, supp_id=1))
frame.add(FixedSupport(supp_coord2, supp_id=2))

# Add hinge
#frame.add(Hinge(hinge_coord))

# Generate frame
frame.generate()

#col1.profile = 'HE 200 A'

# Calculate result
frame.calculate()
#print(frame.nodal_forces[2])





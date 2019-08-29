from src.truss2d import *
from src.frame2d.frame2d import *

import numpy as np

# Simple frame parameters
num_storeys = 1
num_bays = 1
storey_H = 5000
bay_L = 12000

# Simple truss parameters
H0 = storey_H
H1 = 1000
H2 = 1500
H3 = 1000
L1 = bay_L / 2
L2 = bay_L / 2
n = 16

dx = 500

# Create simple frame without beams
frame = Frame2D(simple=[num_storeys, num_bays, storey_H, bay_L],
                supports='fixed',
                beams=False)


# Create truss
truss = Truss2D(simple={'H0': H0,
                        'H1': H1,
                        'H2': H2,
                        'L1': L1,
                        'L2': L2,
                        'n': n,
                        'dx': dx})

frame.add(truss)

for joint in truss.joints.values():
    joint.loc += np.random.randint(-1, 2) * 1e-1
    joint.loc = max(0, joint.loc)
    joint.loc = min(1, joint.loc)


frame.plot(print_text=False)
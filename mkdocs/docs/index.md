# Frame2D

Python library for structural analysis and optimization.

## Getting started

### Installation

Install frame2d using PyPI
```
pip install frame2d
```

### Creating and analysing simple frame

```python
from frame2d import *

# Number of storeys
storeys = 5
# Number of bays
bays = 2
 # Storey height in mm
storey_H = 3400
# Bay length in mm
bay_L = 5000
# Create simple frame using simple keyword
frame = Frame2D(simple=[storeys, bays, storey_H, bay_L], supports='fixed')
# Get all the beams from the frame
beams = [beam for beam in frame.members.values() if beam.mtype =="beam"]
# Iterate through all beams and add line load to them
for beam in beams:
  # Create line load
  # LineLoad(member, force vector, direction)
  ll = LineLoad(beam, [-10, -10] 'y')
  # Add load to frame
  frame.add(ll)
# Generates FEM -model of the frame
frame.generate()
# Runs FEA
frame.calculate()
# Plot bending moment diagram
scaling_factor = 10
frame.bmd(scaling_factor)
```

### Exporting frame to Robot Structural Analysis
Frames can be exported to Robot Structural Analysis by creating a .str -file that contains frame's data.

```python
filename = "example_frame"
# Creates .str -file to current working directory
frame.to_robot(filename, num_frames=3, s=5000)
```

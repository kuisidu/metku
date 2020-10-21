"""
This is an example for creating and analyzing
a simply supported beam
"""

from metku.frame2d import *

# Create empty frame
frame = Frame2D()
# Create beam
coordinates = [[0, 0], [5000, 0]]
beam = SteelBeam(coordinates,
                 profile="ipe 100",
                 material="S355",
                 num_elements=6)
# Add beam to frame
frame.add(beam)
# Create lineload
lineload = LineLoad(beam, [-10, -10], 'y')
# Add load to frame
frame.add(lineload)
# Create supports
sup1 = XYHingedSupport([0, 0])
sup2 = XYHingedSupport([5000, 0])
# Add supports to frame
frame.add(sup1)
frame.add(sup2)
# Generate FE-model of the frame
frame.generate()
# Run linear analysis for frame
frame.calculate()

## Plotting

# Plot bending moment diagram
frame.bmd(scale=50)
# Plot deflected beam
frame.plot_deflection(scale=1)

## Getting force and displacement values

# Print maximum forces
print(f'Maximum normal force: {beam.ned * 1e-3:.2f} kN')
print(f'Maximum shear force: {beam.ved * 1e-3:.2f} kN')
print(f'Maximum bending moment: {beam.med * 1e-6:.2f} kNm')
# Print forces at each node
for node_id, forces in beam.nodal_forces.items():
    N, V, M = forces
    print(f'Forces at node {node_id}: \n'
          f'Normal: {N * 1e-3:.2f} kN '
          f'Shear: {V * 1e-3:.2f} kN '
          f'Moment: {M * 1e-6:.2f} kNm ')

# Print displacements at each node
for node_id, forces in beam.nodal_displacements.items():
    dx, dy, r = forces[0]
    print(f'Displacements at node {node_id}: \n'
          f'X: {dx:.2f} mm '
          f'Y: {dy:.2f} mm '
          f'Rotation: {r:.2f} radians ')


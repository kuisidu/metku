from structures.timber.timber_member import TimberMember
from structures.timber.timber_frame2d import TimberFrame2D
from materials.timber_data import T
from sections.timber.timber_section import TimberSection, GypsumPlasterboardF, WoodenFireProtection, GypsumPlasterboardA
import frame2d.frame2d as f2d
import numpy as np
from math import pi

# fr = TimberFrame2D()
#
# coord1 = [[0, 0], [0, 5000]]
# coord3 = [[0, 5000], [20000, 5000]]
# coord2 = [[20000, 5000], [20000, 0]]
#
# fr.add(f2d.TimberFrameMember(coord1, TimberSection(240, 450), material=T.GL32c, mtype='column', support_nodes_z=0))
# fr.add(f2d.TimberFrameMember(coord2, TimberSection(240, 720), material=T.GL32c, mtype='beam', ldc='mt', num_elements=5,
#                              support_nodes_z=0, edge_load='compression'))
# fr.add(f2d.TimberFrameMember(coord3, TimberSection(240, 450), material=T.GL32c, mtype='column', support_nodes_z=1,
#                              num_elements=5, support_nodes_y=1))
#
# R = 0
# fr.members[0].R = R
# fr.members[1].R = R
# fr.members[2].R = R
#
# fr.add(f2d.LineLoad(fr.members[1], [-15, -15], 'y'))
# fr.add(f2d.LineLoad(fr.members[0], [3, 3], 'x'))
# fr.add(f2d.XYHingedSupport([0, 0]))
# fr.add(f2d.XYHingedSupport([20000, 0]))
# fr.generate()
# fr.calculate()
# fr.bmd(1)
#
# print(fr.f)

print(np.tan(pi/2))

print(np.power(np.array([2, 4]), 2))
print(np.array([2, 4]) ** -0.2)

a = np.array([True, False, True])
b = np.array([0.2, 0.5, 2])
c = b == 0.5
print(c)
print(a * b)

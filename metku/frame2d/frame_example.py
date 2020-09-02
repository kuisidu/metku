# Import dependencies
from fem.frame_fem  import FrameFEM, BeamSection, LineLoad, PointLoad
from fem.elements import EBBeam, EBSemiRigidBeam

fem = FrameFEM()

c0 = [0,0]
c1 = [0,3]
c2 = [0,6]
c3 = [2,6]
c4 = [4,6]
c5 = [8,6]

fem.add_node(*c0)
fem.add_node(*c1)
fem.add_node(*c2)
fem.add_node(*c3)
fem.add_node(*c4)
fem.add_node(*c5)

E = 210e3
nu = 0.3
rho = 7850e-9
fem.add_material(E, nu, rho)

A = 1e3
Iy = 1e6
sect = BeamSection(A, Iy)
fem.add_section(sect)


for i in range(len(fem.nodes)-1):
    n1 = fem.nodes[i]
    n2 = fem.nodes[i+1]
    mat = fem.materials[0]
    ele = EBBeam(n1, n2, sect, mat)
    fem.add_element(ele)

fem.add_support(1, 0, [0, 1, 2], 0)
fem.add_support(1, 4, [1], 0)
fem.add_support(1, 5, [1], 0)

n1 = fem.nodes[fem.nodal_coords.index(c1)]
n2 = fem.nodes[fem.nodal_coords.index(c3)]

pl1 = PointLoad(2, n1, [10, 0, 0], 1)
pl2 = PointLoad(2, n2, [0, -20, 0], 1)

ll1 = LineLoad(2, fem.elements[-1], [0,1], [-5, -5], 1)

fem.add_load(pl1)
fem.add_load(pl2)
fem.add_load(ll1)

fem.add_loadcase(1,2)

fem.nodal_dofs()
fem.linear_statics()

for elem in fem.elements:
    print(elem.bending_moment[0])
print(elem.bending_moment[1])
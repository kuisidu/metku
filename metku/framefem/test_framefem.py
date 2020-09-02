import unittest

from .framefem import FrameFEM, BeamSection, LineLoad, PointLoad
from .elements import EBBeam


class TestFrameFEM(unittest.TestCase):
    def test_adding_material(self):
        fem = FrameFEM()
        fem.add_material(100, 200, 300)
        mat = fem.materials[0]
        self.assertEqual(100, mat.young)
        self.assertEqual(200, mat.nu)
        self.assertEqual(300, mat.density)

    def test_adding_node(self):
        fem = FrameFEM()
        fem.add_node(0, 0)
        node = fem.nodes[0]
        self.assertEqual([0, 0], [node.x, node.y])

    def test_linear_statics(self):
        fem = FrameFEM()
        # Nodes
        fem.add_node(0, 0)
        fem.add_node(0, 1)
        fem.add_node(1, 1)
        fem.add_node(1, 0)
        # Material
        fem.add_material(210e3, 0.3, 7850e-9)
        # Supports
        fem.add_support(1, 0, [0, 1, 2])
        fem.add_support(1, 3, [0, 1, 2])
        # Sections
        sect = BeamSection(1e3, 2e6)
        fem.add_section(sect)
        # Elements
        for nid in range(fem.nnodes() - 1):
            n1 = fem.nodes[nid]
            n2 = fem.nodes[nid + 1]
            mat = fem.materials[0]
            sect = fem.sections[0]
            ele = EBBeam(n1, n2, sect, mat)
            fem.add_element(ele)
        # Loads
        pointload = PointLoad(1, fem.nodes[1], [10, 0, 0], f=1)
        lineload = LineLoad(1, fem.elements[1], [0, 1], [-10, -10], 1)
        fem.add_load(pointload)
        fem.add_load(lineload)
        # Loadcase
        fem.add_loadcase(1, 1)
        # Linear statics
        fem.nodal_dofs()
        fem.linear_statics()


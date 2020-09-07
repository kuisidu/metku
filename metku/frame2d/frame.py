from metku.frame2d.member import *
from metku.frame2d.loads import *
from metku.frame2d.supports import *

import matplotlib.pyplot as plt
import metku.framefem.framefem as ff


class Frame:
    """
    CLASS DEFINITION
    """

    def __init__(self):

        # Initialize empty dicts
        self.members = {}
        self.loads = {}
        self.supports = {}
        self.loadcases = {}


    def generate_fem_model(self, load_id, supp_id):
        """
        Generates new fem model
        :param load_id: loads to include
        :param supp_id: supports to include
        :return: generated frame
        """
        # Initialize frame
        fem_model = ff.FrameFEM()
        # Members
        for mem in self.members.values():

            # for node in mem.nodes.values():
            #     fem_model.add_node(*node.coord)

            for element in mem.elements.values():
                fem_model.add(element)
        # Supports
        for support in self.supports.values():
            if support.sid == supp_id:
                fem_model.add(support)
        # Loads
        for load in self.loads.values():
            if load.load_id == load_id:
                if isinstance(load, PointLoad):
                    pl = ff.PointLoad(load.load_id,
                                      load.node,
                                      load.F,
                                      load.scale)
                    fem_model.add_load(pl)
                elif isinstance(load, LineLoad):
                    pass

        return fem_model

    def add(self, obj, intersect=True):
        """
        Adds object to frame
        Object can be Member, Load, Support or Truss

        :param obj: obj to be added to frame
        """
        if isinstance(obj, FrameMember):
            self.add_member(obj, intersect)

        elif isinstance(obj, PointLoad):
            self.add_pointload(obj)

    def add_member(self, member, intersect=True):
        """
        Adds member to frame
        Adds nodes between intersecting members
        :param member: FrameMember object
        """
        for mem in self.members.values():
            point = member.intersection(mem)
            if len(point) and intersect:
                if list(point) in mem.nodal_coordinates:

                    for nid, node in mem.nodes.items():
                        if np.all(node.coord == point):
                            member.intersecting_nodes[nid] = node
                elif list(point) in member.nodal_coordinates:
                    for nid, node in member.nodes.items():
                        if np.all(node.coord == point):
                            mem.intersecting_nodes[node.nid] = node
                else:
                    nid = member.mem_id * 100 + len(member.nodes)
                    node = FEMNode(nid, *point)
                    mem.intersecting_nodes[nid] = node
                    member.intersecting_nodes[nid] = node
                mem.nels += 1
                member.nels += 1
                mem.update()
                member.update()

        self.members[member.mem_id] = member

    def add_pointload(self, pl):
        """
        Adds point load to frame
        :param pl: PointLoad object
        """
        for mem in self.members.values():
            point = mem.intersection(pl._coord)
            if len(point):
                if list(point) in mem.nodal_coordinates:
                    node = [node for node in mem.nodes.values() if np.all(node.coord == point)][0]
                    pl.node = node
                else:
                    nid = mem.mem_id * 100 + len(mem.nodes)
                    node = FEMNode(nid, *point)
                    mem.nodes[nid] = node
                    pl.node = node



    def calculate(self, load_id=0, supp_id=0):
        """
        Runs linear finite element analysis
        :param load_id: load id
        """
        fem_model = self.generate_fem_model(load_id, supp_id)
        fem_model.add_loadcase(supp_id=supp_id, load_id=load_id)
        fem_model.nodal_dofs()
        fem_model.linear_statics()


    def plot(self, show=True):
        """
        Plots the frame
        :param show: boolean
        """

        for mem in self.members.values():
            mem.plot(show=False)
        if show:
            plt.show()


if __name__ == '__main__':
    frame = Frame()
    mem = FrameMember([0, 0], [1000, 1000], num_elements=1)
    mem2 = FrameMember([-100, 750], [1000, 500], num_elements=1)
    mem3 = FrameMember([100, 100], [-100, 1000], num_elements=3)
    frame.add(mem)
    frame.add(mem2)
    frame.add(mem3)
    for mem in frame.members.values():
        mem.plot_elements(False)
    plt.show()
    fem = frame.generate_fem_model(0,0)
    fem.draw()

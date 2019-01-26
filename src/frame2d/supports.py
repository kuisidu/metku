class Support:
    def __init__(self, coordinate, dofs, supp_id=1):

        self.node = None
        self.node_id = None
        self.coordinate = coordinate
        self.supp_id = supp_id
        self.dofs = dofs

    @property
    def coordinate(self):
        if self.node and list(self.node.x) != self.__coordinate:
            self.__coordinate = list(self.node.x)
        return self.__coordinate

    @coordinate.setter
    def coordinate(self, val):
        self.__coordinate = val
        if self.node:
            self.node.x = np.array(val)

    def add_support(self, fem_model):

        if self.coordinate in fem_model.nodal_coords:
            idx = fem_model.nodal_coords.index(self.coordinate)
            self.node_id = idx
            self.node = fem_model.nodes[idx]
            fem_model.add_support(self.supp_id, self.node_id, self.dofs, val=0.0)
        else:
            pass


class FixedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0, 1, 2], supp_id)


class XHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0], supp_id)


class YHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [1], supp_id)


class XYHingedSupport(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [0, 1], supp_id)


class Hinge(Support):
    def __init__(self, coordinate, supp_id=1):
        super().__init__(coordinate, [], supp_id)
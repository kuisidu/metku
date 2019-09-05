from src.framefem.framefem import Material


class Steel(Material):

    def __init__(self, fy, fu, name):
        super().__init__(210e3, 0.3, 7850e-9)
        self.fy = fy
        self.fu = fu
        self.name = name

    def __repr__(self):
        return self.name


MATERIALS = dict(
    S235=Steel(235, 360, "S235"),
    S275=Steel(275, 430, "S275"),
    S355=Steel(355, 510, "S355"),
    S280GD=Steel(280, 360, "S280GD"),
    S420=Steel(420, 500, "S420"),
)







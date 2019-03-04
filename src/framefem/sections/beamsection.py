class BeamSection(Section):
    """ Class for beam element cross sections

        Parameters:
        -----------
        :param A: cross-sectional area [mm^2]
        :param Iy: second moment of area [mm^4]

        :type A: float
        :type Iy: float

        Variables:
        ----------
        :ivar Iy: second moment of area [mm^4]

        :vartype Iy: float
    """

    def __init__(self, A, Iy):
        Section.__init__(self, A)
        self.Iy = Iy
        """ Second moment of area [mm^4]"""

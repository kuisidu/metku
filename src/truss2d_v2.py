from src.frame2d.frame2d import *


class Truss2D(Frame2D):


    def __init__(self):
        super().__init__()

    def add(self, this):

        super().add(this)




if __name__ == '__main__':

    truss = Truss2D()

    tc = FrameMember([[0,0], [1000,100]])
    truss.add(tc)
    truss.plot()
from src.truss2d import *
from src.frame2d.frame2d import *

if __name__ == '__main__':
    L1 = 5000
    H0 = 5000
    H1 = 1000
    H2 = 1500
    n = 16
    frame = Frame2D(simple=[1, 1, 5000, 5000], supports='fixed', beams=False)

    truss = Truss2D(simple={'H0': H0,
                            'H1': H1,
                            'H2': H2,
                            'L1': L1/2,
                            'n': n,
                            'dx': 500})
    # t1 = TopChord([[0, H0 + H1], [L1/2, H0 + H2]])
    # b1 = BottomChord([[0, H0], [L1, H0]])
    # truss.add(t1)
    # truss.add(b1)
    # w0 = TrussWeb(0, 0.5)
    # w1 = TrussWeb(0.5, 0.5)
    # w2 = TrussWeb(0.25, 0.5)
    # truss.add(w0)
    # truss.add(w1)
    # truss.add(w2)
    # frame.add(truss)
    # frame.generate()
    #frame.plot()
    #frame.f.draw()
    for mem in truss.members.values():
        print(mem.coordinates)

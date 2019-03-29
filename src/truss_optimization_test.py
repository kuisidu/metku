
from truss2d import Truss2D
from frame2d.frame2d import Frame2D, LineLoad, PointLoad, XHingedSupport, XYHingedSupport, FixedSupport
from scipy.optimize import fmin_slsqp, differential_evolution
import random
import numpy as np

H0 = 0
H1 = 1
H2 = 1.5
H3 = H1
L1 = 2
L2 = L1
n = 10
flip = True

# Create truss
truss = Truss2D(simple=[H0, H1, H2,H3, L1,L2, n],
                flip=flip,
                num_elements=2)
# Add supports
truss.add(FixedSupport([0,H0]))
truss.add(FixedSupport([L1+L2, H0]))

# Add loads
truss.add(LineLoad(truss.top_chords[0], [-50, -50], 'y'))
truss.add(LineLoad(truss.top_chords[1], [-50, -50], 'y'))


def costf(X):
    R = constf(X)
    R = np.clip(R, -100, 0)
    R = np.sqrt((R)**2)
    return np.sum(R)

def constf(X):
    """
    for i in range(0, len(X*2), 2):
        prof = 'SHS ' + str(X[i//2]) + "X5"
        if i == 0:
            truss.members['T0'].profile = prof
            truss.members['T1'].profile = prof
        elif i == 2:
            truss.members['B2'].profile = prof
        else:
            mem1 = 'W' + str(i-4)
            mem2 = 'W' + str(i-3)
            truss.members[mem1].profile = prof
            truss.members[mem2].profile = prof
            
            
    
    for i, mem in enumerate(truss.members.values()):
        prof = 'SHS ' + str(X[i]) + "X5"
        mem.profile= prof
    """
    X.sort()
    for i, val in enumerate(X):
        val = round(val, 3)
        truss.joints[i*2].loc = val
        truss.joints[(i*2)+1].loc = 1-val
        
    truss.calculate()
    R = []
    for mem in truss.members.values():
        R.extend(list(mem.r))
        
    return 1-np.asarray(R)


def optimize():

    bounds= [(0.,1)]

    bounds=bounds*6

    differential_evolution(costf,
                            bounds=bounds,
                            maxiter=50,
                            popsize=15,
                            disp=True)




truss.generate()

optimize()
#truss.calculate()
#truss.f.draw()
truss.plot_normal_force()
from optimization.solvers.optsolver import OptSolver




class PSO(OptSolver):

    def __init__(self):

        self.best_X = None

class Particle:

    def __init__(self, X, obj, parent):
        self.X = X
        self.obj = obj
        self.best_X = X

    @property
    def velocity(self):
        pass

    @property
    def best_f(self):
        return self.obj(self.best_X)

    @property
    def swarm_best_X(self):
        return self.parent.best_X

    def move(self):
        """
        Moves particle
        """
        self.X += self.velocity


if __name__ == '__main__':

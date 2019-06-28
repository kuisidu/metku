from optimization.structopt import OptimizationProblem
from abc import ABCMeta, abstractclassmethod

class OptSolver(metaclass=ABCMeta):

    def __init__(self):
        pass

    @abstractclassmethod
    def solve(self, problem):
        pass





class VNS(OptSolver):

    def __intit__(self, step=1):
        self.step = step



    def solve(self, problem):
        print(problem([10]*10))

        return 1, 2



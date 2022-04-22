# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
from frame2d.frame2d import *
from optimization.structopt import *
from truss2d import Truss2D
from optimization.solvers import DiscreteVNS


RHS_PROFILES = ['RHS 40X40X2.0', 'RHS 40X40X2.5', 'RHS 50X50X2.0',
                'RHS 40X40X3.0', 'RHS 60X60X2.0', 'RHS 50X50X2.5',
                'RHS 70X70X2.0', 'RHS 50X50X3.0', 'RHS 60X60X2.5',
                'RHS 80X80X2.0', 'RHS 70X70X2.5', 'RHS 60X60X3.0',
                'RHS 50X50X4.0', 'RHS 80X80X2.5', 'RHS 70X70X3.0',
                'RHS 60X60X4.0', 'RHS 90X90X2.5', 'RHS 80X80X3.0',
                'RHS 100X100X2.5', 'RHS 70X70X4.0', 'RHS 90X90X3.0',
                'RHS 60X60X5.0', 'RHS 100X100X3.0', 'RHS 120X120X2.5',
                'RHS 80X80X4.0', 'RHS 70X70X5.0', 'RHS 90X90X4.0',
                'RHS 140X140X2.5', 'RHS 120X120X3.0', 'RHS 80X80X5.0',
                'RHS 150X150X2.5', 'RHS 100X100X4.0', 'RHS 140X140X3.0',
                'RHS 90X90X5.0', 'RHS 150X150X3.0', 'RHS 120X120X4.0',
                'RHS 100X100X5.0', 'RHS 160X160X3.0', 'RHS 140X140X4.0',
                'RHS 100X100X6.0', 'RHS 120X120X5.0', 'RHS 150X150X4.0',
                'RHS 160X160X4.0', 'RHS 140X140X5.0', 'RHS 120X120X6.0',
                'RHS 180X180X4.0', 'RHS 150X150X5.0', 'RHS 120X120X7.1',
                'RHS 160X160X5.0', 'RHS 200X200X4.0', 'RHS 140X140X6.0',
                'RHS 150X150X6.0', 'RHS 120X120X8.0', 'RHS 180X180X5.0',
                'RHS 160X160X6.0', 'RHS 120X120X8.8', 'RHS 180X180X5.6',
                'RHS 200X200X5.0', 'RHS 150X150X7.1', 'RHS 120X120X10.0',
                'RHS 180X180X6.0', 'RHS 160X160X7.1', 'RHS 220X220X5.0',
                'RHS 150X150X8.0', 'RHS 200X200X6.0', 'RHS 160X160X8.0',
                'RHS 150X150X8.8', 'RHS 180X180X7.1', 'RHS 250X250X5.0',
                'RHS 220X220X6.0', 'RHS 160X160X8.8', 'RHS 150X150X10.0',
                'RHS 180X180X8.0', 'RHS 200X200X7.1', 'RHS 160X160X10.0',
                'RHS 180X180X8.8', 'RHS 250X250X6.0', 'RHS 300X300X5.0',
                'RHS 220X220X7.1', 'RHS 200X200X8.0', 'RHS 260X260X6.0',
                'RHS 180X180X10.0', 'RHS 200X200X8.8', 'RHS 220X220X8.0',
                'RHS 250X250X7.1', 'RHS 300X300X6.0', 'RHS 260X260X7.1',
                'RHS 220X220X8.8', 'RHS 200X200X10.0', 'RHS 250X250X8.0',
                'RHS 180X180X12.5', 'RHS 260X260X8.0', 'RHS 220X220X10.0',
                'RHS 300X300X7.1', 'RHS 250X250X8.8', 'RHS 260X260X8.8',
                'RHS 200X200X12.5', 'RHS 300X300X8.0', 'RHS 250X250X10.0',
                'RHS 400X400X6.0', 'RHS 260X260X10.0', 'RHS 300X300X8.8',
                'RHS 400X400X7.1', 'RHS 250X250X12.5', 'RHS 300X300X10.0',
                'RHS 260X260X12.5', 'RHS 400X400X8.0', 'RHS 400X400X8.8',
                'RHS 300X300X12.5', 'RHS 400X400X10.0', 'RHS 400X400X12.5']



class KTruss24m(OptimizationProblem):
    """
    General class for creating optimization problems that utilize the Frame2D -module

    Finite element analysis is done in the objective function and it's should always
    be calculated before constraints.

    """
    def __init__(self ,name="Test Truss", profiles=RHS_PROFILES):
        super().__init__(name)

        self.profiles = profiles
        self._create_structure()
        self._create_variables()
        self.create_objective()
        self.create_constraints()



    def _create_structure(self):

        truss = Truss2D(simple={'H0': 0,
                                'H1': 1200,
                                'H2': 2400,
                                'L1': 12000,
                                'n': 18})


        truss.add(FixedSupport([0, 1200]))
        truss.add(FixedSupport([24000, 1200]))
        truss.add(LineLoad(truss.top_chords[0], [-23.52, -23.52], 'y'))
        truss.add(LineLoad(truss.top_chords[1], [-23.52, -23.52], 'y'))
        truss.generate()
        truss.calculate()
        self.structure = truss
        truss.plot_normal_force()


    def _create_variables(self):

        # Initialize variables as an empty array
        self.vars = []
        # Create profile variables
        # TODO: symmetry, grouping
        for i, mem in enumerate(self.structure.members.values()):
            name = mem.mem_id
            var = DiscreteVariable(name,
                                profiles=self.profiles,
                                target={"property": "PROFILE",
                                        "objects": [mem]})
            self.vars.append(var)

        # Create geometry variables
        # TODO

    def create_objective(self):

        def objective(X):

            weight = 0
            for x, mem in zip(X, self.structure.members.values()):
                weight += mem.weight
            return weight

        self.obj = objective

    def _create_normal_force_constraints(self):
        """
        Creates normal force constraints
        """
        for i, mem in enumerate(self.structure.members.values()):
            def con_fun(X, i=i):
                return mem.ned / mem.NRd - 1
            constr = NonLinearConstraint(con_fun, name='Normal Force ' + str( i +1))
            self.cons.append(constr)


    def _create_shear_force_constraints(self):
        """
        Creates shear force constraints
        """
        for i, mem in enumerate(self.structure.members.values()):
            def con_fun(X, i=i):
                return mem.ved / mem.VRd - 1
            constr = NonLinearConstraint(con_fun, name='Shear Force ' + str( i +1))
            self.cons.append(constr)


    def _create_bending_moment_constraints(self):
        """
        Creates bending moment constraints about y-axis
        """
        for i, mem in enumerate(self.structure.members.values()):
            def con_fun(X, i=i):
                return mem.med / mem.MRd[0] - 1
            constr = NonLinearConstraint(con_fun, name='Bending Moment ' + str( i +1))
            self.cons.append(constr)

    def create_constraints(self, N=True, V=True, M=True):

        # Initialize constraints as an empty list
        self.cons = []
        # Create normal force constraints
        if N:
            self._create_normal_force_constraints()
        if V:
            self._create_shear_force_constraints()
        if M:
            self._create_bending_moment_constraints()

if __name__ == '__main__':

    truss = KTruss24m()
    solver = DiscreteVNS()
    #solver.solve(truss, maxtime=60)
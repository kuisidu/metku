# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
from metku.optimization.solvers import *
from metku.optimization.problems.truss import TrussProblem
from metku.sections.steel.catalogue import ipe_profiles
from metku.frame2d.frame2d import *

ipe_profiles = list(ipe_profiles.keys())



frame = Frame2D(simple=[3, 2, 3400, 5000], supports='fixed')
beams = []
columns = []
for mem in frame.members.values():
    if mem.mtype == 'beam':
        beams.append(mem)
        frame.add(LineLoad(mem, [-50, -50], 'y'))
    else:
        columns.append(mem)


beam_group = {
    'name': 'Beam Variable',
    'objects': beams,
    'property': 'profile',
    'profiles': ipe_profiles,
    'var_type': 'index'
}

column_group = {
    'name': 'Column Variable',
    'objects': columns,
    'property': 'profile',
    'profiles': ipe_profiles,
    'var_type': 'index'
}

frame.generate()
frame.calculate()


problem = TrussProblem(
    name='Frame Example',
    structure=frame,
    groups=[beam_group, column_group]
)

solver = GA(pop_size=30)
solver.solve(problem, maxiter=50, verb=True)
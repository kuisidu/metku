#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2022. Metku team.
#  All rights reserved.

# @Filename : truss2d_testing
# @Date : 2022
# @Project: metku
# @AUTHOR : Jaakko Huusko
from matplotlib import pyplot as plt

import metku.truss2d as t2d
import metku.frame2d as f2d
import numpy as np

def joint_test():
    n_bays = 1
    n_storeys = 1
    L_bay = 10000
    H_storey = 5000

    simple_frame = [n_storeys, n_bays, H_storey, L_bay]

    frame = f2d.Frame2D(simple=simple_frame, create_beams=False, supports='fixed')

    simple_truss = dict(
        H0=H_storey,
        H1=1500,
        H2=2000,
        L1=L_bay / 2,
        H3=1500,
        dx=1000,
        n=16
    )

    truss = t2d.Truss2D(simple=simple_truss)

    for tc in truss.top_chords:
        truss.add(f2d.LineLoad(tc, [-30, -30], 'y', load_id=f2d.LoadIDs.ULS))
    truss.add(f2d.XYHingedSupport([0, simple_truss['H0'] + simple_truss['H1']]))
    truss.add(f2d.YHingedSupport(
        [truss.L1 + truss.L2, simple_truss['H0'] + simple_truss['H3']]))
    truss.generate()
    truss.calculate(load_id=f2d.LoadIDs.ULS)

    truss.H1 = 1000
    truss.H2 = 1500
    truss.H3 = 1300
    truss.dx = 100



    joint = truss.joints[4]
    w1, w2 = joint.webs.values()

    print("Before: ")
    print('e0', joint.e)
    print('g0', joint.g1)
    print('e1', w1.j1.e)
    print('g1', w1.j1.g1)
    print('e2', w2.j1.e)
    print('g2', w2.j1.g1)

    print("theta1", np.degrees(joint.theta1))
    print("theta2", np.degrees(joint.theta2))
    print("-" * 20)
    val = 50
    print("Changed g0 to:", val)
    print("-" * 20)
    joint.g1 = val
    print("After: ")
    print('e0', joint.e)
    print('g0', joint.g1)
    print('e1', w1.j1.e)
    print('g1', w1.j1.g1)
    print('e2', w2.j1.e)
    print('g2', w2.j1.g1)
    print("theta1", np.degrees(joint.theta1))
    print("theta2", np.degrees(joint.theta2))

    truss.plot(show=False)
    for joint in truss.joints.values():
        joint.plot_joint(show=False)
    plt.show()

def displacement_test():
    simple_truss = dict(
        H0=0,
        H1=1200,
        H2=2600,
        H3=1200,
        L1=12000,
        dx=1600,
        n=16
    )

    truss = t2d.Truss2D(simple=simple_truss)

    for tc in truss.top_chords:
        truss.add(f2d.LineLoad(tc, [-50, -50], 'y', load_id=f2d.LoadIDs.ULS))
    truss.add(f2d.XYHingedSupport([0, simple_truss['H0'] + simple_truss['H1']]))
    truss.add(f2d.YHingedSupport(
        [truss.L1 + truss.L2, simple_truss['H0'] + simple_truss['H3']]))
    truss.generate()
    truss.calculate(load_id=f2d.LoadIDs.ULS)
    print(truss.get_deflection(load_id=f2d.LoadIDs.SLS_Charasteristic))
    truss.plot_deflection()

if __name__ == '__main__':
    displacement_test()
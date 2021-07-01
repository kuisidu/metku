from structures.timber.timber_member import TimberMember
from structures.timber.timber_frame2d import TimberFrame2D
from materials.timber_data import T, Timber
from sections.timber.timber_section import TimberSection, GypsumPlasterboardF, WoodenFireProtection
import frame2d.frame2d as f2d
from sections.timber.timber_section import TimberSection, PulpettiPalkki

import numpy as np


fr = f2d.Frame2D()

coord1 = [[0, 0], [0, 6000]]
coord2 = [[0, 6000], [16000, 6000]]
coord3 = [[16000, 6000], [16000, 0]]
coord4 = [[8000, 0], [8000, 6000]]
coord5 = [[0, 6000], [0, 10000]]
coord6 = [[0, 10000], [16000, 10000]]
coord7 = [[16000, 6000], [16000, 10000]]
coord8 = [[0, 10000], [0, 14000]]
coord9 = [[0, 14000], [16000, 14000]]
coord10 = [[16000, 10000], [16000, 14000]]

fr.add(f2d.TimberFrameMember(coord1, TimberSection(240, 495), material=T.GL30c, mtype='column',
                             num_elements=20, ldc='inst'))
fr.add(f2d.TimberFrameMember(coord2, TimberSection(240, 720), material=T.GL30c, mtype='beam', ldc='mt',
                             num_elements=42, edge_load='compression', Sj1=17e9, Sj2=17e9))
fr.add(f2d.TimberFrameMember(coord3, TimberSection(240, 495), material=T.GL30c, mtype='column',
                             num_elements=20, ldc='inst'))
fr.add(f2d.TimberFrameMember(coord4, TimberSection(240, 360), material=T.C35, mtype='column',
                             num_elements=10))
fr.add(f2d.TimberFrameMember(coord5, TimberSection(240, 450), material=T.GL30c, mtype='column',
                             num_elements=20, ldc='inst'))
fr.add(f2d.TimberFrameMember(coord6, TimberSection(240, 720), material=T.GL30c, mtype='beam', ldc='mt',
                             num_elements=42, edge_load='compression', Sj1=17e9, Sj2=17e9, lateral_support_y=[1/3, 2/3]))
fr.add(f2d.TimberFrameMember(coord7, TimberSection(240, 450), material=T.GL30c, mtype='column',
                             num_elements=20, ldc='inst'))
fr.add(f2d.TimberFrameMember(coord8, TimberSection(240, 450), material=T.GL30c, mtype='column',
                             num_elements=20, ldc='inst'))
fr.add(f2d.TimberFrameMember(coord9, PulpettiPalkki(240, [1000, 400]), material=T.GL30c, mtype='beam', ldc='st',
                             num_elements=42, edge_load='compression', Sj1=17e9, Sj2=17e9))
fr.add(f2d.TimberFrameMember(coord10, TimberSection(240, 450), material=T.GL30c, mtype='column',
                             num_elements=20, ldc='inst'))


fr.add(f2d.LineLoad(fr.members[0], [3, 6], 'x'))
fr.add(f2d.LineLoad(fr.members[1], [-15, -15], 'y'))
fr.add(f2d.LineLoad(fr.members[2], [3, 1], 'x'))
fr.add(f2d.LineLoad(fr.members[4], [6, 10], 'x'))
fr.add(f2d.LineLoad(fr.members[5], [-15, -15], 'y'))
fr.add(f2d.LineLoad(fr.members[6], [3, 5], 'x'))
fr.add(f2d.LineLoad(fr.members[7], [10, 14], 'x'))
#fr.add(f2d.LineLoad(fr.members[8], [-10, -20], 'y'))
fr.add(f2d.LineLoad(fr.members[9], [5, 7], 'x'))
#fr.add(f2d.PointLoad([12000, 6000], [0, -10000, 0]))
#fr.add(f2d.PointLoad([16000, 3000], [-10000, 0, 0]))
#fr.add(f2d.PointLoad([0, 6000], [0, -80000, 140000]))

fr.add_self_weight()

fr.add(f2d.FixedSupport([0, 0]))
fr.add(f2d.FixedSupport([16000, 0]))
fr.add(f2d.FixedSupport([8000, 0]))

fr.generate()
fr.calculate()

# print('ned', fr.members[0].ned)
# print('med', fr.members[0].med)

print('\n1.krs palkin käyttöasteet')
fr.members[1].print_utility_factor()

print('\n1.krs vas.pilarin käyttöasteet')
fr.members[0].print_utility_factor()

print('\n1.krs oik.pilarin käyttöasteet')
fr.members[2].print_utility_factor()

print('\n1.krs keskipilarin käyttöasteet')
fr.members[3].print_utility_factor()

print('\n2.krs palkin käyttöasteet')
fr.members[5].print_utility_factor()

print('\n2.krs vas.pilarin käyttöasteet')
fr.members[4].print_utility_factor()

print('\n2.krs oik.pilarin käyttöasteet')
fr.members[6].print_utility_factor()

print('\n3.krs palkin käyttöasteet')
fr.members[8].print_utility_factor()

print(f'node {fr.members[8].n1.v}')

print('\n3.krs vas.pilarin käyttöasteet')
fr.members[7].print_utility_factor()

print('\n3.krs oik.pilarin käyttöasteet')
fr.members[9].print_utility_factor()


#print(f'Ln {fr.members[1].member.Ln()}')
#print(f'Ln {fr.members[0].member.Ln()}')
# print(f'Ln {fr.members[2].member.Ln()}')
# print(f'lamda {fr.members[2].member.lamda()}')
# print(f'lamda rel {fr.members[2].member.lamda_rel()}')
# print(f'k {fr.members[2].member.k()}')
# print(f'k_c {fr.members[2].member.k_c()}')


#
# for i in range(len(fr.members[1].member.L_c)):
#     print(f'l_ef {fr.members[1].member.l_ef(i)}')

# for n in fr.members[2].nodes.values():
#     print(n.parents)
#     print(fr.members[2].member.sigma_c90d(n))

# ks = list(fr.members[1].nodes.keys())
# for k in ks:
#     print(fr.members[1].to_global(0)[0] == fr.members[1].nodes[k].x and fr.members[1].to_global(0)[1] == fr.members[1].nodes[k].y)

# for i in range(fr.members[1].member.nsect()):
#     #print(fr.members[1].nodes[fr.members[1].nodes.keys()[i]])
#     print(fr.members[1].member.ned[i])


# fr.members[1].check_cross_section()
# print(fr.members[1].r)
# fr.members[1].R = 60
# fr.members[1].check_cross_section()
# print(fr.members[1].r)


#fr.to_robot('test')

fr.bmd(5)
#fr.plot_normal_force(5)
fr.plot_deflection(10)


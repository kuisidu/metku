from structures.timber.timber_member import TimberMember
from structures.timber.timber_frame2d import TimberFrame2D
from materials.timber_data import T, Timber
from sections.timber.timber_section import TimberSection, GypsumPlasterboardF, WoodenFireProtection
import frame2d.frame2d as f2d
from sections.timber.timber_section import TimberSection, PulpettiPalkki, HarjaPalkki, KaarevaHarjaPalkki

import numpy as np

fr = f2d.Frame2D()

coord1 = [[0, 0], [0, 6000]]
coord2 = [[0, 6000], [16000, 6000]]
coord4 = [[16000, 6000], [16000, 0]]

fr.add(f2d.TimberFrameMember(coord1, TimberSection(240, 360), material=T.GL30c, mtype='column', ldc='mt', num_elements=5))
#fr.add(f2d.TimberFrameMember(coord2, PulpettiPalkki(240, [1100, 1100]), material=T.GL30c, mtype='beam', ldc='mt',
#                             num_elements=10, edge_load='compression', lateral_support_z=None))
fr.add(f2d.TimberFrameMember(coord2, PulpettiPalkki(240, [1100, 1100]), material=T.GL30c, mtype='beam', ldc='mt',
                             num_elements=20, edge_load='compression'))
fr.add(f2d.TimberFrameMember(coord4, TimberSection(240, 360), material=T.GL30c, mtype='column', num_elements=5, ldc='mt'))

# fr.add(f2d.SteelFrameMember([[0, 0], [0, 5000]]))
# fr.add(f2d.SteelFrameMember([[0, 5000], [10000, 5000]], mtype='beam', num_elements=18))
# fr.add(f2d.SteelFrameMember([[10000, 5000], [10000, 0]]))

#fr.add(f2d.LineLoad(fr.members[0], [3, 6], 'x'))
fr.add(f2d.LineLoad(fr.members[1], [-10, -15], 'y'))
#fr.add(f2d.PointLoad([0, 6000], [10000, 0, 0]))
#fr.add_self_weight()

fr.add(f2d.FixedSupport([0, 0]))
fr.add(f2d.FixedSupport([16000, 0]))


fr.generate()
fr.calculate()

#print(f'node axial {fr.members[1].elements[5].axial_force}')
print(f'vasemman pilarin ned {fr.members[0].member.ned[-1]}')
print(f'oikean pilarin ned {fr.members[0].member.ned[0]}')

print(f'sig 90 d n1 {fr.members[1].member.sigma_c90d(fr.members[1].n1)}')
print(f'sig 90 d n2 {fr.members[1].member.sigma_c90d(fr.members[1].n2)}')

fr.members[1].print_utility_factor()
# print(f'kmod {fr.members[1].member.kmod}')
# print(fr.members[1].member.k_c())
# print(f'lamda {fr.members[1].member.lamda()}')
# for loc in fr.members[1].member.loc:
#     print(f'k_crit_y {loc}   {fr.members[1].member.k_crit(0, 0, loc)}')
#     print(f'k_crit_z {loc}   {fr.members[1].member.k_crit(0, 1, loc)}')
#     print(f'sigma {loc}   {fr.members[1].member.sigma_md(fr.members[1].member.loc.index(loc))}')
    # print(fr.members[1].member.sigma_md(fr.members[1].member.loc.index(loc))[0] /
    #       (fr.members[1].member.k_crit(0, 1, loc) * fr.members[1].member.fmd) ** 2) + \
    #       fr.members[1].member.sigma_c0d(fr.members[1].member.loc.index(loc)) / (
    #       fr.members[1].member.k_c()[1][0] * fr.members[1].member.fc0d))
# print(f'V {fr.members[1].cross_section.get_H(0.5)}')

#print('(6.41)', fr.members[1].member.check_tapered_beam_bending_tension())
#print('(6.50)', fr.members[1].member.check_perpendicular_to_grain_tension())


# print(f'\npilari vasen')
# fr.members[0].print_utility_factor()
# print(f'kmod {fr.members[0].member.kmod}')
# print(fr.members[0].member.k_c())
# print(f'\npilari oikea')
# fr.members[2].print_utility_factor()
# print(f'kmod {fr.members[2].member.kmod}')
# print(fr.members[2].member.k_c())
#
#
#
# print(fr.members[0].member.Ln())

# print(f' k_c\n{fr.members[1].member.k_c()}')
# print(f' i\n{fr.members[1].member.radius_of_gyration()}')
# print(fr.members[1].member.section_resistance(0))

# for i in range(len(fr.f.loads)):
#     print(fr.f.loads[i].qval)
#     print(fr.f.loads[i].xloc)

fr.bmd(10)
#fr.plot_deflection(10)
#fr.plot()
#fr.plot_loads()
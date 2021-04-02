from structures.timber.timber_member import TimberMember
from structures.timber.timber_frame2d import TimberFrame2D
from materials.timber_data import T, Timber
from sections.timber.timber_section import TimberSection, GypsumPlasterboardF, WoodenFireProtection
import frame2d.frame2d as f2d
from sections.timber.timber_section import TimberSection, PulpettiPalkki, HarjaPalkki, KaarevaHarjaPalkki, KaarevaPalkki

import numpy as np

fr = f2d.Frame2D()

fr.add(f2d.TimberFrameMember([[0, 0], [0, 6000]], TimberSection(174.85, 425.02, material=T.GL30c), num_elements=20,
                             mtype='column', beta=[[0.85], [0.85]], ldc='inst'))
fr.add(f2d.TimberFrameMember([[0, 6000], [16000, 6000]], HarjaPalkki(190.40, 681.23, 3000,
                                                                            material=T.GL30c), num_elements=20,
                             mtype='rigid-beam', lateral_support_z=[0.25, 0.5, 0.75], ldc='mt'))
fr.add(f2d.TimberFrameMember([[16000, 6000], [16000, 0]], TimberSection(174.85, 425.02, material=T.GL30c), num_elements=20,
                             mtype='column', beta=[[0.85], [0.85]], ldc='inst'))

fr.add(f2d.FixedSupport([0,0]))
fr.add(f2d.FixedSupport([16000, 0]))

#fr.members[1].cross_section.plot()

# fr.members[1].add_hinge(0)
# fr.members[1].add_hinge(1)

fr.add(f2d.LineLoad(fr.members[0], [5, 5], 'x', load_id=1))
fr.add(f2d.LineLoad(fr.members[1], [-15, -15], 'y', load_id=2))

#fr.members[1].R = 60

fr.generate()
fr.calculate()


# for i in range(len(fr.members[1].member.loc)):
#     print(fr.members[1].member.tau_tord(i))
#
# for loc in fr.members[1].member.loc:
#     print(fr.members[1].member.I_tor(loc))



# print([md[2] for md in fr.members[1].nodal_displacements.values()])
# print([max(md[2], key=abs) for md in fr.members[1].nodal_displacements.values()])
# print(max([abs(max(md[2], key=abs)) for md in fr.members[1].nodal_displacements.values()]))
# for nd in fr.members[1].nodal_displacements.values():
#     print(nd[2])

from log_to_file import log_to_file
@log_to_file
def pr():
    fr.members[0].print_utility_factor()
    fr.members[1].print_utility_factor()
    fr.members[2].print_utility_factor()
    fr.plot(save=True, show=False)
    fr.plot_deflection(show=False, save=True)
    fr.bmd(10, save=True, show=False)
pr()
#fr.bmd(10)
#fr.plot()

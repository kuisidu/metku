from materials.timber_data import T, Timber
from sections.timber.timber_section import TimberSection, GypsumPlasterboardF, WoodenFireProtection
import frame2d.frame2d as f2d
from sections.timber.timber_section import TimberSection, PulpettiPalkki, HarjaPalkki, KaarevaHarjaPalkki
import eurocodes.en1991.en1991 as en91
import numpy as np


fr = f2d.Frame2D()

coord1 = [[0, 0], [0, 6000]]
coord2 = [[0, 6000], [20000, 6000]]
coord4 = [[20000, 6000], [20000, 0]]

fr.add(f2d.TimberFrameMember(coord1, TimberSection(174.85, 425.02, material=T.GL30c), mtype='column', ldc='inst',
                             num_elements=5, beta=0.85))
fr.add(f2d.TimberFrameMember(coord2, HarjaPalkki(190.40, 681.23, 3000, material=T.GL30c), mtype='rigid-beam', ldc='mt',
                             num_elements=20, edge_load='compression', lateral_support_z=[0.25, 0.5, 0.75]))
fr.add(f2d.TimberFrameMember(coord4, TimberSection(174.85, 425.02, material=T.GL30c), mtype='column', ldc='inst',
                             num_elements=5, beta=0.85))

# fr.members[1].add_hinge(0)
# fr.members[1].add_hinge(1)

fr.add(f2d.FixedSupport([0, 0]))
fr.add(f2d.FixedSupport([20000, 0]))

qw = en91.WindLoad(0.5)
qs = en91.SnowLoad(2.0)
sw = en91.Load(1)
sls = en91.SLSCombiner('B', kk_jako=6)
sls.add([qw, 'x', [0]])
sls.add([qs, 'y', [1]])
sls.add([sw, 'y', [0, 1, 2]])
sls_result = sls.get_result()
uls = en91.ULSCombiner('B', kk_jako=6)
uls.add([qw, 'x', [0]])
uls.add([qs, 'y', [1]])
uls.add([sw, 'y', [0, 1, 2]])
uls_result = uls.get_result()
acc = en91.ACCCombiner('B', kk_jako=6)
acc.add([qw, 'x', [0]])
acc.add([qs, 'y', [1]])
acc.add([sw, 'y', [0, 1, 2]])
acc_result = acc.get_result()

for i, g in enumerate(sls_result[0][0]):
    if g == 0 or g is None:
        continue
    fr.add(f2d.LineLoad(fr.members[i], [g, g], 'x', load_id=1))

for i, g in enumerate(sls_result[0][1]):
    if g == 0 or g is None:
        continue
    fr.add(f2d.LineLoad(fr.members[i], [-g, -g], 'y', load_id=1))

for i, u in enumerate(uls_result[0]):
    if u == 0 or u is None:
        continue
    fr.add(f2d.LineLoad(fr.members[i], [u, u], 'x', load_id=2))

for i, u in enumerate(uls_result[1]):
    if u == 0 or u is None:
        continue
    fr.add(f2d.LineLoad(fr.members[i], [-u, -u], 'y', load_id=2))

for i, u in enumerate(acc_result[0]):
    if u == 0 or u is None:
        continue
    fr.add(f2d.LineLoad(fr.members[i], [u, u], 'x', load_id=3))

for i, u in enumerate(acc_result[1]):
    if u == 0 or u is None:
        continue
    fr.add(f2d.LineLoad(fr.members[i], [-u, -u], 'y', load_id=3))

fr.generate()
fr.calculate()

fr.members[0].print_utility_factor()
fr.members[1].print_utility_factor()
fr.members[2].print_utility_factor()
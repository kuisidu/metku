#from metku.frame2d import *
from metku.raami.raami import *
from metku.sections.steel.RHS import RHS

f = Raami()

L = 5000.0
H = 2000.0
nodes = [[0.0, 0.0], [0.0, H], [L, H], [L, 0.0]]

s = RHS(200, 200, 6)

for node in nodes:
    f.add(FrameNode(node))

f.add(SteelFrameMember(nodes=f.nodes[:2], section=s, mem_type='column', nel=6))
f.add(SteelFrameMember(nodes=f.nodes[1:3], section=s, mem_type='beam', nel=6, hinges=[True, True]))
f.add(SteelFrameMember(nodes=f.nodes[2:4], section=s, mem_type='column', nel=6))

f.add(FixedSupport(f.nodes[0]))
f.add(FixedSupport(f.nodes[-1]))

f.add(LineLoad(f.members[1], [-20, -20], 'y', coord_sys='local', load_id=0, ltype='snow', name='Lumi'))
f.add(LineLoad(f.members[0], [3, 3], 'x', coord_sys='global', load_id=1, ltype='wind', name='Tuuli x+'))

f.add(LoadCombination(comb_id=2, comb_type='ULS', load_cases=list(f.load_cases.values())))

f.load_combs[2].combine_loads()
f.plot()
f.generate_fem()
f.structural_analysis(load_id=2, support_method="REM")
f.bmd(scale=10, load_id=2)

#zelf toegoevoegd:
f.design_members(load_id=2)
f.print_member_utilization()

N = f.members[0].nodal_forces[2]['N']
Vz = f.members[0].nodal_forces[2]['Vz']
My = f.members[0].nodal_forces[2]['My']
print(N, Vz, My)
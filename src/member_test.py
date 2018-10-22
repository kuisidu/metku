
import hollow_sections as hss
import steel_member as member

p = hss.SHS(150,6)
L = 4000

m = member.SteelMember(p,L)

NbRd = m.buckling_strength()

print("NbRd =",format(NbRd[0]*1e-3,"#5.2f"),"kN")
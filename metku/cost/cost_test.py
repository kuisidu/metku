
"""
from workshop import BlastingCost

bc = BlastingCost(length=40,width=10)

bc.costs

print("Blasting")
for key, value in bc.costs.items():
    print("{0} cost = {1:4.2f} €/min".format(key,value))
    
bc.ProductiveTime(4000)

print(bc.cost())
"""

"""
tip_flow = 900000 # [mm3/min]

paints = {"alkyd": {"vol_solid":[0.48,0.45,0.45], "DFT":[80,40,40],"price":[2.8,3.6,3.6],"time":3,"TPP":0.0,"CP":0.0},
          "epoxy": {"vol_solid":[0.50,0.50,0.48], "DFT":[60,60,40],"price":[4.2,4.2,4.9],"time":14,"TPP":0.0,"CP":0.0},
          "polyurethane": {"vol_solid":[0.50,0.50], "DFT":[60,60],"price":[4.2,5.9],"time":27,"TPP":0.0,"CP":0.0},
          "acryl": {"vol_solid":[0.48,0.48,0.40], "DFT":[60,60,40],"price":[6.4,6.4,6.4],"time":7,"TPP":0.0,"CP":0.0},
          }

for paint, films in paints.items():
    print(paint)
    nfilms = len(films["vol_solid"])
    TPP = 0.0
    CP = 0.0
    for i in range(nfilms):
        tloss = 0.01+0.3*films["DFT"][i]*1e-3
        P0 = 1/films["vol_solid"][i]*(films["DFT"][i]*1e-3+tloss)
        TPP += P0
        CP += films["price"][i]*1e-6*P0

    TPP = TPP/tip_flow
    paints[paint]["TPP"] = TPP
    paints[paint]["CP"] = CP
"""

from workshop import AssemblingCostWeld

ac = AssemblingCostWeld(length=15,width=5)
 
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 19:57:54 2019

Column base joint according to EN 1993-1-8

@author: kmela
"""

import math

try:
    from src.eurocodes.en1993 import constants
    from src.eurocodes.en1993 import en1993_1_1
except:
    from eurocodes.en1993 import constants
    from eurocodes.en1993 import en1993_1_1
    from materials import steel_data
    import eurocodes.en1992.en1992_1_1 as en1992
    import eurocodes.en1992.constants as en1992_const
    from sections.steel.ISection import HEA
    from structures.steel.plates import RectPlate


# Anchor bolts
    
HPM_L = {16:{"d":16,"As":157,"L":280,"washer":[40,6],"dh":38,"k":10},
             20:{"d":20,"As":245,"L":350,"washer":[44,6],"dh":46,"k":12},
             24:{"d":24,"As":352,"L":430,"washer":[56,6],"dh":55,"k":13},
             30:{"d":30,"As":561,"L":500,"washer":[65,8],"dh":70,"k":15},
             39:{"d":39,"As":976,"L":700,"washer":[90,10],"dh":90,"k":18},
             }

class ColumnBase:
    
    
    def __init__(self,column_profile, base_plate, concrete="C25/30"):
        """ Constructor
            Input:
                column_profile .. SteelSection class object
                base_plate .. 
        """

        self.column = column_profile
        self.plate = base_plate
        self.anchor_bolt = None
        self.concrete = en1992_const.Concrete(concrete)        
        self.Ac0 = 1.0
        self.Ac1 = 1.0
        self.flange_weld = 0.0
        self.web_weld = 0.0
        self.bolt_rows = []
        self.beta_j = 2.0/3.0
        
    
    def add_bolt_row(self,bolt,z,p):
        """ Input:
            bolt .. bolt type
            z .. position of bolts relative to centroid of column
                 z < 0: left side
                 z > 0: right side
            p .. distance between bolts
        """
    
    @property
    def fcd(self):
        """ Design compression strength of concrete """        
        return self.concrete.fcd()
        #return en1992.fcd(en1992.concrete[self.concrete]["])

    @property
    def fjd(self):
        """ Design compression strength of concrete """                     
        # Strength of the foundation
        return self.beta_j*self.fcd        

    @property
    def tp(self):
        """ Thickness of end plate """
        return self.plate.t
        #return en1992.fcd(en1992.concrete[self.concrete]["])

    @property
    def fyp(self):
        """ Yield strength of end plate """
        return self.plate.material.fy
        #return en1992.fcd(en1992.concrete[self.concrete]["])
        
    @property
    def c(self):
        """ Length 'c' of EN 1993-1-8, Eq. (6.5) """
        return self.tp*math.sqrt(self.fyp/3/self.fjd)
    
    
    def beff(self,c):
        """ Effective width """
        return 2*c+self.column.tf
    
    
    def leff(self,c):
        """ Effective length """
        return 2*c+self.column.b
    
    def FCRd(self):
        """ Compression strength of T-stub 
            EN 1993-1-8, Eq. (6.4)
        """
        c = self.c
        beff = self.beff(c)
        leff = self.leff(c)
        fjd = self.fjd
        
        return fjd*beff*leff

    def Fc_fc_Rd(self):
        """ Column flange and web in compression 
            EN 1993-1-8 6.2.6.7
        """
        McRd = self.column.MRd[0]
        return McRd/(self.column.h-self.column.tf)

if __name__ == '__main__':
    
    col = HEA(240)
    w = col.b + 20
    d = col.h + 140
    plate = RectPlate(width=d,depth=d,thickness=30,material="S355")
        
    
    #conc = en1992.concrete["C25/30"]
    
    col.Ned = -263.22e3
    col.Med = 39.57e6
    col.Ved = 16.70e3
    
    joint = ColumnBase(col,plate,concrete="C25/30")
        
    
    #concrete = en1992_const.Concrete("C25/30")
    #fcd = concrete.fcd
    #beta_j = 2.0/3.0
    # Strength of the foundation
    #fjd = beta_j*fcd
    
    print("fcd = {0:4.2f} MPa".format(joint.fcd))
    print("fjd = {0:4.2f} MPa".format(joint.fjd))
    
    # Weld sizes
    joint.flange_weld = 6
    joint.web_weld = 6
    
    NEd = abs(joint.column.Ned)
    MEd = joint.column.Med
    h = joint.column.h
    tf = joint.column.tf
    NR = 0.5*NEd + MEd/(h-tf)
    NL = 0.5*NEd - MEd/(h-tf)
    
    print("NR = {0:4.2f}; NL ={1:4.2f}".format(NR,NL))
    
    e1 = 50.0
    p = 180.0
    
    if NL < 0:
        zR = 0.5*(h-tf)
        zL = 0.5*h+e1        
        print("zR = {0:4.2f} mm ; zL ={1:4.2f}".format(zR,zL))

    z = zR+zL    
    print("z = {0:4.2f} mm".format(z))
    
    NR_Ed = NEd*zL/z + MEd/z
    NL_Ed = NEd*zR/z - MEd/z
    
    print("NR_Ed = {0:4.2f} kN".format(NR_Ed*1e-3))
    print("NL_Ed = {0:4.2f} kN".format(NL_Ed*1e-3))
    
    c = joint.c
    
    print("c = {0:4.2f} mm".format(c))
    
    FCRd = joint.FCRd()
    Fc_fc_Rd = joint.Fc_fc_Rd()
    
    print("FCRd = {0:4.2f} kN".format(FCRd*1e-3))
    print("Fc_fc_Rd = {0:4.2f} kN".format(Fc_fc_Rd*1e-3))
    
    m = e1-0.8*math.sqrt(2)*joint.flange_weld
    outstand_top = 110
    d_bolt = 24
    
    Fbolt_Rd = 0.9*550*352/1.25
    print("Fbolt_Rd = {0:4.2f} kN".format(Fbolt_Rd*1e-3))
    
    e = outstand_top + 0.5*d_bolt
    ex = c-e1
    
    leff_cp = min([2*math.pi*m,math.pi*m+p,math.pi*m+2*e])
    leff_nc = min([4*m+1.25*e,e+2*m+0.625*ex,0.5*(joint.column.b+2*c),0.5*p+2*m+0.625*ex])
    
    print("m = {0:4.2f} mm".format(m))
    print(leff_cp)
    print(leff_nc)
    
    leff_1 = min(leff_nc,leff_cp)
    Mpl_1_Rd = 0.25*leff_1*joint.tp**2*joint.fyp
    FT_12_Rd = Mpl_1_Rd/m
    
    print("Mpl_1_Rd = {0:4.2f} kNm".format(Mpl_1_Rd*1e-6))
    print("FT_12_Rd = {0:4.2f} kN".format(FT_12_Rd*1e-3))
    
    F_t_wb_Rd = leff_1*joint.column.tw*joint.column.fy
    print("Ft_wb_Rd = {0:4.2f} kN".format(F_t_wb_Rd*1e-3))
    
    
    
    #joint.fcd
    
# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Wed Oct  2 19:57:54 2019

Column base joint according to EN 1993-1-8

@author: kmela
"""

import math

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines


    
from metku.eurocodes.en1993 import constants
from metku.eurocodes.en1993 import en1993_1_1
from metku.materials import steel_data
import metku.eurocodes.en1992.en1992_1_1 as en1992
import metku.eurocodes.en1992.constants as en1992_const
from metku.sections.steel.ISection import HEA
from metku.structures.steel.plates import RectPlate
from metku.materials.steel_data import Steel
from metku.eurocodes.en1993.en1993_1_8.en1993_1_8 import TStubEndPlate, BoltRow, ROW_OUTSIDE_BEAM_TENSION_FLANGE, FIRST_ROW_BELOW_BEAM_TENSION_FLANGE
from metku.eurocodes.en1993.en1993_1_8.en1993_1_8 import ROW_OUTSIDE_BEAM_TENSION_FLANGE, FIRST_ROW_BELOW_BEAM_TENSION_FLANGE, TENSION_ROW, COMPRESSION_ROW, OTHER_INNER_ROW

# Anchor bolts
    
HPM_L = {16:{"d":16,"As":157,"L":280,"washer":[40,6],"dh":38,"k":10, "material": "B500B"},
             20:{"d":20,"As":245,"L":350,"washer":[44,6],"dh":46,"k":12, "material": "B500B"},
             24:{"d":24,"As":352,"L":430,"washer":[56,6],"dh":55,"k":13,"material": "B500B"},
             30:{"d":30,"As":561,"L":500,"washer":[65,8],"dh":70,"k":15, "material": "B500B"},
             39:{"d":39,"As":976,"L":700,"washer":[90,10],"dh":90,"k":18, "material": "B500B"},
             }

class AnchorBolt:
    """ For anchor bolt data """
    
    def __init__(self,bolt_type=HPM_L[24]):
        """ Constructor
            input:
                bolt_type .. from a list of bolts above
                
            attributes:
                d .. diameter of the shank (mm)
                As .. stress area of the threaded shank (mm2)
                L .. length of the bolt (mm)
                washer .. [0] washer diameter (mm)
                          [1] washer thickness (mm)
                dh .. diameter of the head [mm]
                k .. thickness of the head (mm)
                material .. 
        """
        
        self.d = bolt_type["d"]
        self.As = bolt_type["As"]
        self.L = bolt_type["L"]
        self.washer = bolt_type["washer"]
        self.dh = bolt_type["dh"]
        self.k = bolt_type["k"]
        self.material = Steel(bolt_type["material"])

class ColumnBase:
    
    
    def __init__(self,column_profile, plate_edges, plate_thickness, anchor_bolts, bolt_x=None, bolt_p=100, concrete="C25/30"):
        """ Constructor
            Input:
                column_profile .. SteelSection class object
                plate_edges .. list of edge distances
                               [0] .. distance from left column flange to edge
                               [1] .. distance from right column flange to edge
                               [2] .. distance from top edge to column flange
                               [3] .. distance from bottom edge to column flange
                plate_thickness .. thickness of the base plate [mm]
                anchor_bolts
                        
        """

        self.column = column_profile
        self.beam = column_profile
        self._c = plate_edges
        h = column_profile.h
        b = column_profile.b
        hp = h + c[0] + c[1]
        bp = b + c[2] + c[3]
        base_plate = RectPlate(bp,hp,plate_thickness,material="S355")
        self.plate = base_plate
        self.anchor_bolt = anchor_bolts
        self.p = bolt_p
        if bolt_x is None:
            bolt_x = [0,0]
            bolt_x[0] = -0.5*h-0.5*c[0]
            bolt_x[1] = 0.5*h+0.5*c[1]
            
        self.xb = bolt_x
        self.concrete = en1992_const.Concrete(concrete)        
        self.Ac0 = 1.0
        self.Ac1 = 1.0
        self.flange_weld = 5.0
        self.web_weld = 5.0
        self.bolt_rows = []
        self.add_bolt_row(anchor_bolts, bolt_x[0], self.p)
        self.add_bolt_row(anchor_bolts, bolt_x[1], self.p)
        self.beta_j = 2.0/3.0
        
    
    def add_bolt_row(self,bolt,z,p):
        """ Input:
            bolt .. bolt type
            z .. position of bolts relative to centroid of column
                 z < 0: left side
                 z > 0: right side
            p .. distance between bolts
        """
        #(self,bolt,p,z,flange_loc=INNER_ROW, plate_loc=OTHER_INNER_ROW,joint=None, row_type=TENSION_ROW)
        
        if z < -0.5*self.column.h:
            ploc = ROW_OUTSIDE_BEAM_TENSION_FLANGE
            row_type = TENSION_ROW
        elif z < 0:
            ploc = FIRST_ROW_BELOW_BEAM_TENSION_FLANGE
            row_type = TENSION_ROW
        else:
            ploc = OTHER_INNER_ROW
            row_type = COMPRESSION_ROW
        
        new_row = BoltRow(bolt,p,z,plate_loc=ploc,joint=self,row_type=row_type)
        
        self.bolt_rows.add(new_row)
    
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
    
    def bolt_coord(self):
        """ Returns list of tuples that are the xy-coordinates of the bolts """
        coords = []
        
        ypos = 0.5*self.p
        yneg = -ypos
        for xb in self.xb:    
            coords.append((xb,yneg))
            coords.append((xb,ypos))
            
        return coords
    
    def Fc_pl_Rd(self,verb=False):
        """ Compression strength of T-stub 
            EN 1993-1-8, Eq. (6.4)
        """
        c = self.c
        beff = self.beff(c)
        leff = self.leff(c)
        fjd = self.fjd
        
        FRd = fjd*beff*leff
        
        if verb:
            print("Compression resistance of T-stub:")
            print(f"c = {c:4.2f} mm")
            print(f"beff = {beff:4.2f} mm")
            print(f"fjd = {fjd:4.2f} MPa")
            print(f"FC,Rd = {FRd*1e-3:4.2f} kN")
        
        return FRd

    def Fc_fc_Rd(self,verb=False):
        """ Column flange and web in compression 
            EN 1993-1-8 6.2.6.7
        """
        McRd = self.column.MRd
        Fc = McRd/(self.column.h-self.column.tf)
        
        if verb:
            print("Column flange and web in compression")
            print(f"McRd = {McRd*1e-6:4.2f} kNm")
            print(f" h = {self.column.h:4.2f} mm")
            print(f" tf = {self.column.tf:4.2f} mm")
            print(f" Fc_fc_Rd = {Fc*1e-3:4.2f} kN")
        
        return Fc
    
    def FCRd(self,verb=False):
        """ Compression strength of the connection """
        
        FcplRd = self.Fc_pl_Rd(verb)
        FcfcRd = self.Fc_fc_Rd(verb)
        
        FC = min(FcplRd,FcfcRd)
        
        return FC
    
    def Ft_ep_Rd(self,verb=False):
        """ Tension side: bending of base plate.
            EN 1993-1-8, 6.2.6.5. Base plate is treated as end plate.
        """
        
        ex = 0.5*self.plate.h+self.xb[0]
        mx = 0.5*self.column.h+0.8*self.flange_weld*math.sqrt(2) + self.xb[0]
        
        T = TStubEndPlate(self.bolt_rows[0],
                          material=self.plate.material,
                          tf=self.plate.t,
                          emin=ex,
                          m=mx,
                          e=ex,
                          w=self.p,
                          bp=self.plate.b)
        
        F12Rd = T.F_T_12_Rd()
        F3Rd = T.F_T_3_Rd()
        
        
    
    def draw(self, name=""):
        """ Draw the connection """
        fig, ax = plt.subplots(1,2)
        ymax = 150
        
        self.column.draw(axes=ax[1])
        
        # Side view:
        # Origin is located at the bottom of the steel profile center line
        hp, bp, tp = self.plate.h, self.plate.b, self.plate.t
        base_plate = patches.Rectangle((-0.5*hp, -tp), width=hp, height=tp, fill=False, hatch='\\')
        
        base_plate_plan = patches.Rectangle((-0.5*hp, -0.5*bp), width=hp,
                                       height=bp, fill=False, hatch='//')
        
        
        # Draw bolts
        dbolt = self.anchor_bolt.d
        lbolt = self.anchor_bolt.L
        left_bolt = patches.Rectangle((self.xb[0]-0.5*dbolt, -tp-lbolt), width=dbolt,
                                       height=lbolt+tp, fill=True)
        
        right_bolt = patches.Rectangle((self.xb[1]-0.5*dbolt, -tp-lbolt), width=dbolt,
                                       height=lbolt+tp, fill=True)
        """
        bot_flange = patches.Rectangle((-0.5*self.bb, 0), width=self.bb,
                                       height=self.tb, fill=False, hatch='\\')
        web = patches.Rectangle((-0.5*self.tw, self.tb), width=self.tw,
                                height=self.hw, fill=False, hatch='\\')
    
        top_flange = patches.Rectangle((-0.5*self.bt, self.h-self.tt),
                                       width=self.bt, height=self.tt,
                                       fill=False, hatch='\\')
        """
        ax[0].add_patch(base_plate)
        h = self.column.h
        tf = self.column.tf
        ax[0].vlines([-0.5*h,-0.5*h+tf,0.5*h,0.5*h-tf],0,ymax)
        ax[0].vlines(0,0,ymax,linestyles='dashdot')
        ax[0].add_patch(left_bolt)
        ax[0].add_patch(right_bolt)
        ax[0].set_xlim(-0.5*hp, 0.5*hp)
        if self.anchor_bolt is None:
            ymin = -self.plate.t
        else:
            ymin = -self.anchor_bolt.L
        
        ax[0].set_ylim(ymin, ymax)
    
        
        ax[0].set_aspect('equal')
        
        xtick = ax[0].get_xticks()
        #print(xtick)
        
        """ Draw connection from above """
        ax[1].add_patch(base_plate_plan)
        
        # Draw column profile
        
        
        # Draw bolts
        xy_bolts = self.bolt_coord()

        for xy in xy_bolts:
            ax[1].add_patch(patches.Circle(xy,radius=0.5*dbolt,linestyle='solid'))

        
        ax[1].set_xticks(xtick)
        #ax[1].set_xlim(-0.5*hp, 0.5*hp)
        ax[1].set_ylim(-0.5*hp, 0.5*hp)
        ax[1].set_aspect('equal')
        
        return ax
        
if __name__ == '__main__':
    
    col = HEA(240)
    #w = col.b + 20
    #d = col.h + 140
    #plate = RectPlate(width=d,depth=d,thickness=30,material="S355")
        
    
    #conc = en1992.concrete["C25/30"]
    
    col.Ned = -263.22e3
    col.Med[0] = 39.57e6
    col.Ved[1] = 16.70e3
    
    c = [100,100,50,50]
    tp = 40
    bolts = AnchorBolt(HPM_L[24])

    joint = ColumnBase(col,plate_edges = c,\
                       plate_thickness = tp,\
                       anchor_bolts = bolts,\
                       bolt_p = 100,\
                       concrete="C30/37")
        
    ax = joint.draw()
    
    #concrete = en1992_const.Concrete("C25/30")
    #fcd = concrete.fcd
    #beta_j = 2.0/3.0
    # Strength of the foundation
    #fjd = beta_j*fcd
    
    """
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
    
    """
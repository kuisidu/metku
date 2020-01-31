# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 19:41:00 2019

@author: kmela
"""


import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines

from plates import RectPlate
from eurocodes.en1993.en1993_1_8.en1993_1_8 import BoltRow, BoltRowGroup
from eurocodes.en1993.en1993_1_8.en1993_1_8 import END_ROW, INNER_ROW, ROW_OUTSIDE_BEAM_TENSION_FLANGE, FIRST_ROW_BELOW_BEAM_TENSION_FLANGE, OTHER_END_ROW
from eurocodes.en1993.en1993_1_8.en1993_1_8 import TENSION_ROW, SHEAR_ROW, COMBI_ROW
import eurocodes.en1993.en1993_1_8.component_method as cm

from cost.workshop import Workshop
from cost.cost_data import BASIC_STEEL_PRICE, steel_grade_add_on, thickness_add_on
from cost.cost_data import bolt_unit_cost_full_thread as bolt_costs
from cost.cost_data import nut_unit_cost, washer_unit_cost

class EndPlateJoint:
    """ Class for end-plate connections """
    
    
    
    def __init__(self,column,beam,tp,bp,mat_p,etop,ebottom,bolt,y_bolts,e_bolts,
                 bolt_row_pos, groups=[],group_pos=None, row_types=[]):
        """ Contructor
            Input:
                column .. ISection or WISection profile
                beam .. ISection or WISection profile
                tp .. thickness of end plate [mm]
                mat_p .. plate material (e.g. "S355")
                etop .. extension of end plate above tension flange [mm]
                ebottom .. extension of end plate below compression flange
                bolt .. Bolt class object
                y_bolts .. array of vertical positions of bolt rows
                           (origin is located in the centroid of the beam)
                e_bolts .. edge distance of bolts (from centroid of bolt hole)
                            measured to the edge of the end plate [mm]
                bolt_row_pos .. list of dicts notifying the position of the row
                                for column flange in bending and end plate in bending
                                components
                                {"flange": flange_pos, "plate": plate_pos}
                groups .. list of bolt row groups
                group_pos . list of bolt row positions in groups
                    group_pos[i] is a list of length of the number of bolt rows in group 'i'.
                    each item in this list is a bolt row position code, e.g. INNER_ROW, see en1993_1_8.py
        """
        
        self.col = column
        self.beam = beam
        
        # calculate plate height
        hp = etop + ebottom + beam.h        

        
        self.end_plate = RectPlate(bp,hp,tp,material=mat_p)
        self.etop = etop
        self.ebottom = ebottom
        self.ebolts = e_bolts        
        # Make bolt rows
        self.bolt_rows = []
        self.row_groups = []
        
        # Beam flange weld, throat thickness [mm] (to plate)
        self.weld_f = 5
        # Beam web weld, throat thickness [mm] (to plate)
        self.weld_w = 5
        
        # Plate to column weld, throat thickness [mm]
        # Only, if the plate is welded to the column, i.e. no bolts
        self.weld_p = 5
        
        # Distance between centroids of bolts in a row (equal for all rows)
        self.p = bp-2*e_bolts
        
        # Geometrical properties
        
        # Create bolt rows
        for y, row_pos, r_type in zip(y_bolts,bolt_row_pos,row_types):
            """ Determine location of bolt row with respect to
                column flange in bending and
                end plate in bending components
            """
            flange_pos = row_pos["flange"]
            end_plate_pos = row_pos["plate"]
            
            """
            if y > 0.5*beam.h:
                "" Row is above the beam, so it is an end row
                    for column flange and row outside beam tension
                    flange for end plate
                ""
                flange_pos = END_ROW
                end_plate_pos = ROW_OUTSIDE_BEAM_TENSION_FLANGE
                
            else:
                "" Rows below top (tension) flange of the beam:
                    There are several possibilities:
                        - only one bolt row below the beam flange
                        - more than one rows below the beam flange
                        - 
                ""
            """
            self.bolt_rows.append(BoltRow(bolt,self.p,y,flange_pos,end_plate_pos,joint=self,row_type=r_type))
    
        # Create bolt row groups
        for i, group in enumerate(groups):
            self.row_groups.append(BoltRowGroup([self.bolt_rows[g] for g in group],
                                                group_pos[i])
                                   )
        
        self._cost = {'plate':0.0,                                         
                     'welding':0.0,
                     'holes':0.0,
                     'bolts':0.0}
        
    @property
    def nrows(self):
        """ Number of bolt rows """
        return len(self.bolt_rows)
    
    @property
    def tension_rows(self):
        """ Number of rows taking tension """
        n = 0
        
        for row in self.bolt_rows:
            if row.type != SHEAR_ROW:
                n += 1
        
        return n
    
    @property
    def zc(self):
        """ Center of compression: located at the center line of the bottom flange
            of the beam
        """
        return -0.5*(self.beam.h-self.beam.tf)
        
    
    @property
    def bp(self):
        """ Width of end plate """
        return self.end_plate.b
    
    @property
    def tp(self):
        """ Thickness of end plate """
        return self.end_plate.t
    
    def emin(self,end_plate_extension=False):
        """ Distance emin used for T-stubs 
            See EN 1993-1-8, 
            6.2.6.4.1(2) for unstiffened column flange
            6.2.6.4.2(4) for stiffened column flange, and
            6.2.6.5(3) for end plates
            
            They all refer to Fig. 6.8. of EN 1993-1-8
            
        """
        
        if end_plate_extension:
            emin = self.ex()
        
        
    @property
    def beta(self):
        """ Effect of shear in the column web to moment resistance of the joint 
            So far, only one value is given (approximation for one-sided joint)
        """
        return 1.0
        
        
    def V_wp_Rd(self):
        """ Column web shear resistance """
        return cm.col_web_shear(self.col)
    
    def k1(self,z=None):
        """ Stiffnes factor of column web in shear """
        if z is None:
            z = self.zeq()
    
        return cm.col_web_shear_k(self.col,self.beta,z)
            
    
    def Fc_wc_Rd(self):
        """ Column web in compression """
        
        # Calculate compression stre
        sigma_com_Ed = self.col.sigma_com()
        Fc_wc_Rd, beff_wc = cm.col_web_trv_comp(self.col, self.beam, self.end_plate.t, self.ebottom, self.weld_f, self.beta, sigma_com_Ed)
        self.beff_wc = beff_wc
        print("beff_wc_Rd = {0:4.2f} mm".format(beff_wc))
        return Fc_wc_Rd
    
    def beff_c_wc(self):
        
        # Column information
        r_c = self.col.r                                     # Column rounding, VALSSATTU
        t_fc = self.col.tf                                  # Column flange thickness
        t_wc = self.col.tw                                  # Column web thickness
        d_wc = self.col.hw                                  # Column web height
        #A_vc = self.col.Ashear                                  # Shear area of column web
        #f_y_wc = self.col.fy                                # Yield strength of column material
        #E = self.col.column.E                                  # Young's modulus of column material

        # Beam information
        t_fb = self.beam.tf                           # Beam flange thickness

        # Plate information
        t_p = self.tp
        lo_p = self.ebottom
        
        # weld between plate and bottom flange
        a_p = self.weld_f            

        # Width of compressed zone under compressed flange using 45 deg stress distribution
        s_p = min(2.0*t_p, t_p+lo_p-math.sqrt(2.0)*a_p)           # lo_p is overhang of end-plate over lower flange of beam

        # Effective length of compressed web
        s = r_c
        b_eff_c_wc = t_fb+2.0*math.sqrt(2.0)*a_p+5.0*(t_fc+s)+s_p
        
        return b_eff_c_wc
    
    
    def k2(self):
        """ Stiffness factor for column web in compression """
        
        return cm.col_web_trv_comp_k(self.col,self.beff_c_wc())
    
    def Fc_fb_Rd(self):
        """ Beam flange and web in compression """        
        return cm.beam_web_compression(self.beam)
    
    def MjRd(self,verb=False):
        """ Bending moment resistance """
        
        # Calculate column web shear resistance
        # and compression resistances
        V_wp_Rd = self.V_wp_Rd()
        Fc_wc_Rd = self.Fc_wc_Rd()
        Fc_fb_Rd = self.Fc_fb_Rd()
        
        Fcom_Rd = min(V_wp_Rd/self.beta,Fc_wc_Rd,Fc_fb_Rd)
        
        # Sum of tension forces in the bolt rows
        Ft_total = 0.0
        
        # Calculate resistances of bolt rows
        for i, row in enumerate(self.bolt_rows):
            if row.type != SHEAR_ROW:
                if verb:
                    print("  ** ROW {0:1.0f} **".format(i+1))
        
                FtRd = row.tension_resistance(verb)
                
                """ The total force in the bolt rows, including
                    the current bolt row, cannot exceed the
                    compression resistance (including shear of column web panel)
                    of the joint
                """
                row.FtRd = min(FtRd,Fcom_Rd-Ft_total)
                
                Ft_total += row.FtRd
                
                """ Check if the tension force in the row
                    Needs to be reduced due to compression side
                    or shear resistance in the column web
                """
        
        if verb:
            print("--------  ROW RESISTANCES ---------")
            for i, row in enumerate(self.bolt_rows):            
                if row.type != SHEAR_ROW:
                    print("ROW {0:1.0f}: FtRd = {1:4.2f} kN".format(i+1, row.FtRd*1e-3))
            
        # Calculate resistances of bolt row groups 
        if verb:
            print("*** BOLT ROW GROUPS ***")
        
        for i, group in enumerate(self.row_groups):
            if verb:
                print("  * GROUP {0:4.0f} * ".format(i+1))
    
            Ft_group = group.tension_resistance(verb)
            # TODO: REDUCE RESISTANCE OF ROWS BASED ON Ft
            Ft_tot = sum([row.FtRd for row in group.rows])
            
            if verb:
                print("Ft_Group = {0:4.2f} kN".format(Ft_group*1e-3))
                print("Ft_rows = {0:4.2f} kN".format(Ft_tot*1e-3))
            
            if Ft_group < Ft_tot:
                print("Need to reduce strength of rows")
        
        MjRd = 0.0
        for row in self.bolt_rows:
            MjRd += row.FtRd * row.h
        
        if verb:
            print("MjRd = {0:4.2f} kNm".format(MjRd*1e-6))
        
        return MjRd
    
    def hr(self,which_rows='all'):
        """ Returns an array of bolt row distances from 
            Compression centroid
        """
        if which_rows == 'all':
            hr = np.array([row.h for row in self.bolt_rows])
        elif which_rows == 'non-shear':
            hr = np.array([row.h for row in self.bolt_rows if row.type != SHEAR_ROW])
        
        return hr

    def keq(self,verb=False):
        """ Equivalent stiffness of the tension rows """
        keff = np.array([row.keff() for row in self.bolt_rows])
        hr = self.hr()
        
    
        return sum(keff*hr)/self.zeq(verb)
    
    def zeq(self,verb=False):
        """ Equivalent moment arm """
        zeq = 0.0
        
        keff = np.array([row.keff() for row in self.bolt_rows])
        hr = self.hr()
        
        zeq = sum(keff*hr**2)/sum(keff*hr)
        
        return zeq
    
    def zeq_and_keq(self,verb=False):
        """ Return both zeq and keq 
            To avoid double computations
        """
        zeq = 0.0
        
        keff = np.array([row.keff() for row in self.bolt_rows if row.type != SHEAR_ROW])
        hr = self.hr('non-shear')
        
        zeq = sum(keff*hr**2)/sum(keff*hr)
        keq = sum(keff*hr)/zeq
        
        if verb:
            print(" * zeq and keq *")
            for k, h in zip(keff,hr):
                print("keff = {0:4.2f} mm; hr = {1:4.2f} mm".format(k,h))
            
        
        return zeq, keq
        
    def moment_arm(self):
        """ Moment arm for joints with two bolt rows or """
        
        """ hr is an array of distances of bolt rows from
            compression center
        """
        hr = self.hr()
        
        return max(hr)
            

    def k_comp(self,z):
        """ Combined stiffness of compression side components,
            namely
            column web in shear and
            column web in compression
        """
        k1 = self.k1(z)
        k2 = self.k2()
        
        return 1/(1/k1+1/k2)
        

    def Sj_ini(self,verb=False):
        """ Initial rotational stiffness """
        
        if verb:
            print("*** INITIAL ROTATIONAL STIFFNESS ***")
            
        
        E = self.col.material.E
        
        if self.tension_rows == 1:
            """ Only one tension row: moment arm is the distance
                between bottom 
            """
            z = self.moment_arm()
        else:
            z, k_eq = self.zeq_and_keq(verb)
        
        kcomp = self.k_comp(z)
        
        if verb:
            print("   Moment arm: z = {0:4.2f} mm".format(z))
            print("   Stiffness factors:")
            print("     k1 = {0:4.2f} mm".format(self.k1(z)))
            print("     k2 = {0:4.2f} mm".format(self.k2()))
            print("     keq = {0:4.2f} mm".format(k_eq))
            print("     kcomp = {0:4.2f} mm".format(kcomp))
        
        
        
        
        Sj_ini = E*z**2/(1/kcomp+1/k_eq)
        
        
        if verb:
            print("  Sj_ini = {0:4.2f} kNm/rad".format(Sj_ini*1e-6))
        
        return Sj_ini    
        
    def cost(self,workshop,material_cost=BASIC_STEEL_PRICE,verb=False):
        """ Calculate cost of the joint
            units: e
        """
        
        self._cost['plate'] = self.end_plate.cost(workshop,material_cost)
        hole_cut_length = 0.0
        for row in self.bolt_rows:
            """ Costs for bolts and hole forming 
                bolt costs inlucde one nut and two washers
            """            
            self._cost['bolts'] += 2*bolt_costs[int(row.bolt.length)][int(row.bolt.d)]
            self._cost['bolts'] += nut_unit_cost[int(row.bolt.d)]
            self._cost['bolts'] += washer_unit_cost[int(row.bolt.d)]
            """ Hole costs include forming holes both in the plate in the column """
            hole_cut_length += 2*row.bolt.d0*math.pi
           
            #self._cost['holes'] += 2*workshop.cost_centres['cutting'].cost(self.col.tf,row.bolt.d0*math.pi)

        self._cost['holes'] = workshop.cost_centres['cutting'].cost(self.tp,hole_cut_length)

        """ Welding cost """
        flange_weld_length = 2*(self.beam.b + self.beam.tf) - 2*self.beam.r - self.beam.tw + self.beam.r*math.pi
        print(flange_weld_length)
        self._cost['welding'] += workshop.cost_centres['assembly_welding'].cost(self.weld_f,2*flange_weld_length)
        
        web_weld_length = 2*self.beam.hw
        print(web_weld_length)
        self._cost['welding'] += workshop.cost_centres['assembly_welding'].cost(self.weld_w,web_weld_length)
        
        total_cost = 0.0
        for c in self._cost.values():
            total_cost += c
        
        self._total_cost = total_cost
        
        if verb:
            self.cost_distribution()
        
        return total_cost
    
    def cost_distribution(self,pie=False):
        """ Print cost distribution """
        
        if pie:
            pass
        else:
            print(""" Cost distribution of the connection """)
            p_tot = 0.0
            for key, value in self._cost.items():
                p =  value/self._total_cost*100
                p_tot += p
                print("   {0:10s}: {1:5.2f} € [{2:5.2f} %]".format(key,value,p))
            
            print("-----------------------")
            print("   {0:10s}: {1:5.2f} € [{2:5.2f} %]".format("Total",self._total_cost,p_tot))
    
    def info(self,draw=False):
        """ Print relevant information of the joint """
        
        print("Resistance of connection:")
        print(" V_wp_Rd = {0:{fm}} [kN] (Column web panel in shear)".format(self.V_wp_Rd()*1.0e-3, fm='7.2f'))
        print("F_c_wc_Rd = {0:{fm}} [kN] (Column web in transverse compression)".format(self.Fc_wc_Rd()*1.0e-3, fm='7.2f'))
        print("F_c_fb_Rd = {0:{fm}} [kN] (Beam flange and web in compression)".format(self.Fc_fb_Rd()*1.0e-3, fm='7.2f'))
        
        if draw:
            self.draw()
        
        
    def draw(self, name=""):
        """ Draw the connection 
            ax[1] .. front view
            ax[0] .. side view
        """
        
        grey = '0.8'
        fig, ax = plt.subplots(1,2)
        ymax = 150
                        
        t_fc = self.col.tf
        r_c = self.col.r
        h_c = self.col.h
        b_c = self.col.b
        h_b = self.beam.h
        b_b = self.beam.b
        t_fb = self.beam.tf
        t_wb = self.beam.tw
        #d_wb = self.beam.d_w
        a_f = self.weld_f
        r_b = self.beam.r
        h_p = self.end_plate.h
        t_p = self.end_plate.t
        b_p = self.end_plate.b
        up_p = self.etop
        bot_p = self.ebottom
        
        #self.column.draw(axes=ax[1])
        
        # Side view:
        ax_side = ax[0]
        # Origin 
        # Draw column
        #plt.figure("side " + str(self.plot_nro))
        #plt.clf()
        #plt.xlabel("[mm]")
        #plt.ylabel("[mm]")

        axis_lim = max(h_c,h_b,h_p) + 200
        #axis_lim = max(h_c + 500.0, h_b + 500.0, h_p + 200.0)

        # Plot column in figure
        column_y = [-100.0, axis_lim - 100.0]
        
        ax_side.vlines([0,-t_fc,-(h_c-t_fc),-h_c],column_y[0],column_y[1],'k',zorder=2)
        ax_side.vlines([-(t_fc + r_c),-(h_c - (t_fc + r_c))],column_y[0],column_y[1],grey,zorder=1)
        
        # Plot plate in figure
        plate_x = [0.0, t_p, t_p, 0.0]
        plate_y = [0.0, 0.0, h_p, h_p]
        ax_side.plot(plate_x, plate_y, 'k', zorder=2)
        
        # Plot beam in figure
        beam_x = [t_p, axis_lim - (100.0 + h_c + t_p)]
        ax_side.hlines([bot_p,bot_p+t_fb,bot_p+h_b,bot_p+h_b-t_fb],beam_x[0],beam_x[1],'k',zorder=2)
        ax_side.hlines([bot_p+t_fb+r_b,bot_p+h_b-t_fb-r_b],beam_x[0],beam_x[1],grey,zorder=1)

        # Plot welds
        ax_side.plot([t_p, t_p + math.sqrt(2.0)*a_f], [bot_p + h_b + math.sqrt(2.0)*a_f, bot_p+h_b], 'k', zorder=1)
        ax_side.plot([t_p, t_p + math.sqrt(2.0)*a_f], [bot_p + h_b - t_fb - math.sqrt(2.0)*a_f, bot_p + h_b -t_fb], 'k', zorder=1)
        ax_side.plot([t_p, t_p + math.sqrt(2.0)*a_f], [bot_p - math.sqrt(2.0)*a_f, bot_p], 'k', zorder=1)
        ax_side.plot([t_p, t_p + math.sqrt(2.0)*a_f], [bot_p + t_fb + math.sqrt(2.0)*a_f, bot_p + t_fb], 'k', zorder=1)
        
        z0 = self.ebottom + 0.5*h_b
        
        for row in self.bolt_rows:
            h_bolt = row.bolt.head_t
            d_m = row.bolt.head_d
            t_washer = row.bolt.washer_t
            h_nut = row.bolt.nut_t
            L_b = row.bolt.length
            z = row.z + z0
            d = row.bolt.d0
            # Plot bolt hole
            ax_side.plot([-t_fc, t_p], [z + 0.5*d, z + 0.5*d], 'k--', linewidth=0.6)
            ax_side.plot([-t_fc, t_p], [z - 0.5*d, z - 0.5*d], 'k--', linewidth=0.6)
        
            # Plot bolt head
            bolt_head = patches.Rectangle((t_p, z - 0.5*d_m), width=h_bolt, height=d_m,edgecolor='k', facecolor="w", zorder=2)
            ax_side.add_patch(bolt_head)
            
            nut_head = patches.Rectangle((-t_fc - t_washer - h_nut, z - 0.5*d_m), width=h_nut, height=d_m, edgecolor='k', facecolor="w", zorder=2)
            ax_side.add_patch(nut_head)
        
            # Plot center line
            ax_side.plot([-0.5 * L_b, 0.5 * L_b], [z, z], 'b-.', linewidth=0.6, zorder=3)
        
        ax_side.set_aspect('equal', 'datalim')
        
        
        ax_side.set_axis_off()
        
        """ Draw front view """
        ax_front = ax[1]
        
        #axis_lim = max(b_c + 200.0, h_b + 200.0, h_p + 200.0)

        # Draw column
        ax_front.vlines([-0.5*b_c, -0.5*b_c,0.5*b_c, 0.5*b_c], -100.0, axis_lim - 100.0, 'k', linewidth=1.0)

        # Draw plate
        plate_front = patches.Rectangle((-0.5*b_p, 0.0), width=b_p, height=h_p,\
                                        linewidth=1.0, edgecolor='k', facecolor="w", zorder=2)
        ax_front.add_patch(plate_front)
        
        # Draw beam profile
        self.beam.draw(axes=ax_front,origin=[0,self.ebottom+0.5*h_b])
        
        # Draw bolt rows
        for row in self.bolt_rows:
            w = 0.5*row.p
            
            h_bolt = row.bolt.head_t
            d_m = row.bolt.head_d
            d_washer = row.bolt.washer_d
            h_nut = row.bolt.nut_t
            L_b = row.bolt.length
            z = row.z + z0
            d = row.bolt.d0
            
            # Plot washers
            ax_front.add_patch(patches.Circle((-w, z), 0.5*d_washer, edgecolor='k', facecolor="w", zorder=4))
            ax_front.add_patch(patches.Circle((w, z), 0.5*d_washer, edgecolor='k', facecolor="w", zorder=4))
            # Plot bolts            
            ax_front.add_patch(patches.Circle((-w, z), 0.5*d_m, edgecolor='k', facecolor="w", zorder=4))
            ax_front.add_patch(patches.Circle((w, z), 0.5*d_m, edgecolor='k', facecolor="w", zorder=4))
            ax_front.add_patch(patches.Circle((-w, z), 0.5*d, linestyle='--', linewidth=0.6, edgecolor='k',
                                                    facecolor="w", zorder=4))
            ax_front.add_patch(patches.Circle((w, z), 0.5*d, linestyle='--', linewidth=0.6, edgecolor='k',
                                                    facecolor="w", zorder=4))
            
            ax_front.vlines([-w,w],z - 0.6 * d_m, z + 0.6 * d_m,'b', linewidth=0.6, zorder=5)
            ax_front.hlines(z,-w - 0.6 * d_m, -w + 0.6 * d_m,'b', linewidth=0.6, zorder=5)
            ax_front.hlines(z,w - 0.6 * d_m, w + 0.6 * d_m,'b', linewidth=0.6, zorder=5)
        
            #ax_front.plot([-w, -w], [z - 0.6 * d_m, z + 0.6 * d_m], 'b-.', linewidth=0.6, zorder=5)
            #ax_front.plot([-0.5 * w - 0.6 * d_m, -0.5 * w + 0.6 * d_m], [z, z], 'b-.', linewidth=0.6, zorder=5)
            #ax_front.plot([0.5 * w, 0.5 * w], [z - 0.6 * d_m, z + 0.6 * d_m], 'b-.', linewidth=0.6, zorder=5)
            #ax_front.plot([0.5 * w - 0.6 * d_m, 0.5 * w + 0.6 * d_m], [z, z], 'b-.', linewidth=0.6, zorder=5)
        
        ax_front.set_aspect('equal', 'datalim')
        
        ax_front.set_axis_off()
        
def example_1():
    bolt = Bolt(24,10.9)
    bolt.washer_t = 0.0
    col = HEA(340,fy=235)
    beam = IPE(500,fy=235)
    
    hp = 50+85+80+320+65
    bp = 240
    ebottom = hp-85-beam.h
    
    y0 = ebottom + 0.5*beam.h
    y_bolts = []
    y_bolts.append(hp-50-y0)
    y_bolts.append(hp-50-85-y0)
    y_bolts.append(hp-50-85-80-y0)
    y_bolts.append(hp-50-85-80-320-y0)
    
    bolt_row_pos = [{"flange": END_ROW, "plate": ROW_OUTSIDE_BEAM_TENSION_FLANGE},
                    {"flange": INNER_ROW, "plate": FIRST_ROW_BELOW_BEAM_TENSION_FLANGE},
                    {"flange": END_ROW, "plate": OTHER_END_ROW},
                    {"flange": END_ROW, "plate": OTHER_END_ROW},
                    ]
    bolt_row_groups = [[0,1],[0,1,2],[1,2]]
    
    row_types = [TENSION_ROW,TENSION_ROW,TENSION_ROW,SHEAR_ROW]
    
    # POSITION OF ROWS IN GROUPS
    # For each row, position with respect to column flange and end plate are
    # required.
    group_pos = [[{"flange":END_ROW, "plate": ROW_OUTSIDE_BEAM_TENSION_FLANGE},
                  {"flange":END_ROW, "plate": FIRST_ROW_BELOW_BEAM_TENSION_FLANGE}
                 ],
                 [{"flange":END_ROW, "plate": ROW_OUTSIDE_BEAM_TENSION_FLANGE},
                  {"flange":INNER_ROW, "plate": FIRST_ROW_BELOW_BEAM_TENSION_FLANGE},
                  {"flange":END_ROW, "plate": OTHER_END_ROW}
                 ],
                 [{"flange":END_ROW, "plate": FIRST_ROW_BELOW_BEAM_TENSION_FLANGE},
                  {"flange":END_ROW, "plate": OTHER_END_ROW}
                 ]
                 ]
    conn = EndPlateJoint(col,beam,tp=15,bp=bp,mat_p="S235",etop=85,\
                         ebottom=ebottom,bolt=bolt,y_bolts=y_bolts,\
                         e_bolts=60,\
                         bolt_row_pos=bolt_row_pos,
                         groups=bolt_row_groups,
                         group_pos=group_pos,
                         row_types=row_types)
        
    conn.weld_f = 8
    
    return conn

def example_diaz():
    d = 12
    bolt = Bolt(d,8.8)
    col = HEB(140,fy=355)
    beam = IPE(200,fy=355)
    
    aw = 5.0
    af = 3.0
    
    bp = col.b
    #tp = 12
    tp = 7.15
    ebottom = math.floor(math.sqrt(2)*af + tp)
    
    e = 28.8
    #ex = 22.52
    #px = 79.32
    ex = 37.61
    px = 58.01
    
    hp = beam.h + 4*d + ebottom
    
    """ Generate bolt rows
        y is the vertical position of a bolt row. The
        origin is located in the centroid of the beam
    """
    y0 = ebottom + 0.5*beam.h
    y_bolts = []
    #y_bolts.append(hp-ex-y0)
    #y_bolts.append(hp-ex-px-y0)
    y_bolts.append(0.5*beam.h + 4*d-ex)
    y_bolts.append(y_bolts[0]-px)
    
    bolt_row_pos = [{"flange": END_ROW, "plate": ROW_OUTSIDE_BEAM_TENSION_FLANGE},
                    {"flange": INNER_ROW, "plate": FIRST_ROW_BELOW_BEAM_TENSION_FLANGE},
                    ]
    bolt_row_groups = [[0,1]]
    
    row_types = [TENSION_ROW,TENSION_ROW]
    
    # POSITION OF ROWS IN GROUPS
    # For each row, position with respect to column flange and end plate are
    # required.
    group_pos = [[{"flange":END_ROW, "plate": ROW_OUTSIDE_BEAM_TENSION_FLANGE},
                  {"flange":END_ROW, "plate": FIRST_ROW_BELOW_BEAM_TENSION_FLANGE}
                 ]]
    conn = EndPlateJoint(col,beam,tp=tp,bp=bp,mat_p="S355",etop=4*d,\
                         ebottom=ebottom,bolt=bolt,y_bolts=y_bolts,\
                         e_bolts=e,\
                         bolt_row_pos=bolt_row_pos,
                         groups=bolt_row_groups,
                         group_pos=group_pos,
                         row_types=row_types)
    
    conn.weld_f = af
    conn.weld_w = aw
    
    return conn

if __name__ == '__main__':
    
    from sections.steel.ISection import HEA, HEB, IPE
    from eurocodes.en1993.en1993_1_8.en1993_1_8 import Bolt
    
    conn = example_diaz()
    
    ws = Workshop()
    
    conn.info(draw=False)
    
    conn.MjRd(True)
    conn.Sj_ini(True)
    conn.cost(ws,verb=True)
    """
    for row in conn.bolt_rows:
        print(row.stiffness_factors)
    """
    """
        for group in row.groups:
            group.Tstub_col_flange.leff_row(row)
            group.Tstub_end_plate.leff_row(row)
    """
    """
    for group in conn.row_groups:
        print(group.Tstub_col_flange._leff_rows)
        print(group.Tstub_end_plate._leff_rows)
    """
    
    #conn.row_groups[0].Tstub_col_flange.leff_nc_row(r)
    #conn.row_groups[0].Tstub_col_flange.leff_cp_row(r)
        
        
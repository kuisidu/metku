# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 19:41:00 2019

@author: kmela
"""


import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines

from plates import RectPlate
from eurocodes.en1993.en1993_1_8.en1993_1_8 import BoltRow
import eurocodes.en1993.en1993_1_8.component_method as cm

class EndPlateJoint:
    """ Class for end-plate connections """
    
    
    
    def __init__(self,column,beam,tp,bp,mat_p,etop,ebottom,bolt,y_bolts,e_bolts):
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
                e_bolts .. edge distance of bolts (from centroid) [mm[]]
        """
        
        self.col = column
        self.beam = beam
        
        # calculate plate height
        hp = etop + ebottom + beam.h        

        
        self.end_plate = RectPlate(bp,hp,tp,material="S355")
        self.etop = etop
        self.ebottom = ebottom
        self.ebolts = e_bolts        
        # Make bolt rows
        self.bolt_rows = []
        
        # Beam flange weld, throat thickness [mm] (to plate)
        self.weld_f = 5
        # Beam web weld, throat thickness [mm] (to plate)
        self.weld_w = 5
        
        # Plate to column weld, throat thickness [mm]
        self.weld_p = 5
        
        self.p = bp-2*e_bolts
        print(bp-2*e_bolts)
        for y in y_bolts:
            self.bolt_rows.append(BoltRow(bolt,self.p,y))
            
    @property
    def beta(self):
        """ Effect of shear in the column web to moment resistance of the joint 
            So far, only one value is given (approximation for one-sided joint)
        """
        return 1.0
        
        
    def V_wp_Rd(self):
        """ Column web shear resistance """
        return cm.col_web_shear(self.col)
    
    def Fc_wc_Rd(self):
        """ Column web in compression """
        
        # Calculate compression stre
        sigma_com_Ed = self.col.sigma_com()
        Fc_wc_Rd, beff_wc_Rd = cm.col_web_trv_comp(self.col, self.beam, self.end_plate.t, self.ebottom, self.weld_p, self.beta, sigma_com_Ed)
        return Fc_wc_Rd
    
    def Fc_fb_Rd(self):
        """ Beam flange and web in compression """        
        return cm.beam_web_compression(self.beam)
    
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
        beam.draw(axes=ax_front,origin=[0,ebottom+0.5*h_b])
        
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
        """
        
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
        
         Draw connection from above 
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
        """

if __name__ == '__main__':
    
    from ISection import HEA, IPE
    from eurocodes.en1993.en1993_1_8.en1993_1_8 import Bolt
    
    bolt = Bolt(24,10.9)
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
    
    conn = EndPlateJoint(col,beam,tp=15,bp=bp,mat_p="S235",etop=85,\
                         ebottom=ebottom,bolt=bolt,y_bolts=y_bolts,\
                         e_bolts=60)
    
    conn.info(draw=False)
        
        
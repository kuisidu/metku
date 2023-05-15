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
from metku.eurocodes.en1993.en1993_1_8.en1993_1_8 import TENSION_ROW, COMPRESSION_ROW, SHEAR_ROW, COMBI_ROW
from metku.eurocodes.en1993.en1993_1_8.en1993_1_8 import TStubEndPlate, BoltRow

ROW_OUTSIDE_COLUMN_FLANGE = "Row outside column flange"
ROW_INSIDE_COLUMN_FLANGE = "Row inside column flange"


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
        
class AnchorBoltRow:
    """ Class for anchor bolt rows 
        
        Methods and functionality:
            - compute tension resistance of each component
            - compute stiffness of each component
        
        The components are:
            - end plate in bending
            - column flange in bending
            - column web in tension
            - beam web in tension
    
    """
    
    def __init__(self,bolt,p,z,loc=ROW_OUTSIDE_COLUMN_FLANGE,joint=None, row_type=TENSION_ROW):
        """ Constructor
            input:
                bolt .. Bolt class object defining the bolt
                p .. array of distances between adjacent bolt centroids
                z .. horizontal position of bolt row (origin is at the centroid of the column section)
                loc .. position of the row with respect to column flange                
                joint .. the joint object to which the bolt row belongs
                row_type .. TENSION_ROW (default)
        """
        
        self.type = row_type
        self.bolt = bolt
        #self.p = p
        self.z = z
        self.loc = loc        
        
        """ Location of bolt row in a group. This varies,
            if the row belongs to several groups        
        """
        self.group_loc_flange = None
        self.group_loc_plate = None
        
        self.joint = joint
        
        # Edge distance
        self._e1 = 0.5*self.joint.column.h+self.joint.etop - z
        #self.y = y
        """ p can be a single number, in which case there are two bolts
            in the row
            
            if p is an array, the number of bolts is length of p + 1
        """
        if isinstance(p,float) or isinstance(p,int):
            self.bolts = 2
        else:
            self.bolts = len(p)+1
        
        """ Initialize effective lengths corresponding to
            column flange in bending (leff_flange) and
            end plate in bending (leff_plate), when
            the row is treated as an individual row
            
            first element is for Mode 1 and second for Mode 2
        """
        #self.Tstub_col_flange = TStubColumnFlange([bolt,bolt],flange_loc,self.joint.col.material,tf=10,e=20,emin=20,m=10)
        self.Tstub_col_flange = TStubColumnFlange(self,self.joint.col.material,tf=10,e=20,emin=20,m=10)
        self.Tstub_end_plate = TStubEndPlate(self,self.joint.end_plate.material,tf=10,emin=20,m=10,e=20,w=20,bp=self.joint.bp)
        #self.leff_flange = [0.0,0.0]
        #self.leff_plate = [0.0,0.0]
        
        
        
        self.FtRd = -1.0
    
        # Dictionary of stiffness factors
        self.stiffness_factors = {"col_web": np.inf,
                                  "col_flange": np.inf,
                                  "plate":np.inf,
                                  "bolt":1.6*self.bolt.As/self.Lb}
        
        # List of row groups that contain this row
        # this list appended when bolt row groups are created.
        self.groups = []
    
        # INITIALIZATION: CALCULATE VARIOUS GEOMETRIC QUANTITIES
        
        """ Effective length of column flange """
        #self._leff_nc_flange = self.Tstub_col_flange.leff_nc()
        #self._leff_cp_flange = self.Tstub_col_flange.leff_cp()
        
        """ Effective length of end plate """
        #self._leff_nc_flange = self.Tstub_end_plate.leff_nc()
        #self._leff_cp_flange = self.Tstub_end_plate.leff_cp()
    
    def __repr__(self):
        
        return "Bolt row"
    
    
    def info(self):
        """ Prints info """
        
        print("Bolt row:")
        print("Location (flange): {0:s}".format(self.loc_col_flange))
        print("Location (end plate): {0:s}".format(self.loc_end_plate))
        
    @property
    def bolt_size(self):
        return self.bolt.d
    
    @bolt_size.setter
    def bolt_size(self,val):
        self.bolt.size = val
        
    @property
    def z_top(self):
        """ Distance of the bolt row from the top of the end plate """
        return self._z_top
    
    @z_top.setter
    def z_top(self,val):
        """ Distance of the bolt row from the top of the end plate """
        self._z_top = val
        self.z = 0.5*self.joint.beam.h+self.joint.etop - val
        
    @property
    def Lb(self):
        """ Ruuvirivin venymäpituus """
        return self.joint.tp + self.joint.col.tf + 2*self.bolt.washer_t + 0.5*(self.bolt.head_t + self.bolt.nut_t)
    
    @property
    def h(self):
       """ Distance from the center of compression of the joint """
       return self.z-self.joint.zc 
    
    @property
    def emin_flange(self):
        """ Distance 'emin' used for T-stubs for the component
            'column flange in bending'
        """
        """ If the plate is narrower than column flange,
            take 'emin' to be the distance to plate edge.
                
            Otherwise, use the distance to column flange edge.
        """
        if self.joint.bp < self.joint.col.b:
            emin = self.joint.ebolts
        else:
            emin = self.joint.ebolts - 0.5*(self.joint.bp-self.joint.col.b)
        
        return emin
    
    @property
    def emin_plate(self):
        """ Distance 'emin' used for T-stubs for the component
            'end plate in bending'
        """
        if self.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            """ For the row in the plate extension, use the
                distance to plate end, ex.
            """
            emin = self.ex
        else:
            """ If the plate is narrower than column flange,
                take 'emin' to be the distance to plate edge.
                
                Otheriwse, use the distance to column flange edge.
            """
            if self.joint.bp < self.joint.col.b:
                emin = self.joint.ebolts
            else:
                emin = self.joint.ebolts - 0.5*(self.joint.bp-self.joint.col.b)
        
        return emin
    
    @property
    def ex(self):
        """ Distance from bolt row center to top edge of the end plate
        """
        
        return 0.5*self.joint.beam.h + self.joint.etop - self.z
    
    @property
    def mx(self):
        """ Distance from bolt row center to top of the top flange
            of the beam. This is used for bolt rows outside the
            beam tension flange.
        """
        
        return self.z - 0.5*self.joint.beam.h - 0.8*math.sqrt(2)*self.joint.weld_f
    
    @property
    def w(self):
        """ Horizontal distance between bolt holes """
        return self.joint.p        
        
    @property
    def col_flange_m(self):
        """ Distance from column web to bolt hole center """
        
        return 0.5*(self.w-self.joint.col.tw)-0.8*self.joint.col.r
    
    @property
    def col_flange_e(self):
        """ Distance of column flange edge from the bolt hole center """
        if self.joint.bp < self.joint.col.b:
            e = self.joint.ebolts + 0.5*(self.joint.col.b-self.joint.bp)
        else:
            e = self.joint.ebolts - 0.5*(self.joint.bp-self.joint.col.b)
        
        return e
    
    @property
    def end_plate_m(self):
        """ Distance from beam web to bolt hole center """
        if self.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            m = self.mx
        else:
            m = 0.5*(self.w-self.joint.beam.tw)-0.8*math.sqrt(2)*self.joint.weld_w 
        return m
    
    @property
    def end_plate_m2(self):
        """ Distance from beam flange to bolt hole center """
        
        return 0.5*self.joint.beam.h-self.joint.beam.tf-0.8*math.sqrt(2)*self.joint.weld_f-self.z
    
    @property
    def end_plate_e(self):
        """ Distance of end plate edge from the bolt hole center """
        
        return self.joint.ebolts
    
    @property
    def top_flange_stiffener_m2(self):
        """ Distance of bolt row from top flange stiffener """
        zbeam = 0.5*self.joint.beam.h
        if self.z > zbeam:
            """ Row is above the stiffener """
            zweld = zbeam - 0.5 * self.joint.beam.tf + \
                    + 0.5 * self.joint.stiffeners['col_top_flange'].khp + \
                    + 0.8 * math.sqrt(2) * self.joint.stiffeners['col_top_flange'].weld_size
            return self.z-zweld
        else:
            """ Row is below the stiffener """
            zweld = zbeam - 0.5 * self.joint.beam.tf + \
                    - 0.5 * self.joint.stiffeners['col_top_flange'].khp + \
                    - 0.8 * math.sqrt(2) * self.joint.stiffeners['col_top_flange'].weld_size
            
            return zweld-self.z
    
    def leff_flange(self):
        """ Effective length of column flange in bending component """
        return self.Tstub_col_flange.leff()
    
    def leff_plate(self):
        """ Effective length of end plate in bending component """
        return self.Tstub_end_plate.leff()
    
    def leff_flange_group(self):
        """ Effective length of column flange in bending for the row
            when considered to be in a group
            
            IDEA:
                for all groups that have the row:
                    1) get row location in that group
                    2) calculate effective width for the row
                        considered as in a group and with the
                        obtained location.
                    3) add the calculated effective width to the
                        list of effective widths
                    4) return minimum
    
            QUESTION: how to know, which groups have the row
                in them.
        """
        leff = 1e5
        """ The for loop is over all the groups that involve this bolt row
            (the groups are stored in 'groups' attribute)
        """
        for group in self.groups:
            """ The '_leff_rows' attribute of the column flange T-stub
                for the group is a dict with BoltRow objects as keys. Using
                'self' will get the effective lengths of this bolt row
                appearing in the given group. Then, take the minimum of those values.
            """
            #leff_g = min(group.Tstub_col_flange._leff_rows[self].values())
            leff_g = min(group.Tstub_col_flange.leff_cp_row(self),group.Tstub_col_flange.leff_nc_row(self))
            if leff_g > 0.0:
                """ If the group does not have column flange in bending
                    component, then leff_g = 0 and that group can be skipped.
                    Otherwise, update the effective length.
                """
                leff = min(leff,leff_g)        
            
        return leff
        
    
    def leff_plate_group(self):
        """ Effective length of end plate in bending component
            for the row when considered to be in a group
        """
        leff = 1e5
        for group in self.groups:
            #leff_g = min(group.Tstub_end_plate._leff_rows[self].values())
            
            # If the group does not have the end plate in bending component,
            # return a large value.
            
            if not group.Tstub_end_plate is None:
                leff_g = min(group.Tstub_end_plate.leff_cp_row(self),group.Tstub_end_plate.leff_nc_row(self))                
            else:
                leff_g = leff
            if leff_g > 0.0:
                leff = min(leff,leff_g)
            
        return leff
    
    def k3(self,verb=False):
        """ Stiffness factor for column web in tension component """
        
        #beff_row = self.leff_flange()
        # NOTE! Here, the _leff attribute of the T-stub corresponds to the 
        # effective lengths calculated during bending resistance evaluation.
        # This can be a problem if the bending resistance has not been
        # calculated before rotational stiffness, or when calculating derivatives
        # in optimization using finite differences.
        #
        # beff_row = min(self.Tstub_col_flange._leff)
        
        #self.Tstub_col_flange.tf = self.joint.col.tf
        #self.Tstub_col_flange.emin = self.emin_flange
        #self.Tstub_col_flange.m = self.col_flange_m
        #self.Tstub_col_flange.e = self.col_flange_e
        
        beff_row = min(self.Tstub_col_flange.leff_cp,self.Tstub_col_flange.leff_nc)
        
        #beff_group = 1e5
        beff_group = self.leff_flange_group()
        
        beff = min(beff_row,beff_group)
        
        if verb:
            print(f'beff_row = {beff_row:4.2f} mm')
            print(f'beff_group = {beff_group:4.2f} mm')
        
        return cm.column_web_tension_k(self.joint.col,beff,verb)
    
    
    def k4(self,verb=False):
        """ Stiffness factor for column flange in bending component """
        
        """ This way picks that effective length that corresponds to the
            governing failure mode.
        beff_row = self.leff_flange()
        """
        
        """ This way takes the smallest effective length for the row """
        #beff_row = min(self.Tstub_col_flange._leff)
        #self.Tstub_col_flange.tf = self.joint.col.tf
        #self.Tstub_col_flange.emin = self.emin_flange
        #self.Tstub_col_flange.m = self.col_flange_m
        #self.Tstub_col_flange.e = self.col_flange_e
        
        beff_row = min(self.Tstub_col_flange.leff_cp,self.Tstub_col_flange.leff_nc)
        #beff_group = 1e5
        beff_group = self.leff_flange_group()
        leff = min(beff_row,beff_group)        
        
        
        m = self.Tstub_col_flange.m
        return cm.column_flange_bending_k(self.joint.col,leff,m,verb)
        
    def k5(self,verb=False):
        """ Stiffness factor for end plate in bending component """
        
        #beff_row = self.leff_plate()
        #beff_row = min(self.Tstub_end_plate._leff)
        #self.Tstub_end_plate.m2 = self.end_plate_m2
        
        #self.Tstub_end_plate.bp = self.joint.bp
        #self.Tstub_end_plate.tf = self.joint.tp
        #self.Tstub_end_plate.w = self.p
        
        #if self.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
        #    self.Tstub_end_plate.emin = self.ex
        #    self.Tstub_end_plate.e = self.end_plate_e
        #    self.Tstub_end_plate.m = self.mx
        #else:
        #    self.Tstub_end_plate.emin = self.emin_plate
        #    self.Tstub_end_plate.e = self.end_plate_e
        #    self.Tstub_end_plate.m = self.end_plate_m
        #    self.Tstub_end_plate.m2 = self.end_plate_m2
        
        beff_row = min(self.Tstub_end_plate.leff_cp,self.Tstub_end_plate.leff_nc)
        
        
        #beff_group = 1e5
        beff_group = self.leff_plate_group()
        
        #print(f'beff (row) = {beff_row:4.2f}')
        #print(f'beff (group) = {beff_group:4.2f}')
        
        leff = min(beff_row,beff_group)                
        
        m = self.Tstub_end_plate.m
        
        #print(beff_row,beff_group,m)
        return cm.end_plate_bending_k(self.joint.tp,leff,m,verb)
    
    def k10(self,verb=False):
        """ Stiffness factor for bolts in tension """
        
        k10 = 1.6*self.bolt.As/self.Lb
        
        return k10
        
    
    def column_web_in_tension(self,verb=False):
        """ Resistance of component """
        #self.Tstub_col_flange.tf = self.joint.col.tf
        #self.Tstub_col_flange.emin = self.emin_flange
        #self.Tstub_col_flange.m = self.col_flange_m
        #self.Tstub_col_flange.e = self.col_flange_e
        
        """ The effective length of the component
            corresponds to the effective length of the column flange.
            This is obtained only after the resistance of the column
            flange has been determined, because the governing failure
            mode also defines the effective width
        """
        #l_eff = self.Tstub_col_flange.leff(self.loc_col_flange)
        
        #l_eff = [0,0]
        #l_eff[0] = self.Tstub_col_flange.leff_1()
        #l_eff[1] = self.Tstub_col_flange.leff_2()
        
        l_eff = self.Tstub_col_flange.leff()
        #print(self.Tstub_col_flange.leff)
        
        F_t_wc_Rd = cm.column_web_tension(self.joint.col, l_eff, self.joint.beta,verb)
        
        if verb:
            print("COlUMN WEB IN TENSION:")
            print("l_eff  = {0:4.2f} mm".format(l_eff))
            #print("l_eff (Mode 2) = {0:4.2f} mm".format(self.Tstub_col_flange._leff[1]))
            print("Ft_wc_Rd = {0:4.2f} kN".format(F_t_wc_Rd*1e-3))
        
        return F_t_wc_Rd
    
    def column_flange_in_bending(self,verb=False):
        """ Resistance of component """
        
        """ Set properties for the T-stub
            This is basically for optimization purposes,
            because the dimensions will change
        """
        #self.Tstub_col_flange.tf = self.joint.col.tf
        #self.Tstub_col_flange.emin = self.emin_flange
        #self.Tstub_col_flange.m = self.col_flange_m
        #self.Tstub_col_flange.e = self.col_flange_e
        if self.joint.stiffeners['col_top_flange'] is not None:
            self.Tstub_col_flange.m2 = self.top_flange_stiffener_m2
        """ Determine effective length 
            The 'leff_1' and 'leff_2' methods store the
            effective length of each failure mode to the T-stub,
            so it can be
        """
        
        """ Calculate the effective length of the T-stub
            The effective length will be stored in the T-stub object.
        """
        l_eff = [0.0,0.0]
        l_eff[0] = self.Tstub_col_flange.leff_1
        l_eff[1] = self.Tstub_col_flange.leff_2
        
        if verb:
            print("COLUMN FLANGE IN BENDING:")
            print("l_eff (Mode 1) = {0:4.2f} mm".format(l_eff[0]))
            print("l_eff (Mode 2) = {0:4.2f} mm".format(l_eff[1]))
        
        """ The resistance is based on the T-stub """
        F_t_cf_Rd = cm.column_flange_bending(self.Tstub_col_flange,verb)
        
        if verb:            
            print("Ft_cf_Rd = {0:4.2f} kN".format(F_t_cf_Rd*1e-3))
        
        return F_t_cf_Rd
    
    def end_plate_in_bending(self,verb=False):
        """ Resistance of component """
        
        """
        self.Tstub_end_plate.bp = self.joint.bp
        self.Tstub_end_plate.tf = self.joint.tp
        self.Tstub_end_plate.w = self.p
        
        if self.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            self.Tstub_end_plate.emin = self.ex
            self.Tstub_end_plate.e = self.end_plate_e
            self.Tstub_end_plate.m = self.mx
        else:
            self.Tstub_end_plate.emin = self.emin_plate
            self.Tstub_end_plate.e = self.end_plate_e
            self.Tstub_end_plate.m = self.end_plate_m
            self.Tstub_end_plate.m2 = self.end_plate_m2
        """
        
        l_eff = [0.0,0.0]
        l_eff[0] = self.Tstub_end_plate.leff_1
        l_eff[1] = self.Tstub_end_plate.leff_2
        
        if verb:
            print("END PLATE IN BENDING:")
            self.Tstub_end_plate.info()
            print("l_eff (Mode 1) = {0:4.2f} mm".format(l_eff[0]))
            print("l_eff (Mode 2) = {0:4.2f} mm".format(l_eff[1]))
        
        F_t_ep_Rd = cm.end_plate_bending(self.Tstub_end_plate,verb)
        
        return F_t_ep_Rd
    
    def beam_web_in_tension(self,verb=False):
        """ Resistance of component """
        """
        self.Tstub_end_plate.tf = self.joint.tp
        self.Tstub_end_plate.w = self.p
        
        if self.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            self.Tstub_end_plate.emin = self.ex
            self.Tstub_end_plate.e = self.ex
            self.Tstub_end_plate.m = self.mx
        else:
            self.Tstub_end_plate.emin = self.emin_plate
            self.Tstub_end_plate.e = self.end_plate_e
            self.Tstub_end_plate.m = self.end_plate_m
            self.Tstub_end_plate.m2 = self.end_plate_m2
        """
        
        if self.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            F_t_bw_Rd = 1e8
        else:        
            l_eff = self.Tstub_end_plate.leff()
            #l_eff[0] = self.Tstub_end_plate.leff_1(self.loc_end_plate)
            #l_eff[1] = self.Tstub_end_plate.leff_2(self.loc_end_plate)
        
            F_t_bw_Rd = cm.beam_web_tension(self.joint.beam, l_eff)
        
        if verb:
            print("BEAM WEB IN TENSION:")
            if self.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
                print("Row outside beam tension flange: no beam web in tension.")
            else:
                print("l_eff = {0:4.2f} mm".format(l_eff))
                print("Ft_bf_Rd = {0:4.2f} kN".format(F_t_bw_Rd*1e-3))
        
        return F_t_bw_Rd
    
    def tension_resistance(self,verb=False):
        """ Tension resistance of individual row, not considered
            as a group.
            
            Evaluate first the components
            'column flange in bending' and
            'end plate in bending', because
            they define the effective lengths also used in
            'column web in tension' and
            'beam web in tension' components
        """
        
        Ft_i_Rd = np.array([0.0,0.0,0.0,0.0])
        components = ['Column flange in bending',
                      'End plate in bending',
                      'Column web in tension',
                      'Beam web in tension'
                      ]
        
        Ft_i_Rd[0] = self.column_flange_in_bending(verb)
        Ft_i_Rd[1] = self.end_plate_in_bending(verb)
        Ft_i_Rd[2] = self.column_web_in_tension(verb)
        Ft_i_Rd[3] = self.beam_web_in_tension(verb)
            
        min_comp = np.argmin(Ft_i_Rd)
        
        self.FtRd = Ft_i_Rd[min_comp]
        
        if verb:
            print("Ft_Rd = {0:4.2f} kN ({1:s})".format(self.FtRd*1e-3,components[min_comp]))
        
        return self.FtRd
    
    def keff(self):
        """ Effective stiffness of bolt row 
            EN 1993-1-8, Eq. (6.30)
            
            Requires that the stiffness factors have been evaluated
            and stored in 'stffiness_factors' dict.
        """
        
        # Inverse of keff
        keff_inv = 0.0
        
        self.stiffness_factors['col_web'] = self.k3()
        self.stiffness_factors['col_flange'] = self.k4()
        self.stiffness_factors['plate'] = self.k5()
        self.stiffness_factors['bolt'] = self.k10()
        
        #print(self.stiffness_factors)
        
        for k in self.stiffness_factors.values():
            keff_inv += 1/k
        
        return 1/keff_inv

class ColumnBase:
    
    
    def __init__(self,column_profile, plate_edges, plate_thickness, anchor_bolts, bolt_z=None, bolt_w=100, concrete="C25/30"):
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
        self._c = plate_edges
        h = column_profile.h
        b = column_profile.b
        hp = h + self._c[0] + self._c[1]
        bp = b + self._c[2] + self._c[3]
        base_plate = RectPlate(bp,hp,plate_thickness,material="S355")
        self.plate = base_plate
        self.anchor_bolt = anchor_bolts
        self.w = bolt_w
        
        # By default, place anchor bolt rows halfway at outstand part of the base plate
        if bolt_z is None:
            bolt_z = [-0.5*h-0.5*self._c[0], 0.5*h+0.5*self._c[1]]
        
        self.zb = bolt_z
        self.concrete = en1992_const.Concrete(concrete)        
        self.Ac0 = 1.0
        self.Ac1 = 1.0
        self.flange_weld = 5.0
        self.web_weld = 5.0
        self.bolt_rows = []
        self.add_bolt_row(anchor_bolts, bolt_z[0], self.w)
        self.add_bolt_row(anchor_bolts, bolt_z[1], self.w)
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
            loc = ROW_OUTSIDE_COLUMN_FLANGE
            row_type = TENSION_ROW
        elif z < 0:
            loc = ROW_INSIDE_COLUMN_FLANGE
            row_type = TENSION_ROW
        else:
            loc = ROW_OUTSIDE_COLUMN_FLANGE
            row_type = COMPRESSION_ROW
        
        new_row = AnchorBoltRow(bolt,p,z,loc=loc,joint=self,row_type=row_type)
        
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
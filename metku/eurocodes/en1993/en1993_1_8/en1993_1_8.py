# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-8 Rules for connections

@author: kmela
"""

import math
import numpy as np

from materials.steel_data import Steel
from eurocodes.en1993.constants import gammaM0, gammaM2, gammaM5
from cost.cost_data import bolt_unit_cost_part_thread, bolt_unit_cost_full_thread
import eurocodes.en1993.en1993_1_8.component_method as cm

# Bolt materials [MPa], SFS EN 1993-1-8, table 3.1
mat_bolt = {4.6: {"f_yb": 240.0, "f_ub": 400.0},
            4.8: {"f_yb": 320.0, "f_ub": 400.0},
            5.6: {"f_yb": 300.0, "f_ub": 500.0},
            5.8: {"f_yb": 400.0, "f_ub": 500.0},
            6.8: {"f_yb": 480.0, "f_ub": 600.0},
            8.8: {"f_yb": 640.0, "f_ub": 800.0},
            10.9: {"f_yb": 900.0, "f_ub": 1000.0}}

# Standard bolt sizes
bolt_sizes = [12, 16, 20, 24, 30, 36]

# Standard bolt lengths
bolt_lengths = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 90.0, 100.0, 110.0,
                120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0]

""" Bolt sizes, stress area [mm^2]
    As .. stress area
    k .. head thickness
    s .. avainväli
    e .. ulkomitta

"""
bolt_size = {12: {"A_s": 84.3, "k": 7.5, "s": 18, "e": 20.03},
             16: {"A_s": 157.0, "k": 10.0, "s": 24, "e": 26.75},
             20: {"A_s": 245.0, "k": 12.5, "s": 30, "e": 33.53},
             24: {"A_s": 353.0, "k": 15, "s": 36, "e": 39.98},
             30: {"A_s": 561.0, "k": 18.7, "s": 46, "e": 50.85},
             36: {"A_s": 817.0, "k": 22.5, "s": 55, "e": 60.79},
             }

""" EN ISO 4032
    m .. nut thickness
    s .. avainväli
    e .. ulkohalkaisija
"""
nut_size = {12: {"m": 10.8, "s": 18, "e": 20.03},
             16: {"m": 14.8, "s": 24, "e": 26.75},
             20: {"m": 18.0, "s": 30, "e": 32.95},
             24: {"m": 21.5, "s": 36, "e": 39.55},
             30: {"m": 25.6, "s": 46, "e": 50.85},
             36: {"m": 31.0, "s": 55, "e": 60.79},
             }

""" EN ISO 7089
    h .. thickness
    d2 .. outside diameter
    d1 .. clearance hole
"""
washer_size = {12: {"h": 2.5, "d1": 13, "d2": 24.0},
             16: {"h": 3.0, "d1": 17, "d2": 30.0},
             20: {"h": 3.0, "d1": 21, "d2": 37.0},
             24: {"h": 4.0, "d1": 25, "d2": 44.0},
             30: {"h": 4.0, "d1": 31, "d2": 56.0},
             36: {"h": 5.0, "d1": 37, "d2": 66.0},
             }

correlation_coefficient = {"S235": 0.8,
                           "S275": 0.85,
                           "S355": 0.9,
                           "S420": 1.0,
                           "S460": 1.0}

""" Bolt row types """
TENSION_ROW = 'tension row'
SHEAR_ROW = 'shear row'
COMBI_ROW = 'tension and shear row'
COMPRESSION_ROW = 'compression row'

""" Bolt row positions 

    1. Unstiffened column flange
"""
INNER_ROW = 'Inner Row'
END_ROW = 'End Row'
""" 2. stiffened column flange """
ADJACENT_TO_STIFFENER_ROW = "Adjacent to stiffener"
OTHER_INNER_ROW = "Other inner row"
OTHER_END_ROW = "Other end row"
END_ROW_ADJACENT_TO_STIFFENER = "End row adjacent to stiffener"
""" 3. end plate """
ROW_OUTSIDE_BEAM_TENSION_FLANGE = "Outside beam tension flange"
FIRST_ROW_BELOW_BEAM_TENSION_FLANGE = "First row below beam tension flange"


""" CHAPTER 3: Bolts """
class Bolt:
    """ Class for bolts """
    
    def __init__(self,size,bolt_class=8.8,length=50.0):
        """ Constructor
             input: size .. bolt size
                    bolt_class .. bolt class
                    length .. length of bolt
        """
        self.d = size # diameter of unthreaded shank
        # This is the minimum length
        self._length = length
        
        
        bolt_type = 'Full thread'
        bolt_available = False
        for key, value in bolt_unit_cost_part_thread.items():
            if key >= length:
                if value[size] != 'N/A':
                    length = key                    
                    bolt_available = True
                    bolt_type = 'Partial thread'
                    break
        
        """ If bolt is not available in this length,
            try fully threaded bolts
        """
        if not bolt_available:
            for key, value in bolt_unit_cost_full_thread.items():
                if key >= length:
                    if value[size] != 'N/A':
                        length = key                    
                        bolt_available = True                        
                        break
        
        if not bolt_available:
            print("Warning: Bolt of size {0:2.0f} is not available with minimum length of {1:4.2f} mm".format(size,length))
        #else:
        #    print("Bolt length {0:4.0f} mm".format(length))
        
        # bolt hole: EN 1090-2 (2018): Table 11
        if size < 16:
            self.d0 = size+1
        elif size < 27:
            self.d0 = size+2
        else:
            self.d0 = size+3
        
        self.A = math.pi*size**2/4
        self.As = bolt_size[size]["A_s"]
        self.head_t = bolt_size[size]["k"]
        self.head_d = bolt_size[size]["s"]
        
        self.nut_t = nut_size[size]["m"]
        self.nut_e = nut_size[size]["e"]
        
        self.washer_t = washer_size[size]["h"]
        self.washer_d = washer_size[size]["d2"]
        
        self.length = length
        self.fyb = mat_bolt[bolt_class]["f_yb"]
        self.fub = mat_bolt[bolt_class]["f_ub"]
        self.bolt_class = bolt_class
        self.bolt_type = bolt_type
        
        
    @property
    def dw(self):
        """ Maximum diameter covere by washer, bolt head or bolt nut """
        return max(self.washer_d,self.head_d,self.nut_e)
    
    @property
    def size(self):
        return self.d
    
    @size.setter
    def size(self,val):
        """ If the size of the bolt is changed, the following must be done:
            1. set d
            2. set d0
            3. recalculate A, As, head_t, nut_t, nut_e, washer_t, washer_d
            4. check the length
        """
        self.d = val
        
        if val < 16:
            self.d0 = val+1
        elif val < 27:
            self.d0 = val+2
        else:
            self.d0 = val+3
        
        self.A = math.pi*val**2/4
        self.As = bolt_size[val]["A_s"]
        self.head_t = bolt_size[val]["k"]
        self.head_d = bolt_size[val]["s"]
        
        self.nut_t = nut_size[val]["m"]
        self.nut_e = nut_size[val]["e"]
        
        self.washer_t = washer_size[val]["h"]
        self.washer_d = washer_size[val]["d2"]
        
        bolt_available = False
        for key, value in bolt_unit_cost_part_thread.items():
            if key >= self._length:
                if value[val] != 'N/A':
                    length = key                    
                    bolt_available = True
                    bolt_type = 'Partial thread'
                    break
        
        """ If bolt is not available in this length,
            try fully threaded bolts
        """
        if not bolt_available:
            for key, value in bolt_unit_cost_full_thread.items():
                if key >= self._length:
                    if value[val] != 'N/A':
                        length = key                    
                        bolt_available = True
                        bolt_type = 'Full thread'
                        break
        
        self.length = length
        self.bolt_type = bolt_type
    
    def shear_resistance(self,threads_in_plane=True,verb=False):
        """ EN 1993-1-8, Table 3.4
        Input:
            fub .. ultimate strength of bolt [MPa]
            A .. cross-sectional area [mm^2]
            threads_at_plane .. True, if shear plane passes through threaded
                                portion of the bolt
                                False, otherwise
        Output:
            FvRd .. bolt shear resistance [N]
        """
        
        if threads_in_plane:
            As = self.As
            if self.bolt_class in {10.9,6.8,5.8,4.8}:
                av = 0.5
            elif self.bolt_class in {4.6,5.6,8.8}:
                av = 0.6
        else:
            av = 0.6
            As = self.A
            
    
        
        FvRd = av*self.fub*As/gammaM2
        
        if verb:
            print("Shear resistance:")
            if threads_in_plane:
                print("Bolt threads in shear plane")
            else:
                print("Bolt threads not in shear plane")
                
            print("Area: As = {0:4.2f}".format(As))            
            print("alpha_v = {0:4.2f}".format(av))            
            print("fub = {0:4.2f}".format(self.fub))            
            print("FvRd = {0:4.2f} kN".format(FvRd*1e-3))    
        
        
        return FvRd
    
    def tension_resistance(self,k2=0.9):
        """ EN 1993-1-8, Table 3.4
            Input:
                k2 = 0.63 for countersunk bolts
                   = 0.9 otherwise (default)
            Output:
                FtRd .. bolt tension resistance [N]
        """
        FtRd = k2*self.fub*self.As/gammaM2
        return FtRd
    
 
    def bearing_resistance(self,fu,t,e,p,pos_perp="edge",pos_load="edge",\
                           verb=False):
        """ EN 1993-1-8, Table 3.4
            Input:
                t .. thickness of plate [mm]
                pos_perp .. position of bolt in direction perpendicular to
                            load transfer ("edge" or "inner")
                pos_load .. position of bolt in direction of load transfer
                            ("edge" or "inner")
                e .. array for edge distances e[0] = e1, e[1] = e2
                p .. array for bolt row distances p[0] = p1, p[1] = p2
                fu .. ultimate strength of plate [MPa]
            Output:
                FbRd .. bearing resistance [N]
        """
        e1 = e[0]
        e2 = e[1]
        p1 = p[0]
        p2 = p[1]
        
        if pos_load == "edge":
            ad = e1/3/self.d0
        else:
            ad = p1/3/self.d0 - 0.25    
    
        if pos_perp == "edge":
            k1 = min(2.8*e2/self.d0-1.7,2.5)
        else:
            k1 = min(1.4*p2/self.d0-1.7,2.5)
            
        
        ab = min(ad,self.fub/fu,1.0)
        
        FbRd = k1*ab*fu*self.d*t/gammaM2
        
        if verb:
            print("Bearing resistance:")
            print("Edge distances: e1 = {0:4.2f}, e2 = {1:4.2f}".format(e1,e2))
            print("Bolt row distances: p1 = {0:4.2f}, p2 = {1:4.2f}".format(p1,p2))
            print("alpha_d = {0:4.2f}".format(ad))
            print("alpha_b = {0:4.2f}".format(ab))
            print("k1 = {0:4.2f}".format(k1))
            print("d = {0:4.2f}".format(self.d))
            print("t = {0:4.2f}".format(t))
            print("t = {0:4.2f}".format(t))
            print("FbRd = {0:4.2f} kN".format(FbRd*1e-3))
        
        return FbRd

class BoltRow:
    """ Class for bolt rows 
        
        Methods and functionality:
            - compute tension resistance of each component
            - compute stiffness of each component
        
        The components are:
            - end plate in bending
            - column flange in bending
            - column web in tension
            - beam web in tension
    
    """
    
    def __init__(self,bolt,p,z,flange_loc=INNER_ROW, plate_loc=OTHER_INNER_ROW,joint=None, row_type=TENSION_ROW):
        """ Constructor
            input:
                bolt .. Bolt class object defining the bolt
                p .. array of distances between adjacent bolt centroids
                z .. vertical position of bolt row
                flange_loc .. position of the row with respect to column flange in bending
                plate_loc .. position of the row with respect to end plate in bending
                joint .. the joint object to which the bolt row belongs
        """
        
        self.type = row_type
        self.bolt = bolt
        #self.p = p
        self.z = z
        self.loc_col_flange = flange_loc
        self.loc_end_plate = plate_loc
        
        """ Location of bolt row in a group. This varies,
            if the row belongs to several groups        
        """
        self.group_loc_flange = None
        self.group_loc_plate = None
        
        self.joint = joint
        
        self._z_top = 0.5*self.joint.beam.h+self.joint.etop - z
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
        
class BoltRowGroup:
    """ Class for group of bolt rows 
        Bolt row group consists of the following data:
            1) Set of bolt rows
            2) Location of each bolt row in terms of column flange
                and end plate. This data can be in the form of a
                list of dictionaries:
                    [{'flange':LOCATION, "plate": LOCATION}]
            NOTE! For some groups, not all components are active. For example,
            a group consisting of two rows, above and below the beam tension flange,
            has the component 'column flange in bending' but not end plate in bending,
            because the column flange prevents yield lines between the rows on opposite
            sides of the beam flange.
        
        Each bolt row has a T-stub for column flange in bending and for end plate in bending.
        These T-stubs can be used for bolt row groups as well.
    
    """
    
    def __init__(self,bolt_rows,row_loc):
        """ Constructor
            input:
                bolt_rows .. array of BoltRow objects
                row_loc .. location of rows
                
            attributes:
                rows .. list of BoltRow objects
                nrows .. number of bolt rows (MAYBE THIS IS NOT NEEDED!)
                row_loc .. list of bolt row location codes.
                p .. distance between adjacent bolt rows
        """
        
        self.rows = bolt_rows
        self.nrows = len(bolt_rows)
        self.row_loc = row_loc
        #self.p = []
        # Tension resistance of the group
        self.FtRd = 0.0
        
        # Effective length of the bolt row group
        self.leff_flange = [0.0,0.0]
        self.leff_plate = [0.0,0.0]
        
        """ Evaluate distance between adjacent groups
            it is assumed that the rows as inserted in ordered list
        """
        #for i in range(len(bolt_rows)-1):
        #    self.p.append(bolt_rows[i].z-bolt_rows[i+1].z)

        
        self.Tstub_col_flange = TStubGroupColumnFlange(self)
        
        # If the group has two bolt rows and the first bolt row is
        # outside the beam tension flange, then the group does not have the
        # component "End plate in bending".
        #print('Bolt row group created:')
        #print(f'Number of rows: {len(bolt_rows)}')
        #print(f"First row position: {row_loc[0]['plate']}")
        if len(bolt_rows) == 2 and row_loc[0]['plate'] == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            self.Tstub_end_plate = None
        else:
            self.Tstub_end_plate = TStubGroupEndPlate(self)

        # Link this group to the bolt rows
        for row in bolt_rows:
            row.groups.append(self)
    
    @property
    def p(self):
        """ Returns an array of vertical distances between consecutive
            rows.
        """
        return [self.rows[i].z-self.rows[i+1].z for i in range(len(self.rows)-1)]
    def info(self):
        """ Prints info """
        
        print("Bolt row group:")
        for i, row in enumerate(self.row_loc):
            print(" Row {0:1.0f}".format(i))
            print("Location (flange): {0:s}".format(row['flange']))
            print("Location (end plate): {0:s}".format(row['plate']))

    def column_web_in_tension(self,verb=False):
        """ Resistance of component 
        
        """
        
        """ The effective length of the component
            corresponds to the sum of the effective lengths of the bolt rows
            for column flange.
            
            This is obtained only after the resistance of the column
            flange has been determined, because the governing failure
            mode also defines the effective width
        """
        
        l_eff = self.Tstub_col_flange.leff()
        
        F_t_wc_Rd = cm.column_web_tension(self.rows[0].joint.col, l_eff, self.rows[0].joint.beta)
        
        if verb:
            print("COlUMN WEB IN TENSION:")
            print("l_eff  = {0:4.2f} mm".format(l_eff))
            #print("l_eff (Mode 2) = {0:4.2f} mm".format(self.Tstub_col_flange._leff[1]))
            print("Ft_wc_Rd = {0:4.2f} kN".format(F_t_wc_Rd*1e-3))
        
        return F_t_wc_Rd
        
    def column_flange_in_bending(self,verb=False):
        """ Resistance of component """
        
        """ Determine effective length as the sum of effective lengths of the bolt rows
            leff_nc = sum_r leff_nc(r)
            leff_cp = sum_r leff_cp(r)
            
            leff_1 = min(leff_nc,leff_cp)
            leff_2 = leff_nc
            
        """
        if verb:
            print("COLUMN FLANGE IN BENDING:")
        
        leff_nc = self.Tstub_col_flange.leff_nc
        leff_cp = self.Tstub_col_flange.leff_cp
        self.leff_flange[0] = min(leff_nc,leff_cp)
        self.leff_flange[1] = leff_nc
        
        
        if verb:
            print("l_eff (Mode 1) = {0:4.2f} mm".format(self.leff_flange[0]))
            print("l_eff (Mode 2) = {0:4.2f} mm".format(self.leff_flange[1]))
        
        F_t_cf_Rd = cm.column_flange_bending(self.Tstub_col_flange,verb)
        
        if verb:            
            print("Ft_cf_Rd = {0:4.2f} kN".format(F_t_cf_Rd*1e-3))
        
        return F_t_cf_Rd
    
    def end_plate_in_bending(self,verb=False):
        """ Resistance of component 
            It is assumed that the dimensions of the T-stubs
            of individual rows have been set earlier
        """
        
        leff_nc = self.Tstub_end_plate.leff_nc
        leff_cp = self.Tstub_end_plate.leff_cp
        self.leff_plate[0] = min(leff_nc,leff_cp)
        self.leff_plate[1] = leff_nc
        
        
        if verb:
            print("END PLATE IN BENDING:")            
            print("l_eff (Mode 1) = {0:4.2f} mm".format(self.leff_plate[0]))
            print("l_eff (Mode 2) = {0:4.2f} mm".format(self.leff_plate[1]))
        
        F_t_ep_Rd = cm.end_plate_bending(self.Tstub_end_plate,verb)
        
        return F_t_ep_Rd
    
    def beam_web_in_tension(self,verb=False):
        """ Resistance of component """
        l_eff = self.Tstub_end_plate.leff()        
        
        F_t_bw_Rd = cm.beam_web_tension(self.rows[0].joint.beam, l_eff)
        
        if verb:
            print("BEAM WEB IN TENSION:")
            print("l_eff = {0:4.2f} mm".format(l_eff))
            print("Ft_bf_Rd = {0:4.2f} kN".format(F_t_bw_Rd*1e-3))
        
        return F_t_bw_Rd
    
    def tension_resistance(self,verb=False):
        """ Tension resistance of the bolt row group.
        """
        Ft_i_Rd = np.array([0.0,0.0,0.0,0.0])
        components = ['Column flange in bending',
                      'End plate in bending',
                      'Column web in tension',
                      'Beam web in tension'
                      ]
        
        if self.row_loc[0]['flange'] == END_ROW_ADJACENT_TO_STIFFENER:
            """ In this case, the first row is on the opposite side
                of the stiffener than the other rows. This implies that
                the resistance related to column flange in bending of this
                group is infinite, i.e. the component need not be considered.
            """
            Ft_i_Rd[0] = np.inf
            Ft_i_Rd[2] = np.inf
        elif self.row_loc[0]['flange'] == ADJACENT_TO_STIFFENER_ROW and self.row_loc[0]['plate'] == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            Ft_i_Rd[0] = np.inf
            Ft_i_Rd[2] = np.inf
        else:
            Ft_i_Rd[0] = self.column_flange_in_bending(verb)
            Ft_i_Rd[2] = self.column_web_in_tension(verb)
        
        """ If the group has two rows and if the first row
            is above the tension flange of the beam,
            then end plate in bending group component can be
            neglected as well as beam web in tension
        """

        #if self.nrows == 2 and self.row_loc[0]['plate'] == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
        if self.row_loc[0]['plate'] == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            #print("First row outside tension flange.")
            Ft_i_Rd[1] = np.inf
            Ft_i_Rd[3] = np.inf
        else:
            Ft_i_Rd[1] = self.end_plate_in_bending(verb)        
            Ft_i_Rd[3] = self.beam_web_in_tension(verb)
            
        min_comp = np.argmin(Ft_i_Rd)
        
        #print(min_comp,Ft_i_Rd,Ft_i_Rd[min_comp])
        
        self.FtRd = Ft_i_Rd[min_comp]
        
        if verb:
            print("Ft_Rd = {0:4.2f} kN ({1:s})".format(self.FtRd*1e-3,components[min_comp]))
        
        return self.FtRd
    
    
class TStub:
    """ Class for T-stubs 
        Abstract class that serves as a basis for the following components:
            
         - column flange in bending; 
         - end-plate in bending; 
         - flange cleat in bending; 
         - base plate in bending under tension
    
    """
    
    #def __init__(self,bolts,row_loc,material,tf,emin,m):
    def __init__(self,bolt_row,material,tf,emin,m):
        """ Constructor
            Input:
                bolt_row .. BoltRow object
                row_loc .. location of the row (e.g. INNER_ROW)
                material .. Steel class object, for the material properties
                tf .. thicknes of the flange
                n .. edge distance from the centroid of bolt hole
                m .. distance from center of bolt hole towards web of the T-stub
                
        """
        
        self.bolts = bolt_row
        
        #self.bolts = bolts
        #self.row_loc = row_loc
        self.mat = material
        #self.tf = tf
        #self.emin = emin
        #self.m = m
        self._leff = np.array([-1.0,-1.0])
        self._leff_group = {'cp':-1.0,'nc':-1.0}
        
        self.failure_mode = -1
    
    def info(self,name=None):
        """ Prints relevant information of the T-stub """
        
        if name is None:
            print("T-stub:")
        else:
            print("Tstub (" + name + "): ")
    
        print("tf = {0:4.2f}".format(self.tf))
        print("emin = {0:4.2f}".format(self.emin))
        print("m = {0:4.2f}".format(self.m))
        
    @property
    def fy(self):
        return self.mat.fy
    
    @property
    def emin(self):
        """ Distance from bolt row center to the edge of the
            T-stub flange.
        """
        pass
    
    @property
    def tf(self):
        """ Thickness of T-stub flange """
        pass
    
    @property
    def m(self):
        """ Distance from bolt hole center to the web of the
            T-stub, taking into account weld and rounded corner of
            the profile.
        """
        pass
    
    @property
    def m2(self):
        """ Distance from bolt hole center to the flange weld,
            see Fig. 6.11. of EN 1993-1-8:2005
        """
        pass
    
    @property
    def n(self):
        return min(self.emin,1.25*self.m)
            
    def leff(self):
        """ Returns the effective width corresponding
            to the governing failure mode
            
        """
        
        """ If no failure mode has been detected,
            calculate resistance of the T-stub
        """
        if self.failure_mode == -1:
            self.FT_Rd()
         
        if self.failure_mode == 0:
            """ Complete mechanism governs """
            leff = self._leff[0]
        elif self.failure_mode == 1:
            """ Flange yield with bold failure governs """
            leff = self._leff[1]
        else:
            """ Bolt tension resistance governs 
                NOTE! This is not clear!
            """
            leff = self._leff[1]
        
        return leff
    
    @property
    def leff_1(self):
        """ Effective length for Mode 1 """
            
        self._leff[0] = min(self.leff_cp,self.leff_nc)
        
        return self._leff[0]
    
    @property
    def leff_2(self):
        """ Effective length for Mode 2 """
        
        self._leff[1] = self.leff_nc
        
        return self._leff[1]
    
    @property
    def leff_nc(self):
        """ Non-circular yield line patterns """
        pass
    
    @property
    def leff_cp(self):
        """ Circular yield line patterns """
        pass
    
    def leff_nc_as_group(self,position=None):
        pass
    
    def leff_cp_as_group(self,position=None):
        pass
    
    @property
    def MplRd_1(self):
        """ Plastic moment of mode 1 """
        
        return 0.25*self.leff_1*self.tf**2*self.fy/gammaM0
    
    @property
    def MplRd_2(self):
        """ Plastic moment of mode 2 """
        
        return 0.25*self.leff_2*self.tf**2*self.fy/gammaM0

    def F_T_1_Rd(self):
        """ Tension resistance of mode 1 """
        return 4*self.MplRd_1/self.m
    
    def F_T_2_Rd(self):
        """ Tension resistance of mode 2 """
        return (2*self.MplRd_2 + self.n*self.F_T_3_Rd())/(self.m+self.n)
    
    def F_T_12_Rd(self):
        """ Tension resistance without prying forces """
        return 2*self.MplRd_1/self.m
    
    def F_T_3_Rd(self):
        """ Tension resistance of mode 3 
            Bolt failure
        """
        #FRd = 0.0
        
        return self.bolts.bolts*self.bolts.bolt.tension_resistance()
        
        """
        for bolt in self.bolts:
            FRd += bolt.tension_resistance()
            
        return FRd
        """
    def FT_Rd(self,verb=False):
        """ Design tension resistance """
        
        FT_Rd = np.array([0.0,0.0,0.0])
        
        FT_Rd[0] = self.F_T_1_Rd()
        FT_Rd[1] = self.F_T_2_Rd()
        FT_Rd[2] = self.F_T_3_Rd()
        
        if verb:
            for i, FRd in enumerate(FT_Rd):
                print("FT_{0:1.0f}_Rd = {1:4.2f} kN".format(i+1,FRd*1e-3))
        
        fmode = np.argmin(FT_Rd)
        self.failure_mode = fmode
        
        return FT_Rd[fmode]

class TStubColumnFlange(TStub):
    """ T-stub for column flange """
    
    #def __init__(self,bolts,row_loc,material,tf,e,emin,m,e1=math.inf,p=0):
    def __init__(self,bolt_row,material,tf,e,emin,m,e1=math.inf,p=0):
        """ Constructor
            
            Attributes:
                e1 .. distance of end bolt row from column end
                e .. distance of column flange edge from bolt hole center
                p .. distance to the next bolt row
        """
        
        #super().__init__(bolts,row_loc,material,tf,emin,m)
        super().__init__(bolt_row,material,tf,emin,m)
        self.e1 = e1
        #self.e = e
        self.p = p
        #self.m2 = 0
    
    @property
    def tf(self):
        """ T-stub flange thickness is the thickness of the column flange """
        return self.bolts.joint.col.tf
    
    @property
    def emin(self):
        """ See Fig. 6.8 of EN 1993-1-8:2005 """   
        return self.bolts.emin_flange
    
    @property
    def m(self):
        """ See Fig. 6.8 of EN 1993-1-8:2005 """
        return self.bolts.col_flange_m
    
    @property
    def m2(self):
        """ See Fig. 6.11 of EN 1993-1-8:2005 
            NOTE: This value is relevant only if the bolt row
                  is the first below a transverse column web stiffener
        """
        return self.bolts.top_flange_stiffener_m2
    
    @property
    def e(self):
        """ Distance between bolt hole center and edge of column flange """
        return self.bolts.col_flange_e
    
    @property
    def leff_cp(self):
        """ circular pattern """
        
        """ Treat row individually """
        if self.bolts.loc_col_flange == INNER_ROW or self.bolts.loc_col_flange == OTHER_INNER_ROW:
            leff = 2*math.pi*self.m
        elif self.bolts.loc_col_flange == END_ROW or self.bolts.loc_col_flange == OTHER_END_ROW:
            leff = min(2*math.pi*self.m,math.pi*self.m+2*self.e1)
        elif self.bolts.loc_col_flange == ADJACENT_TO_STIFFENER_ROW:
            leff = 2*math.pi*self.m
        elif self.bolts.loc_col_flange == END_ROW_ADJACENT_TO_STIFFENER:
            leff = min(2*math.pi*self.m,math.pi*self.m+2*self.e1)
            
        return leff

    @property
    def leff_nc(self):        
        """ non-circular pattern """        
        
        if self.bolts.loc_col_flange == INNER_ROW or self.bolts.loc_col_flange == OTHER_INNER_ROW:
            leff = 4*self.m + 1.25*self.e            
        elif self.bolts.loc_col_flange == END_ROW or self.bolts.loc_col_flange == OTHER_END_ROW:
            leff = min(4*self.m + 1.25*self.e,2*self.m + 0.625*self.e + self.e1)
        elif self.bolts.loc_col_flange == ADJACENT_TO_STIFFENER_ROW:
            try:
                lambda1 = self.m/(self.m+self.e)
                lambda2 = self.m2/(self.m+self.e)
                alpha = cm.find_alpha(lambda1,lambda2)
                #print(lambda1,lambda2,alpha)
            except:
                alpha = cm.new_alpha(self.e, self.m, self.m2)                            
                #print(alpha)
            leff = alpha*self.m
        elif self.bolts.loc_col_flange == END_ROW_ADJACENT_TO_STIFFENER:
            try:
                lambda1 = self.m/(self.m+self.e)
                lambda2 = self.m2/(self.m+self.e)
                alpha = cm.find_alpha(lambda1,lambda2)
                #print(lambda1,lambda2,alpha)
            except:
                alpha = cm.new_alpha(self.e, self.m, self.m2)                            
                #print(alpha)

            leff = self.e1 + alpha*self.m - (2*self.m + 0.625*self.e)
        else:
            print("Warning: unidentified row position.")

        return leff

    def leff_cp_as_group(self,position=None):
        """ circular pattern for bolt row considered as part of a group """
        """ Treat row as in a bolt row group """
        if position is None:
            position = self.bolts.group_loc_flange
        
        
        if position == INNER_ROW or position == OTHER_INNER_ROW:
            leff = 2*self.p
        elif position == END_ROW or position == OTHER_END_ROW:                        
            leff = min(math.pi*self.m+self.p,2*self.e1+self.p) 
        elif position == ADJACENT_TO_STIFFENER_ROW:
            leff = math.pi*self.m+self.p        
            
        self._leff_group['cp'] = leff
        
        return leff
    
    def leff_nc_as_group(self,position=None):
        """ non-circular pattern """
        if position is None:
            position = self.bolts.group_loc_flange
        
        #print(position)
        
        if position == INNER_ROW or position == OTHER_INNER_ROW:
            leff = self.p
        elif position == END_ROW or position == OTHER_END_ROW:
            leff = min(2*self.m+0.625*self.e+0.5*self.p, self.e1+0.5*self.p)    
        elif position == ADJACENT_TO_STIFFENER_ROW:
            try:
                lambda1 = self.m/(self.m+self.e)
                lambda2 = self.m2/(self.m+self.e)
                alpha = cm.find_alpha(lambda1,lambda2)
                #print(lambda1,lambda2,alpha)
            except:
                alpha = cm.new_alpha(self.e, self.m, self.m2)                            
                #print(alpha)
            print(self.p,alpha,self.m,self.e)
            leff = 0.5*self.p + alpha*self.m - (2*self.m + 0.625*self.e)
            
            
        self._leff_group['nc'] = leff
        
        return leff
        
class TStubEndPlate(TStub):
    """ T-stub for end plate in bending """
    
    #def __init__(self,bolts,row_loc,material,tf,emin,m,e,w,bp,p=0):
    def __init__(self,bolt_row,material,tf,emin,m,e,w,bp,p=0):
        """ Constructor
            
            Attributes:
                tf .. end plate thickness
                emin .. 'ex' for end plate extension and value from Fig. 6.8 of EN 1993-1-8 otherwise
                m .. distance from beam web to bolt hole center
                e .. distancef from bolt hole center to plate edge
                w .. distance between bolts in the row
                p .. distance to the next bolt row
                bp .. width of end plate
        """
        
        super().__init__(bolt_row,material,tf,emin,m)
        #self.e = e
        self.p = p
        #self.w = w
        #self.m2 = 0
        #self.bp = bp
    
    def info(self):
        """ Prints relevant info """
        
        super().info("end plate")
        
        print("e = {0:4.2f}".format(self.e))
        print("w = {0:4.2f}".format(self.w))
    
    @property
    def bp(self):
        """ Width of end plate """
        return self.bolts.joint.bp
    
    @property
    def tf(self):
        """ T-stub flange thickness is the thickness of the end plate """
        return self.bolts.joint.tp
    
    
    @property
    def emin(self):
        """ See Fig. 6.10 of EN 1993-1-8:2005 """
        if self.bolts.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            emin = self.bolts.ex
        else:
            emin = self.bolts.end_plate_e
            
        return emin
    
    @property
    def m(self):
        """ See Fig. 6.8 of EN 1993-1-8:2005 """
        return self.bolts.end_plate_m
    
    @property
    def m2(self):
        """ See Fig. 6.11 of EN 1993-1-8:2005 
            NOTE: This value is relevant only if the bolt row
                  is the first below a transverse column web stiffener
        """
        return self.bolts.end_plate_m2
    
    @property
    def e(self):
        """ Distance between bolt hole center and edge of column flange """
        return self.bolts.end_plate_e
    
    @property
    def w(self):
        """ Horizontal distance between bolt holes """
        return self.bolts.w
    
    @property
    def leff_cp(self):
        """ circular pattern """
        
        if self.bolts.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            """ self.m is the 'mx' and
                self.e is the 'ex' 
                    
                of Table 6.6. of EN 1993-1-8
            """
            leff = min(2*math.pi*self.m,math.pi*self.m+self.w,math.pi*self.m+2*self.e)
        elif self.bolts.loc_end_plate == FIRST_ROW_BELOW_BEAM_TENSION_FLANGE:
            leff = 2*math.pi*self.m
        elif self.bolts.loc_end_plate == OTHER_INNER_ROW:
            leff = 2*math.pi*self.m
        elif self.bolts.loc_end_plate == OTHER_END_ROW:
            leff = 2*math.pi*self.m
            
        return leff

    @property
    def leff_nc(self):        
        """ non-circular pattern """
        
        if self.bolts.loc_end_plate == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            leff = min(4*self.m + 1.25*self.emin,
                       self.e+2*self.m+0.625*self.emin,
                       0.5*self.bp,
                       0.5*self.w+2*self.m+0.625*self.emin)
        elif self.bolts.loc_end_plate == FIRST_ROW_BELOW_BEAM_TENSION_FLANGE:
            lambda1 = self.m/(self.m+self.e)
            lambda2 = self.m2/(self.m+self.e)
            
            #print(lambda1,lambda2)
            #alpha = cm.par_alfa(lambda1, lambda2)        
            
            alpha = cm.new_alpha(self.e, self.m, self.m2)
            #print(alpha)
            
            leff = alpha*self.m
        elif self.bolts.loc_end_plate == OTHER_INNER_ROW:           
            leff = 4*self.m+1.25*self.e
        elif self.bolts.loc_end_plate == OTHER_END_ROW:
            leff = 4*self.m + 1.25*self.e

        return leff
    
    def leff_cp_as_group(self,position=None):
        
        if position is None:
            position = self.bolts.group_loc_plate
        
        if position == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            leff = 0.0
        elif position == FIRST_ROW_BELOW_BEAM_TENSION_FLANGE:
            leff = math.pi*self.m+self.p
        elif position == OTHER_INNER_ROW:
            leff = 2*self.p
        elif position == OTHER_END_ROW:
            leff = math.pi*self.m + self.p
        
        return leff
    
    def leff_nc_as_group(self,position=None):
        
        if position is None:
            position = self.bolts.group_loc_plate
        
        if position == ROW_OUTSIDE_BEAM_TENSION_FLANGE:
            leff = self.p
        elif position == FIRST_ROW_BELOW_BEAM_TENSION_FLANGE:
            lambda1 = self.m/(self.m+self.e)
            lambda2 = self.m2/(self.m+self.e)   
            
            #print(lambda1,lambda2)
            #alpha = cm.par_alfa(lambda1, lambda2)
            alpha = cm.new_alpha(self.e, self.m, self.m2)
            leff = 0.5*self.p + alpha*self.m -(2*self.m+0.625*self.e)                
            #print(alpha,self.p,self.m,self.e,leff)
        elif position == OTHER_INNER_ROW:
            leff = self.p
        elif position == OTHER_END_ROW:
            leff = 2*self.m+0.625*self.e + 0.5*self.p
        
        return leff
        
class TStubGroup:
    """ Class for T-stubs for bolt row groups """
    
    def __init__(self,row_group):
        """ Constructor
            input:
                row_group .. BoltRowGroup object that constitutes the T-stub
        """
        
        self.group = row_group
        
        self.failure_mode = -1
        
        self._leff = np.array([-1.0,-1.0])
        
        # Storage for effective widths for bolt rows in the group
        self._leff_rows = {}
        
        for row in self.group.rows:
            """ In the dictionary, the bolt row objects themselves
                act as the keys.
            """
            self._leff_rows[row] = {'leff_nc':0.0,
                                    'leff_cp':0.0}
    
    @property
    def fy(self):
        pass
    
    @property
    def tf(self):
        """ Thickness of the flange of the T-stub """
        pass
    
    def n(self):
        pass
    
    def m(self):
        pass
        
    def leff(self):
        """ Returns the effective width of the T-stub
        """
        
        """ If no failure mode has been detected,
            calculate resistance of the T-stub
        """
        if self.failure_mode == -1:
            self.FT_Rd()
         
        if self.failure_mode == 0:
            """ Complete mechanism governs """
            leff = self._leff[0]
        elif self.failure_mode == 1:
            """ Flange yield with bold failure governs """
            leff = self._leff[1]
        else:
            """ Bolt tension resistance governs 
                NOTE! This is not clear!
            """
            leff = self._leff[1]
        
        return leff
    
    @property
    def leff_1(self):
        """ Effective length for Mode 1 """
            
        self._leff[0] = min(self.leff_cp,self.leff_nc)
        
        return self._leff[0]
    
    @property
    def leff_2(self):
        """ Effective length for Mode 2 """
        
        self._leff[1] = self.leff_nc
        
        return self._leff[1]
    
    @property
    def leff_nc(self):
        """ Non-circular yield line patterns """
        pass
    
    @property
    def leff_cp(self):
        """ Circular yield line patterns """
        pass
    
    def leff_nc_as_group(self):
        pass
    
    def leff_cp_as_group(self):
        pass
    
    def leff_nc_row(self,row):
        pass
    
    def leff_cp_row(self,row):
        pass
    
    def leff_row(self,row):
       """ Calculates effective lengths for a given bolt row """
       self.leff_nc_row(row)
       self.leff_cp_row(row)
    
    @property
    def MplRd_1(self):
        """ Plastic moment of mode 1 """
        
        return 0.25*self.leff_1*self.tf**2*self.fy/gammaM0
    
    @property
    def MplRd_2(self):
        """ Plastic moment of mode 2 """
        
        return 0.25*self.leff_2*self.tf**2*self.fy/gammaM0

    def F_T_1_Rd(self):
        """ Tension resistance of mode 1 """
        return 4*self.MplRd_1/self.m
    
    def F_T_2_Rd(self):
        """ Tension resistance of mode 2 """
        return (2*self.MplRd_2 + self.n*self.F_T_3_Rd())/(self.m+self.n)
        
    def F_T_3_Rd(self):
        """ Tension resistance of mode 3 
            Bolt failure
        """
        FRd = 0.0
        
        for row in self.group.rows:
            FRd += row.bolts*row.bolt.tension_resistance()
            
        return FRd
    
    def FT_Rd(self,verb):
        """ Design tension resistance """
        
        FT_Rd = np.array([0.0,0.0,0.0])
        
        FT_Rd[0] = self.F_T_1_Rd()
        FT_Rd[1] = self.F_T_2_Rd()
        FT_Rd[2] = self.F_T_3_Rd()
        
        if verb:
            for i, FRd in enumerate(FT_Rd):
                print("FT_{0:1.0f}_Rd = {1:4.2f} kN".format(i+1,FRd*1e-3))
        
        fmode = np.argmin(FT_Rd)
        self.failure_mode = fmode
        
        return FT_Rd[fmode]
    
class TStubGroupColumnFlange(TStubGroup):
    """ Class for T-stub of a bolt row group for the
        component 'column flange in bending'
    """
    
    def __init__(self,row_group):
        """ Constructor """
        
        super().__init__(row_group)        
        
    
    @property
    def fy(self):
        return self.group.rows[0].Tstub_col_flange.mat.fy
    
    @property
    def tf(self):
        return self.group.rows[0].Tstub_col_flange.tf
        
    @property
    def n(self):
        
        return self.group.rows[0].Tstub_col_flange.n
    
    @property
    def m(self):
        
        return self.group.rows[0].Tstub_col_flange.m
    
    @property
    def leff_nc(self):
        """ Non-circular yield lines for the group """
        
        leff_nc = 0.0
        
        for i, row in enumerate(self.group.rows):
            
            position = self.group.row_loc[i]['flange']
            if position == END_ROW or position == ADJACENT_TO_STIFFENER_ROW:
                if i == 0:
                    row.Tstub_col_flange.p = self.group.p[i]
                else:
                    row.Tstub_col_flange.p = self.group.p[i-1]
            elif position == INNER_ROW:
                row.Tstub_col_flange.p = 0.5*(self.group.p[i-1]+self.group.p[i])
            
            #print("p = ",row.Tstub_col_flange.p)
            
            self._leff_rows[row]['leff_nc'] = row.Tstub_col_flange.leff_nc_as_group(position)
            leff_nc += self._leff_rows[row]['leff_nc']
            #leff_nc += row.Tstub_col_flange.leff_nc_as_group(position)
            #print(leff_nc)
        
        return leff_nc
    
    @property
    def leff_cp(self):
        
        leff_cp = 0.0
        
        for i, row in enumerate(self.group.rows):
            position = self.group.row_loc[i]['flange']
            
            if position == END_ROW or position == ADJACENT_TO_STIFFENER_ROW:
                if i == 0:
                    row.Tstub_col_flange.p = self.group.p[i]
                else:
                    row.Tstub_col_flange.p = self.group.p[i-1]
            elif position == INNER_ROW:
                row.Tstub_col_flange.p = 0.5*(self.group.p[i-1]+self.group.p[i])
        
            self._leff_rows[row]['leff_cp'] = row.Tstub_col_flange.leff_cp_as_group(position)
            leff_cp += self._leff_rows[row]['leff_cp']
            #leff_cp += row.Tstub_col_flange.leff_cp_as_group(position)
        
        return leff_cp
    
    def leff_nc_row(self,row):
        """ Calculates non-circular effective width for
            the 'row' (BoltRow class object belonging to the group)
        """
        try:
            i = self.group.rows.index(row)
            
            position = self.group.row_loc[i]['flange']
            if position == END_ROW:
                if i == 0:
                    row.Tstub_col_flange.p = self.group.p[i]
                else:
                    row.Tstub_col_flange.p = self.group.p[i-1]
            elif position == INNER_ROW:
                row.Tstub_col_flange.p = 0.5*(self.group.p[i-1]+self.group.p[i])
            
            self._leff_rows[row]['leff_nc'] = row.Tstub_col_flange.leff_nc_as_group(position)
            
            return self._leff_rows[row]['leff_nc']
        except:
            print("leff_nc_row: row not found in group.")
    
    
    def leff_cp_row(self,row):
        """ Calculates circular effective width for
            the 'row' (BoltRow class object belonging to the group)
        """
        try:
            i = self.group.rows.index(row)
            
            position = self.group.row_loc[i]['flange']
            if position == END_ROW:
                if i == 0:
                    row.Tstub_col_flange.p = self.group.p[i]
                else:
                    row.Tstub_col_flange.p = self.group.p[i-1]
            elif position == INNER_ROW:
                row.Tstub_col_flange.p = 0.5*(self.group.p[i-1]+self.group.p[i])
            
            self._leff_rows[row]['leff_cp'] = row.Tstub_col_flange.leff_cp_as_group(position)
            
            return self._leff_rows[row]['leff_cp']
        except:
            print("leff_nc_row: row not found in group.")
         
    
class TStubGroupEndPlate(TStubGroup):
    """ Class for T-stub of a bolt row group for the
        component 'end plate in bending'
    """
    
    def __init__(self,row_group):
        """ Constructor """
        
        super().__init__(row_group)
    
    @property
    def fy(self):
        return self.group.rows[0].Tstub_end_plate.mat.fy
    
    @property
    def tf(self):
        return self.group.rows[0].Tstub_end_plate.tf
        
    @property
    def n(self):
        
        return self.group.rows[0].Tstub_end_plate.n
    
    @property
    def m(self):
        
        return self.group.rows[0].Tstub_end_plate.m
    
    @property
    def leff_nc(self):
    
        leff_nc = 0.0
        
        for i, row in enumerate(self.group.rows):
            position = self.group.row_loc[i]['plate']
            if position == FIRST_ROW_BELOW_BEAM_TENSION_FLANGE:
                """ First row below tension flange is always
                    the first in any group
                """
                row.Tstub_end_plate.p = self.group.p[i]
            elif position == OTHER_END_ROW:
                """ Other inner row can be the first in a group
                    or the last
                """
                if i == 0:
                    """ Row is the first in a group """
                    row.Tstub_end_plate.p = self.group.p[i]
                else:
                    """ Row is the last in a group """
                    row.Tstub_end_plate.p = self.group.p[i-1]
                    
            elif position == OTHER_INNER_ROW:
                row.Tstub_end_plate.p = 0.5*(self.group.p[i-1]+self.group.p[i])
                
            self._leff_rows[row]['leff_nc'] = row.Tstub_end_plate.leff_nc_as_group(position)
            leff_nc += self._leff_rows[row]['leff_nc']
            #leff_nc += row.Tstub_end_plate.leff_nc_as_group(position)
            
        
        return leff_nc
    
    @property
    def leff_cp(self):
        
        leff_cp = 0.0
        
        for i, row in enumerate(self.group.rows):
            position = self.group.row_loc[i]['plate']
            
            if position == END_ROW:
                if i == 0:
                    row.Tstub_end_plate.p = self.group.p[i]
                else:
                    row.Tstub_end_plate.p = self.group.p[i-1]
            elif position == INNER_ROW:
                row.Tstub_end_plate.p = 0.5*(self.group.p[i-1]+self.group.p[i])
            
            self._leff_rows[row]['leff_cp'] = row.Tstub_end_plate.leff_cp_as_group(position)
            leff_cp += self._leff_rows[row]['leff_cp']
            #leff_cp += row.Tstub_end_plate.leff_cp_as_group(position)
        
        return leff_cp
    
    def leff_nc_row(self,row):
        """ Calculates non-circular effective width for
            the 'row' (BoltRow class object belonging to the group)
        """
        try:
            i = self.group.rows.index(row)
            
            position = self.group.row_loc[i]['plate']
            
            if position == FIRST_ROW_BELOW_BEAM_TENSION_FLANGE:
                """ First row below tension flange is always
                    the first in any group
                """
                row.Tstub_end_plate.p = self.group.p[i]
            elif position == OTHER_END_ROW:
                """ Other inner row can be the first in a group
                    or the last
                """
                if i == 0:
                    """ Row is the first in a group """
                    row.Tstub_end_plate.p = self.group.p[i]
                else:
                    """ Row is the last in a group """
                    row.Tstub_end_plate.p = self.group.p[i-1]
                    
            elif position == OTHER_INNER_ROW:
                row.Tstub_end_plate.p = 0.5*(self.group.p[i-1]+self.group.p[i])
            
            
            self._leff_rows[row]['leff_nc'] = row.Tstub_end_plate.leff_nc_as_group(position)
            
            return self._leff_rows[row]['leff_nc']
        except:
            print("leff_nc_row: row not found in group.")
    
    
    def leff_cp_row(self,row):
        """ Calculates circular effective width for
            the 'row' (BoltRow class object belonging to the group)
        """
        try:
            i = self.group.rows.index(row)
            
            position = self.group.row_loc[i]['plate']
            if position == FIRST_ROW_BELOW_BEAM_TENSION_FLANGE:
                """ First row below tension flange is always
                    the first in any group
                """
                row.Tstub_end_plate.p = self.group.p[i]
            elif position == OTHER_END_ROW:
                """ Other inner row can be the first in a group
                    or the last
                """
                if i == 0:
                    """ Row is the first in a group """
                    row.Tstub_end_plate.p = self.group.p[i]
                else:
                    """ Row is the last in a group """
                    row.Tstub_end_plate.p = self.group.p[i-1]
                    
            elif position == OTHER_INNER_ROW:
                row.Tstub_end_plate.p = 0.5*(self.group.p[i-1]+self.group.p[i])
            
            self._leff_rows[row]['leff_cp'] = row.Tstub_end_plate.leff_cp_as_group(position)
            
            return self._leff_rows[row]['leff_cp']
        except:
            print("leff_nc_row: row not found in group.")

    
def bolt_shear_resistance(fub,A,bolt_class=8.8,threads_in_plane=False):
    """ EN 1993-1-8, Table 3.4
        Input:
            fub .. ultimate strength of bolt [MPa]
            A .. cross-sectional area [mm^2]
            threads_at_plane .. True, if shear plane passes through threaded
                                portion of the bolt
                                False, otherwise
        Output:
            FvRd .. bolt shear resistance [N]
    """
    
    if threads_in_plane:
        av = 0.6
    else:
        if bolt_class in {10.9,6.8,5.8,4.8}:
            av = 0.5
        elif bolt_class in {4.6,5.6,8.8}:
            av = 0.6
    
    FvRd = av*fub*A/gammaM2
    
    
    
    return FvRd

def bolt_bearing_resistance(fub,fu,d,t,e,p,d0,pos_perp="edge",pos_load="edge"):
    """ EN 1993-1-8, Table 3.4
        Input:
            fub .. ultimate strength of bolt [MPa]
            fu .. ultimate strength of plate [MPa]
            d .. diameter of bolt hole [mm]
            t .. thickness of plate [mm]
            pos_perp .. position of bolt in direction perpendicular to
                        load transfer ("edge" or "innter")
            pos_load .. position of bolt in direction of load transfer
                        ("edge" or "inner")
        Output:
            FbRd .. bearing resistance [N]
    """
    e1 = e[0]
    e2 = e[1]
    p1 = p[0]
    p2 = p[1]
    
    if pos_load == "edge":
        ad = e1/3/d0
    else:
        ad = p1/3/d0 - 0.25    

    if pos_perp == "edge":
        k1 = min(2.8*e2/d0-1.7,2.5)
    else:
        k1 = min(1.4*p2/d0-1.7,2.5)
        
    
    ab = min(ad,fub/fu,1.0)
    
    FbRd = k1*ab*fu*d*t    
    
    return FbRd

def bolt_tension_resistance(fub,As,k2=0.9):
    """ EN 1993-1-8, Table 3.4
        Input:
            fub .. ultimate strength of bolt [MPa]
            As .. stress area of bolt [mm^2]
            k2 = 0.63 for countersunk bolts
               = 0.9 otherwise (default)
        Output:
            FtRd .. bolt tension resistance [N]
    """
    FtRd = k2*fub*As/gammaM2
    return FtRd

def bolt_punching_shear_resistance(dm,tp,fu):
    """ EN 1993-1-8, Table 3.4
        Input:
            dm .. mean of the across points and across flats 
                dimensions of the bolt head or the nut, whichever is smaller
            tp .. thickness of the plate under the bolt or nut
        Output:
            BpRd .. punching shear resistance [N]
    """
    BpRd = 0.6*math.pi*dm*tp*fu/gammaM2
    return BpRd

def shear_and_tension_resistance(FvEd,FvRd,FtEd,FtRd):
    """ EN 1993-1-8, Table 3.4
        Input:
            FvEd .. shear force in bolt
            FvRd .. shear resistance of bolt
            FtEd .. tension force in bolt
            FtRd .. tension resistance of bolt
        Output:
            U .. utilization ratio
    """
    U = FvEd/FvRd + FtEd/FtRd/1.4
    return U

def block_tearing(fy,fu,Ant,Anv,concentric_load=True,verb=False):
    """ Check block tearing resistance 
        fy .. yield strength of the plate [MPa]
        fu .. ultimate strength of the plate [MPa]
        Ant .. net area subject to tension [mm2]
        Anv .. net area subject to shear [mm2]
    """
    
    if concentric_load:
        Veff = fu*Ant/gammaM2 + fy/math.sqrt(3)*Anv/gammaM0
    else:
        Veff = 0.5*fu*Ant/gammaM2 + fy/math.sqrt(3)*Anv/gammaM0
    
    if verb:
        print("Block tearing:")
        print("fy = {0:4.2f} MPa, fu = {1:4.2f} MPa".format(fy,fu))
        print("Ant = {0:4.2f} mm^2, Anv = {1:4.2f} mm^2".format(Ant,Anv))
        print("Veff,Rd = {0:4.2f} kN".format(Veff*1e-3))
    
    return Veff

def stiffness_ratio(MjRatio,psi=2.7):
    """ Stiffness ratio (mu) of Eq. (6.28) """
    
    if MjRatio <= 2.0/3:
        mu = 1
    else:
        mu = (1.5*MjRatio)**psi
    
    return mu

def full_strength_weld_throat(mat=Steel("S355")):
    """ Throat thickness of double-sided fillet weld,
        when there is tension in the welded plate
    """
    
    beta = correlation_coefficient[mat.name]
    fy = mat.fy
    fu = mat.fu
    
    return beta/math.sqrt(2)*gammaM2/gammaM0*fy/fu

def full_strength_weld_tube(mat=Steel("S355")):
    """ Throat thickness of single-sided fillet weld,
        when there is tension in the welded plate
    """
    
    beta = correlation_coefficient[mat.name]
    fy = mat.fy
    fu = mat.fu
    
    return 2*beta/math.sqrt(2)*gammaM2/gammaM0*fy/fu
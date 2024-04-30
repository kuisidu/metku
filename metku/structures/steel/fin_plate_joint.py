#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2024. Metku team.
#  All rights reserved.



import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import os

from metku.sections.steel.ISection import HEA, HEB, IPE
from metku.structures.steel.plates import RectPlate, RectPlateWithHoles
from metku.eurocodes.en1993.constants import gammaM0, gammaM1, gammaM2
from metku.eurocodes.en1993.en1993_1_1 import buckling_reduction_factor, buckling_curve
from metku.eurocodes.en1993.en1993_1_8.en1993_1_8 import Bolt, BoltRow, BoltRowGroup, stiffness_ratio, full_strength_weld_throat
from metku.eurocodes.en1993.en1993_1_8.en1993_1_8 import nut_size, washer_size
import metku.eurocodes.en1993.en1993_1_8.component_method as cm

from metku.cost.workshop import Workshop
from metku.cost.cost_data import BASIC_STEEL_PRICE, steel_grade_add_on, thickness_add_on
from metku.cost.cost_data import bolt_unit_cost_full_thread as bolt_costs_ft
from metku.cost.cost_data import bolt_unit_cost_part_thread as bolt_costs_pt
from metku.cost.cost_data import nut_unit_cost, washer_unit_cost


class FinPlateJoint:
    """ Class for fin plate joints """
    def __init__(self, beam, plate, bolt=Bolt(20,8.8), gh=10,workshop=None):
        """
            beam .. beam profile (ISection)
            plate .. RectPlateWithHoles object
            bolt .. Bolt class object for bolts
            gh .. gap between beam end and support
            
            -- TWO of these three should be given --
            hp .. height of fin plate
            e1 .. vertical edge distance
            p1 .. vertical distance between bolt centers
            
            
            bp .. width of fin plate            
            z .. horizontal distance from support to bolt group center
            e2 .. horizontal edge distance
            p2 .. horizontal distance between bolt centers
            
        """
        self.beam = beam
        self.plate = plate
        # Vertical position of the center of the plate relative to
        # centerline of the beam.
        self.plate_pos_y = 0.0
        self.bolt = bolt
        self.plate.d0 = bolt.d0
        self.gh = gh
        self.workshop = workshop
        
        self.weld_size = 5.0 # throat thickness of a fillet weld, mm
        
        self._cost = {'plate':0.0,                                         
                     'welding':0.0,
                     'beam holes':0.0,
                     'bolts':0.0}
    
        self.vertical_symmetry = True
    
    @property
    def e1(self):
        """ Vertical edge distance to fin plate edge """
        if self.vertical_symmetry:
            # Here, we require that the vertical edge distances are the same
            # from the top and from the bottom of the plate
            self.plate.x0[1] = 0.5*self.hp - 0.5*(self.plate.n1-1)*self.p1
        return self.plate.eTop()
    
    @property
    def e2(self):
        """ Horizontal edge distance to fin plate edge """
        return self.plate.eRight()
    
    @property
    def e1b(self):
        """ Vertical edge distance to beam end edge """
        return 0.5*self.beam.h-(self.plate_pos_y+0.5*self.hp-self.e1)
    
    @property
    def e2b(self):
        """ Horizontal edge distance to beam end edge """
        return self.plate.eLeft()-self.gh
        
    @property
    def z(self):
        """ Horizontal distance from support to bolt group center """
        return self.plate.eLeft() + 0.5*(self.plate.n2-1)*self.plate.px
    
    @property
    def bp(self):
        return self.plate.b
    
    @bp.setter
    def bp(self,value):
        self.plate.b = value
    
    @property
    def hp(self):
        return self.plate.h

    @hp.setter 
    def hp(self,value):
        self.plate.h = value
    
    @property
    def tp(self):
        return self.plate.t
    
    @tp.setter
    def tp(self,value):
        self.plate.t = value
        
    @property
    def fyp(self):
        """ Yield strength of fin plate """
        return self.plate.fy
    
    @property
    def fup(self):
        """ Ultimate strength of fin plate """
        return self.plate.fu
    
    @property
    def d0(self):
        """ Bolt hole diameter """
        return self.bolt.d0
        
    @property
    def p1(self):
        """ Vertical distance between bolt hole centers """
        return self.plate.py

    @p1.setter
    def p1(self,value):
        """ Vertical distance between bolt hole centers """
        self.plate.py = value

    @property
    def p2(self):
        """ Horizontal distance between bolt hole centers """
        return self.plate.px
    
    @p2.setter
    def p2(self,value):
        """ Vertical distance between bolt hole centers """
        self.plate.px = value
    
    @property
    def I(self):
        """ Moment of inertia of bolt row group with two vertical rows """
        return 0.5*self.plate.n1*self.p2**2 + 1/6*self.plate.n1*(self.plate.n1**2-1)*self.p1**2
    
    def bolt_shear(self,verb=False):
        """ Shear resistance of bolts """
        n = self.plate.n
        z = self.z
        p1 = self.p1
        
        FvRd = self.bolt.shear_resistance()
        
        if self.plate.n2 == 1:
            C = n/np.sqrt(1+(6*z/(n+1)/p1)**2)
        else:
            p2 = self.p2
            Ibolts = self.I
            n1 = self.plate.n1
            C = 1/np.sqrt((0.5*z*p2/Ibolts+1/n)**2 + (0.5*z*p1/Ibolts*(n1-1))**2)
        
        VRd1 = C*FvRd
        
        if verb:
            print(f'VRd1 = {VRd1*1e-3:4.2f} kN')
        
        return VRd1
    
    def fin_plate_bearing(self,verb=False):
        """ Bearing resistance of fin plate """
        n1 = self.plate.n1
        if self.plate.n2 == 1:
            p2 = 1e5
            a = 0
            b = 6*self.z/n1/(n1+1)/self.p1  
        else:
            p2 = self.p2
            a = 0.5*self.z*p2/self.I
            b = 0.5*self.z*self.p1/self.I*(n1-1)
        
        ever = [self.e1,self.e2]        
        pver = [self.p1,p2]
        # NOTE! This bearing resistance distinguishes between edge and other bolts.
        # ECCS guide takes the minimum edge distance term using e and p regardless of
        # bolt position.
        FbverRd = self.bolt.bearing_resistance(self.fup, self.tp, ever, pver)
        
        
        ehor = [self.e2,self.e1]
        phor = [p2,self.p1]
        FbhorRd = self.bolt.bearing_resistance(self.fup, self.tp, ehor, phor)
        
                
        n = self.plate.n
        VRd2 = n/np.sqrt(((1+n*a)/FbverRd)**2+(n*b/FbhorRd)**2)
        
        if verb:
            print('Fin plate bearing resistance:')
            print(f'Fb,ver,Rd = {FbverRd*1e-3:4.2f} kN')
            print(f'Fb,hor,Rd = {FbhorRd*1e-3:4.2f} kN')
            print(f'     alpha = {a:4.2f}')
            print(f'     beta  = {b:4.2f}')
            print(f'     VRd2 = {VRd2*1e-3:4.2f} kN')
        
        return VRd2
    
    def fin_plate_shear_gross(self,verb=False):
        """ Shear resistance of fin plate """
        
        VRd3 = self.hp*self.tp*self.fyp/1.27/np.sqrt(3)/gammaM0

        if verb:
            print(f'Fin plate shear resistance, gross section: VRd3 = {VRd3*1e-3:4.2f} kN')
        
        return VRd3
    
    def fin_plate_shear_net(self,verb=False):
        """ Shear resistance of fin plate, net section """
        
        Avnet = self.tp*(self.hp-self.plate.n1*self.d0)
        VRd4 = Avnet*self.fup/np.sqrt(3)/gammaM2

        if verb:
            print(f'Fin plate shear resistance, net section: VRd4 = {VRd4*1e-3:4.2f} kN')
        
        return VRd4

    def fin_plate_block_tearing(self,verb=False):
        """ Shear resistance of fin plate, block tearing """
        
        Ant = self.tp*(self.e2+self.p2-(2*self.plate.n2-1)*0.5*self.d0)
        Anv = self.tp*(self.hp-self.e1-(self.plate.n1-0.5)*self.d0)
        VRd5 = 0.5*self.fup*Ant/gammaM2 + self.fyp*Anv/np.sqrt(3)/gammaM0

        if verb:
            print(f'Fin plate shear resistance, block tearing: VRd5 = {VRd5*1e-3:4.2f} kN')
        
        return VRd5
    
    def fin_plate_bending(self,verb=False):
        """ Fin plate in bending resistance """
        
        Wel = self.tp*self.hp**2/6
        VRd6 = self.fyp*Wel/self.z/gammaM0
    
        if verb:
            print(f'Fin plate in bending: VRd6 = {VRd6*1e-3:4.2f} kN')
            if self.hp > 2.73*self.z:
                print(f'Condition hp = {self.hp:4.2f} mm >= 2.73*z = {2.73*self.z:4.2f} mm is satisfied, plate bending not relevant.')
            else:
                print('Condition hp >= 2.73*z is not satisfied, plate bending is relevant.')            
                print(f'hp = {self.hp:4.2f} mm, 2.73*z = {2.73*self.z:4.2f} mm')
        
        return VRd6
    
    def fin_plate_ltb(self,verb=False):
        """ Lateral-torsional buckling of fin plate """
        
        VpRd = self.fyp*self.tp*self.hp**2/6/self.z
        
        if self.z > self.tp/0.15:      
            if verb:
                print('Fin plate LTB need not be checked.')
            C = 1.0
        else:
            if verb:
                print('Fin plate LTB needs to be checked.')
            l1 = np.pi*np.sqrt(self.plate.E/self.fyp)
            lLT = 2.8/l1*np.sqrt(self.z*self.hp/1.5/(self.tp**2))
            
            # !TODO: Check how chi_LT should really be determined!
            chiLT = buckling_reduction_factor(lLT,buckling_curve['d'])
            C =  min(chiLT/0.6/gammaM1,1/gammaM0)
            
        VbRd = C*VpRd
        
        if verb:
            print(f'Fin plate LTB resistance: VRd7 = {VbRd*1e-3:4.2f} kN')
            print(f'                          VpRd = {VpRd*1e-3:4.2f} kN')
            print(f'                            l1 = {l1:4.2f}')
            print(f'                           lLT = {lLT:4.2f}')
            print(f'                             C = {C:4.2f}')
        
        return VbRd
    
    def beam_web_bearing(self,verb=False):
        """ Beam web bearing resistance """
        
        n1 = self.plate.n1
        if self.plate.n2 == 1:
            p2 = 1e5
            a = 0
            b = 6*self.z/n1/(n1+1)/self.p1  
        else:
            p2 = self.p2
            a = 0.5*self.z*p2/self.I
            b = 0.5*self.z*self.p1/self.I*(n1-1)
        
        ever = [self.e1,self.e2b]        
        pver = [self.p1,p2]
        # NOTE! This bearing resistance distinguishes between edge and other bolts.
        # ECCS guide takes the minimum edge distance term using e and p regardless of
        # bolt position.
        FbverRd = self.bolt.bearing_resistance(self.beam.fu, self.beam.tw, ever, pver)
        
        
        ehor = [self.e2b,self.e1]
        phor = [p2,self.p1]
        FbhorRd = self.bolt.bearing_resistance(self.beam.fu, self.beam.tw, ehor, phor)
        
                
        n = self.plate.n
        VRd8 = n/np.sqrt(((1+n*a)/FbverRd)**2+(n*b/FbhorRd)**2)
        
        if verb:
            print('Beam web bearing resistance:')
            print(f'Fb,ver,Rd = {FbverRd*1e-3:4.2f} kN')
            print(f'Fb,hor,Rd = {FbhorRd*1e-3:4.2f} kN')
            print(f'     alpha = {a:4.2f}')
            print(f'     beta  = {b:4.2f}')
            print(f'     VRd8 = {VRd8*1e-3:4.2f} kN')
        
        return VRd8
    
    def beam_web_shear_gross(self,verb=False):
        """ Shear resistance of beam web """
        
        VRd9 = self.beam.VzRd

        if verb:
            print(f'Beam shear resistance, gross section: VRd9 = {VRd9*1e-3:4.2f} kN')
        
        return VRd9
    
    def beam_web_shear_net(self,verb=False):
        """ Shear resistance of beam web, net section """
        
        Abvnet = self.beam.Ashear[1] - self.plate.n1*self.d0*self.beam.tw
        VRd10 = Abvnet*self.beam.fu/np.sqrt(3)/gammaM2

        if verb:
            print(f'Beam web shear resistance, net section: VRd10 = {VRd10*1e-3:4.2f} kN')
        
        return VRd10
    
    def beam_web_block_tearing(self,verb=False):
        """ Shear resistance of beam web, block tearing """
        
        Ant = self.beam.tw*(self.e2b+self.p2-(2*self.plate.n2-1)*0.5*self.d0)
        Anv = self.beam.tw*(self.e1b+(self.plate.n1-1)*self.p1-(self.plate.n1-0.5)*self.d0)
        VRd11 = 0.5*self.beam.fu*Ant/gammaM2 + self.beam.fy*Anv/np.sqrt(3)/gammaM0

        if verb:
            print(f'Beam shear resistance, block tearing: VRd11 = {VRd11*1e-3:4.2f} kN')
        
        return VRd11
    
    def draw(self,axes=None):
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        h = self.beam.h
        tf = self.beam.tf
        bp = self.bp
        hp = self.hp
        gh = self.gh
        dx = 100
        dy = 20
        
        ax.vlines([0],-0.5*h-dy,0.5*h+dy,'k')
        
        ax.hlines([0.5*h,0.5*h-tf,-0.5*h,-0.5*h+tf],gh,gh+bp+dx,'k')
        ax.vlines([gh],0.5*h-tf,0.5*h,'k')
        ax.vlines([gh],-0.5*h,-0.5*h+tf,'k')
        
        # Origin of the coordinate system is at the centerline of the beam
        # and at the face of the support.
        x0 = (0,self.plate_pos_y-0.5*hp)
        self.plate.draw(x0,ax)
        
        ax.vlines([gh],0.5*hp,0.5*h-tf,'k')
        ax.vlines([gh],-0.5*h+tf,-0.5*hp,'k')
    
        plt.rcParams['hatch.color'] = 'k'
        hatch_rect = patches.Rectangle((-30,-0.5*h-dy), 30, h+2*dy, linewidth=0,hatch='/',edgecolor='k',facecolor='white')
        ax.add_patch(hatch_rect)
    
    def check_bolt_position(self,verb=False):
        """ Checks the bolt positioning rules of EN 1993-1-8 """
        res = True
        
        if self.e1 < 1.2*self.d0:
            if verb:
                print(f'Condition e1 >= 1.2*d0 is not satisfied: e1 = {self.e1:4.1g}, 1.2*d0 = {1.2*self.d0:4.1g}')
            res = False
        elif self.e2 < 1.2*self.d0:
            res = False
        elif self.e2b < 1.2*self.d0:
            res = False
        elif self.p1 < 2.2*self.d0:
            res = False
        elif self.p2 > 0 and self.p2 < 2.4*self.d0:
            res = False
        
        if res and verb:
            print('All bolt positioning conditions are OK!')
            
        return res
    
    def cost(self,workshop=None,material_cost=BASIC_STEEL_PRICE,verb=False):
        """ Calculate cost of the joint
            units: euro
            
            Total costs consist of:
                plate: including material, cutting, hole forming, painting
                welding: fin plate is welded to the column by two-sided fillet weld.
                beam holes: holes are drilled to the beam web
                bolts: material costs of bolts
                
            !TODO:
                Check beam hole forming, is it cutting, drilling or what
                Check welding, is it assembly welding or part welding.
        """
        
        for key, value in self._cost.items():
            self._cost[key] = 0.0
        
        if workshop is not None:
            self.workshop = workshop        
        
        plate_cost = self.workshop.steel_price["plates"]
        
        self._cost['plate'] = self.plate.cost(self.workshop,plate_cost)
        
        # Bolt costs: bolt assemply includes
        # 1) Bolt
        # 2) Nut
        # 3) One washer
        if self.bolt.bolt_type == "Full thread":
            self._cost['bolts'] += self.plate.n * bolt_costs_ft[int(self.bolt.length)][int(self.bolt.d)]
        else:
            self._cost['bolts'] += self.plate.n * bolt_costs_pt[int(self.bolt.length)][int(self.bolt.d)]
            
        self._cost['bolts'] += self.plate.n * nut_unit_cost[int(self.bolt.d)]
        self._cost['bolts'] += self.plate.n * washer_unit_cost[int(self.bolt.d)]
        
        hole_cut_length = self.plate.n * self.d0*math.pi
        
        self._cost['beam holes'] = self.workshop.cost_centres['cutting'].cost(self.beam.tw,hole_cut_length)

        """ Welding cost """
        plate_weld_length = 2*self.hp
        #print(flange_weld_length)
        self._cost['welding'] += self.workshop.cost_centres['assembly_welding'].cost(self.weld_size,plate_weld_length)
        
            
        total_cost = 0.0
        
        for c in self._cost.values():            
            total_cost += c
        
        self._total_cost = total_cost
        
        #self.cost_distribution()
        
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
    

def example():
    
    beam = IPE(300,fy=235)
    hp = 230
    e1 = 45
    n1 = 3
    n2 = 1
    b = Bolt(20,8.8)
    
    plate = RectPlateWithHoles(width=110, depth=230, thickness=10.0, d0=22, x0=[60,45],py=70,n2=1,material="S235")
    
    conn = FinPlateJoint(beam, plate, bolt=b)
    
    verb = True
    
    conn.check_bolt_position(verb)
    
    conn.bolt_shear(verb)
    conn.fin_plate_bearing(verb)
    conn.fin_plate_shear_gross(verb)
    conn.fin_plate_shear_net(verb)
    conn.fin_plate_block_tearing(verb)
    conn.fin_plate_bending(verb)
    conn.fin_plate_ltb(verb)
    conn.beam_web_bearing(verb)
    conn.beam_web_shear_gross(verb)
    conn.beam_web_shear_net(verb)
    conn.beam_web_block_tearing(verb)

    conn.workshop = Workshop()
    conn.cost(verb=verb)
    conn.draw()

    return conn


class FinPlateConnection:
    def __init__(self, plate_thickness, plate_width, plate_length, bolt_diameter, number_of_bolts):
        self.plate_thickness = plate_thickness
        self.plate_width = plate_width
        self.plate_length = plate_length
        self.bolt_diameter = bolt_diameter
        self.number_of_bolts = number_of_bolts

    def calculate_connection_details(self):
        # Implement the calculations to determine the connection details

        
        # just plot the plate dimensions
        pass



    def plot_connection_diagram(self):
        fig, ax = plt.subplots()
        ax.set_aspect('equal', adjustable='box')

        # Plot the fin plate
        ax.plot([0, self.plate_length], [0, 0], 'k-', linewidth=2)  # Bottom edge
        ax.plot([0, 0], [0, self.plate_width], 'k-', linewidth=2)   # Left edge
        ax.plot([self.plate_length, self.plate_length], [0, self.plate_width], 'k-', linewidth=2)  # Right edge
        ax.plot([0, self.plate_length], [self.plate_width, self.plate_width], 'k-', linewidth=2)  # Top edge

        # Plot the bolts (for simplicity, let's assume they're evenly spaced)
        bolt_spacing = self.plate_length / (self.number_of_bolts - 1)
        bolt_positions = [i * bolt_spacing for i in range(self.number_of_bolts)]
        for pos in bolt_positions:
            ax.plot([pos], [self.plate_width / 2], 'ro', markersize=8)  # Bolt position

        ax.set_xlabel('Length (mm)')
        ax.set_ylabel('Width (mm)')
        ax.set_title('Fin Plate Connection Diagram')

        plt.grid(True)
        plt.show()



# Example usage:
if __name__ == "__main__":
    plate_thickness = 10  # in mm
    plate_width = 150     # in mm
    plate_length = 300    # in mm
    bolt_diameter = 20    # in mm
    number_of_bolts = 3

    #connection = FinPlateConnection(plate_thickness, plate_width, plate_length, bolt_diameter, number_of_bolts)
    #connection.plot_connection_diagram()

    conn = example()    



















# Design Checks to be done in Fin plate connections

# Shear capacity of the bolt


#Tension capacity of the bolt


#Bearing capacity of the bolt


#Block shear capacity of the bolt


#Combined shear and tension capacity


#Weld check
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 23:12:31 2018

Rectangular hollow sections

@author: kmela
"""

import math
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines


from eurocodes.en1993 import en1993_1_1
from sections.steel.steel_section import SteelSection


class RHS(SteelSection):
    """ Rectangular hollow sections """

    def __init__(self, H, B, T, fy=355):
        """ Rectangular hollow sections
            
            H -- height
            B -- width
            T -- wall thickness
            fy -- yield strength
        """
        R = RHS_outer_radius(T)
        A = RHS_area(H, B, T, R)
        I = [0.0, 0.0]
        I[0] = RHS_second_moment(H, B, T, R)
        I[1] = RHS_second_moment(B, H, T, R)
        Ashear = [RHS_shear_area(A,B,H),RHS_shear_area(A, H, B)]
        
        Au = RHS_paint_area(H, B, R)
        Wel = [0.0, 0.0]
        Wel[0] = RHS_elastic_modulus(H, I[0])
        Wel[1] = RHS_elastic_modulus(B, I[1])
        Wpl = [0.0, 0.0]
        Wpl[0] = RHS_plastic_modulus(H, B, R, T)
        Wpl[1] = RHS_plastic_modulus(B, H, R, T)

        super().__init__(fy,
                         A,
                         I,
                         Au,
                         Wpl,
                         Wel,
                         Ashear)

        self.H = H
        self.B = B
        self.T = T
        self.Iw = 0.0
        self.It = self.torsional_constant()
        self.imp_factor = [en1993_1_1.buckling_curve["c"],
                           en1993_1_1.buckling_curve["c"]]
        self.imp_factor_LT_gen = en1993_1_1.buckling_curve["c"]  # TODO
        self.imp_factor_LT = en1993_1_1.buckling_curve["c"]  # TODO
        
        
    

    @property
    def h(self):
        return self.H
    
    @property
    def b(self):
        return self.B

    @property
    def t(self):
        return self.T

    @property
    def r(self):
        return self.R
    
    def __repr__(self):
        return f"{type(self).__name__} {self.H:.0f}X{self.B:.0f}X{self.T:.1f}"
    
    def robot(self):
        
        return f"RRHS {self.H:.0f}x{self.B:.0f}x{self.T:.1f}"

    def info(self,latex=False):
        """ Prints a bunch of section properties """
        
        if latex:
            print("*** " + self.__repr__() + " ***")
            print(" $h = {0:.0f}{1:s}$".format(self.H,'\\,\\unit{mm}'))
            print(" $b = {0:.0f}{1:s}$".format(self.B,'\\,\\unit{mm}'))
            print(" $t = {0:.1f}{1:s}$".format(self.t,'\\,\\unit{mm}'))            
            print(" $r = {0:.1f}{1:s}$".format(self.r,'\\,\\unit{mm}'))
            print(" $A = {0:.2f}{1:s}$".format(self.A*1e-2,'\\ee{2}\\,\\squnit{mm}'))
            print(" $I_y = {0:.2f}{1:s}$".format(self.I[0]*1e-4,'\\ee{4}\\,\\quunit{mm}'))
            print(" $I_z = {0:.2f}{1:s}$".format(self.I[1]*1e-4,'\\ee{4}\\,\\quunit{mm}'))
            print(" $W_el,y = {0:.2f}{1:s}$".format(self.Wel[0]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))
            print(" $W_el,z = {0:.2f}{1:s}$".format(self.Wel[1]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))
            print(" $W_pl,y = {0:.2f}{1:s}$".format(self.Wpl[0]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))
            print(" $W_pl,z = {0:.2f}{1:s}$".format(self.Wpl[1]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))        
            print(" $Ashear = {0:.2f}{1:s}$".format(self.Ashear*1e-2,'\\ee{2}\\,\\squnit{mm}'))
        else:
            print("*** " + self.__repr__() + " ***")
            print(" h = {0:.0f} mm".format(self.h))
            print(" b = {0:.0f} mm".format(self.b))            
            print(" t = {0:.1f} mm".format(self.t))
            print(" r = {0:.1f} mm".format(self.r))
            print(" A = {0:.2f} (100 mm2)".format(self.A*1e-2))
            print(" Iy = {0:.2f} (10^4 mm4)".format(self.I[0]*1e-4))
            print(" Iz = {0:.2f} (10^4 mm4)".format(self.I[1]*1e-4))
            print(" Wel,y = {0:.2f} (10^3 mm4)".format(self.Wel[0]*1e-3))
            print(" Wel,z = {0:.2f} (10^3 mm4)".format(self.Wel[1]*1e-3))
            print(" Wpl,y = {0:.2f} (10^3 mm4)".format(self.Wpl[0]*1e-3))
            print(" Wpl,z = {0:.2f} (10^3 mm4)".format(self.Wpl[1]*1e-3))
            print(" Ashear = {0:.2f} (100 mm2)".format(self.Ashear*1e-2))

    def update_properties(self):
        """
        Updates cross-sectional properties
        """

        try:
            # Using values from catalog, all catalog values
            # have maximum of 1 decimal
            # if (self.H * 10).is_integer() and (self.B * 10).is_integer() and (self.T * 10).is_integer():
            self.R = RHS_outer_radius(self.T)
            # else:
            # self.R = 2 * self.T
            self.A = RHS_area(self.H, self.B, self.T, self.R)
            I = [0.0, 0.0]
            I[0] = RHS_second_moment(self.H, self.B, self.T, self.R)
            I[1] = RHS_second_moment(self.B, self.H, self.T, self.R)
            self.I = np.asarray(I)
            self.Ashear = [RHS_shear_area(self.A, self.B, self.B),RHS_shear_area(self.A, self.H, self.B)]
            self.Au = RHS_paint_area(self.H, self.B, self.R)
            Wel = [0.0, 0.0]
            Wel[0] = RHS_elastic_modulus(self.H, I[0])
            Wel[1] = RHS_elastic_modulus(self.B, I[1])
            self.Wel = np.asarray(Wel)
            Wpl = [0.0, 0.0]
            Wpl[0] = RHS_plastic_modulus(self.H, self.B, self.R, self.T)
            Wpl[1] = RHS_plastic_modulus(self.B, self.H, self.R, self.T)
            self.Wpl = np.asarray(Wpl)
            self.It = self.torsional_constant()

        except AttributeError:
            pass

    def __setattr__(self, key, value):

        dim_vars = ['H', 'B', 'T']
        super().__setattr__(key, value)
        # Whenever dimensions change, update cross-sectional properties
        if key in dim_vars:
            self.update_properties()

    @property
    def c_web(self):
        return self.H - 2 * self.R

    @property
    def c_flange(self):
        return self.B - 2 * self.R

    def torsional_constant(self):
        """
        Calculated according to EN 10219-2
        
        """

        # Rc .. corner radius at center line
        Rc = self.R - 0.5 * self.T
        h = 2*((self.B-self.T) + (self.H-self.T)) - 2*Rc*(4-math.pi)
        Ah = (self.B-self.T) * (self.H-self.T) - Rc**2*(4-math.pi)
        K = 2*Ah*self.T/h
        It = (self.T**3*h/3 + 2*K*Ah)

        return It


    def flange_class(self, verb=False):
        """ Determine class of compressed flange """
        cf = self.B - 2 * self.R
        rf = cf / self.T
        cFlange = en1993_1_1.internal_part_in_compression(rf, self.eps)

        return cFlange

    def web_class_comp(self, verb=False):
        """ Determine class of compressed web """
        cw = self.H - 2 * self.R
        rw = cw / self.T
        cWeb = en1993_1_1.internal_part_in_compression(rw, self.eps)

        return cWeb

    def web_class_bend(self, verb=False):
        """ Determine class of web in bending """
        cw = self.H - 2 * self.R
        rw = cw / self.T
        cWeb = en1993_1_1.internal_part_in_bending(rw, self.eps)

        return cWeb

    def web_class_comp_bend(self, Ned, verb=False):

        return 2

    def moment_axial_force_interact(self, UN, MRd):
        """ Interaction rule for section resistance for combined
            axial force and bending
            
            input: UN .. NEd/NRd
                  MRd .. moment resistance
        """
        aw = min((self.A - 2 * self.B * self.T) / self.A, 0.5)

        if UN > 1.0:
            MNRd = 0.0
        else:
            MNRd = MRd * min((1 - UN) / (1 - 0.5 * aw), 1)

        return MNRd
    
    def draw(self,axes=None,theta=0):
        """ Draw the profile 
            input:
                axes .. matplotlib.axes.Axes object. If this is given,
                        the profile will be drawn to that Axes
                theta .. rotation (optional)
        """
        
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        if theta is not 0:
            phi = np.radians(theta)
            c, s = np.cos(phi), np.sin(phi)
            R = np.array(((c,-s), (s, c)))
        
        lw = 1.5
        col = 'b'
        
        h = self.H
        b = self.B
        t = self.T
        r = self.R
        
        # Draw webs
        y_web = [-0.5*h+r,0.5*h-r]                            
        ax.vlines([-0.5*b,-0.5*b+t,0.5*b-t,0.5*b],y_web[0],y_web[1],colors='b',linewidth=lw)
        
        # Draw flanges
        x_flange = [-0.5*b+r,0.5*b-r]
        y_top = [0.5*h,0.5*h-t]
        y_bottom = [-0.5*h,-0.5*h+t]
        y_flanges = y_top + y_bottom
        
        ax.hlines(y_flanges,x_flange[0],x_flange[1],colors='b',linewidth=lw)
        
        # Draw corners
        # Upper left corner
        x_up_left = (-0.5*b+r,0.5*h-r)
        
        up_left_corner_in = patches.Arc(x_up_left,2*(r-t),2*(r-t),0,90,180,color=col,linewidth=lw)
        up_left_corner_out = patches.Arc(x_up_left,2*r,2*r,0,90,180,color=col,linewidth=lw)
        
        ax.add_patch(up_left_corner_in)
        ax.add_patch(up_left_corner_out)
        
        x_bot_left = (-0.5*b+r,-0.5*h+r)
        
        bot_left_corner_in = patches.Arc(x_bot_left,2*(r-t),2*(r-t),0,180,270,color=col,linewidth=lw)
        bot_left_corner_out = patches.Arc(x_bot_left,2*r,2*r,0,180,270,color=col,linewidth=lw)
        
        ax.add_patch(bot_left_corner_in)
        ax.add_patch(bot_left_corner_out)
                            
        # Upper left corner
        x_up_right = (0.5*b-r,0.5*h-r)
        
        up_right_corner_in = patches.Arc(x_up_right,2*(r-t),2*(r-t),0,0,90,color=col,linewidth=lw)
        up_right_corner_out = patches.Arc(x_up_right,2*r,2*r,0,0,90,color=col,linewidth=lw)
        
        ax.add_patch(up_right_corner_in)
        ax.add_patch(up_right_corner_out)
        
        
        x_bot_right = (0.5*b-r,-0.5*h+r)
        
        bot_right_corner_in = patches.Arc(x_bot_right,2*(r-t),2*(r-t),0,270,360,color=col,linewidth=lw)
        bot_right_corner_out = patches.Arc(x_bot_right,2*r,2*r,0,270,360,color=col,linewidth=lw)
        
        ax.add_patch(bot_right_corner_in)
        ax.add_patch(bot_right_corner_out)
        
        
        a = 0.1
        
        ax.set_xlim(-(0.5+a)*b, (0.5+a)*b)
        ax.set_ylim(-(0.5+a)*h, (0.5+a)*h)
        
        ax.set_aspect('equal')
    
    def abaqus(self,filename,matname='Steel',setname='Set-1'):
        """ Writes section data for exporting to abaqus """
        unit_conv = 1e-3
        
        # Calculate coordinates for abaqus cross-section
        W = self.B*unit_conv
        H = self.H*unit_conv
        T = self.T*unit_conv
        NomT = self.T*unit_conv
        R = self.R*unit_conv
        
        RMid=R-T/2
            
        angle=math.radians(22.5)
        
        X=[]
        Y=[]
        
        X.append(0.0)
        Y.append((H-T)/2)
        
        X.append(W/2-R-2*T)
        Y.append((H-T)/2)  
        
        X.append(W/2-R)
        Y.append((H-T)/2)
        
        X.append(W/2-R+RMid*math.sin(angle))
        Y.append((H-T)/2-(RMid-RMid*math.cos(angle)))
        
        X.append(W/2-R+RMid*math.sin(2*angle))
        Y.append((H-T)/2-(RMid-RMid*math.cos(2*angle))) 
        
        X.append(W/2-R+RMid*math.sin(3*angle))
        Y.append((H-T)/2-(RMid-RMid*math.cos(3*angle))) 
        
        X.append((W-T)/2)
        Y.append(H/2-R) 
        
        X.append((W-T)/2)
        Y.append(H/2-R-2*T) 
        
        X.append((W-T)/2)
        Y.append(0.0)
        

        X += X[7::-1]
        Y += [-y for y in Y[7::-1]]
        
        X += [-x for x in X[1:9]]
        Y += [-y for y in Y[1:9]]
        
        X += [-x for x in X[7:0:-1]]
        Y += Y[7:0:-1]
        
        X.append(X[0])
        Y.append(Y[0])
        
        # Write coordinates to file
        MatName = matname + "_RedMat"
        nx = len(X)
        with open(filename,'w') as file:
            #file.write(f"** Member {mem_name}\n")
            file.write("** Cross-section dimensions \n")
            file.write(f"** WIDTH = {W}\n")
            file.write(f"** HEIGHT = {H}\n")
            file.write(f"** THICKNESS = {T}\n")
            file.write(f"** CORNER RADIUS = {R}\n**\n")
            file.write(f"*Beam Section, elset={setname}, material={MatName}, section=ARBITRARY\n")
            file.write(f'{nx-1}, {X[0]:8.6f},{Y[0]:8.6f},{X[1]:8.6f},{Y[1]:8.6f},{T:8.6f}\n')                 
                   
            for x, y in zip(X[2:],Y[2:]):
                file.write(f"{x:8.6f},{y:8.6f},{T:8.6f}\n")
        
            file.write("0.,1.,0.\n")
        

class SHS(RHS):
    """ Square hollow sections
        Subclass of RHS
    """

    def __init__(self, H, T, fy=355):
        RHS.__init__(self, H, H, T, fy)

    def __repr__(self):
        return f"{type(self).__name__} {self.H:.0f}X{self.B:.0f}X{self.T:.1f}"


# Functions for computing RHS section properties
def RHS_area(H, B, T, R):
    """ cross-sectional area """
    A = 2 * (H - 2 * R) * T + 2 * (B - 2 * R) * T + math.pi * (
                R ** 2 - (R - T) ** 2)
    return A


def RHS_paint_area(H, B, R):
    """ circumference """
    Au = 2 * (H + B) - 8 * R + 2 * math.pi * R
    return Au


def RHS_shear_area(A, H, B):
    Ashear = A * H / (H + B)
    return Ashear


def RHS_second_moment(H, B, T, R):
    """ Second moment of area """
    I1 = 2 * (1 / 12 * T * (H - 2 * R) ** 3 + 1 / 12 * (B - 2 * R) * T ** 3 + (
                0.5 * (H - T)) ** 2 * (B - 2 * R) * T)
    IR = 4 * ((math.pi / 16 - 4 / 9 / math.pi) * R ** 4 + (
                0.5 * H - R + 4 * R / 3 / math.pi) ** 2 * math.pi * R ** 2 / 4)
    IRneg = 4 * ((math.pi / 16 - 4 / 9 / math.pi) * (R - T) ** 4 + (
                0.5 * H - R + 4 * (R - T) / 3 / math.pi) ** 2 * math.pi * (
                             R - T) ** 2 / 4)
    I1 = I1 + IR - IRneg

    return I1


def RHS_elastic_modulus(H, I1):
    Wel = I1 / (0.5 * H)

    return Wel


def RHS_plastic_modulus(H, B, R, T):
    Wp = 2 * ((0.5 * H - R) * T * (0.5 * H - R)) + (B - 2 * R) * T * (H - T)
    WR = 2 * math.pi / 4 * R ** 2 * 2 * (0.5 * H - R + 4 * R / 3 / math.pi)
    Wneg = 2 * math.pi / 4 * (R - T) ** 2 * 2 * (
                0.5 * H - R + 4 * (R - T) / 3 / math.pi)
    Wp = Wp + WR - Wneg

    return Wp


def RHS_outer_radius(T):
    if T <= 6.0:
        R = 2 * T
    elif T <= 10.0:
        R = 2.5 * T
    else:
        R = 3 * T
    return R

if __name__ == "__main__":
    
    p = RHS(250,150,10)
    #p.draw()

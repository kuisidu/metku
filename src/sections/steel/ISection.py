# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 23:12:31 2018

Rectangular hollow sections

@author: kmela
"""
import math

from src.eurocodes.en1993 import en1993_1_1, constants
from src.sections.steel.steel_section import SteelSection
from src.sections.steel.catalogue import profile


class ISection(SteelSection):
    """ Base class for I -sections """

    def __init__(self, H, B, tf, tw, r, fy=355):
        """
            H -- height
            B -- width
            T -- wall thickness
            fy -- yield strength
        """
        self.h = H
        self.b = B
        self.tf = tf
        self.tw = tw
        self.r = r
        self.fy = fy

        A, Au, Ashear, I, It, Iw, Wpl, Wel = cross_section_properties(B, H, tf, tw, r)
        self.Iw = Iw
        self.It = It
        super().__init__(fy, A, I, Au, Wpl, Wel, Ashear)

        """ Determine buckling curve: EN 1993-1-1, Table 6.2 """
        if H / B > 1.2:
            if tf <= 40:
                if fy <= 420:
                    buck_y = "a"
                    buck_z = "b"
                else:
                    buck_y = "a0"
                    buck_z = "a0"
            else:
                if fy <= 420:
                    buck_y = "b"
                    buck_z = "c"
                else:
                    buck_y = "a"
                    buck_z = "a"
        else:
            if tf <= 100:
                if fy <= 420:
                    buck_y = "b"
                    buck_z = "c"
                else:
                    buck_y = "a"
                    buck_z = "a"
            else:
                if fy <= 420:
                    buck_y = "d"
                    buck_z = "d"
                else:
                    buck_y = "c"
                    buck_z = "c"

        self.imp_factor = [en1993_1_1.buckling_curve[buck_y],
                           en1993_1_1.buckling_curve[buck_z]]

    @property
    def hw(self):
        """ Straight part of the web """
        return self.h - 2 * self.tf - 2 * self.r

    @property
    def cf(self):
        """ Straight part of the flange """
        return 0.5 * (self.b - self.tw) - self.r

    def flange_class(self,verb=False):
        """ Determine class of compressed flange """
        # cf = 0.5*(self.b - self.tw) - self.r
        rf = self.cf / self.tf
        cFlange = en1993_1_1.outstand_part_in_compression(rf, self.eps)

        if verb:
            print("Flange classification (outstand element):")
            print("cf = {0:4.2f}, tf = {1:4.2f}".format(self.cf,self.tf))
            print("cf/tf = {0:4.2f}".format(rf))
            print("Flange class = {0}".format(cFlange))
            

        return cFlange

    def web_class_comp(self,verb=False):
        """ Determine class of compressed web """
        rw = self.hw / self.tw
        cWeb = en1993_1_1.internal_part_in_compression(rw, self.eps)

        if verb:
            print("Web classification (internal part in compression):")
            print("hw = {0:4.2f}, tw = {1:4.2f}".format(self.hw,self.tw))
            print("hw/tw = {0:4.2f}".format(rw))
            print("Web class = {0}".format(cWeb))

        return cWeb

    def web_class_bend(self,verb=False):
        """ Determine class of web in bending """
        rw = self.hw / self.tw
        cWeb = en1993_1_1.internal_part_in_bending(rw, self.eps)

        if verb:
            print("Web classification (internal part in bending):")
            print("hw = {0:4.2f}, tw = {1:4.2f}".format(self.hw,self.tw))
            print("hw/tw = {0:4.2f}".format(rw))
            print("Web class = {0}".format(cWeb))

        return cWeb

    def web_class_comp_bend(self, Ned, verb=False):
        """ Determine class of web in combined bending and compression 
            TODO
        """
        # cw = self.h-2*self.tf-2*self.r
        rw = self.hw / self.tw
        Ac = 0.5 * (self.A + Ned / (self.fy / constants.gammaM0))
        a = Ac / self.A
        p = -1

        cWeb = en1993_1_1.internal_part_comp_bend(rw, self.eps, a, p)
        
        if verb:
            print("Web classification (internal part in compression and bending):")
            print("hw = {0:4.2f}, tw = {1:4.2f}".format(self.hw,self.tw))
            print("hw/tw = {0:4.2f}".format(rw))
            print("Web class = {0}".format(cWeb))

        return cWeb

    def moment_axial_force_interact(self, UN, MRd):
        """ Interaction rule for section resistance for combined
            axial force and bending
            
            input: UN .. NEd/NRd
                  MRd .. moment resistance                 
        """
        aw = min((self.A - 2 * self.b * self.tf) / self.A, 0.5)

        if UN > 1.0:
            MNRd = 0.0
        else:
            MNRd = MRd * min((1 - UN) / (1 - 0.5 * aw), 1)

        return MNRd


class IPE(ISection):
    """ IPE sections
        Subclass of ISection
    """

    def __init__(self, height=100, fy=355, catalogue=False):
        name = 'IPE ' + str(height)
        # H,B,tf,tb,r,fy=355
        if catalogue or name in profile.keys():
            try:
                H = profile[name]['h']
                B = profile[name]['b']
                tf = profile[name]['t_f']
                tw = profile[name]['t_w']
                r = profile[name]['r']
            except KeyError:
                print('Error: No profile named ' + name + ' in the catalogue.')
        else:
            H = height
            K = [0.0256, 3.2831, 0.0155, 2.4921, 0.3486, 30.178, 0.0155, 2.4921]
            tf = K[0] * H + K[1]
            tw = K[2] * H + K[3]
            B = K[4] * H + K[5]
            r = K[6] * H + K[7]

        ISection.__init__(self, H, B, tf, tw, r, fy=fy)


class HEA(ISection):
    """ European wide flange sections
        Subclass of ISection        
    """

    def __init__(self, height=100, fy=355, catalogue=False):
        name = 'HE ' + str(height) + ' ' + 'A'
        # HEIGHT =! TO PROFILE NUMBER

        # H,B,tf,tb,r,fy=355
        if catalogue:
            H = profile[name]['h']
            B = profile[name]['b']
            tf = profile[name]['t_f']
            tw = profile[name]['t_w']
            r = profile[name]['r']
        else:
            H = height
            K = [0.0294, 5.7651, 0.014, 4.2949, 1.0323, 2.5368, 0.0781, 2.9206]
            tf = K[0] * H + K[1]
            tw = K[2] * H + K[3]

            if H <= 270:
                B = K[4] * H + K[5]
                r = K[6] * H + K[7]
            else:
                B = 300
                if H <= 690:
                    r = 27
                else:
                    r = 30

        ISection.__init__(self, H, B, tf, tw, r, fy=fy)


class HEAA(ISection):
    """ European wide flange sections
        Subclass of ISection        
    """

    def __init__(self, height=100, fy=355, catalogue=False):
        name = 'HE ' + str(height) + ' ' + 'AA'
        # HEIGHT =! TO PROFILE NUMBER

        # H,B,tf,tb,r,fy=355
        if catalogue:
            H = profile[name]['h']
            B = profile[name]['b']
            tf = profile[name]['t_f']
            tw = profile[name]['t_w']
            r = profile[name]['r']
        else:
            H = height
            K = [0.0179, 4.9342, 0.0144, 3.2122, 1.0403, 6.2902, 0.0786, 3.221]
            tf = K[0] * H + K[1]
            tw = K[2] * H + K[3]

            if H <= 264:
                B = K[4] * H + K[5]
                r = K[6] * H + K[7]
            else:
                B = 300
                if H <= 670:
                    r = 27
                else:
                    r = 30

        ISection.__init__(self, H, B, tf, tw, r, fy=fy)


class HEB(ISection):
    """ European wide flange sections
        Subclass of ISection        
    """
    def __init__(self,height=100,fy=355,catalogue=True):
        name = 'HE ' + str(height) + ' ' + 'B'
        # H,B,tf,tb,r,fy=355
        if name in profile.keys():
            H=profile[name]['h']
            B=profile[name]['b']
            tf = profile[name]['t_f']
            tw = profile[name]['t_w']
            r = profile[name]['r']
        else:
            H = height
            K = [0.0308, 9.5736, 0.0147, 6.2014, 1, 0, 0.0755, 2.7636]
            tf = K[0] * H + K[1]
            tw = K[2] * H + K[3]

            if H <= 280:
                B = K[4] * H + K[5]
                r = K[6] * H + K[7]
            else:
                B = 300
                if H <= 700:
                    r = 27
                else:
                    r = 30

        ISection.__init__(self, H, B, tf, tw, r, fy=fy)


class HEC(ISection):
    """ HEC sections
        Subclass of ISection
    """

    def __init__(self, height=100, fy=355, catalogue=False):
        name = 'HE ' + str(height) + ' ' + 'C'
        # H,B,tf,tb,r,fy=355
        if catalogue:
            H = profile[name]['h']
            B = profile[name]['b']
            tf = profile[name]['t_f']
            tw = profile[name]['t_w']
            r = profile[name]['r']
        else:
            H = height
            K = [0.0669, 6.7987, 0.0317, 5.4267, 0.9256, 5.593, 0.0747, 2.0905]
            tf = K[0] * H + K[1]
            tw = K[2] * H + K[3]
            B = K[4] * H + K[5]
            r = K[6] * H + K[7]

        ISection.__init__(self, H, B, tf, tw, r, fy=fy)


class HEM(ISection):
    """ European wide flange sections
        Subclass of ISection        
    """

    def __init__(self, height=100, fy=355, catalogue=False):
        name = 'HE ' + str(height) + ' ' + 'M'
        # HEIGHT =! TO PROFILE NUMBER

        # H,B,tf,tb,r,fy=355
        if catalogue:
            H = profile[name]['h']
            B = profile[name]['b']
            tf = profile[name]['t_f']
            tw = profile[name]['t_w']
            r = profile[name]['r']
        else:
            H = height
            K = [0.0824, 8.5346, 0.0392, 6.8259, 0.9344, -2.8966, 0.0732, 1.5707]

            if H <= 340:
                tf = K[0] * H + K[1]
                tw = K[2] * H + K[3]
                B = K[4] * H + K[5]
                r = K[6] * H + K[7]
            else:
                tf = 40
                tw = 21
                B = -0.0112 * H + 312.41
                if H <= 716:
                    r = 27
                else:
                    r = 30

        ISection.__init__(self, H, B, tf, tw, r, fy=fy)


class CustomISection(ISection):
    def __init__(self, H, B, tf, tw, r, fy=355):
        self.H = H
        self.B = B
        self.tf = tf
        self.tw = tw
        self.r = r

        ISection.__init__(self, H, B, tf, tw, r, fy=fy)

    def cross_section_properties(self):
        args = cross_section_properties(self.B, self.H, self.tf, self.tw, self.r)
        self.A, self.Au, self.Ashear, self.I, self.It, self.Iw, self.Wpl, self.Wel = args


def cross_section_properties(b, h, tf, tw, r):
    """ Cross-section properties        
        Calculation according to http://sections.arcelormittal.com/fileadmin/redaction/4-Library/
        1-Sales_programme_Brochures/Sales_programme/Sections_MB-ArcelorMittal_FR_EN_DE-V2017-3.pdf        
        page 193 ->
    """
    # Area of section [mm^2]
    A = 2.0 * tf * b + (h - 2.0 * tf) * tw + (4.0 - math.pi) * r ** 2.0

    # Perimeter of cross-section
    Au = 4 * (b - 2 * r) + 2 * (h - 2 * tw) + 2 * math.pi * r

    # Shear area [mm^2]    
    Ashear = A - 2.0 * b * tf + (tw + 2.0 * r) * tf

    # Second moment of area, I_y is stronger direction [mm^4]
    I = [0.0, 0.0]
    I[0] = 1.0 / 12.0 * (b * h ** 3.0 - (b - tw) * (h - 2.0 * tf) ** 3.0) + 0.03 * r ** 4 + 0.2146 * r ** 2.0 * (
                                                                                                                h - 2.0 * tf - 0.4468 * r) ** 2.0
    I[1] = 1.0 / 12.0 * (2.0 * tf * b ** 3.0 + (h - 2.0 * tf) * tw ** 3.0) + 0.03 * r ** 4.0 + 0.2146 * r ** 2.0 * (
                                                                                                                   tw + 0.4468 * r) ** 2.0

    # Torsional constant [mm^4]
    It = 2.0 / 3.0 * (b - 0.63 * tf) * tf ** 3.0 + 1.0 / 3.0 * (h - 2.0 * tf) * tw ** 3.0
    + 2.0 * (tw / tf) * (0.145 + 0.1 * r / tf) * (((r + tw / 2.0) ** 2.0 + (r + tf) ** 2.0 - r ** 2.0) / (
    2.0 * r + tf)) ** 4

    # Warping constant [mm^6]
    Iw = (tf * b ** 3.0) / 24.0 * (h - tf) ** 2.0

    # Plastic section modulus [mm^3]
    Wpl = [0.0, 0.0]
    Wpl[0] = (tw * h ** 2.0) / 4.0 + (b - tw) * (h - tf) * tf + (4.0 - math.pi) / 2.0 * r ** 2.0 * (h - 2.0 * tf) + (
                                                                                                                    3.0 * math.pi - 10.0) / 3.0 * r ** 3.0
    Wpl[1] = (b ** 2.0 * tf) / 2.0 + (h - 2.0 * tf) / 4.0 * tw ** 2.0 + r ** 3.0 * (10.0 / 3.0 - math.pi) + (
                                                                                                            2.0 - math.pi / 2.0) * tw * r ** 2.0

    # Elastic section modulus [mm^3]
    Wel = [0.0, 0.0]
    Wel[0] = I[0] / (0.5 * h)
    Wel[1] = I[1] / (0.5 * b)

    # Mass per unit length [kg/mm]
    # g = A*rho

    # Variable epsilon used in design
    # eps = math.sqrt(235.0/fy)

    return A, Au, Ashear, I, It, Iw, Wpl, Wel

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 23:12:31 2018

Hot-rolled I- and H-section

@author: kmela
"""
import math
try:
    from src.eurocodes.en1993 import en1993_1_1, constants
    from src.sections.steel.steel_section import SteelSection
    from src.sections.steel.catalogue import profile
except:
    from eurocodes.en1993 import en1993_1_1, constants
    from sections.steel.steel_section import SteelSection
    from sections.steel.catalogue import profile
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.lines as mlines
    import matplotlib.path as mpath
    import numpy as np

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
        
        #self.fy = fy

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
        
        
        if H/B <= 2.0:
            self.imp_factor_LT_gen = en1993_1_1.buckling_curve["a"]
            self.imp_factor_LT = en1993_1_1.buckling_curve["b"]
        else:
            self.imp_factor_LT_gen = en1993_1_1.buckling_curve["b"]
            self.imp_factor_LT = en1993_1_1.buckling_curve["c"]


    def __repr__(self):
        return f"{type(self).__name__} {self.h:.0f}"
    
    def info(self,latex=False):
        """ Prints a bunch of section properties """
        
        if latex:
            print("*** " + self.__repr__() + " ***")
            print(" $h = {0:.0f}{1:s}$".format(self.h,'\\,\\unit{mm}'))
            print(" $t_w = {0:.1f}{1:s}$".format(self.tw,'\\,\\unit{mm}'))
            print(" $b = {0:.0f}{1:s}$".format(self.b,'\\,\\unit{mm}'))
            print(" $t_f = {0:.1f}{1:s}$".format(self.tf,'\\,\\unit{mm}'))
            print(" $r = {0:.0f}{1:s}$".format(self.r,'\\,\\unit{mm}'))
            print(" $A = {0:.2f}{1:s}$".format(self.A*1e-2,'\\ee{2}\\,\\squnit{mm}'))
            print(" $I_y = {0:.2f}{1:s}$".format(self.I[0]*1e-4,'\\ee{4}\\,\\quunit{mm}'))
            print(" $I_z = {0:.2f}{1:s}$".format(self.I[1]*1e-4,'\\ee{4}\\,\\quunit{mm}'))
            print(" $W_el,y = {0:.2f}{1:s}$".format(self.Wel[0]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))
            print(" $W_el,z = {0:.2f}{1:s}$".format(self.Wel[1]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))
            print(" $W_pl,y = {0:.2f}{1:s}$".format(self.Wpl[0]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))
            print(" $W_pl,z = {0:.2f}{1:s}$".format(self.Wpl[1]*1e-3,'\\ee{3}\\,\\cuunit{mm}'))
            print(" $I_t = {0:.2f}{1:s}$".format(self.It*1e-4,'\\ee{4}\\,\\quunit{mm}'))
            print(" $I_w = {0:.2f}{1:s}$".format(self.Iw*1e-12,'\\ee{12}\\,\\text{mm}^{6}'))
        else:
            print("*** " + self.__repr__() + " ***")
            print(" h = {0:.0f} mm".format(self.h))
            print(" tw = {0:.1f} mm".format(self.tw))
            print(" b = {0:.0f} mm".format(self.b))
            print(" tf = {0:.1f} mm".format(self.tf))
            print(" r = {0:.0f} mm".format(self.r))
            print(" A = {0:.2f} (100 mm2)".format(self.A*1e-2))
            print(" Iy = {0:.2f} (10^4 mm4)".format(self.I[0]*1e-4))
            print(" Iz = {0:.2f} (10^4 mm4)".format(self.I[1]*1e-4))
            print(" Wel,y = {0:.2f} (10^3 mm4)".format(self.Wel[0]*1e-3))
            print(" Wel,z = {0:.2f} (10^3 mm4)".format(self.Wel[1]*1e-3))
            print(" Wpl,y = {0:.2f} (10^3 mm4)".format(self.Wpl[0]*1e-3))
            print(" Wpl,z = {0:.2f} (10^3 mm4)".format(self.Wpl[1]*1e-3))
            print(" It = {0:.2f} (10^4 mm4)".format(self.It*1e-4))
            print(" Iw = {0:.2f} (10^12 mm6)".format(self.Iw*1e-12))
        

    @property
    def hw(self):
        """ Straight part of the web """
        return self.h - 2 * self.tf - 2 * self.r

    @property
    def cf(self):
        """ Straight part of the flange """
        return 0.5 * (self.b - self.tw) - self.r
    
    @property
    def box_section_factor(self):
        """ Section factor for fire, considering the profile
            as a box
        """
        return 2*(self.b+self.h)/self.A*1e3
    
    
    def shadow_effect(self):
        """ Shadow effect factor ksh, from EN 1993-1-2 """
        
        return 0.9*self.box_section_factor/self.section_factor

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

    def draw(self,axes=None,origin=[0,0],theta=0):
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
        
        if theta != 0:
            phi = np.radians(theta)
            c, s = np.cos(phi), np.sin(phi)
            R = np.array(((c,-s), (s, c)))
        
        lw = 1.5
        
        h = self.h
        b = self.b
        tf = self.tf
        tw = self.tw
        r = self.r
        hw = self.hw
                
        xdata = [-0.5*tw-r,-0.5*b,-0.5*b,0.5*b,0.5*b,0.5*tw+r]
        ydata = [-0.5*h+tf, -0.5*h+tf, -0.5*h,-0.5*h,-0.5*h+tf,-0.5*h+tf]
        
        if origin[0] != 0:
            xdata = [x+origin[0] for x in xdata]
        
        if origin[1] != 0:
            ydata = [y+origin[1] for y in ydata]
        
        if theta != 0:
            xrot = [R.dot(np.array([x,y])) for x, y in zip(xdata,ydata)]
            xdata = [x[0] for x in xrot]
            ydata = [x[1] for x in xrot]
        
        bflange = mlines.Line2D(xdata,ydata,linewidth=lw,color='k')
                
        ax.add_line(bflange)
        
        # Bottom chamfers
        if theta != 0:
            x0 = tuple(R.dot(np.array([0.5*tw+r,-0.5*self.hw])))
        else:
            x0 = (0.5*tw+r+origin[0],-0.5*self.hw+origin[1])
                
        x1 = (-(0.5*tw+r)+origin[0],x0[1])        
        
        
        br_corner = patches.Arc(x0,2*r,2*r,0,180+theta,270+theta,linewidth=lw,zorder=4,color='k')
        bl_corner = patches.Arc(x1,2*r,2*r,0,270+theta,360+theta,linewidth=lw,zorder=4,color='k')
        
        ax.add_patch(br_corner)
        ax.add_patch(bl_corner)
        
    
        #ax.plot(0.5*tw+r,-0.5*self.hw,'bo')
        
        top_x_data = xdata
        top_y_data = [0.5*h-tf, 0.5*h-tf, 0.5*h,0.5*h,0.5*h-tf,0.5*h-tf]        
        
        if origin[1] != 0:
            top_y_data = [y+origin[1] for y in top_y_data]
        
        tflange =  mlines.Line2D(top_x_data,top_y_data,linewidth=lw,color='k')
        ax.add_line(tflange)
        
        
        
        # Top chamfers
        tr_corner = patches.Arc((0.5*tw+r+origin[0],0.5*hw+origin[1]),2*r,2*r,0,90,180,zorder=4,color='k',linewidth=lw)
        tl_corner = patches.Arc((-0.5*tw-r+origin[0],0.5*hw+origin[1]),2*r,2*r,0,0,90,zorder=4,color='k',linewidth=lw)
        
        ax.add_patch(tr_corner)
        ax.add_patch(tl_corner)
        
        # Draw web
        xweb_l = [-0.5*tw+origin[0],-0.5*tw+origin[0]]
        yweb = [-0.5*hw+origin[1],0.5*hw+origin[1]]
        
        if theta is not 0:
            web_rot = [R.dot(np.array([x,y])) for x, y in zip(xweb_l,yweb)]
            xweb_l = [x[0] for x in web_rot]
            yweb_l = [x[1] for x in web_rot]
        else:
            yweb_l = yweb
        
        lweb = mlines.Line2D(xweb_l,yweb_l,linewidth=lw,color='k')
        ax.add_line(lweb)
        
        xweb_r = [0.5*tw+origin[0],0.5*tw+origin[0]]
        
        rweb = mlines.Line2D(xweb_r,yweb,linewidth=lw,color='k')
        ax.add_line(rweb)
        
        """
        verts = [(0,0),
                 (1,0),
                 (0,1)]
        
        codes = [mpath.Path.]
        
        Iprof_path = mpath.Path(verts,codes)
        
        
        Iprof = patches.PathPatch(Iprof_path,facecolor="none",lw=1.5)
        ax.add_patch(Iprof)
        """
        
        #ax.add_patch(bot_flange)
        #ax.add_patch(web)
        #ax.add_patch(top_flange)
        #ax.set_xlim(-0.5*b, 0.5*b)
        #ax.set_ylim(-0.5*h, 0.5*h)
        
        if axes is  None:
            ax.set_xlim(-b, b)
            ax.set_ylim(-h, h)
        
            ax.set_aspect('equal')
            

class IPE(ISection):
    """ IPE sections
        Subclass of ISection
    """

    def __init__(self, height=100, fy=355, catalogue=True):
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

    def __init__(self, height=100, fy=355, catalogue=True):        
        name = 'HE ' + str(height) + ' ' + 'A'
        # HEIGHT =! TO PROFILE NUMBER

        self.name = name
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

    def __repr__(self):
        return self.name
        #return "HEA " + str(self.b-(self.b%10))

class HEAA(ISection):
    """ European wide flange sections
        Subclass of ISection        
    """

    def __init__(self, height=100, fy=355, catalogue=True):
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

    def __init__(self, height=100, fy=355, catalogue=True):
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

    def __init__(self, height=100, fy=355, catalogue=True):
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
    # Taken from Arcelor Mittal, with source SCI P385 Design of Members in Torsion
    #It = 2.0 / 3.0 * (b - 0.63 * tf) *tf ** 3.0 + 1.0 / 3.0 * (h - 2.0 * tf) * tw ** 3.0 + 2.0 * (tw / tf) * (0.145 + 0.1 * r / tf) * (((r + tw / 2.0) ** 2.0 + (r + tf) ** 2.0 - r ** 2.0) / (2.0 * r + tf)) ** 4
    a1 = -0.042+0.2204*tw/tf + 0.1355*r/tf - 0.0865*r*tw/tf**2 - 0.0725*tw**2/tf**2
    D1 = ((tf+r)**2+(r+0.25*tw)*tw)/(2*r+tf)
    It = 2.0/3.0*b*tf**3 + (h-2*tf)*tw**3/3.0 + 2*a1*D1**4 - 0.420*tf**4
    
    """
    It = 2.0/3.0*(b - 0.63*tf)*tf**3.0 + \
         + 1.0/3.0*(h-2.0*tf)*tw**3.0 + \
         + (0.042 + 0.2204*(tw / tf)+0.1355*(r/tf)-0.0865*(r*tw/tf**2) -0.0725*(tw/tf)**2) * (((r + tw/2.0)**2.0 + (r + tf)**2.0 - r**2.0)/(2.0*r + tf))**4.0
    """
    #print(tf,tw,b,h,It)

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


if __name__ == '__main__':
    
    p = IPE(400)
    p.info(latex=True)
    #p.draw(theta = 0)
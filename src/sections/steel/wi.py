# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:42 2019

Welded I-sections

@author: kmela
"""
import math
from eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
from sections.steel.steel_section import SteelSection

import matplotlib.pyplot as plt
import matplotlib.patches as patches

class WISection(SteelSection):
    """ Welded I-sections """

    def __init__(self, h, tw, b, tf, fy=355, weld_throat=5):
        """
            Mono-symmetric welded I-section
            
            h -- height of the profile
            b -- width of the flanges
                b[0] .. top flange
                b[1] .. bottom flange
            tw .. thickness of web
            tf .. tf thickness of flanges
                tf[0] .. top flange
                tf[1] .. bottom flange
            fy -- yield strength
            weld_throat .. throat thickness of the weld between
                           flanges and web (fillet weld assumed)
        """
        self.h = h
        self.b = b
        self.tf = tf
        self.tw = tw
        self.weld_throat = weld_throat
        
        self.fy = fy
        
        A = self.area()
        Au = self.paint_area()
        I = self.second_moment()
        Ashear = self.shear_area()
        Wel = self.section_modulus()
        Wpl = self.plastic_section_modulus()        
        
        super().__init__(fy, A, I, Au, Wpl, Wel, Ashear)

        """ Determine buckling curve: EN 1993-1-1, Table 6.2 """        
        self.imp_factor = [en1993_1_1.buckling_curve["b"],
                           en1993_1_1.buckling_curve["b"]]
            
    def __getattribute__(self,name):
        """ override the attribute access for those attributes
            that depend on the section dimensions and that are
            constant for rolled sections
        """
        if name=="A":            
            return self.area()
        elif name=="I":
            return self.second_moment()
        elif name=="Wel":
            return self.section_modulus()
        elif name=="Wpl":
            return self.plastic_section_modulus()
        elif name=="Ashear":
            return self.shear_area()
        elif name=="Au":
            return self.paint_area()
        else:
            return super().__getattribute__(name)
    
    
    @property
    def hw(self):
        """ Web height """
        return self.h - self.tf[0] - self.tf[1]
    
    @property
    def Aw(self):
        """ Area of web plate """
        return self.hw*self.tw
    
    @property
    def zw(self):
        """ Distance of centroid of web from bottom of section """
        return self.tb + 0.5*self.hw

    @property
    def tb(self):
        """ Thickness of bottom flange """
        return self.tf[1]
    
    @property
    def bb(self):
        """ Width of bottom flange """
        return self.b[1]
    
    @property
    def Ab(self):
        """ Area of bottom flange """
        return self.bb*self.tb
    
    @property
    def zb(self):
        """ Distance of centroid of bottom flange from bottom of section """
        return 0.5*self.tb
    
    @property
    def tt(self):
        """ Thickness of top flange """
        return self.tf[0]
    
    @property
    def bt(self):
        """ Width of top flange """
        return self.b[0]
    
    @property
    def At(self):
        """ Area of top flange """
        return self.bt*self.tt
    
    @property
    def zt(self):
        """ Distance of centroid of top flange from bottom of section """
        return self.h -0.5*self.tt
    
    @property
    def cf_top(self):
        """ Straight part of the top flange """
        return 0.5*(self.b[0]-self.tw) - math.sqrt(2)*self.weld_throat

    @property
    def cf_bottom(self):
        """ Straight part of the bottom flange """
        return 0.5*(self.b[1]-self.tw) - math.sqrt(2)*self.weld_throat

    @property
    def cw(self):
        """ Straight part of the web """
        return self.hw - 2*math.sqrt(2)*self.weld_throat
   
    def flange_class(self,verb=False):
        """ Determine class of compressed flange """
        # cf = 0.5*(self.b - self.tw) - self.r
        rf = self.cf_top / self.tt
        cFlange = en1993_1_1.outstand_part_in_compression(rf, self.eps)

        if verb:
            print("Flange classification (outstand element):")
            print("cf = {0:4.2f}, tf = {1:4.2f}".format(self.cf_top,self.tt))
            print("cf/tf = {0:4.2f}".format(rf))
            print("Flange class = {0}".format(cFlange))

        return cFlange

    def web_class_comp(self,verb=False):
        """ Determine class of compressed web """
        rw = self.cw / self.tw
        cWeb = en1993_1_1.internal_part_in_compression(rw, self.eps)

        if verb:
            print("Web classification (internal part in compression):")
            print("hw = {0:4.2f}, tw = {1:4.2f}".format(self.cw,self.tw))
            print("hw/tw = {0:4.2f}".format(rw))
            print("Web class = {0}".format(cWeb))

        return cWeb

    def web_class_bend(self,verb=False):
        """ Determine class of web in bending """
        zel = self.elastic_neutral_axis()
        dpl = self.plastic_neutral_axis()
        rw = self.cw / self.tw
        psi = zel/(self.h-zel)
        alpha = dpl/self.hw
        cWeb = en1993_1_1.internal_part_comp_bend(rw, self.eps,alpha,psi)

        if verb:
            print("Web classification (internal part in compression and bending):")
            print("hw = {0:4.2f}, tw = {1:4.2f}".format(self.hw,self.tw))
            print("hw/tw = {0:4.2f}".format(rw))
            print("alpha = {0:4.2f}".format(alpha))
            print("psi = {0:4.2f}".format(psi))
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
            TODO
        """
        aw = min((self.A - 2 * self.b * self.tf) / self.A, 0.5)

        if UN > 1.0:
            MNRd = 0.0
        else:
            MNRd = MRd * min((1 - UN) / (1 - 0.5 * aw), 1)

        return MNRd

    def area(self):
        """ Cross-sectional area """
        return self.Aw + self.At + self.Ab
    
    def shear_area(self):
        """ Shear area """
        return en1993_1_5.shear_eta(self.fy)*self.hw*self.tw
    
    def paint_area(self):
        """ Area to be painted (circumference of the section) """
        return 2*self.tt + self.bt + 2*self.hw + 2*self.tb + 2*self.bb - 2*self.tw
        
    
    def elastic_neutral_axis(self):
        """ Elastic neutral axis 
            measured from the bottom of the section
        """
        zel = (self.zb*self.Ab + self.zw*self.Aw + self.zt*self.At)/self.area()
        return zel
    
    def plastic_neutral_axis(self):
        """ Plastic neutral axis
            measured downwards from the bottom of the top flange
        """
        # with respect to y-axis
        dpl = 0.5*(self.Ab + self.Aw - self.At)/self.tw
        return dpl
    
    def second_moment(self):
        """ Second moment of area with respect to both principal axes"""
        zel = self.elastic_neutral_axis()
                
        # with respect to z axis
        Iy = 1/12*self.bb*self.tb**3 + (zel-self.zb)**2*self.Ab + \
                1/12*self.tw*self.h**3 + (zel-self.zw)**2*self.Aw + \
                1/12*self.bt*self.tt**3 + (zel-self.zt)**2*self.At

        # with respect to z axis        
        Iz = 1/12*self.tb*self.bb**3 + 1/12*self.tt*self.bt**3 + \
                1/12*self.hw*self.tw**3

        return [Iy,Iz]
    
    def section_modulus(self):
        """ Elastic section modulus for bending """
        
        I = self.second_moment()
        zel = self.elastic_neutral_axis()
        
        # distance of neutral axis from the top of the section
        ztop = self.h - zel
        
        Wel =[0,0]
        Wel[0] = I[0]/max(zel,ztop)
        
        # with respect to z-axis
        Wel[1] = I[1]/(0.5*max(self.b))
        
        return Wel
    
    def plastic_section_modulus(self):
        """ Plastic section modulus for bending """
        
        """ distance of plastic neutral axis from the bottom of the top flange """
        Wpl = [0,0]
    
        # with respect to y-axis
        dpl = self.plastic_neutral_axis()

        
        if dpl < 0:
            # neutral axis is in the top flange
            # TODO
            Wpl[0] = 0
        elif dpl > self.hw:
            # neutral axis is in the bottom flange
            # TODO
            Wpl[0] = 0        
        else:
            # neutral axis is in the web
            rt = dpl + 0.5*self.tt
            rw1 = 0.5*(dpl+self.tt)
            rw2 = 0.5*(self.hw-dpl)
            rb = self.hw-dpl+0.5*self.tb
            
            Wpl[0] = self.At*rt + 2*(dpl+self.tt)*self.tw*rw1 + \
                    2*(self.hw-dpl)*self.tw*rw2 + self.Ab*rb
            
        # with respect to z-axis
        Wpl[1] = 1/4*self.Ab*self.bb + 1/4*self.At*self.bt + \
                    1/4*self.Aw*(self.bt+self.tw)
        return Wpl
    
    def effective_flange(self,flange="Top"):
        """ Effective flange, assumed to be in pure compression """
        
        ksigma = en1993_1_5.buckling_factor_outstand(psi=1.0)
                
        if flange == "Top":
            c = self.cf_top
            t = self.tt
        else:
            c = self.cf_bottom
            t = self.tb

        lp = en1993_1_5.lambda_p(c,t,self.eps,ksigma)
        rho = en1993_1_5.reduction_factor_outstand(lp)

        return rho                
    
    def draw(self):
        """ Draw the profile """
        fig, ax = plt.subplots(1)
        
        # create section parts
        bot_flange = patches.Rectangle((-0.5*self.bb,0),width=self.bb,height=self.tb,\
                                       fill=False, hatch='\\')
        web = patches.Rectangle((-0.5*self.tw,self.tb),width=self.tw,height=self.hw,\
                                 fill=False, hatch='\\')
    
        top_flange = patches.Rectangle((-0.5*self.bt,self.h-self.tt),\
                                       width=self.bt,height=self.tt,\
                                       fill=False, hatch='\\')
        ax.add_patch(bot_flange)
        ax.add_patch(web)
        ax.add_patch(top_flange)
        ax.set_xlim(-0.5*max(self.b),0.5*max(self.b))
        ax.set_ylim(0,self.h)
        
        zel = self.elastic_neutral_axis()
        dpl = self.plastic_neutral_axis()        
        
        print(zel)
        print(dpl)
        
        plt.plot(0.0,zel,'or')
        plt.plot(0.0,self.tb+self.hw-dpl,'db')
        
        ax.set_aspect('equal')
        
        plt.show()
        
def test_sym():
    """ symmetric I-beam """
    p = WISection(500,8,[250,250],[12,12])
    p.draw()

def test_non_sym():
    """ symmetric I-beam """
    print("Create section")
    p = WISection(500,8,[300,250],[10,15])
    print("Done.")
    print(p.I)
    # p.draw()


if __name__ == '__main__':
    
    test_non_sym()
    
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 22:17 2019

WQ sections (welded box beam)

@author: kmela
"""

from eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
from sections.steel.steel_section import SteelSection

import matplotlib.pyplot as plt
import matplotlib.patches as patches

class WQSection(SteelSection):
    """ WQ-sections """

    def __init__(self, h, tw, b, tf, fy=355):
        """
            h -- height from top of profile to top of bottom flange
            b -- b[0] = width of top flange (between webs)
                 b[1] = width of bottom flange
            tw .. thickness of web
            tf .. tf[0] thickness of top flange
                  tf[1] thickness of bottom flange
            fy -- yield strength
        """
        self.h = h
        self.b = b
        self.tf = tf
        self.tw = tw
        
        #self.fy = fy
        
        A = self.area()
        Au = self.paint_area()
        I = self.second_moment()
        #Ashear = self.shear_area()
        Wel = self.section_modulus()
        Wpl = self.plastic_section_modulus()        
        
        
        super().__init__(fy, A, I, Au, Wpl, Wel, Ashear=1.0)

        self.Ashear = self.shear_area()

        """ Determine buckling curve: EN 1993-1-1, Table 6.2 """        
        self.imp_factor = [en1993_1_1.buckling_curve["b"],
                           en1993_1_1.buckling_curve["b"]]
        
    @property
    def height(self):
        """ Total height of the profile """
        return self.h+self.tb
    
    @property
    def hw(self):
        """ Part of the web between flanges"""
        return self.h - self.tf[0]
    @property
    def Aw(self):
        """ Area of one web plate """
        return self.h*self.tw
    
    @property
    def zw(self):
        """ Distance of centroid of web from bottom of section """
        return self.tb + 0.5*self.h

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
        return self.tb + self.h -0.5*self.tt
    
    @property
    def cf(self):
        """ Straight part of the flange """
        return self.b[0]

    @property
    def lip(self):
        """ length of outstand part of bottom flange (lip) """
        return 0.5*(self.bb-self.bt) - self.tw

    def flange_class(self,verb=False):
        """ Determine class of compressed flange """
        # cf = 0.5*(self.b - self.tw) - self.r
        rf = self.bt / self.tt
        cFlange = en1993_1_1.internal_part_in_compression(rf, self.eps)

        if verb:
            print("Flange classification (internal element):")
            print("cf = {0:4.2f}, tf = {1:4.2f}".format(self.bt,self.tt))
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
        zel = self.elastic_neutral_axis()
        dpl = self.plastic_neutral_axis()
        rw = self.hw / self.tw
        psi = zel/(self.height-zel)
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
        return 2*self.Aw + self.At + self.Ab
    
    def shear_area(self):
        """ Shear area """
        return 2*en1993_1_5.shear_eta(self.fy)*self.hw*self.tw
    
    def paint_area(self):
        """ Area to be painted (circumference of the section) """
        return self.bb+2*self.tb + 2*self.lip + 2*(self.h+self.tw) + self.bt
        
    
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
        dpl = 0.5*(self.hw-self.tt)-0.25*(self.At-self.Ab)/self.tw
        return dpl
    
    def second_moment(self):
        """ Second moment of area with respect to both principal axes"""
        zel = self.elastic_neutral_axis()
                
        # with respect to z axis
        Iy = 1/12*self.bb*self.tb**3 + (zel-self.zb)**2*self.Ab + \
                1/6*self.tw*self.h**3 + (zel-self.zw)**2*2*self.Aw + \
                1/12*self.bt*self.tt**3 + (zel-self.zt)**2*self.At

        # with respect to z axis
        # yw .. distance of web centroid from z axis
        yw = 0.5*(self.bt+self.tw)
        Iz = 1/12*self.tb*self.bb**3 + 1/12*self.tt*self.bt**3 + \
                1/6*self.h*self.tw**3 + yw**2*2*self.Aw

        return [Iy,Iz]
    
    def section_modulus(self):
        """ Elastic section modulus for bending """
        
        I = self.second_moment()
        zel = self.elastic_neutral_axis()
        
        # distance of neutral axis from the top of the section
        ztop = self.h+self.bb - zel
        
        Wel =[0,0]
        Wel[0] = I[0]/max(zel,ztop)
        Wel[1] = I[1]/(0.5*self.bb)
        
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
                    self.Aw*(self.bt+self.tw)
        return Wpl
    
    def draw(self,axes=None,origin=[0,0],plot_centroid=False,
             plot_plastic_centroid=False):
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
        
        # create section parts
        bot_flange = patches.Rectangle((0,0),width=self.bb,height=self.tb,\
                                       fill=False, hatch='\\')
        web1 = patches.Rectangle((self.lip,self.tb),width=self.tw,height=self.h,\
                                 fill=False, hatch='\\')
        web2 = patches.Rectangle((self.bb-self.lip-self.tw,self.tb),\
                                 width=self.tw,height=self.h,\
                                 fill=False, hatch='\\')
        top_flange = patches.Rectangle((self.lip+self.tw,self.hw+self.tb),\
                                       width=self.bt,height=self.tt,\
                                       fill=False, hatch='\\')
        ax.add_patch(bot_flange)
        ax.add_patch(web1)
        ax.add_patch(web2)
        ax.add_patch(top_flange)
        ax.set_xlim((0,self.bb))
        ax.set_ylim((0,self.tb+self.h))
        
        if plot_centroid:
            zel = self.elastic_neutral_axis()
            ax.plot(0.5*self.bb,zel,'or')
        
        if plot_plastic_centroid:
            dpl = self.plastic_neutral_axis()        
        
            ax.plot(0.5*self.bb,self.tb+self.hw-dpl,'db')
        
        plt.show()
        # bottom flange coordinates
        #bf_xy = [[0,0],[self.bb,0.0],[self.bb,self.tb],[0,self.tb]]
        
        return ax
        

if __name__ == '__main__':
    p = WQSection(390, 10, [200,400], [20,25])
    #print(p.elastic_neutral_axis())
    p.web_class_bend(verb=True)
    p.flange_class(verb=True)
    p.draw()
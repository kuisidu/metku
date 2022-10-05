# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Fri May 17 10:42 2019

Welded box sections

@author: kmela
"""
import math
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches


from metku.eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
from metku.sections.steel.steel_section import SteelSection



class WBSection(SteelSection):
    """ Welded box sections """

    def __init__(self, h, tw, b, tf, c, fy=355, weld_throat=5):
        """
            Mono-symmetric welded I-section
            
            h -- height of the profile
            b -- width of the flanges
            tw .. thickness of web
            tf .. tf thickness of flanges                
            c .. length of the lip
            fy -- yield strength
            weld_throat .. throat thickness of the weld between
                           flanges and web (fillet weld assumed)
        """
        self.h = h
        self.b = b
        self.tf = tf
        self.tw = tw
        self.weld_throat = weld_throat
        self.lip = c        
        #self.fy = fy
        
        # Reduction factors for effective width of parts
        self.rho_web = 1.0
        self.rho_top = 1.0
        self.rho_bottom = 1.0
        
        A = self.area()
        Au = self.paint_area()
        I = self.second_moment()
        Ashear = self.shear_area(fy)
        Wel = self.section_modulus()
        Wpl = self.plastic_section_modulus()

        self.It = self.torsional_constant()
        self.Iw = self.warping_constant()

        super().__init__(fy, A, I, Au, Wpl, Wel, Ashear)

        """ Determine buckling curve: EN 1993-1-1, Table 6.2 """ 
        if max(tf) <= 40:
            self.imp_factor = [en1993_1_1.buckling_curve["b"],
                               en1993_1_1.buckling_curve["c"]]
        else:
            self.imp_factor = [en1993_1_1.buckling_curve["c"],
                               en1993_1_1.buckling_curve["d"]]

        if h/min(b) <= 2.0:
            self.imp_factor_LT_gen = en1993_1_1.buckling_curve["c"]
            self.imp_factor_LT = en1993_1_1.buckling_curve["c"]
        else:
            self.imp_factor_LT_gen = en1993_1_1.buckling_curve["d"]
            self.imp_factor_LT = en1993_1_1.buckling_curve["d"]
            
    def __getattribute__(self, name):
        """ override the attribute access for those attributes
            that depend on the section dimensions and that are
            constant for rolled sections
        """
        if name == "A":
            if self.area() < 0:
                print(self.h, self.b, self.tt, self.tb, self.tw)
                
            return self.area()
        elif name == "I":
            I = self.second_moment()
            if I[1] < 0:
                print(self.h, self.b, self.tt, self.tb, self.tw)
            return self.second_moment()
        elif name == "Wel":
            return self.section_modulus()
        elif name == "Wpl":
            return self.plastic_section_modulus()
        elif name == "Ashear":
            return self.shear_area()
        elif name == "Au":
            return self.paint_area()
        elif name == "It":
            return self.torsional_constant()
        elif name == "Iw":
            return self.warping_constant()
        elif name == "zs":
            return self.shear_centre()
        elif name == "zj":
            return self.wagners_factor()
        elif name == "b":
            return np.array([self.bt, self.bb])
        elif name == "tf":
            return np.array([self.tt, self.tb])
        else:
            return super().__getattribute__(name)

    def __setattr__(self, key, value):
        if key == "b":        
            self.bt = value
            self.bb = value
        elif key == "tf":
            self.tt = value
            self.tb = value
        #elif key == "bt":
        #    self.b[0] = value
        #elif key == "bb":
        #    self.b[1] = value
        else:
            super().__setattr__(key, value)

    @property
    def hw(self):
        """ Web height """
        return self.h - 2*self.tf
    
    @property
    def Aw(self):
        """ Area of web plate """
        return self.hw * self.tw
    
    @property
    def Izw(self):
        """ Second moment of area of web with respect to its own major axis """
        return self.hw * self.tw**3/12

    @property
    def Iyw(self):
        """ Second moment of area of web with respect to its own major axis """
        return self.tw * self.hw ** 3 / 12
    
    @property
    def zw(self):
        """ Distance of centroid of web from bottom of section """
        return self.tb + 0.5 * self.hw

    # @property
    # def tb(self):
    #     """ Thickness of bottom flange """
    #     return self.tf[1]
    #
    # @property
    # def bb(self):
    #     """ Width of bottom flange """
    #     return self.b[1]
    
    @property
    def Af(self):
        """ Area of one flange """
        return self.b * self.tf
    
    @property
    def Izb(self):
        """ Second moment of area of bottom flange with respect to its own
            major axis
        """
        return self.tb * self.bb**3/12

    @property
    def Iyb(self):
        """ Second moment of area of bottom flange with respect to its own
            major axis
        """
        return self.bb * self.tb ** 3 / 12
    
    @property
    def zb(self):
        """ Distance of centroid of bottom flange from bottom of section """
        return 0.5 * self.tb
        
    @property
    def Izt(self):
        """ Second moment of area of top flange with respect to its own major
            axis
        """
        return self.tt * self.bt**3/12

    @property
    def Iyt(self):
        """ Second moment of area of top flange with respect to its own major
            axis
        """
        return self.bt * self.tt ** 3 / 12
    
    @property
    def zt(self):
        """ Distance of centroid of top flange from bottom of section """
        return self.h - 0.5 * self.tt
    
    @property
    def cf_top(self):
        """ Straight part of the top flange """
        return 0.5*(self.bt-self.tw) - math.sqrt(2)*self.weld_throat

    @property
    def cf_bottom(self):
        """ Straight part of the bottom flange """
        return 0.5*(self.bb-self.tw) - math.sqrt(2)*self.weld_throat

    @property
    def cw(self):
        """ Straight part of the web """
        return self.hw - 2*math.sqrt(2)*self.weld_throat
   
    @property
    def MfRd(self):
        """ Plastic bending resistanc of flanges only
            Commentary and Worked Examples to EN 1993-1-5, pp.70, Fig. 5.11
        """
        
        # Distance between center lines of flanges
        hf = self.h-0.5*(self.tb+self.tt)
        Af1 = self.At
        Af2 = self.Ab
        return min(hf*Af1*self.fy,hf*Af2*self.fy)
    
    @property
    def NtopRd(self):
        """ Axial resistance of top flange
        """
        return self.At*self.fy
    
    @property
    def NbottomRd(self):
        """ Axial resistance of bottom flange
        """
        return self.Ab*self.fy
    
    
    def Aeff(self):
        """ Effective area of the section """
        
    
    def shear_centre(self):
        """ Position of the shear centre S from centroid G
            Source: Maquoi et al. (2003)
            Lateral torsional buckling in steel and composite beams,
            Book 2, 7.1.2 Mono-symmetrical I cross-sections
        """
        
        """ Position of the centroids of the different plates
            with respect to centroid of the cross-section
        """
        zG = self.elastic_neutral_axis()
        zft = self.zt-zG
        zfb = self.zb-zG
        zwS = self.zw-zG
        
        return (self.Izt*zft + self.Izb*zfb + self.Izw*zwS)/self.I[1]

    def wagners_factor(self):
        """ Source: Maquoi et al. (2003)
            Lateral torsional buckling in steel and composite beams,
            Book 2, 7.1.2 Mono-symmetrical I cross-sections
        """

        zf1 = self.zt
        zf2 = self.zb
        zw = self.zw
        Izf1 = self.Izt
        Izf2 = self.Izb
        Izw = self.Izw
        Iyf1 = self.Iyt
        Iyf2 = self.Iyb
        Iyw = self.Iyw
        Af1 = self.At
        Af2 = self.Ab
        Aw = self.Aw

        zs = self.zs

        zj = zs - (1 / (2 * self.I[0])) * (
                zf1 * (Izf1 + Af1 * zf1 ** 2 + 3 * Iyf1)
                + (zf2 * (Izf2 + Af2 * zf2 ** 2 + 3 * Iyf2))
                + (zw * (Izw + Aw * zw ** 2 + 3 * Iyw)))

        return zj
    
    def flange_class(self, verb=False):
        """ Determine class of compressed flange """
        # cf = 0.5*(self.b - self.tw) - self.r
        rf = self.cf_top / self.tt
        cFlange = en1993_1_1.outstand_part_in_compression(rf, self.eps)

        if verb:
            print("Flange classification (outstand element):")
            print("cf = {0:4.2f}, tf = {1:4.2f}".format(self.cf_top, self.tt))
            print("cf/tf = {0:4.2f}".format(rf))
            print("Flange class = {0}".format(cFlange))

        return cFlange

    def web_class_comp(self, verb=False):
        """ Determine class of compressed web """
        rw = self.cw / self.tw
        cWeb = en1993_1_1.internal_part_in_compression(rw, self.eps)

        if verb:
            print("Web classification (internal part in compression):")
            print("cw = {0:4.2f}, tw = {1:4.2f}".format(self.cw, self.tw))
            print("cw/tw = {0:4.2f}".format(rw))
            print("Web class = {0}".format(cWeb))

        return cWeb

    def web_class_bend(self, verb=False):
        """ Determine class of web in bending """
        zel = self.elastic_neutral_axis()
        dpl = self.plastic_neutral_axis()
        rw = self.cw / self.tw
        psi = (-zel+self.tb)/(self.h-zel-self.tt)

        if dpl <= self.tb + math.sqrt(2) * self.weld_throat:
            alpha = 1
        elif dpl >= self.h - self.tt - math.sqrt(2) * self.weld_throat:
            alpha = 0  # TODO tarkista tämä
        else:
            alpha = (self.h - self.tt - math.sqrt(2) * self.weld_throat - dpl)\
                    / (self.h - self.tt - self.tb - 2 * math.sqrt(2))

        cWeb = en1993_1_1.internal_part_comp_bend(rw, self.eps, alpha, psi)

        if verb:
            print("Web classification (internal part in compression and "
                  "bending):")
            print("cw = {0:4.2f}, tw = {1:4.2f}".format(self.cw, self.tw))
            print("cw/tw = {0:4.2f}".format(rw))
            print("alpha = {0:4.2f}".format(alpha))
            print("psi = {0:4.2f}".format(psi))
            print("Web class = {0}".format(cWeb))

        return cWeb

    def web_class_comp_bend(self, Ned, verb=False):
        """ Determine class of web in combined bending and compression 
            TODO tarkista tämä
        """
        # cw = self.h-2*self.tf-2*self.r
        rw = self.hw / self.tw
        Ac = 0.5 * (self.A + Ned / (self.fy / constants.gammaM0))
        a = Ac / self.A
        p = -1

        cWeb = en1993_1_1.internal_part_comp_bend(rw, self.eps, a, p)

        if verb:
            print("Web classification (internal part in compression and "
                  "bending):")
            print("hw = {0:4.2f}, tw = {1:4.2f}".format(self.hw, self.tw))
            print("hw/tw = {0:4.2f}".format(rw))
            print("Web class = {0}".format(cWeb))

        return cWeb

    def moment_axial_force_interact(self, UN, MRd):
        """ Interaction rule for section resistance for combined
            axial force and bending
            
            input: UN .. NEd/NRd
                  MRd .. moment resistance           
            TODO tarkista tämä
        """
        aw = min((self.A - (self.bt * self.tt - self.bb * self.tb))
                 / self.A, 0.5)

        if UN > 1.0:
            MNRd = 0.0
        else:
            MNRd = MRd * min((1 - UN) / (1 - 0.5 * aw), 1)

        return MNRd

    def area(self):
        """ Cross-sectional area """
        return self.Aw + 2*self.Af
    
    def shear_area(self,fy=None):
        """ Shear area """
        # if self.hw < 1e-6 or self.tw < 1e-3:
        #     print(self.hw, self.tw)
        if fy is None:
            fy = self.fy
        
        return en1993_1_5.shear_eta(fy)*self.hw*self.tw
    
    def paint_area(self):
        """ Area to be painted (circumference of the section) """
        return 2*self.tt + self.bt + 2*self.hw + 2*self.tb + 2*self.bb \
               - 2*self.tw

    def elastic_neutral_axis(self):
        """ Elastic neutral axis 
            measured from the bottom of the section
        """
        zel = (self.zb*self.Ab + self.zw*self.Aw + self.zt*self.At)/self.area()
        return zel
    
    def plastic_neutral_axis(self):
        """ Plastic neutral axis
            measured from the bottom of the section
        """
        # measured downwards from the bottom of the top flange
        # dpl = 0.5*(self.Ab + self.Aw - self.At)/self.tw

        # with respect to y-axis
        if self.At >= self.A / 2:
            dpl = self.h - self.A / (2 * self.bt)
        elif self.Ab >= self.A / 2:
            dpl = self.A / (2 * self.bb)
        else:
            dpl = (self.A / 2 - self.Ab) / self.tw + self.tb

        #  print(dpl)
        return dpl
    
    def second_moment(self):
        """ Second moment of area with respect to both principal axes"""
        zel = self.elastic_neutral_axis()
                
        # with respect to z axis
        Iy = 1/12*self.bb*self.tb**3 + (zel-self.zb)**2*self.Ab + \
                1/12*self.tw*self.hw**3 + (zel-self.zw)**2*self.Aw + \
                1/12*self.bt*self.tt**3 + (zel-self.zt)**2*self.At

        # with respect to z axis        
        Iz = 1/12*self.tb*self.bb**3 + 1/12*self.tt*self.bt**3 + \
                1/12*self.hw*self.tw**3

        return [Iy, Iz]
    
    def section_modulus(self):
        """ Elastic section modulus for bending """
        
        I = self.second_moment()
        zel = self.elastic_neutral_axis()
        
        # distance of neutral axis from the top of the section
        ztop = self.h - zel
        
        Wel = [0, 0]
        Wel[0] = I[0]/max(zel, ztop)
        
        # with respect to z-axis
        Wel[1] = I[1]/(0.5*max(self.b))
        
        return np.asarray(Wel)
    
    def plastic_section_modulus(self):
        """ Plastic section modulus for bending """
        
        """ distance of plastic neutral axis from the bottom of the top 
        flange """

        Wpl = [0, 0]
    
        # with respect to y-axis
        dpl = self.plastic_neutral_axis()

        if self.At >= self.A / 2:
            # neutral axis is in the top flange
            Wpl[0] = (self.bt * pow(self.h - dpl, 2) / 2) \
                     + (self.bt * pow(dpl - self.hw - self.tb, 2) / 2) \
                     + (self.Aw * (dpl - self.hw/2-self.tb)) \
                     + (self.Ab * (dpl - self.tb / 2))

        elif self.Ab >= self.A / 2:
            # neutral axis is in the bottom flange
            Wpl[0] = (self.At * (self.h - dpl - self.tt / 2)) \
                     + (self.Aw * (self.h - dpl - self.tt - self.hw / 2)) \
                     + (self.bb * pow(self.h - dpl -
                                      self.tt - self.hw, 2) / 2) \
                     + (self.bb * pow(dpl, 2) / 2)

        else:
            # neutral axis is in the web
            Wpl[0] = (self.At * (self.h - dpl - self.tt / 2)) \
                     + (self.tw * pow(self.h - dpl - self.tt, 2) / 2) \
                     + (self.tw * pow(dpl - self.tb, 2) / 2) \
                     + (self.Ab * (dpl - self.tb / 2))
            
        # with respect to z-axis
        Wpl[1] = (1/4 * self.Ab * self.bb) + (1/4 * self.At * self.bt) \
                 + (1/4 * self.Aw * self.tw)

        return np.asarray(Wpl)

    def torsional_constant(self):
        """ Torsional constant """

        It = 1/3 * (self.bt * pow(self.tt, 3) + self.bb *
                    pow(self.tb, 3) + self.hw * pow(self.tw, 3))

        # Tarkempi kaava (hidastaa laskentaa) ja lisätermit lähellä arvoa 1
        # It = 1/3 * (self.bt * pow(self.tt, 3) * (1 - 0.63 * (
        #         self.tt / self.bt) * (1 - pow(self.tt, 4) / (12 * pow(
        #             self.bt, 4)))) + self.bb * pow(self.tb, 3) * (1 - 0.63 * (
        #                 self.tb / self.bb) * (1 - pow(self.tb, 4) / (12 * pow(
        #                     self.bb, 4)))) + self.hw * pow(self.tw, 3))

        return It

    def warping_constant(self):
        """ Warping constant on the weaker axis """

        Izf1 = (self.tt * pow(self.bt, 3)) / 12
        Iw = Izf1 * (1 - Izf1 / self.I[1]) * \
             pow(self.h - self.tt / 2 - self.tb / 2, 2)

        return Iw

    def effective_flange(self, flange="Top"):
        """ Effective flange, assumed to be in pure compression """
        
        ksigma = en1993_1_5.buckling_factor_outstand(psi=1.0)
                
        if flange == "Top":
            c = self.cf_top
            t = self.tt
        else:
            c = self.cf_bottom
            t = self.tb

        lp = en1993_1_5.lambda_p(c, t, self.eps, ksigma)
        rho = en1993_1_5.reduction_factor_outstand(lp)

        return rho
    
    def effective_web(self,psi=None,verb=False):
        """ Effective web of the cross-section """
        
        cw = self.cw
        
        """ If stress ratio is not provided, calculate it """
        if psi == None:
            zel = self.elastic_neutral_axis()
            psi = (-zel+self.tb)/(self.h-zel-self.tt)
        
        ksigma = en1993_1_5.buckling_factor_internal(psi)
        lp = en1993_1_5.lambda_p(cw, self.tw, self.eps, ksigma)
        rho = en1993_1_5.reduction_factor_internal(lp,psi)
        
        beff, be = en1993_1_5.effective_width_internal(cw,psi,rho)
        
        if verb:
            print("Effective Web:")
            print("psi = {0:4.2f}".format(psi))
            print("k_sigma = {0:4.3f}".format(ksigma))
            print("lambda_p = {0:4.3f}".format(lp))
            print("rho = {0:4.3f}".format(rho))
            print("beff = {0:4.3f}".format(beff))
            print("be1 = {0:4.3f}, be2 = {1:4.3f} ".format(be[0],be[1]))
        
        

    
    def shear_buckling_resistance(self,stiffener_spacing=10e5,end_post="non-rigid"):
        """ Shear buckling resistance according to EN 1993-1-5, Clause 5 
            input:
                stiffener_spacing .. spacing between adjacent stiffeners
                end_post .. "non-rigid" or "rigid"
            
        """
    
        hw = self.hw
        tw = self.tw
        b = self.bt
        a = stiffener_spacing
        rw = hw/tw
        fyw = self.fy
        eta = en1993_1_5.shear_eta(fyw)
    
        if rw < 72*self.eps/eta:
            """ No need for shear buckling:
                return shear resistance of section
            """
            VbRd = self.VRd
        else:           
            """ Contribution from the web """
            tau_cr = en1993_1_5.tau_crit(hw,a,tw,b)            
            slend_w = en1993_1_5.shear_buckling_slenderness(fyw,tau_cr)
            chi_w = en1993_1_5.shear_buckling_reduction_factor(slend_w,eta,end_post)
            Vbw_Rd = en1993_1_5.shear_buckling_web(chi_w,fyw,hw,t)
            
            """ Contribution from the flange """
            NEd = abs(self.Ned)
            
            """ If axial force is present, the resistance MfRd must
                be reduced by the factor (1-rN). (EN 1993-1-5, 5.4(2))
            """
            if NEd > 1e-6:
                rN = NEd/((self.At+self.Ab)*fy)
            else:
                rN = 0
            
            rM = abs(self.Med)/((1-rN)*self.MfRd)
            
            if rM < 1.0:
                if self.NtopRd <= self.NbottomRd:
                    bf = self.bt
                    tf = self.tt
                else:
                    bf = self.bb
                    tf = self.tb
                
                bf = min(bf,30*e*tf+tw)
                
                Vbf_Rd = en1993_1_5.shear_buckling_flanges(bf,tf,fy,a,hw,t,fy,rM)
            
    
    def draw(self, name=""):
        """ Draw the profile """
        fig, ax = plt.subplots(1)
        
        # create section parts
        bot_flange = patches.Rectangle((-0.5*self.bb, 0), width=self.bb,
                                       height=self.tb, fill=False, hatch='\\')
        web = patches.Rectangle((-0.5*self.tw, self.tb), width=self.tw,
                                height=self.hw, fill=False, hatch='\\')
    
        top_flange = patches.Rectangle((-0.5*self.bt, self.h-self.tt),
                                       width=self.bt, height=self.tt,
                                       fill=False, hatch='\\')
        ax.add_patch(bot_flange)
        ax.add_patch(web)
        ax.add_patch(top_flange)
        ax.set_xlim(-0.5*max(self.b), 0.5*max(self.b))
        ax.set_ylim(0, self.h)
        
        zel = self.elastic_neutral_axis()
        dpl = self.plastic_neutral_axis()        
        
        print(zel)
        print(dpl)
        
        plt.plot(0.0, zel, 'or')
        plt.plot(0.0, dpl, 'db')
        
        ax.set_aspect('equal')
        
        #  plt.show()

        plt.savefig('testi' + 'plot.png')


def test_sym():
    """ symmetric I-beam """
    p = WISection(500, 8, [250, 250], [12, 12])
    print(p.MRd)
    print("Done.")
    # p.draw()


def test_non_sym():
    """ symmetric I-beam """
    print("Create section")
    p = WISection(500, 8, [300, 250], [10, 15])
    print(p.MRd)
    print("Done.")
    # p.draw()


if __name__ == '__main__':

    from copy import deepcopy
    
    #test_non_sym()
    #from frame2d.frame2d import *
    #frame = Frame2D()
    #col = SteelColumn([[0, 0], [0, 5000]])
    #frame.add(col)
    
    #p = WISection(1594, 7.0, [300.0, 300.0], [22.0, 22.0],fy=235,weld_throat=3.5)    
    p = WISection(1000, 6.0, [300.0, 300.0], [15.0, 15.0],fy=355,weld_throat=4)    
    p.Med = 2.228e6    
    p.section_class(verb=True)
    p.effective_web(verb=True)
    p.draw()
    
    #s = p
    
    s = WISection(p.h,p.tw,p.b,p.tf,p.fy,p.weld_throat)
    s.h = 1200
    
    print(s.h)
    print(p.h)
    
    #profile = "WI 200-5-200-8-100-5"
    # frame.plot()
    #print(col.cross_section.b)
    #print(col.cross_section.tf)
    #col.cross_section.b = 150
    #col.cross_section.tf = 6
    #print(col.cross_section.b)
    #print(col.cross_section.tf)
    # col.cross_section.draw()


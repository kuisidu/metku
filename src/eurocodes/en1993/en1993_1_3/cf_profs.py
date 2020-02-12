# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 16:06:13 2018

Cold-formed profiles

@author: kmela
"""

from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from materials.steel_data import Steel
from open_prof import OpenProf
from eurocodes.en1993.en1993_1_3 import en1993_1_3
from eurocodes.en1993 import en1993_1_5

class ZSection:
    """ Class for unstiffened Z-sections """
    
    def __init__(self,t_nom=2.0,h=120,a=60,b=60,ca=20,cb=20,r=2.0,
                 material="S350GD",t_coat=0.04,
                 MEd=0.0, NEd=0.0, VEd=0.0):
        """ Constructor 
            input:
                t_nom .. nominal wall thickness
                h .. total height
                a .. width of top flange
                b .. width of bottom flange                
                ca .. length of top flange edge stiffener (=0 if no stiffener)
                cb .. length of bottom flange edge stiffener (=0 if no stiffener)
                r .. rounding
                material .. 
                t_coat = 0.04 .. coating [mm]        
        """
        self.t_nom = t_nom
        self._h = h
        self._a = a
        self._b = b
        self._ca = ca
        self._cb = cb
        self.material = Steel(material)
        self.t_coat = t_coat
        self.t_core = t_nom-t_coat
        
        self.Ned = NEd
        self.Med = MEd
        self.Ved = VEd
        
        """ Reduction factors for effective cross-section """
        self._rho_t = 1.0
        self._rho_t_lip = 1.0
        self._rho_w = 1.0
        self._rho_b = 1.0
        self._rho_b_lip = 1.0
        
        if r == 0:
            """ According to EN 10162, Table A.2;
                S320GD+Z and S350GD+Z
            """
            self.r = 1.5*t_nom
        else:
            self.r = r
            
        """ Calculation models for gross cross section and
            effective cross section.
        """
        self.prof = None
        self.prof_eff = None
        
        """ Calculation models for effective top and bottom
            edge stiffeners
        """
        self.edge_stiff = {'top':None, 'bottom':None}
        
    
    def issymmetric(self):
        """ Cross-section is considered symmetric if flanges have
            equal widths and edge stiffeners have equal lengths
        """
        return self.a == self.b and self.ca == self.cb
    
    @property
    def fya(self):
        """ Average yield strength """        
        return en1993_1_3.fya(self.fyb,self.fu,self.t,Ag=self.Agross(),n=4)
        
    @property
    def fyb(self):
        """ Basic yield strength """
        return self.material.fy
    
    @property
    def fu(self):
        """ Ultimate strengt """
        return self.material.fu
    
    @property
    def t(self):
        """ Design thickness of the profile """
        return self.t_core
    
    @property
    def h(self):
        """ Design value of the height of the profile """
        return self._h
        #return self.design_value('_h')
    
    @property
    def a(self):
        """ Design value of the width of the top flange """
        return self._a
        #return self.design_value('_a')
    
    @property
    def b(self):
        """ Design value of the width of the bottom flange """
        return self._b
        #return self.design_value('_b')
    
    @property
    def ca(self):
        """ Design value of the top flange edge stiffener """
        return self._ca
        #return self.design_value('_ca')
    
    @property
    def cb(self):
        """ Design value of the bottom flange edge stiffener """
        return self._cb
        #return self.design_value('_cb')
    
    @property
    def bp_w(self):
        """ Notional width of the web """
        return en1993_1_3.notional_width(self.hw,[self.rm,self.rm],[90,90])
    
    @property
    def bp_t(self):
        """ Notional width of the top flange """
        if self.ca > 0.0:
            bp = en1993_1_3.notional_width(self.c_top,[self.rm,self.rm],[90,90])
        else:
            bp = en1993_1_3.notional_width(self.c_top,[self.rm],[90])
        return bp
    
    @property
    def bp_b(self):
        """ Notional width of the bottom flange """
        if self.cb > 0.0:
            bp = en1993_1_3.notional_width(self.c_bottom,[self.rm,self.rm],[90,90])
        else:
            bp = en1993_1_3.notional_width(self.c_bottom,[self.rm],[90])
        return bp
    
    @property
    def bp_t_lip(self):
        """ Notional width of the top flange lip"""
        if self.ca > 0.0:
            bp = en1993_1_3.notional_width(self.ca_straight,[self.rm],[90])
        else:
            bp = 0.0
        return bp
    
    @property
    def bp_b_lip(self):
        """ Notional width of the bottom flange lip """
        if self.cb > 0.0:
            bp = en1993_1_3.notional_width(self.cb_straight,[self.rm],[90])
        else:
            bp = 0.0
        return bp
    
    def design_value(self,attr):
        """ Calculates the design value of a dimension 
            This means that the coating is reduced from the dimension.
        """
        nom_val = getattr(self,attr)        
        return nom_val-self.t_coat
        
    def create_model(self):
        """ Creates calculation model according to 
            EN 1993-1-3:2006, Annex C 
        
            Origin is located at the middle of the web
        """
        bp_web = self.bp_w
        bp_bottom = self.bp_b
        bp_top = self.bp_t
        bp_bottom_lip = self.bp_b_lip
        bp_top_lip = self.bp_t_lip
        web_nodes = [[0.0,-0.5*bp_web],[0.0,0.5*bp_web]]
        if bp_bottom_lip > 0.0:
            bottom_nodes = [[-bp_bottom,web_nodes[0][1]+bp_bottom_lip],
                            [-bp_bottom,web_nodes[0][1]]
                            ]
        else:
            bottom_nodes = [[-bp_bottom,web_nodes[0][1]]]
            
        if bp_top_lip > 0.0:
            top_nodes = [[bp_top,web_nodes[1][1]],
                             [bp_top,web_nodes[1][1]-bp_top_lip]]
        else:
            top_nodes = [[bp_top,web_nodes[1][1]]]
        
        nodes = bottom_nodes + web_nodes + top_nodes
        
        self.prof = OpenProf(nodes,self.t)
                            
    
    @property
    def hw(self):
        """ Straight part of the web """
        return self.h - 2*(self.r+self.t)
    
    @property
    def hw_center(self):
        """ Height of the profile from center line of the top flange
            to the center line of the bottom flange.
        """
        return self.h - self.t_nom
    
    @property
    def c_top(self):
        """ Straight part of the top flange """
        if self.ca > 0:
            c = self.a-2*(self.r+self.t)
        else:
            c = self.a-(self.r+self.t)
        return c
    
    @property
    def c_bottom(self):
        """ Straight part of the bottom flange """
        if self.cb > 0:
            c = self.b-2*(self.r+self.t)
        else:
            c = self.b-(self.r+self.t)
        return c
    
    @property
    def ca_straight(self):
        """ Straight part of the top lip """
        if self.ca > 0:
            c = self.ca-self.t-self.r
        else:
            c = 0.0
        return c
    
    @property
    def cb_straight(self):
        """ Straight part of the bottom lip """
        if self.cb > 0:
            c = self.cb-self.t-self.r
        else:
            c = 0.0
        
        return c
    
    @property
    def rm(self):
        """ Radius to the center line """
        return self.r + 0.5*self.t
    
    def center_line_length(self):    
        straight_parts = self.cb_straight + self.c_bottom + self.hw + self.c_top + self.ca_straight 
        
        n = 2
        if self.ca > 0.0:
            n += 1
        
        if self.cb > 0.0:
            n += 1
        
        return straight_parts +n*np.pi*self.rm/2
    
    def Agross(self):
        """ Gross cross-sectional area """
        
        return self.t_nom * self.center_line_length()
    
    def effective_width_flange(self,psi,flange='top',verb=False):
        """ Effective width of a flange """
        if flange == "top":
            bp = self.bp_t
            edge_stiffener = abs(self.ca) > 1e-10
        else:
            bp = self.bp_b
            edge_stiffener = abs(self.cb) > 1e-10
        
        if edge_stiffener:
            ksigma = en1993_1_5.buckling_factor_internal(psi)
            lambda_p = en1993_1_5.lambda_p(bp,self.t,self.material.eps(),ksigma)
            rho = en1993_1_5.reduction_factor_internal(lambda_p, psi)
            beff, be = en1993_1_5.effective_width_internal(bp,psi,rho)
        else:
            ksigma = en1993_1_5.buckling_factor_outstand(psi)
            lambda_p = en1993_1_5.lambda_p(bp,self.t,self.material.eps(),ksigma)
            rho = en1993_1_5.reduction_factor_outstand(lambda_p)
            beff = en1993_1_5.effective_width_outstand(bp,psi,rho)
            
        if flange == 'top':
            self._rho_t = rho
        else:
            self._rho_b = rho
        
        if verb:
            print("** Effective {0:s} flange **".format(flange))
            print("  Notional width: bp = {0:4.2f}".format(bp))
            print("  Buckling factor: ksigma = {0:4.3f}".format(ksigma))
            print("  Slenderness: lambda_p = {0:4.3f}".format(lambda_p))
            print("  Reduction factor: rho = {0:4.3f}".format(rho))
            print("  Effective width: beff = {0:4.3f} mm".format(beff))
            
    def effective_width_edge_stiffener(self,flange='top',verb=False):
        """ Effective width of an edge stiffener """
        if flange == "top":
            bp = self.bp_t
            bpc = self.bp_t_lip
        else:
            bp = self.bp_b
            bpc = self.bt_b_lip
        
        # NOTE: Buckling factor according to EN 1993-1-3, not EN 1993-1-5!
        ksigma = en1993_1_3.buckling_factor_edge_stiffener(bpc/bp)
        lambda_p = en1993_1_5.lambda_p(bpc,self.t,self.material.eps(),ksigma)
        rho = en1993_1_5.reduction_factor_outstand(lambda_p)
        ceff = en1993_1_5.effective_width_outstand(bpc,psi=1.0,rho=rho)
            
        if flange == 'top':
            self._rho_t = rho
        else:
            self._rho_b = rho
        
        if verb:
            print("** Effective {0:s} edge stiffener **".format(flange))
            print("  Notional width of flange: bp = {0:4.2f}".format(bp))
            print("  Notional width of stiffener: bp,c = {0:4.2f}".format(bpc))
            print("  Buckling factor: ksigma = {0:4.3f}".format(ksigma))
            print("  Slenderness: lambda_p = {0:4.3f}".format(lambda_p))
            print("  Reduction factor: rho = {0:4.3f}".format(rho))
            print("  Effective width: ceff = {0:4.3f} mm".format(ceff))
        
        
    
    def effective_width_web(self,psi,verb=False):
        """ Effective width of web """
        bp = self.bp_w
        ksigma = en1993_1_5.buckling_factor_internal(psi)
        lambda_p = en1993_1_5.lambda_p(bp,self.t,self.material.eps(),ksigma)
        rho = en1993_1_5.reduction_factor_internal(lambda_p, psi)
        beff, be = en1993_1_5.effective_width_internal(bp,psi,rho)
        
        if verb:
            print("** Effective web **")
            print("  Notional width: bp = {0:4.2f} mm".format(bp))
            print("  Stress ratio: psi = {0:4.2f}".format(psi))
            print("  Buckling factor: ksigma = {0:4.3f}".format(ksigma))
            print("  Slenderness: lambda_p = {0:4.3f}".format(lambda_p))
            print("  Reduction factor: rho = {0:4.3f}".format(rho))
            print("  Effective width: beff = {0:4.3f} mm".format(beff))
    
    
    def area(self):
        """ Gross cross-sectional area """
        
        r = [self.r,self.r]
        phi = [90,90]
        
        bp = [self.bp_b, self.bp_w, self.bp_t]
        
        if self.ca > 0.0:
            r.append(self.r)
            phi.append(90)
            bp.append(self.bp_t_lip)
        
        if self.cb > 0.0:
            r.append(self.r)
            phi.append(90)
            bp.append(self.bp_b_lip)
        
        
        d = en1993_1_3.rounding_factor(r,phi,bp)
        
        return self.prof.area()*(1-d)
    
    def centroid(self):
        
        ygc, zgc, A = self.prof.centroid()
        
        return ygc, zgc
    
    def centroid_eff(self):
        """ Centroid of the effective cross-section """
        
        ygc, zgc, A = self.prof_eff.centroid()
        
        return ygc, zgc
    
    def second_moment(self):
        """ second moment of area for the gross cross section 
            with respect to 
        """
    
    def effective_section(self,verb=False):
        """ Returns effective section, including local buckling
            and distortional buckling.
            
            The effective section contains:
                _rho_t .. effective top flange
                _rho_w .. effective web
                _rho_b .. effective bottom flange
                _rho_t_lip .. effective top edge stiffener
                _rho_t_bottom .. effective bottom edge stiffener
                
                _chid_t .. reduction factor for top edge stiffener for
                            distortional buckling
                _chid_b .. reduction factotr for bottom edge stiffener
                            for distortional buckling
        """
        
        if self.Med > 0.0:
            """ Positive bending:
                top flange in compression
                bottom flange in tension
            """
            self.effective_width_flange(psi=1.0,flange='top',verb=verb)
            """ Determine stress ratio for the web, using effecive flange """
            self.create_effective_model()
            psi = self.stress_ratio()
            
            self.effective_width_web(psi,verb=verb)
            
    def stress_ratio(self):
        """ Calculates the ratio of the top flange stress and
            bottom flange stress
        """
        
        if self.prof_eff is None:
            ygc, zgc = self.centroid()
            print("No effective model created")
        else:
            ygc, zgc = self.centroid_eff()
        
        """ Height of the compressed part of the profile
            NOTE:
                zgc is determined using the simplified computational model,
                where the height of the profile is smaller than 'h', namely
                the thickness of the profile is reduced as notional widths
                are employed.
                
                Maybe the height of the computational profile should be
                used instead of 'h'
        """
        
        # Use hw_center to go from center line to center line.
        hc = 0.5*self.hw_center - zgc
        
        print("Compressed part of the profile: {0:4.2f} mm".format(hc))
        
        return -(self.hw_center/hc-1)

    def edge_stiffness(self,flange='top'):
        """ EN 1993-1-3, Eq. (5.10b)
            Stiffness provided by the profile to the edge stiffener
        """
        
        C = 0.25*self.material.E*self.t**3/(1-self.material.nu**2)        
        hw = self.hw_center
        
        """ TODO:
            Implement b1 and b2 (distance from far edge of the flange to
            the centroid of effective stiffener)
        """
        
        if flange == 'top':
            f1 = 'top'
            f2 = 'bottom'
        else:
            f1 = 'bottom'
            f2 = 'top'
        
        if self.Ned < 0.0:
            if self.issymmetric():
                kf = 1.0
            else:
                kf = self.edge_stiff[f2].A/self.edge_stiff[f1].A
        elif self.Med > 0.0 and flange == 'top':            
            kf = 0.0
        elif self.Med < 0.0 and flange == 'bottom':
            kf = 0.0
        
        return C/(b1**2*hw + b1**3 + 0.5*b1*b2*hw*kf)
        
        
    def draw(self):
        """ Draw the profile 
            
        """        
        
        fig, ax = plt.subplots(1)        
        
        lw = 1.5    # Line width of the outer line
        lw_c = 1.25 # Line width of the center line
        col = 'b'
        
        h = self.h
        hweb = self.hw
        b_top = self.a
        b_bot = self.b
        ca = self.ca
        cb = self.cb
        tnom = self.t_nom
        r = self.r
        
        # Draw webs
        y_web = [-0.5*hweb,0.5*hweb]                            
        ax.vlines([-0.5*tnom,0.5*tnom],y_web[0],y_web[1],colors='k',linewidth=lw)
        ax.vlines([0.0],y_web[0],y_web[1],colors='b',linewidth=lw_c,linestyle='-.')
        
        
        # Draw top flange
        ctop = self.c_top # Straight part of the top flange
        x_top = [0.5*tnom+r,0.5*tnom+r+ctop]
        y_top = [0.5*h,0.5*h-tnom]
        
        ax.hlines(y_top,x_top[0],x_top[1],colors='k',linewidth=lw)
        ax.hlines([0.5*(h-tnom)],x_top[0],x_top[1],colors='b',linewidth=lw_c,linestyle='-.')
        
        # Corner from web to top flange
        x_web = (0.5*tnom+r,0.5*h-tnom-r)
        
        up_corner_in = patches.Arc(x_web,2*r,2*r,0,90,180,color='k',linewidth=lw)
        up_corner_out = patches.Arc(x_web,2*(r+tnom),2*(r+tnom),0,90,180,color='k',linewidth=lw)
        up_corner_mid = patches.Arc(x_web,2*r+tnom,2*r+tnom,0,90,180,color='b',linewidth=lw_c,linestyle='-.')
        
        ax.add_patch(up_corner_in)
        ax.add_patch(up_corner_out)
        ax.add_patch(up_corner_mid)
        
        # Draw top lip
        if self.ca > 0.0:
            y_top_lip = [0.5*hweb-(ca-tnom-r),0.5*hweb]                            
            ax.vlines([b_top-1.5*tnom,b_top-0.5*tnom],y_top_lip[0],y_top_lip[1],colors='k',linewidth=lw)
            ax.vlines([b_top-tnom],y_top_lip[0],y_top_lip[1],colors='b',linewidth=lw_c,linestyle='-.')
        
            # Corner of top lip
            x_top_right = (b_top-1.5*tnom-r,0.5*h-tnom-r)
        
            top_right_corner_in = patches.Arc(x_top_right,2*r,2*r,0,0,90,color='k',linewidth=lw)
            top_right_corner_out = patches.Arc(x_top_right,2*(r+tnom),2*(r+tnom),0,0,90,color='k',linewidth=lw)
        
            ax.add_patch(top_right_corner_in)
            ax.add_patch(top_right_corner_out)
        
            ax.hlines(0.5*hweb-(ca-tnom-r),b_top-1.5*tnom,b_bot-0.5*tnom,colors='k',linewidth=lw)
        
        # Bottom flange
        cbottom = self.c_bottom # Straight part of the top flange
        x_bot = [-(0.5*tnom+r),-(0.5*tnom+r+cbottom)]
        y_bot = [-(0.5*h-tnom),-0.5*h]
        
        ax.hlines(y_bot,x_bot[0],x_bot[1],colors='k',linewidth=lw)
        ax.hlines([-0.5*(h-tnom)],x_bot[0],x_bot[1],colors='b',linewidth=lw_c,linestyle='-.')
        
        # Corner from web to bottom flange
        x_web = (-(0.5*tnom+r),-(0.5*h-tnom-r))
        
        bot_corner_in = patches.Arc(x_web,2*r,2*r,0,270,360,color='k',linewidth=lw)
        bot_corner_out = patches.Arc(x_web,2*(r+tnom),2*(r+tnom),0,270,360,color='k',linewidth=lw)
        bot_corner_mid = patches.Arc(x_web,2*r+tnom,2*r+tnom,0,270,360,color='b',linewidth=lw_c,linestyle='-.')
        
        ax.add_patch(bot_corner_in)
        ax.add_patch(bot_corner_out)
        ax.add_patch(bot_corner_mid)
        
        # Left lip      
        if self.cb > 0.0:
            y_bot_lip = [-(0.5*hweb-(cb-tnom-r)),-(0.5*hweb)]                            
            ax.vlines([-(b_bot-1.5*tnom),-(b_bot-0.5*tnom)],y_bot_lip[0],y_bot_lip[1],colors='k',linewidth=lw)
            ax.vlines([-(b_bot-tnom)],y_bot_lip[0],y_bot_lip[1],colors='b',linewidth=lw_c,linestyle='-.')
            
            # Corner of the lip
            x_bot_left = (-(b_bot-1.5*tnom-r),-(0.5*h-tnom-r))
            
            bot_left_corner_in = patches.Arc(x_bot_left,2*r,2*r,0,180,270,color='k',linewidth=lw)
            bot_left_corner_out = patches.Arc(x_bot_left,2*(r+tnom),2*(r+tnom),0,180,270,color='k',linewidth=lw)
            bot_left_corner_mid = patches.Arc(x_bot_left,2*r+tnom,2*r+tnom,0,180,270,color='b',linewidth=lw_c,linestyle='-.')
            
            ax.add_patch(bot_left_corner_in)
            ax.add_patch(bot_left_corner_out)
            ax.add_patch(bot_left_corner_mid)
            
            ax.hlines(-(0.5*hweb-(cb-tnom-r)),-b_bot+0.5*tnom,-b_bot+1.5*tnom,colors='k',linewidth=lw)
        
        
        
        ax.set_aspect('equal')
        
class CSection(ZSection):
    """ C section
        This profile has many things in common with the Z section,
        so it inherits the basic functionalities
    """
    
    def __init__(self,t_nom=2.0,h=120,a=60,b=60,ca=20,cb=20,r=2.0,material="S350GD",t_coat=0.04):
        
        super().__init__(t_nom,h,a,b,ca,cb,r,material,t_coat)
    
    def create_model(self):
        """ Creates calculation model according to 
            EN 1993-1-3:2006, Annex C 
        
            Origin is located at the middle of the web
        """
        bp_web = self.bp_w
        bp_bottom = self.bp_b
        bp_top = self.bp_t
        bp_bottom_lip = self.bp_b_lip
        bp_top_lip = self.bp_t_lip
        web_nodes = [[0.0,-0.5*bp_web],[0.0,0.5*bp_web]]
        if bp_bottom_lip > 0.0:
            bottom_nodes = [[bp_bottom,web_nodes[0][1]+bp_bottom_lip],
                            [bp_bottom,web_nodes[0][1]]
                            ]
        else:
            """ No bottom edge stiffener """
            bottom_nodes = [[bp_bottom,web_nodes[0][1]]]
        
        if bp_top_lip > 0.0:
            top_nodes = [[bp_top,web_nodes[1][1]],
                             [bp_top,web_nodes[1][1]-bp_top_lip]]
        else:
            top_nodes = [[bp_top,web_nodes[1][1]]]
            
        nodes = bottom_nodes + web_nodes + top_nodes
        
        self.prof = OpenProf(nodes,self.t)
    
    def create_effective_model(self):
        """ Create calculation model for the effective cross-section 
            This is done based on the model of the gross cross-section
        """
    
        if self.prof is None:
            self.create_model()
        
        prof_eff = deepcopy(self.prof)
        
        if self._rho_t < 1.0:
            """ Effective top flange
                Top flange nodes are the last nodes in the model
            """
            if self.ca > 0.0:
                """ The effective top flange is an internal element """
                pass
            else:
                """ The effective top flange is an outstand element """
                prof_eff.nodes[-1].y -= (1-self._rho_t)*self.bp_t
                
        
        self.prof_eff = prof_eff
            
    

    def draw(self):
        """ Draw the profile 
            
        """        
        
        fig, ax = plt.subplots(1)        
        
        lw = 1.5    # Line width of the outer line
        lw_c = 1.25 # Line width of the center line
        col = 'b'
        
        h = self.h
        hweb = self.hw
        b_top = self.a
        b_bot = self.b
        ca = self.ca
        cb = self.cb
        tnom = self.t_nom
        r = self.r
        
        # Draw webs
        y_web = [-0.5*hweb,0.5*hweb]                            
        ax.vlines([-0.5*tnom,0.5*tnom],y_web[0],y_web[1],colors='k',linewidth=lw)
        ax.vlines([0.0],y_web[0],y_web[1],colors='b',linewidth=lw_c,linestyle='-.')
        
        
        # Draw top flange
        ctop = self.c_top # Straight part of the top flange
        x_top = [0.5*tnom+r,0.5*tnom+r+ctop]
        y_top = [0.5*h,0.5*h-tnom]
        
        ax.hlines(y_top,x_top[0],x_top[1],colors='k',linewidth=lw)
        ax.hlines([0.5*(h-tnom)],x_top[0],x_top[1],colors='b',linewidth=lw_c,linestyle='-.')
        
        # Corner from web to top flange
        x_web = (0.5*tnom+r,0.5*h-tnom-r)
        
        up_corner_in = patches.Arc(x_web,2*r,2*r,0,90,180,color='k',linewidth=lw)
        up_corner_out = patches.Arc(x_web,2*(r+tnom),2*(r+tnom),0,90,180,color='k',linewidth=lw)
        up_corner_mid = patches.Arc(x_web,2*r+tnom,2*r+tnom,0,90,180,color='b',linewidth=lw_c,linestyle='-.')
        
        ax.add_patch(up_corner_in)
        ax.add_patch(up_corner_out)
        ax.add_patch(up_corner_mid)
        
        # Draw top lip
        if self.ca > 0.0:
            y_top_lip = [0.5*hweb-(ca-tnom-r),0.5*hweb]                            
            ax.vlines([b_top-1.5*tnom,b_top-0.5*tnom],y_top_lip[0],y_top_lip[1],colors='k',linewidth=lw)
            ax.vlines([b_top-tnom],y_top_lip[0],y_top_lip[1],colors='b',linewidth=lw_c,linestyle='-.')
            
            # Corner of top lip
            x_top_right = (b_top-1.5*tnom-r,0.5*h-tnom-r)
            
            top_right_corner_in = patches.Arc(x_top_right,2*r,2*r,0,0,90,color='k',linewidth=lw)
            top_right_corner_out = patches.Arc(x_top_right,2*(r+tnom),2*(r+tnom),0,0,90,color='k',linewidth=lw)
            
            ax.add_patch(top_right_corner_in)
            ax.add_patch(top_right_corner_out)
            
            ax.hlines(0.5*hweb-(ca-tnom-r),b_top-1.5*tnom,b_bot-0.5*tnom,colors='k',linewidth=lw)
        
        # Bottom flange
        cbottom = self.c_bottom # Straight part of the top flange
        x_bot = [0.5*tnom+r,0.5*tnom+r+cbottom]
        y_bot = [-(0.5*h-tnom),-0.5*h]
        
        ax.hlines(y_bot,x_bot[0],x_bot[1],colors='k',linewidth=lw)
        ax.hlines([-0.5*(h-tnom)],x_bot[0],x_bot[1],colors='b',linewidth=lw_c,linestyle='-.')
        
        # Corner from web to bottom flange
        x_web = (0.5*tnom+r,-(0.5*h-tnom-r))
        
        bot_corner_in = patches.Arc(x_web,2*r,2*r,0,180,270,color='k',linewidth=lw)
        bot_corner_out = patches.Arc(x_web,2*(r+tnom),2*(r+tnom),0,180,270,color='k',linewidth=lw)
        bot_corner_mid = patches.Arc(x_web,2*r+tnom,2*r+tnom,0,180,270,color='b',linewidth=lw_c,linestyle='-.')
        
        ax.add_patch(bot_corner_in)
        ax.add_patch(bot_corner_out)
        ax.add_patch(bot_corner_mid)
        
        # Left lip        
        if self.cb > 0.0:
            y_bot_lip = [-(0.5*hweb-(cb-tnom-r)),-(0.5*hweb)]                            
            ax.vlines([b_bot-1.5*tnom,b_bot-0.5*tnom],y_bot_lip[0],y_bot_lip[1],colors='k',linewidth=lw)
            ax.vlines([b_bot-tnom],y_bot_lip[0],y_bot_lip[1],colors='b',linewidth=lw_c,linestyle='-.')
            
            # Corner of the lip
            x_bot_left = (b_bot-1.5*tnom-r,-(0.5*h-tnom-r))
            
            bot_left_corner_in = patches.Arc(x_bot_left,2*r,2*r,0,270,360,color='k',linewidth=lw)
            bot_left_corner_out = patches.Arc(x_bot_left,2*(r+tnom),2*(r+tnom),0,270,360,color='k',linewidth=lw)
            bot_left_corner_mid = patches.Arc(x_bot_left,2*r+tnom,2*r+tnom,0,270,360,color='b',linewidth=lw_c,linestyle='-.')
            
            ax.add_patch(bot_left_corner_in)
            ax.add_patch(bot_left_corner_out)
            ax.add_patch(bot_left_corner_mid)
            
            ax.hlines(-(0.5*hweb-(cb-tnom-r)),b_bot-0.5*tnom,b_bot-1.5*tnom,colors='k',linewidth=lw)
                
        ax.set_aspect('equal')

class ZProf:
    """ ZProf Class for open Z-sections
        Main purpose is to calculate cross-sectional properties as
        described in EN 1993-1-3, Annex C.    
    """
    """    
        properties (Depent = true)
            hWeb % Web height
            bTop % Width of Top flange
            bBottom % Width of Bottom flange
            gr % length needed for notional widths
            rm % corner rounding to the midline
        
    """
    def __init__(self,Tnom,H,A,B,CA,CB,R=0,fyb=350,fub=420,Tcoat=0.04):
        """ Constructor 
            input:
                Tnom .. nominal wall thickness
                H .. total height
                A .. width of top flange
                B .. width of bottom flange
                CA .. length of top fold
                CB .. length of bottom fold
                R .. rounding
                fyb = 350 .. basic yield strength [MPa]
                fub = 420 .. basic ultimate strength [MPa]
                Tcoat = 0.04 .. coating [mm]        
        """
        
        """
            %rounding_method = 1; % corner node at intersection of midlines
            %rounding_method = 2; % corner by two elements
            rounding_method = 3; % corner by three elements
        """
        self.tnom = Tnom
        self.h = H
        self.a = A
        self.b = B
        self.ca = CA
        self.cb = CB            
        self.fyb = fy
        self.fub = fub
        self.tcoat = Tcoat
        self.tcore = Tnom-Tcoat
            
        if R == 0:
            """ According to EN 10162, Table A.2;
                S320GD+Z and S350GD+Z
            """
            self.r = 1.5*T
        else:
            self.r = R
            
        # store profile data here (instance of OpenProf class)
        self.prof = None
            
        """ Create nodes.
            a) The origin is in the intersection of
               midlines of bottom flange and web.
            b) corner roundings are taken into account by the "German
               approach"
               The nodes are numbered starting from the tip of
               the fold of the bottom flange. If CB = 0, there is
               no fold
            switch rounding_method
                case 1
                    y = [B-Tnom; 0; 0; -(A-Tnom)];
                    z = [0; 0; H-Tnom; H-Tnom];
                    if nargin > 4 && ~isempty(CA) && CA > 0
                        y = [y; -(A-Tnom)];
                        z = [z; H-Tnom - (CA-0.5*Tnom)];
                    
                    if nargin > 5 && ~isempty(CB) && CB > 0
                        y = [B-Tnom;y];
                        z = [CB-0.5*Tnom;z];
                    
                case 2
                case 3                    
                    rm = R+0.5*Tnom;
                    %  bottom flange edge stiffener
                    % (assume 90 degree angle)
                    if CB > 0                      
                        yBottomLip = B-Tnom;                      
                        y = [yBottomLip*[1;1;1];yBottomLip-0.7*rm;yBottomLip-rm];                      
                        z = [CB-0.5*Tnom;rm;0.7*rm;0;0];
                    else
                        y = B-0.5*Tnom;
                        z = 0;
                    
                    % bottom flange corner
                    y = [y;rm;0.7*rm;0;0];
                    z = [z;0;0;0.7*rm;rm];
                    % web;
                    hw = H-Tnom;
                    y = [y;0;0];
                    z = [z;hw-rm;hw-0.7*rm];
                    % top flange corner
                    y = [y;-0.7*rm;-rm];
                    z = [z;hw;hw];
                    % top flange                    
                    if CA > 0
                        yTopLip = A-Tnom;
                        y = [y;-(yTopLip-rm);-(yTopLip-0.7*rm);-yTopLip*[1;1;1]];
                        z = [z;hw;hw;hw-0.7*rm;hw-rm;hw-(CA-0.5*Tnom)];
                    else
                        y = [y;-(A-0.5*Tnom)];
                        z = [z;hw];                       
                    
                    %size(y)
                    %size(z)
            """        
            
        """
        def self = set.Tnom(self,val)
            self.Tnom = val;
            self.T = val-self.Tcoat;
            self.t(:) = self.T;
        
        
        def self = set.A(self,val)
            self.A = val;
            if self.CA > 0
                rm = self.rm;
                y0 = val-self.Tnom;
                self.y(-4:) = -[y0-rm,y0-0.7*rm,y0,y0,y0];
            else
                self.y() = -(val-0.5*self.Tnom);
        """
        
                
        def hWeb(self):
            """ Web heigth """
            return self.h-self.tnom
        
        
        def bTop(self):
            """ Width of the top flange (center line) """
            if self.ca == 0:
                r = self.a-0.5*self.tnom
            else:
                r = self.a-self.tnom
            return r
                    
        def bBottom(self):
            """ Width of the bottom flange (center line) """
            if self.cb == 0:
                r = self.b-0.5*self.tnom
            else:
                r = self.b-self.tnom
            return r
            
        def rounding_center(self):
            """ Corner rounding to the center line """
            return self.r+0.5*self.tnom
                
        def gr(self):
            """ Length gr from EN 1993-1-3
            Rm = self.rounding_center()
            a = 90
            gr = Rm*(tand(0.5*a)-sind(0.5*a))
            return gr
                
        def center_line_length(self):
            return self.ca+self.a+self.h+self.b+self.cb-8*(self.r+self.tnom)+2*math.pi()*(self.r+0.5*self.tnom)
                
        def average_yield(self):
            k = 7  # coefficient for roll forming
            n = 4  # number of 90 degree bs with r <= 5t
            fave = 0.5*(self.fyb+self.fub)
            fya = min(self.fyb + (self.fub-self.fyb)*k*n*self.t**2/self.Area,fave)
            
        
        
        def stCom = StressTopFlange(self)
            gc = self.Centroid;
            zgc = gc(2);
            zb1 = self.hWeb-zgc;
            r = min(zb1/(self.hWeb-zb1),1.0);
            stCom = self.fyb*r;
        
        
        % computes the notional width of a flange
        % input:
        %        C -- length of edge fold
        %  Bflange -- flange width
        def bp = NotionalWidthTopFlange(self)             
            bp = self.NotionalWidthFlange(self.A,self.CA,self.Tnom,self.gr);
        
        
        def bp = NotionalWidthBottomFlange(self)
            bp = self.NotionalWidthFlange(self.B,self.CB,self.T,self.gr);
        
        
        def bcp = NotionalWidthTopLip(self) 
            if self.CA > 1e-6
                bcp = self.CA-0.5*self.T-self.gr;
            else
                bcp = 0;
            
        
        
        def bcp = NotionalWidthBottomLip(self) 
            if self.CB > 1e-6
                bcp = self.CB-0.5*self.T-self.gr;
            else
                bcp = 0;
            
        
        
        def bp = NotionalWidthWeb(self)             
            bp = self.hWeb-2*self.gr;            
        
        
        
        def [beff,be,rho] = EffectiveFlange(self,flange,e)
            if nargin < 2 || isempty(flange)
                flange = 'top';
            
            
            switch flange
                case {'top',1}   
                    bp = self.NotionalWidthTopFlange;
                    C = self.CA;
%                     C = self.CA;
%                     Bflange = self.A;
                case {'bottom',0}
                    bp = self.NotionalWidthBottomFlange;
                    C = self.CB;
%                     C = self.CB;
%                     Bflange = self.B;
            
            
           % bp = self.NotionalWidthFlange(C,Bflange);
                        
            r = bp/self.T;    
            if nargin < 3
                e = Epsilon(self.fyb);
            
                                                    
            if abs(C) < 1e-6
                % no edge stiffener
                % assume constant stress in top flange:
                % this implies that stress ratio is 1.
                stressRatio = 1.0;
                ksigma = BucklingFactorOutstand(stressRatio);
                slred = PlateSlerness(r,e,ksigma);
                rho = OutstandReduction(slred);
                beff = rho*bp;
                be = beff;
            else
                % assume constant stress in top flange:
                % this implies that stress ratio is 1.
                stressRatio = 1.0;
                ksigma = BucklingFactorInternal(stressRatio);
                slred = PlateSlerness(r,e,ksigma);
                rho = InternalReduction(slred,stressRatio);
                beff = rho*bp;
                be = 0.5*beff*[1;1];
            
                        
            %self.DrawEffectiveFlange(flange,be(1),bp,stressRatio)
            %self.DrawEffectiveFlange(flange,be(1),bp,rho);            
        
        
        def [ceff,rho,bcp] = EffectiveLip(self,flange,e)
            verb = 0;
           
            if nargin < 2 || isempty(flange)
                flange = 'top';
            
            
            switch flange
                case {'top',1}               
                    bp = self.NotionalWidthTopFlange;
                    C = self.CA;
                    %Bflange = self.A;
                case {'bottom',0}
                    bp = self.NotionalWidthBottomFlange;
                    C = self.CB;
                    %Bflange = self.B;
            
            % unstiffened top flange:
            % plane element without stiffener
            %bp = self.NotionalWidthFlange(C,Bflange);
            bcp = self.NotionalWidthLip(C);
            r = bcp/self.T;                
            if nargin < 3
                e = Epsilon(self.fyb);
            
                                    
            rb = bcp/bp;
            ksigma = BucklingFactorLip(rb);            
            slred = PlateSlerness(r,e,ksigma);
            rho = OutstandReduction(slred);
            ceff = rho*bcp;
            
            if verb
                fprintf(1,'Tehollinen reunajäykiste\n');
                fprintf(1,'bpc = %6.3f\n',bcp);
                fprintf(1,'rc = %6.3f\n',r);
                fprintf(1,'bcp/bp = %6.3f\n',rb);
                fprintf(1,'ksigma = %6.3f\n',ksigma);
                fprintf(1,'slred = %6.3f\n',slred);
                fprintf(1,'rho = %6.3f\n',rho);
                fprintf(1,'ceff = %6.3f\n',ceff);
            
        
        
        def bcp = NotionalWidthLip(self,C)            
            bcp = C-0.5*self.Tnom-self.gr;
        
        
        def d = ReductionDelta(self)
            % Compute the delta to be used in reduction of 
            % cross-sectional properties taking into account
            % rounded corners
            bp(1) = self.NotionalWidthTopLip;
            bp(2) = self.NotionalWidthTopFlange;            
            bp(3) = self.NotionalWidthWeb;
            bp(4) = self.NotionalWidthBottomFlange;
            bp(5) = self.NotionalWidthBottomLip;
            d = 0.43*4*self.R/sum(bp);
        
        
        % Calculates the required shear stiffness to
        % stabilize the purlin
        def Sreq = RequiredShearStiffness(self,L)
            E = 210; % GPa
            G = 81;  % GPa
            Iw = self.WarpingConstant;
            It = self.TorsionConstant;
            I2 = self.SecondMomentArea;
            Iz = I2(2);
            h = self.H;
            Sreq = (E*Iw*pi^2/L^2 + G*It + 0.25*E*Iz*pi^2/L^2*h^2)*70/h^2;
        
        
        % Calculate effective cross section in bing
        def EffectiveSectionBing(self)
            
            % Start with effective top flange
            [beff,be,rhoF] = self.EffectiveFlange('top')
            
            % Get effective lip
            [ceff,rhoL,bcp] = self.EffectiveLip('top')
            
            if rhoF < 1.0
                % here, two points should added to mark
                % the non-effective part of the flange
            else
                % if the entire flange is effective,
                % one point should be added to the flange
            
        
        
    
    
    methods (Static = true)
        def bp = NotionalWidthFlange(Bflange,C,T,gr)
            if abs(C) < 1e-6
                b = Bflange-0.5*T;                
                bp = b-gr;
            else
                b = Bflange-T;
                bp = b-2*gr;
               
        
    """
    
def Cexample():
    
    c = CSection(t_nom=2.0,h=150.0,a=47,b=41,ca=0.0,cb=0.0,r=3.0,t_coat=0.04)
    
    c.Med = 10.0
    #c.effective_width_flange(psi=1.0, flange='top',verb=True)
    c.effective_section(True)
    #c.create_model()
    #c.create_effective_model()
    
    #c.prof_eff.draw()

    #ygc, zgc = c.centroid_eff()
    #print(ygc,zgc)

    return c

if __name__ == "__main__":
    
    c = Cexample()
    
    """
    z = ZSection(ca=0.0,cb=0.0)
    c = CSection(1.48,97.3,37.0,37.0,12.5,12.5,0.85,t_coat=0.0)
    c.material.fy = 505.0
    #c.draw()
    z.create_model()
    z.prof.draw()
    c.create_model()
    #c.prof.draw()
    c.effective_width_lip(1.0,'top',verb=True)
    c.effective_width_lip(1.0,'bottom',verb=True)
    c.effective_width_flange(1.0,'top',verb=True)
    c.effective_width_flange(1.0,'bottom',verb=True)
    c.effective_width_web(1.0,verb=True)
    """
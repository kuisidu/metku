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
from eurocodes.en1993.en1993_1_3.open_prof import OpenProf
from eurocodes.en1993.en1993_1_3 import en1993_1_3
from eurocodes.en1993 import en1993_1_5
from eurocodes.en1993.constants import gammaM0, gammaM1

weckman_c_z = {100:{'a':42,'b':48,'c':18,'r':4,'t':[1,1.2,1.5,2],'width':[216,216,212,210]},
             120:{'a':42,'b':48,'c':18,'r':4,'t':[1,1.2,1.5,2],'width':[238,236,234,230]},
             150:{'a':42,'b':48,'c':18,'r':4,'t':[1.2,1.5,2],'width':[266,264,260]},
             200:{'a':42,'b':48,'c':18,'r':4,'t':[1.5,2,2.5],'width':[312,310,305]},
             250:{'a':65,'b':73,'c':21.5,'r':4,'t':[1.5,2,2.5],'width':[419,415,411]}
             }

def ruukki_z(h,t):
    
    if h <= 150:
        c = 18.0
        r = 2.0
        if t == 1.0:
            a = 45.0
            b = 39.0
        elif t == 1.2:
            a = 45.4
            b = 39.4
        elif t == 1.5:
            a = 46.0
            b = 40.0
        elif t == 2.0:
            a = 47.0
            b = 41.0
    elif h <= 250:
        c = 26.0
        r = 2.0
        if t == 1.5:
            a = 70.0
            b = 62.0
        elif t == 2.0:
            a = 71.0
            b = 63.0
        elif t == 2.5:
            a = 72.0
            b = 64.0
        elif t == 3.0:
            a = 73.0
            b = 65.0
    elif h == 300:
        c = 26.0
        r = 1.0
        if t == 1.5:
            a = 89.0
            b = 81.0
        elif t == 2.0:
            a = 90.0
            b = 82.0
        elif t == 2.5:
            a = 91.0
            b = 83.0
        elif t == 3.0:
            a = 92.0
            b = 84.0
    elif h == 350:
        c = 30.0
        r = 2.0
        if t == 2.0:
            a = 90.0
            b = 82.0
        elif t == 2.5:
            a = 91.0
            b = 83.0
        elif t == 3.0:
            a = 92
            b = 84
            
    z = ZSection(t_nom=t,h=h,a=a,b=b,ca=c,cb=c,r=r)
    
    return z
        
def corner_points(xc,r,phi=[0,np.pi],n=6):
    """ Creates corner points
        :param xc: center point
        :param r: radius of the circle
        :param phi: [start_angle,end_angle] (radians)
        :param n: number of division points
    """
    angles = np.linspace(phi[0],phi[1],n)
    return [xc + r*np.array([np.cos(angle),np.sin(angle)]) for angle in angles]


class ZSection:
    """ Class for unstiffened Z-sections """
    
    def __init__(self,t_nom=2.0,h=120,a=60,b=60,ca=20,cb=20,r=2.0,
                 material="S350GD",t_coat=0.04,
                 MEd=0.0, NEd=0.0, VEd=0.0,corners=True):
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
        
        """ Flag: if True, then rounded corners are taken in to account in cross-section
            modelling. if False, then cross-section is modelled using centerline intersection
            points.
        """
        self.rounded_corners = corners
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
        
        self._points = []
        
        """ Reduction factors for effective cross-section """
        self._lambda_t = 0.1
        self._rho_t = 1.0
        self._lambda_t_lip = 0.1
        self._rho_t_lip = 1.0
        self._lambda_w = 0.1
        self._rho_w = 1.0
        self._psi_w = 1.0 # Stress ratio for web
        self._lambda_b = 0.1
        self._rho_b = 1.0
        self._lambda_b_lip = 0.1
        self._rho_b_lip = 1.0
        """ Reduction factors for distortional buckling """
        self._chi_d_t = 1.0
        self._chi_d_b = 1.0
        
        if r == 0:
            """ According to EN 10162, Table A.2;
                S320GD+Z and S350GD+Z
            """
            self.r = 1.5*t_nom
        else:
            self.r = r
            
        
        if corners is False: 
            if self.r/self.t <= 5 and self.r/self.bp_t <= 0.1 and self.r/self.bp_b <= 0.1:
                self.rounded_corners = False
            else:
                self.rounded_corners = True
        
        """ Calculation models for gross cross section and
            effective cross section.
        """
        
        self.parts = {'top_flange':{'nodes':[],'segments':[]},
                      'web':{'nodes':[],'segments':[]},
                      'bottom_flange':{'nodes':[],'segments':[]},                      
                      'top_lip':{'nodes':[],'segments':[]},
                      'bottom_lip':{'nodes':[],'segments':[]},
                      }
        
        self.prof = None
        self.create_points(n=5)
        self.create_model()
        self.prof_eff = None
        
        """ Calculation models for effective top and bottom
            edge stiffeners
        """
        self.edge_stiff = {'top':None, 'bottom':None}
        
        self.free_flange = {'top':None, 'bottom':None}
    
    
    def issymmetric(self):
        """ Cross-section is considered symmetric if flanges have
            equal widths and edge stiffeners have equal lengths
        """
        return self.a == self.b and self.ca == self.cb
    
    @property
    def npoints(self):
        """ Number of points in the section """
        
        return len(self._points)
    
    @property
    def Ro(self):
        """ Outer radius """
        return self.r + self.t_nom
    
    @property
    def Rm(self):
        """ Mid radius """
        return self.r + 0.5*self.t_nom
    
    @property
    def E(self):
        """ Young's modulus """
        return self.material.E
    
    @property
    def G(self):
        """ Shear modulus """
        return self.material.G
    
    @property
    def A(self):
        """ Gross area of cross section """
        return self.area()

    @property
    def Aeff(self):
        self.prof_eff = deepcopy(self.prof)
        self.effective_section('compression')
        return self.prof_eff.area()
    
    @property
    def Iw(self):
        """ Warping constant of the gross cross-section """
        return self.prof.warping_constant()
    
    @property
    def It(self):
        """ Torsion constant """
        return self.prof.torsion_constant()
    
    @property
    def Iy(self):
        """ Second moment of area of gross cross-section with respect
            to the major axis
        """
        return self.prof.Iy()
    
    @property
    def Iz(self):
        """ Second moment of area of gross cross-section with respect
            to the minor axis
        """
        return self.prof.Iz()
    
    @property
    def Iyz(self):
        """ Product moment of the gross cross-section """
        return self.prof.product_moment_gc()
    
    def Ifz(self,flange='bottom'):
        """ Second moment of area of the free flange """
        return self.free_flange[flange].Iz()
    
    @property
    def NRd(self):
        return self.axial_force_resistance()
    
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
        if self.rounded_corners:
            return en1993_1_3.notional_width(self.hw,[self.rm,self.rm],[90,90])
        else:
            return self.hw_center
    
    @property
    def bp_t(self):
        """ Notional width of the top flange """
        if self.rounded_corners:
            if self.ca > 0.0:            
                    bp = en1993_1_3.notional_width(self.c_top,[self.rm,self.rm],[90,90])
            else:
                    bp = en1993_1_3.notional_width(self.c_top,[self.rm],[90])            
        else:
            bp = self.bt_center
        return bp
    
    @property
    def bp_b(self):
        """ Notional width of the bottom flange """
        if self.cb > 0.0:
            if self.rounded_corners:
                bp = en1993_1_3.notional_width(self.c_bottom,[self.rm,self.rm],[90,90])
            else:
                bp = self.bb_center
        else:
            if self.rounded_corners:
                bp = en1993_1_3.notional_width(self.c_bottom,[self.rm],[90])
            else:
                bp = self.bb_center
        return bp
    
    @property
    def bp_t_lip(self):
        """ Notional width of the top flange lip"""
        if self.rounded_corners:
            if self.ca > 0.0:
                bp = en1993_1_3.notional_width(self.ca_straight,[self.rm],[90])
            else:
                bp = 0.0
        else:
            bp = self.bt_lip_center
        return bp
    
    @property
    def bp_b_lip(self):
        """ Notional width of the bottom flange lip """
        if self.rounded_corners:
            if self.cb > 0.0:
                bp = en1993_1_3.notional_width(self.cb_straight,[self.rm],[90])
            else:
                bp = 0.0
        else:
            bp = self.bb_lip_center
        return bp
    
    @property
    def shear_center(self):
        """ Coordinates of the shear center """
        y_sc, z_sc = self.prof.shear_center()
        return [y_sc, z_sc]
    
    
    
    @property
    def gr(self):
        """ Dimension gr in the rounded corner """
        return en1993_1_3.gr(self.r,self.t_nom,45)
    
    def add(self,this,part=None):
        """ Adds stuff to the section """
        
        if isinstance(this,list) or isinstance(this,np.ndarray):
            """ Item to be added is a coordinate point """
            self._points.append(this)
            if part is not None:
                self.parts[part]["nodes"].append(self.npoints-1)
    
    def design_value(self,attr):
        """ Calculates the design value of a dimension 
            This means that the coating is reduced from the dimension.
        """
        nom_val = getattr(self,attr)        
        return nom_val-self.t_coat
    
    def weight(self):
        """ Weight per unit length (kg/mm) """
        
        return self.material.rho * self.Agross()
    
    def self_weight(self):
        """ Self-weight (kN/m) """
        
        return 9.81*self.weight()
    
    def create_points(self,n=5):
        """ Creates points for the calculation model 
            :param n: number of divisions in the corner areas
            
            Number starts from the end of the flange of the extension
            and goes counter-clockwise
        """
        # Number of points; used as an index
        npt = 0
        
        self.corner_pts = n
        h = self.h
        hweb = self.hw
        
        tnom = self.t_nom
        rm = self.Rm
        ro = self.Ro
        r = self.r
        
        bp_web = self.hw_center
        bp_bottom = self.bb_center
        bp_top = self.bt_center
        bp_bottom_lip = self.bb_lip_center
        bp_top_lip = self.bt_lip_center
        web_nodes = [[0.0,-0.5*bp_web],[0.0,0.5*bp_web]]
       
        if bp_bottom_lip > 0.0:
            bottom_lip_nodes = [[-bp_bottom,web_nodes[0][1]+bp_bottom_lip],
                                [-bp_bottom,web_nodes[0][1]+bp_bottom_lip-0.5*self.cb_straight]]
            for node in bottom_lip_nodes:
                self.add(node,part='bottom_lip')
            
            # Corner points
            # Lip-to-flange corner, bottom
            x_bottom_lip = (-bp_bottom+rm,web_nodes[0][1]+rm)
        
            newPts = corner_points(x_bottom_lip,rm,phi=[np.pi,1.5*np.pi],n=n)
        
            for newPt in newPts[:-1]:            
                self.add(newPt)
            
            # The second and third nodes are for helping to cope
            # with the effective cross section.
            bottom_nodes = [[-bp_bottom+rm,web_nodes[0][1]],
                            [-bp_bottom+rm+1/3*self.c_bottom,web_nodes[0][1]],
                            [-bp_bottom+rm+2/3*self.c_bottom,web_nodes[0][1]],
                            [-bp_bottom+rm+self.c_bottom,web_nodes[0][1]]
                            ]            
        else:
            bottom_nodes = [[-bp_bottom,web_nodes[0][1]],
                            [-bp_bottom+0.5*self.c_bottom,web_nodes[0][1]],
                            [-bp_bottom+self.c_bottom,web_nodes[0][1]]
                            ]
        
        for node in bottom_nodes:
            self.add(node,part='bottom_flange')
        
        
        # Corner points
        # flange-to-web corner, bottom
        x_bottom_web = (-rm,web_nodes[0][1]+rm)
        
        newPts = corner_points(x_bottom_web,rm,phi=[1.5*np.pi,2.0*np.pi],n=n)
        
        for newPt in newPts[1:-1]:            
            self.add(newPt)
        
        # Web points
        self.add([0,-0.5*bp_web+rm],part='web')
        self.add([0,-0.5*bp_web+rm+0.25*hweb],part='web')
        self.add([0,-0.5*bp_web+rm+0.75*hweb],part='web')
        self.add([0,0.5*bp_web-rm],part='web')
        
        # Corner points
        # web-to-flange corner, top
        x_top_web = (+rm,web_nodes[1][1]-rm)
        
        newPts = corner_points(x_top_web,rm,phi=[np.pi,0.5*np.pi],n=n)
        
        for newPt in newPts[1:-1]:            
            self.add(newPt)
        
        # Top flange
        top_nodes = [[rm,web_nodes[1][1]]]
        
        if bp_top_lip > 0.0:
            top_nodes += [[0.25*bp_top,web_nodes[1][1]],
                          [0.75*bp_top,web_nodes[1][1]],
                          [bp_top-rm,web_nodes[1][1]]]
            
            # Corner
            # web-to-flange corner, top
            x_top_lip = (bp_top-rm,web_nodes[1][1]-rm)
        
            for node in top_nodes:
                self.add(node,part='top_flange')
        
            newPts = corner_points(x_top_lip,rm,phi=[0.5*np.pi,0],n=n)
        
            for newPt in newPts[1:-1]:            
                self.add(newPt)
                
            # Top lip nodes
            
            self.add([bp_top,web_nodes[1][1]-rm],'top_lip')
            self.add([bp_top,web_nodes[1][1]-0.5*bp_top_lip],'top_lip')
            self.add([bp_top,web_nodes[1][1]-bp_top_lip],'top_lip')
        else:
            top_nodes += [[0.5*bp_top,web_nodes[1][1]],[bp_top,web_nodes[1][1]]]
        
            for node in top_nodes:
                self.add(node,part='top_flange')
                
    
    def create_model(self):
        """ Creates calculation model according to 
            EN 1993-1-3:2006, Annex C 
        
            Origin is located at the middle of the web
        """
        
        """
        bp_web = self.bp_w
        bp_bottom = self.bp_b
        bp_top = self.bp_t
        bp_bottom_lip = self.bp_b_lip
        bp_top_lip = self.bp_t_lip
        web_nodes = [[0.0,-0.5*bp_web],[0.0,0.5*bp_web]]
        """
        
        self.prof = OpenProf(self._points,self.t)
        
        """ Set nodes and segments to parts """
        for i, node in enumerate(self.prof.nodes):
            x = node.coord
        
        """
        bp_web = self.hw_center
        bp_bottom = self.bb_center
        bp_top = self.bt_center
        bp_bottom_lip = self.bb_lip_center
        bp_top_lip = self.bt_lip_center
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
        """
        
    def create_effective_model(self):
        """ Create calculation model for the effective cross-section 
            This is done based on the model of the gross cross-section
        """
    
        if self.prof is None:
            self.create_model()
                    
        prof_eff = deepcopy(self.prof)
        
        if self._rho_w < 1.0:
            """ Effective web """            
            if self._psi_w < 0.0:
                hc = self.hw_center/(1-self._psi_w)
                beff = self._rho_w * hc
                be1 = 0.4*beff
                be2 = 0.6*beff                
            elif self._psi_w < 1.0:
                hc = self.hw_center
                beff = self._rho_w * self.hw_center
                be1 = 2*beff/(5-self._psi_w)
                be2 = beff-be1            
            else:
                hc = self.hw_center
                beff = self._rho_w * self.hw_center
                be1 = 0.5*beff
                be2 = 0.5*beff
        
            bneg = hc-beff
            # Part of the web in tension
            ht = self.hw_center-hc
        
            if self.cb > 0:
                nweb = 2
            else:
                nweb = 1
            
            nw_top = prof_eff.parts['web']['nodes'][2]
            nw_bot = prof_eff.parts['web']['nodes'][1]
            
            
            prof_eff.nodes[nw_top].z = prof_eff.nodes[nw_top+1].z + self.gr - be1
            prof_eff.nodes[nw_bot].z = prof_eff.nodes[nw_top].z - bneg
            
            prof_eff.segments[nw_bot].t = 0.0
            
            #s = (be2+ht)/self.hw_center
            #prof_eff.split_segment(nweb,s)
            #s2 = bneg/(be1+bneg)
            #prof_eff.split_segment(nweb+1,s2)
            #prof_eff.segments[nweb+1].t = 0.0
        
        if self.ca > 0.0:            
            """ Edge stiffener is present """
            # Number of segments in the section            
            #nseg = prof_eff.n-1
            
            # Move the node of the edge stiffener up
            #prof_eff.nodes[-1].z += (1-self._rho_t_lip)*self.bt_lip_center
                    
            # Split the top flange into two parts
            #s = 0.5*self._rho_t            
            #prof_eff.split_segment(nseg-2,s)
                        
            #prof_eff.draw(node_labels=True)
            if self._rho_t < 1.0:
                s2 = (1-self._rho_t)*self.bt_center/prof_eff.segments[-2].length()
                prof_eff.split_segment(nseg-1,s2)
                prof_eff.segments[nseg-1].t = 0.0
            prof_eff.segments[-2].t *= self._chi_d_t
            prof_eff.segments[-1].t *= self._chi_d_t
                        
        else:
            """ The effective top flange is an outstand element """
            #prof_eff.nodes[-2].y -= (1-self._rho_t)*self.bp_t + self.gr
            prof_eff.nodes[-2].y = self._rho_t*self.bp_t + self.gr
            prof_eff.segments[-1].t = 0.0

        if self.cb > 0.0:
            """ Edge stiffener is present in the bottom flange """
            
            # Move the node of the edge stiffener down
            prof_eff.nodes[0].z -= (1-self._rho_b_lip)*self.bb_lip_center
            
            # Split the bottom flange into two parts
            s = 0.5*self._rho_b
            prof_eff.split_segment(1,s)
                        
            """ If there is a non-effective part of the bottom flange
                remove that too.
            """
            if self._rho_b < 1.0:
                s2 = (1-self._rho_b)*self.bb_center/prof_eff.segments[2].length()
                prof_eff.split_segment(2,s2)
                prof_eff.segments[2].khp = 0.0
            prof_eff.segments[0].t *= self._chi_d_b
            prof_eff.segments[1].t *= self._chi_d_b

        self.prof_eff = prof_eff
                            
    
    def create_edge_stiffener(self,flange='top'):
        """ Creates a model for the edge stiffener of a flange """
        if flange == 'top':
            n1 = self.parts['top_flange']['nodes'][-2]
            n2 = self.parts['top_lip']['nodes'][-1]
            nodes = self.prof.nodes[n1:n2+1]
            
            #corner = list(self.prof.nodes[-2].coord)
            #y_edge = corner[0]-0.5*self.bt_center*self._rho_t
            #z_lip  = corner[1]-self.bt_lip_center*self._rho_t_lip
            #nodes = [[y_edge,corner[1]],corner,[corner[0],z_lip]]
            self.edge_stiff[flange] = OpenProf(nodes,self.t)
            
            self.edge_stiff[flange].draw()
        else:
            corner = list(self.prof.nodes[1].coord)            
            y_edge = corner[0]+0.5*self.bb_center*self._rho_b
            z_lip  = corner[1]+self.bb_lip_center*self._rho_b_lip
            nodes = [[y_edge,corner[1]],corner,[corner[0],z_lip]]
            self.edge_stiff[flange] = OpenProf(nodes,self.t)
    
    def create_free_flange(self,flange='bottom'):
        """ Creates a model for the free flange
            This includes the bottom flange with edge stiffener and
            1/5 of the height of the web
        """
        if flange == 'top':            
            corner = list(self.prof.nodes[-2].coord)
            y_edge = corner[0]-0.5*self.bt_center*self._rho_t
            z_lip  = corner[1]-self.bt_lip_center*self._rho_t_lip
            nodes = [[y_edge,corner[1]],corner,[corner[0],z_lip]]
            self.edge_stiff[flange] = OpenProf(nodes,self.t)
        else:
            nodes = []
            for i in range(0,3):
                nodes.append(list(self.prof.nodes[i].coord))
                        
            web_node = [0,nodes[-1][1]]
            web_node[1] += 0.2*self.hw 

            nodes.append(web_node)
            
            self.free_flange[flange] = OpenProf(nodes,self.t)
        
        
    
    @property
    def hw(self):
        """ Straight part of the web """
        return self.h - 2*(self.r+self.t_nom)
    
    @property
    def hw_center(self):
        """ Height of the profile from center line of the top flange
            to the center line of the bottom flange.
        """
        return self.h - self.t_nom
    
    @property
    def bt_center(self):
        """ width of the top flange either from the edge to the
            center line of the web or from center line of edge stiffener
            to the center line of the web
        """
        if self.ca == 0:
            bt = self.a - 0.5*self.t_nom
        else:
            bt = self.a - self.t_nom
        return bt
    
    @property
    def bb_center(self):
        """ width of the bottom flange either from the edge to the
            center line of the web or from center line of edge stiffener
            to the center line of the web
        """
        if self.cb == 0:
            bb = self.b - 0.5*self.t_nom
        else:
            bb = self.b - self.t_nom
        return bb
    
    @property
    def bt_lip_center(self):
        """ width of the top edge stiffener to the 
            to the center line of the top flange
        """
        if self.ca == 0:
            return 0
        else:
            return self.ca-0.5*self.t_nom
        
    @property
    def bb_lip_center(self):
        """ width of the bottom edge stiffener to the 
            to the center line of the bottom flange
        """
        if self.cb == 0:
            return 0
        else:
            return self.cb-0.5*self.t_nom
    
    @property
    def c_top(self):
        """ Straight part of the top flange """
        if self.ca > 0:
            c = self.a-2*(self.r+self.t_nom)
        else:
            c = self.a-(self.r+self.t_nom)
        return c
    
    @property
    def c_bottom(self):
        """ Straight part of the bottom flange """
        if self.cb > 0:
            c = self.b-2*(self.r+self.t_nom)
        else:
            c = self.b-(self.r+self.t_nom)
        return c
    
    @property
    def ca_straight(self):
        """ Straight part of the top lip """
        if self.ca > 0:
            c = self.ca-self.t_nom-self.r
        else:
            c = 0.0
        return c
    
    @property
    def cb_straight(self):
        """ Straight part of the bottom lip """
        if self.cb > 0:
            c = self.cb-self.t_nom-self.r
        else:
            c = 0.0
        
        return c
    
    @property
    def rm(self):
        """ Radius to the center line """
        return self.r + 0.5*self.t_nom
    
    def flange_lip(self,flange='top'):
        
        if flange == 'top':
            return self.ca
        else:
            return self.cb
    
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
        
        return self.t * self.center_line_length()
    
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
            self._lambda_t = lambda_p
        else:
            self._rho_b = rho
            self._lambda_b = lambda_p
        
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
            bpc = self.bp_b_lip
        
        # NOTE: Buckling factor according to EN 1993-1-3, not EN 1993-1-5!
        ksigma = en1993_1_3.buckling_factor_edge_stiffener(bpc/bp)
        lambda_p = en1993_1_5.lambda_p(bpc,self.t,self.material.eps(),ksigma)
        rho = en1993_1_5.reduction_factor_outstand(lambda_p)
        ceff = en1993_1_5.effective_width_outstand(bpc,psi=1.0,rho=rho)
            
        if flange == 'top':
            self._rho_t_lip = rho
            self._lambda_t_lip = lambda_p
        else:
            self._rho_b_lip = rho
            self._lambda_b_lip = lambda_p
            
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
            print("  Stress ratio: psi = {0:4.3f}".format(psi))
            print("  Buckling factor: ksigma = {0:4.3f}".format(ksigma))
            print("  Slenderness: lambda_p = {0:4.3f}".format(lambda_p))
            print("  Reduction factor: rho = {0:4.3f}".format(rho))
            print("  Effective width: beff = {0:4.3f} mm".format(beff))
    
        self._rho_w = rho
        self._psi_w = psi
        self._lambda_w = lambda_p
    
    def distortional_buckling(self,flange='top',verb=False):
        """ Distortional buckling of the edge stiffene 
            1. Generate edge stiffener effective part, As, Is
            2. Calculate stiffness K
            3. Calculate critical stress
            4. Calculate slenderness
            5. Calculate reduction factor for thickness reduction.
        
        """        
        K = self.edge_stiffness(flange)
        Is = self.edge_stiff[flange].second_moments_gc()[0]
        As = self.edge_stiff[flange].area()
        E = self.material.E
        
        scr = en1993_1_3.distortional_buckling_stress(Is,As,E,K)
        lambda_d = np.sqrt(self.fyb/scr)
        chi_d = en1993_1_3.distortional_buckling_reduction_factor(lambda_d)
        if verb:
            print("Distortional buckling:")
            print("  K = {0:4.3f} N/mm2".format(K))
            print("  s_cr = {0:4.3f} N/mm2".format(scr))
            print("  lambda_d = {0:4.3f}".format(lambda_d))
            print("  chi_d = {0:4.3f}".format(chi_d))
        
        if flange == 'top':
            self._chi_d_t = chi_d
            self._lambda_t = lambda_d
        else:
            self._chi_d_b = chi_d
            self._lambda_b = lambda_d
            
        return chi_d
    
    def area(self,delta=False):
        """ Gross cross-sectional area 
        
            :param : delta .. Flag for stating whether or not
                    the Eq. (5.1) or EN 1993-1-3 is used or not.
        """
        
        if delta:
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
        else:
            d = 0
        
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
    
    def effective_section(self,load='compression',verb=False):
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
                _chid_b .. reduction factor for bottom edge stiffener
                            for distortional buckling
        """
        
        if load == 'bending_y_pos':
            """ Positive bending:
                top flange in compression
                bottom flange in tension
            """
            self.effective_width_flange(psi=1.0,flange='top',verb=verb)
            
            if self.ca > 0:
                self.effective_width_edge_stiffener(flange='top',verb=verb)                
                self.create_edge_stiffener(flange='top')
                self.distortional_buckling(flange='top',verb=verb)
            
            
            """ Determine stress ratio for the web, using effecive flange """
            self.create_effective_model()
            
            self.prof_eff.draw(coord_axes=True)
            psi = self.stress_ratio()            
                        
            if verb:
                print("Psi = {0:4.3f}".format(psi))
            
            self.effective_width_web(psi,verb=verb)
            self.create_effective_model()
        elif load == 'compression':
            """ Compression: all parts are checked """
            for flange in ['top','bottom']:
                self.effective_width_flange(psi=1.0,flange=flange,verb=verb)
                if self.flange_lip(flange) > 0:
                    self.effective_width_edge_stiffener(flange,verb)
                    self.create_edge_stiffener(flange)
                
                    
                    
            for flange in ['top','bottom']:
                if self.flange_lip(flange) > 0:
                    self.distortional_buckling(flange=flange,verb=verb)
            
            
            self.effective_width_web(psi=1.0,verb=verb)            
            self.create_effective_model()
            
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
        
        print("Centroid: ygc = {0:4.3f}, zgc = {1:4.3f}".format(ygc,zgc))
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
                kf = self.edge_stiff[f2].area()/self.edge_stiff[f1].area()
        elif self.Med > 0.0 and flange == 'top':            
            kf = 0.0
        elif self.Med < 0.0 and flange == 'bottom':
            kf = 0.0
        
        b1 = abs(self.edge_stiff[f1].centroid()[0])
        #print('b1 = {0:4.3f}'.format(b1))
        if kf != 0.0:
            b2 = abs(self.edge_stiff[f2].centroid()[0])
        else:
            b2 = 0.0
        
        return C/(b1**2*hw + b1**3 + 0.5*b1*b2*hw*kf)
    
    def kh(self,load='gravity',a=0.0):
        """ Equivalent lateral load for purlins """
    
        xs = self.shear_center
        e = xs[0]
        h = self.h
        gs = 0.5*h-xs[1]
    
        kh0 = self.Iyz/self.Iy * gs/h        
        
        if load == 'gravity':
            kh = kh0+e/h
            if e < 0:
                """ shear center is left of qEd:
                    load is to the left, i.e. in the negative
                    direction
                """
                kh *= -1.0                
        elif load == 'uplift':
            aS = a-e
            kh = kh0+aS/h
            
            if aS/h >= kh0:
                kh *= -1.0                

        return kh
    
    def KB(self,load='gravity',afp=0.0):
        """ Stiffness of purlin-sheeting connection """
        
        if afp == 0.0:
            """ If no fastener position is given, use half of the width
                of the top flange
            """
            afp = 0.5*self.a
        
        if self.kh(load) > 0:
            bmod = afp
        else:
            bmod = 2*afp + self.b
            
        hd = self.hw_center
                
        return 0.25*self.material.E*self.t**3/(self.hw_center**2*(1-self.material.nu**2)*(hd+bmod))
        
    def K(self,load='gravity',afp=0.0,CD=0.0):
        """ Stiffness of the purlin-sheeting connection """
        KB = self.KB(load,afp)
        
        if CD != 0.0:
            KA_inv = self.hw_center**2/CD
        else:
            ba = self.a
            tnom = 1.0
            position = 'pos'
            A = 10
            bT = 50
            bR = 150
            
            CD = en1993_1_3.CD_A(ba,tnom,position,A,bT,bR,load)
            print("C_D = {0:4.3f}".format(CD))
            KA_inv = self.hw_center**2/CD
        
        print("KB = {0:4.3f}".format(KB))
        print("KA = {0:4.3f}".format(1/KA_inv))
        
        return 1/(1/KB+KA_inv)

    def axial_force_resistance(self,verb=False):
        if self.Ned >= 0:
            NRd = self.fya*self.A/gammaM0
        else:
            Aeff = self.Aeff
            Ag = self.A
            if Aeff < Ag:
                NRd = self.fyb*Aeff/gammaM0
            else:
                lambda_ratio = self._lambda_w/0.673
                
                if self.ca > 0.0:
                    top_lambda_ratio = self._lambda_t/0.65
                else:
                    top_lambda_ratio = self._lambda_t/0.673
                
                if self.cb > 0.0:
                    bot_lambda_ratio = self._lambda_b/0.65
                else:
                    bot_lambda_ratio = self._lambda_b/0.673
                    
                lambda_ratio = max(lambda_ratio,top_lambda_ratio,bot_lambda_ratio)
                NcRd = Ag*(self.fyb+4*(self.fya-self.fyb)*(1-lambda_ratio))
                
                NRd = min(NcRd,self.fya*Ag)/gammaM0

        if verb:
            print("NRd = {0:4.2f} kN".format(NRd*1e-3))

        return NRd

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
    
    def kh0(self):
        """ Coefficient kh0 from Fig. 10.3 """
        sc = self.shear_center
        gs = 0.5*self.h-sc[1]
        kh0 = self.Iyz/self.Iy*gs/self.h
        
        return kh0, sc[0]
        
    def kh(self,load='down'):
        """ Coefficient kh from Fig. 10.3 """
        
        kh0, y_sc = self.kh0()
        
        if load == 'down':
            kh = kh0 + y_sc/self.h
        else:
            f = 0.5*self.bt_center - y_sc
            kh = kh0 - f/self.h
            
        return kh
    
    def create_points(self,n=5):
        """ Creates points for the calculation model 
            :param n: number of divisions in the corner areas
            
        """
        # Number of points; used as an index
        npt = 0
        
        self.corner_pts = n
        h = self.h
        hweb = self.hw
        
        tnom = self.t_nom
        rm = self.Rm
        ro = self.Ro
        r = self.r
        
        bp_web = self.hw_center
        bp_bottom = self.bb_center
        bp_top = self.bt_center
        bp_bottom_lip = self.bb_lip_center
        bp_top_lip = self.bt_lip_center
        web_nodes = [[0.0,-0.5*bp_web],[0.0,0.5*bp_web]]
        
        if bp_bottom_lip > 0.0:
            bottom_lip_nodes = [[bp_bottom,web_nodes[0][1]+bp_bottom_lip],
                                [bp_bottom,web_nodes[0][1]+bp_bottom_lip-0.5*self.cb_straight]]
            for node in bottom_lip_nodes:
                self.add(node,part='bottom_lip')
            
            # Corner points
            # Lip-to-flange corner, bottom
            x_bottom_lip = (bp_bottom-rm,web_nodes[0][1]+rm)
        
            newPts = corner_points(x_bottom_lip,rm,phi=[0,-0.5*np.pi],n=n)
        
            for newPt in newPts[:-1]:            
                self.add(newPt)
            
            # The second and third nodes are for helping to cope
            # with the effective cross section.
            bottom_nodes = [[bp_bottom-rm,web_nodes[0][1]],
                            [bp_bottom-rm-1/3*self.c_bottom,web_nodes[0][1]],
                            [bp_bottom-rm-2/3*self.c_bottom,web_nodes[0][1]],
                            [bp_bottom-rm-self.c_bottom,web_nodes[0][1]]
                            ]            
        else:
            bottom_nodes = [[bp_bottom,web_nodes[0][1]],
                            [bp_bottom-0.5*self.c_bottom,web_nodes[0][1]],
                            [bp_bottom-self.c_bottom,web_nodes[0][1]]
                            ]
        
        for node in bottom_nodes:
            self.add(node,part='bottom_flange')
        
        
        # Corner points
        # flange-to-web corner, bottom
        x_bottom_web = (rm,web_nodes[0][1]+rm)
        
        newPts = corner_points(x_bottom_web,rm,phi=[1.5*np.pi,np.pi],n=n)
        
        for newPt in newPts[1:-1]:            
            self.add(newPt)
        
        # Web points
        self.add([0,-0.5*bp_web+rm],part='web')
        self.add([0,-0.5*bp_web+rm+0.25*hweb],part='web')
        self.add([0,-0.5*bp_web+rm+0.75*hweb],part='web')
        self.add([0,0.5*bp_web-rm],part='web')
        
        # Corner points
        # web-to-flange corner, top
        x_top_web = (+rm,web_nodes[1][1]-rm)
        
        newPts = corner_points(x_top_web,rm,phi=[np.pi,0.5*np.pi],n=n)
        
        for newPt in newPts[1:-1]:            
            self.add(newPt)
        
        # Top flange
        top_nodes = [[rm,web_nodes[1][1]]]
        
        if bp_top_lip > 0.0:
            top_nodes += [[0.25*bp_top,web_nodes[1][1]],
                          [0.75*bp_top,web_nodes[1][1]],
                          [bp_top-rm,web_nodes[1][1]]]
            
            # Corner
            # web-to-flange corner, top
            x_top_lip = (bp_top-rm,web_nodes[1][1]-rm)
        
            for node in top_nodes:
                self.add(node,part='top_flange')
        
            newPts = corner_points(x_top_lip,rm,phi=[0.5*np.pi,0],n=n)
        
            for newPt in newPts[1:-1]:            
                self.add(newPt)
                
            # Top lip nodes
            
            self.add([bp_top,web_nodes[1][1]-rm],'top_lip')
            self.add([bp_top,web_nodes[1][1]-0.5*bp_top_lip],'top_lip')
            self.add([bp_top,web_nodes[1][1]-bp_top_lip],'top_lip')
        else:
            top_nodes += [[0.5*bp_top,web_nodes[1][1]],[bp_top,web_nodes[1][1]]]
        
            for node in top_nodes:
                self.add(node,part='top_flange')
                
    
    #def create_model(self):
        """ Creates calculation model according to 
            EN 1993-1-3:2006, Annex C 
        
            Origin is located at the middle of the web
        """
        """
        bp_web = self.bp_w
        bp_bottom = self.bp_b
        bp_top = self.bp_t
        bp_bottom_lip = self.bp_b_lip
        bp_top_lip = self.bp_t_lip
        """
        
        """
        bp_web = self.hw_center
        bp_bottom = self.bb_center
        bp_top = self.bt_center
        bp_bottom_lip = self.bb_lip_center
        bp_top_lip = self.bt_lip_center
        web_nodes = [[0.0,-0.5*bp_web],[0.0,0.5*bp_web]]
        if bp_bottom_lip > 0.0:
            bottom_nodes = [[bp_bottom,web_nodes[0][1]+bp_bottom_lip],
                            [bp_bottom,web_nodes[0][1]]
                            ]
        else:
            # No bottom edge stiffener
            bottom_nodes = [[bp_bottom,web_nodes[0][1]]]
        
        if bp_top_lip > 0.0:
            top_nodes = [[bp_top,web_nodes[1][1]],
                             [bp_top,web_nodes[1][1]-bp_top_lip]]
        else:
            top_nodes = [[bp_top,web_nodes[1][1]]]
            
        nodes = bottom_nodes + web_nodes + top_nodes
        
        self.prof = OpenProf(nodes,self.t)
        """
    
    def create_effective_model(self):
        """ Create calculation model for the effective cross-section 
            This is done based on the model of the gross cross-section
        """
    
        if self.prof is None:
            self.create_model()
                    
        prof_eff = deepcopy(self.prof)
                
        if self._rho_w < 1.0:
            """ Effective web """            
            if self._psi_w < 0.0:
                hc = self.bp_w/(1-self._psi_w)
                beff = self._rho_w * hc
                be1 = 0.4*beff
                be2 = 0.6*beff                
            elif self._psi_w < 1.0:
                hc = self.bp_w
                beff = self._rho_w * hc
                be1 = 2*beff/(5-self._psi_w)
                be2 = beff-be1            
            else:
                hc = self.bp_w
                beff = self._rho_w * hc
                be1 = 0.5*beff
                be2 = 0.5*beff
        
            bneg = hc-beff
            # Part of the web in tension
            #ht = self.hw_center-hc
        
            if self.cb > 0:
                nweb = 2
            else:
                nweb = 1
            
            """
            s = (be2+ht)/self.hw_center
            prof_eff.split_segment(nweb,s)
            s2 = bneg/(be1+bneg)
            prof_eff.split_segment(nweb+1,s2)
            prof_eff.segments[nweb+1].t = 0.0
            """
            nw_top = self.parts['web']['nodes'][2]
            nw_bot = self.parts['web']['nodes'][1]
            
            print(be1,be2, bneg)
                    
            prof_eff.nodes[nw_top].z = prof_eff.nodes[nw_top+1].z + self.rm*np.sin(np.radians(45)) - be1
            prof_eff.nodes[nw_bot].z = prof_eff.nodes[nw_top].z - bneg
            
            prof_eff.segments[nw_bot].t = 0.0
        
        if self.ca > 0.0:            
            """ Edge stiffener is present """
            # Number of segments in the section
            
            ntop = self.parts['top_flange']['nodes'][-2]
            
            # Reduce the thickness of edge stiffener, if needed
            for segment in prof_eff.segments[ntop:]:
                segment.t *= self._chi_d_t
                
            # If the top flange is not entirely effective,
            # make the thickness of the non-effective segment zero.
            if self._rho_t < 1.0:
                prof_eff.segments[ntop-1].t = 0
            
                              
            
            """
            nseg = prof_eff.n-1
            # Move the node of the edge stiffener up
            prof_eff.nodes[-1].z += (1-self._rho_t_lip)*self.bt_lip_center
                    
            # Split the top flange into two parts
            s = 0.5*self._rho_t            
            prof_eff.split_segment(nseg-2,s)
                        
            #prof_eff.draw(node_labels=True)
            if self._rho_t < 1.0:
                s2 = (1-self._rho_t)*self.bt_center/prof_eff.segments[-2].length()
                prof_eff.split_segment(nseg-1,s2)
                prof_eff.segments[nseg-1].khp = 0.0
                
            prof_eff.segments[-2].t *= self._chi_d_t
            prof_eff.segments[-1].t *= self._chi_d_t
            """          
        else:
            """ The effective top flange is an outstand element """
            #prof_eff.nodes[-1].y -= (1-self._rho_t)*self.bp_t
            if self._rho_t < 1.0:
                prof_eff.nodes[-2].y = self._rho_t*self.bp_t + self.gr
                prof_eff.segments[-1].t = 0.0

        
        if self.cb > 0.0:
            """ Edge stiffener is present in the bottom flange """
            
            nbot = self.parts['bottom_flange']['nodes'][1]
            
            # Reduce the thickness of edge stiffener, if needed
            for segment in prof_eff.segments[:nbot]:
                segment.t *= self._chi_d_b
            
            # If the bottom flange is not entirely effective,
            # make the thickness of the non-effective segment zero.
            if self._rho_b < 1.0:
                prof_eff.segments[nbot].t = 0
            
            """
            # Move the node of the edge stiffener down
            prof_eff.nodes[0].z -= (1-self._rho_b_lip)*self.bb_lip_center
            
            # Split the bottom flange into two parts
            s = 0.5*self._rho_b
            prof_eff.split_segment(1,s)
            """            
            """ If there is a non-effective part of the bottom flange
                remove that too.
            """
            """
            if self._rho_b < 1.0:
                s2 = (1-self._rho_b)*self.bb_center/prof_eff.segments[2].length()
                prof_eff.split_segment(2,s2)
                prof_eff.segments[2].khp = 0.0
            prof_eff.segments[0].t *= self._chi_d_b
            prof_eff.segments[1].t *= self._chi_d_b        
            """
        else:
            
            """ The effective bottom flange is an outstand element """
            #prof_eff.nodes[-1].y -= (1-self._rho_t)*self.bp_t
            if self._rho_b < 1.0:
                prof_eff.nodes[1].y = self._rho_b*self.bp_b + self.gr
                
            

        self.prof_eff = prof_eff
            
    def create_edge_stiffener(self,flange='top'):
        """ Creates a model for the edge stiffener of a flange """
        if flange == 'top':            
            n1 = self.parts['top_flange']['nodes'][-2]
            n2 = self.parts['top_lip']['nodes'][-1]
            nodes = self.prof.nodes[n1:n2+1]
            
            nodes[0].y = self.bp_t+self.gr-0.5*self.bp_t*self._rho_t
            
            if self._rho_t_lip < 1.0:            
                nodes[-2].z = nodes[-1].z + (1-self._rho_t_lip)*self.bp_t_lip
            
            #corner = list(self.prof.nodes[-2].coord)
            #y_edge = corner[0]-0.5*self.bt_center*self._rho_t
            #z_lip  = corner[1]-self.bt_lip_center*self._rho_t_lip
            #nodes = [[y_edge,corner[1]],corner,[corner[0],z_lip]]
            self.edge_stiff[flange] = OpenProf(nodes,self.t)
            #self.edge_stiff[flange].draw()
            """
            self.edge_stiff[flange].draw()
            corner = list(self.prof.nodes[-2].coord)
            y_edge = corner[0]-0.5*self.bt_center*self._rho_t
            z_lip  = corner[1]-self.bt_lip_center*self._rho_t_lip
            nodes = [[y_edge,corner[1]],corner,[corner[0],z_lip]]
            self.edge_stiff[flange] = OpenProf(nodes,self.t)
            """
        else:
            n1 = self.parts['bottom_lip']['nodes'][0]
            n2 = self.parts['bottom_flange']['nodes'][1]            
            nodes = self.prof.nodes[n1:n2+1]
            
            nodes[-1].y = self.bp_b+self.gr-0.5*self.bp_b*self._rho_b
            
            if self._rho_b_lip < 1.0:
                nodes[1].z = nodes[0].z + (1-self._rho_b_lip)*self.bp_b_lip
            
            self.edge_stiff[flange] = OpenProf(nodes,self.t)
            """
            corner = list(self.prof.nodes[1].coord)            
            y_edge = corner[0]-0.5*self.bb_center*self._rho_b
            z_lip  = corner[1]+self.bb_lip_center*self._rho_b_lip
            nodes = [[y_edge,corner[1]],corner,[corner[0],z_lip]]
            self.edge_stiff[flange] = OpenProf(nodes,self.t)
            """
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
    
def Cexample(t=2.0,ca=0,cb=0):
    
    c = CSection(t_nom=t,h=150.0,a=47,b=41,ca=ca,cb=cb,r=3.0,t_coat=0.04)
    
    #c.Med = 10.0
    c.Ned = -10.0
    #c.effective_width_flange(psi=1.0, flange='top',verb=True)
    c.effective_section('compression',True)
    #c.create_model()
    #c.prof.draw(node_labels=True)
    #c.create_effective_model()
    
    #c.prof_eff.draw(seg_labels=True,node_labels=True)

    #ygc, zgc = c.centroid_eff()
    #print(ygc,zgc)

    return c

def Cdesign(t=2.0,ca=0,cb=0):
    
    c = CSection(t_nom=t,h=150.0,a=47,b=47,ca=ca,cb=cb,r=3.0,t_coat=0.04)
    
    #c.Med = 10.0
    c.Ned = -10.0
    #c.effective_width_flange(psi=1.0, flange='top',verb=True)
    #c.effective_section('compression',True)

    print("Axial resistance: {0:4.2f} kN".format(c.NRd*1e-3))
    
    c.prof.draw(coord_axes=True)
    plt.grid(True)
    c.prof_eff.draw(coord_axes=True)
    plt.grid(True)
    
    return c

if __name__ == "__main__":
    
    #c = Cexample(t=1.0,ca=16,cb=16)
    #c = Cexample(t=1.5,ca=18,cb=18)
    c = Cdesign(t=1.5,ca=18,cb=18)
    
    #z = ZSection(t_nom=1.5,h=200,a=74,b=66,ca=21.2,cb=21.2,r=3,
    #             material="S350GD",t_coat=0.04)
    
    #z = ruukki_z(300,2.5)
    
    #z.prof.draw(node_labels=True)
    
    #print(z.Agross()*1e-2)
    
    #c = CSection(1.0,150,47.0,41.0,16,16,3)
    #c.prof.draw(node_labels=False)
    
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
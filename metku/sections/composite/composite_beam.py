# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 19:53:07 2019

Composite beam sections

@author: kmela
"""

import math
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

try:
    from metku.eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
    
except:
    from eurocodes.en1993 import en1993_1_1, en1993_1_5, constants
    from eurocodes.en1994 import en1994_1_1
    

class ConcreteSlab:
    """ Simple concrete slab with rectangular cross-section
        
    
    """
    
    def __init__(self,thickness,material,reinforcement=None):
        """ Constructor
            input:
                thickness [mm] thickness of slab
                material .. Concrete class material
        """
        
        self.hc = thickness
        self.material = material
        
    @property
    def h(self):
        """ Total height of the slab """
        return self.hc
    
    @property
    def h0(self):
        """ Notional size of slab thickness, for drying (EN 1992-1-1, Annex B) """
        return 2*self.hc
        
    def weight(self):
        """ Weight of the slab """
        return (self.material.density+1)*self.hc*1e-3
        
    def draw(self,axes=None,y0=0,width=100):
        """ Draw slab """
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        
        
        slab = patches.Rectangle((-0.5*width, y0), width=width,
                                       height=self.hc, hatch='\\')
        
        ax.add_patch(slab)
        
        
    
class CompositeSlab:
    """ Simple composite slab with rectangular cross-section
        
    
    """
    
    def __init__(self,h_total,concrete,h_steel):
        """ Constructor
            input:
                thickness [mm] thickness of slab
                material .. Concrete class material
        """
        
        self.hc = h_total-h_steel
        self.material = concrete
        self.ha = h_steel
        
    @property
    def h(self):
        """ Total height of the slab """
        return self.hc + self.ha
    
    @property
    def h0(self):
        """ Notional size of slab thickness, for drying (EN 1992-1-1, Annex B) """
        return 2*self.hc
    
    def concrete_weight(self):
        """ Weight of the concrete """
        return (self.material.density+1)*self.h*1e-3
    
    def draw(self,axes=None,y0=0,width=100):
        """ Draw slab """
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        
        
        slab_steel = patches.Rectangle((-0.5*width, y0), width=width,
                                       height=self.ha, fill=False)
        
        slab_concrete = patches.Rectangle((-0.5*width, y0+self.ha), width=width,
                                       height=self.hc, hatch='\\')
        ax.add_patch(slab_steel)
        ax.add_patch(slab_concrete)
        

class CompositeIBeam:
    """ Class for composite beams with an I section as
        steel part and a slab (composite or concrete) on top
    """
    
    def __init__(self,length,steel_part,slab,studs=None,RH=50):
        """ Constructor
            Input:
                length .. length of the beam
                steel_part: ISection class profile or WI profile steel section
                slab: concrete or composite slab
                studs: shear studs
                RH: relative humidity of concrete
                
            parameters
                length .. length of the beam
                steel .. steel profile
                slab .. slab
                studs .. shear connector, of class ShearStud
                st .. stud spacing
                sf .. transverse reinforcement spacing
                
        """
        
        self.length = length
        self.steel = steel_part
        self.slab = slab
        self.slab.material.RH = RH
        self.studs = studs
        self.st = 150
        self.sf = 150
        self.b0 = 0                 
        self.t0 = 28 # Time for creep (days)         
        
    def creep_coefficient(self,t0=None):
        
        if t0 is None:
            t0 = self.t0
            
        return self.slab.material.creep_coefficient(self.slab.h0,1e8,t0)
    
    @property
    def n0(self):
        """ Ratio of Young's modulae of steel and concrete (short term) """
        return self.Ea/self.Ecm
    
    @property
    def nC(self):
        """ Ratio of Young's modulae of steel and concrete (creep)"""
        return self.Ea/self.Ecc
    
    @property
    def nS(self):
        """ Ratio of Young's modulae of steel and concrete (shrinkage)"""
        return self.Ea/self.Ecs
        
    @property 
    def beff(self):
        """ Effective width of concrete flange """
        return self.b0 + 2.0*self.length/8.0
        
    @property
    def yc(self):
        """ Distance of centroid of concrete slab from the top
            of the composite section
        """
        return 0.5*self.slab.hc
    
    @property
    def ya(self):
        """ Distance of centroid of steel part from the top
            of the composite section
        """
        return self.slab.h + 0.5*self.steel.h
    
    @property
    def Aa(self):
        """ Area of steel part """
        return self.steel.A
    
    @property
    def Ia(self):
        """ Second moment of area of steel part with respect to
            the major axis at centroid of steel part
        """
        return self.steel.I[0]
    
    @property
    def Ea(self):
        """ Young's modulus of steel part """
        return self.steel.E
    
    @property
    def EAa(self):
        """ Axial stiffness of steel part """
        return self.Ea*self.Aa
    
    @property
    def EIa(self):
        """ Bending stiffness of steel part 
            with respect to the neutral axis of steel part
        """
        return self.Ea*self.Ia
    
    @property
    def Ecm(self):
        """ Young's modulus of concrete """
        return self.slab.material.Ecm
    
    @property
    def Ecc(self):
        """ Young's modulus of concrete with creep effects """        
        phi = self.creep_coefficient()        
        return self.Ecm/(1+1.1*phi)
    
    @property
    def Ecs(self):
        """ Young's modulus of concrete with shrinkage effects """        
        phi = self.creep_coefficient(t0=1)        
        return self.Ecm/(1+0.55*phi)
    
    @property
    def EIc(self):
        """ Bending stiffness of concrete part 
            with respect to the neutral axis of concrete part
        """
        return self.Ecm*self.Ic
    
    @property
    def Ac(self):
        """ Area of concrete part """
        return self.slab.hc*self.beff
    
    @property
    def Ic(self):
        """ Second moment of area of concrete part with respect to
            the major axis at centroid of concrete part
        """
        return 1/12*self.beff*self.slab.hc**3
    
    @property
    def EAc(self):
        """ Axial stiffness of concrete """
        return self.Ecm*self.Ac
        
    @property
    def EAcc(self):
        """ Axial stiffness of concrete in creep"""
        return self.Ecc*self.Ac
    
    @property
    def EIcc(self):
        """ Bending stiffness of concrete part with creep effects
            with respect to the neutral axis of concrete part
        """
        return self.Ecc*self.Ic
    
    def centroid(self,long_term=None):
        """ Centroid of the composite section 
            Assume full interaction
        """        
        
        if long_term is None:
            yN = (self.EAa*self.ya + self.EAc*self.yc)/self.EAcom()
        elif long_term == "creep":
            yN = (self.EAa*self.ya + self.EAcc*self.yc)/self.EAcom("creep")
                        
        return yN
    
    def EAcom(self,long_term=None):
        """ axial stiffness composite section """
        if long_term is None:
            EAcom = self.EAa + self.EAc
        elif long_term == "creep":            
            EAcom = self.EAa + self.EAcc
        
        return EAcom


    def ya_com(self,long_term=None):
        """ distance of centroid of steel part from the neutral
            axis of the composite section
        """
        
        return abs(self.ya-self.centroid(long_term))
    
    def yc_com(self,long_term=None):
        """ distance of centroid of concrete part from the neutral
            axis of the composite section
        """
        
        return abs(self.yc-self.centroid(long_term))

    def EIcom(self,long_term=None):
        """ sending stiffness of composite section """        
        if long_term is None:
            EIcom = self.Ea*(self.Ia+self.ya_com()**2*self.Aa) + self.Ecm*(self.Ic+self.yc_com()**2*self.Ac)
        elif long_term == "creep":
            EIcom = self.Ea*(self.Ia+self.ya_com(long_term)**2*self.Aa) + self.Ecc*(self.Ic+self.yc_com(long_term)**2*self.Ac)
        
        return EIcom

    def distance_centroids(self):
        """ distance between centroids of steel and concrete parts """
        return self.ya - self.yc
    
    def composite_coefficient(self):
        """ Coefficint alpha that indicates the degree of composite action 
            TRY/BY, pp. 24        
        """
        a = self.distance_centroids()
        return a**2*self.EAa*self.EAc/(self.EAa+self.EAc)/(self.EIa + self.EIc)
    
    def Ra(self):
        """ Full resultant of steel part """
        
        return self.steel.fy*self.Aa
    
    def Rc(self):
        """ Full resultant of concrete part """
        
        return 0.85*self.slab.material.fcd()*self.Ac
    
    def Rf(self):
        """ Resultant of top flange of the steel beam """
        
        Af = self.steel.tf * self.steel.b        
        
        return self.steel.fy*Af
    
    def Rw(self):
        """ Resultant of web of the steel beam """
        
        Aw = self.steel.hw * self.steel.tw        
        
        return self.steel.fy*Aw
    
    def plastic_neutral_axis(self):
        """ Location of the plastic neutral axis
            from the top of the concrete slab
        """
        Ra = self.Ra()
        Rc = self.Rc()
        
        hc = self.slab.hc
        
        if Ra < Rc:
            ec0 = Ra/Rc*hc
        elif Rc < Ra-2*self.Rf():
            hl = self.slab.h
            ec0 = 0.5*(Ra-Rc)/self.steel.b/self.steel.fy + self.slab.h
        else:
            Rf = self.Rf()
            Rw = self.Rw()
            hw = self.steel.hw
            hl = self.slab.h
            tf = self.steel.tf
            
            ec0 = 0.5*(Ra-Rc-2*Rf)/Rw*hw+hl+tf
            
        return ec0
        
    
    def MplRd(self,verb=False):
        """ Plastic moment resistance based on full shear connection """
        
        Ra = self.Ra()
        Rc = self.Rc()
        
        ei = self.distance_centroids()
        hc = self.slab.hc
        if Ra < Rc:
            """ Neutral axis in concrete part """            
            ec0 = Ra/Rc*hc
            MplRd = Ra*(ei+0.5*hc-0.5*ec0)
            if verb:
                print("Plastic neutral axis in concrete:")
                print("Ra = {0:4.2f} kN".format(Ra*1e-3))
                print("Rc = {0:4.2f} kN".format(Rc*1e-3))
                print("ec0 = {0:4.2f} [mm]".format(ec0))
        elif Rc < Ra-2*self.Rf():
            """ Neutral axis in top flange of steel beam """
            hl = self.slab.h
            ec0 = 0.5*(Ra-Rc)/self.steel.b/self.steel.fy + self.slab.h
            MplRd = Ra*(ei-hl+0.5*hc) + Rc*(hl-0.5*hc) - 0.25*(Ra-Rc)**2/self.Rf()*self.steel.tf
        else:
            """ Neutral axis in the web of the steel part """
            Rf = self.Rf()
            Rw = self.Rw()
            hw = self.steel.hw
            hl = self.slab.h
            tf = self.steel.tf
            
            ec0 = 0.5*(Ra-Rc-2*Rf)/Rw*hw+hl+tf
            
            MplRd = Ra*(ei-hl-tf+0.5*hc) + Rc*(hl+tf-0.5*hc) + Rf*tf - 0.25*(Ra-Rc-2*Rf)**2/Rw*hw
                    
        if verb:
            print("MplRd = {0:4.2f} [kNm]".format(MplRd*1e-6))
        
        return MplRd
    
    def VRd(self):
        """ Shear resistance
            This equals shear resistance of the steel part
        """
        
        return self.steel.shear_force_resistance()
    
    def draw(self):
        """ Draw the profile             
        """
        
        fig, ax = plt.subplots(1)

        beff = self.beff

        self.steel.draw(axes=ax)        
                    
        
        
        self.slab.draw(axes=ax,y0=0.5*self.steel.h,width=beff)    
            
        yel = 0.5*self.steel.h + self.slab.h-self.centroid()
        ypl = 0.5*self.steel.h + self.slab.h-self.plastic_neutral_axis()      
        
                        
        plt.plot(0.0, yel, 'or')
        plt.plot(0.0, ypl, 'db')
    
        """
        tflange =  mlines.Line2D(top_x_data,top_y_data,linewidth=lw)
        ax.add_line(tflange)
        
        
        # Top chamfers
        tr_corner = patches.Arc((0.5*tw+r,0.5*hw),2*r,2*r,0,90,180)
        tl_corner = patches.Arc((-0.5*tw-r,0.5*hw),2*r,2*r,0,0,90)
        
        ax.add_patch(tr_corner)
        ax.add_patch(tl_corner)
        
        # Draw web
        xweb_l = [-0.5*tw,-0.5*tw]
        yweb = [-0.5*hw,0.5*hw]
        
        if theta is not 0:
            web_rot = [R.dot(np.array([x,y])) for x, y in zip(xweb_l,yweb)]
            xweb_l = [x[0] for x in web_rot]
            yweb_l = [x[1] for x in web_rot]
        else:
            yweb_l = yweb
        
        lweb = mlines.Line2D(xweb_l,yweb_l,linewidth=lw)
        ax.add_line(lweb)
        
        xweb_r = [0.5*tw,0.5*tw]
        
        rweb = mlines.Line2D(xweb_r,yweb,linewidth=lw)
        ax.add_line(rweb)
        """
        
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
         
        ax.set_xlim(-0.5*beff, 0.5*beff)
        ax.set_ylim(-0.5*self.steel.h, 0.5*self.steel.h+self.slab.h)
        
        ax.set_aspect('equal')


if __name__ == '__main__':
    
        
    from sections.steel.ISection import IPE
    from eurocodes.en1992 import en1992_1_1, constants


    """    
    steel = IPE(300)
    #slab = ConcreteSlab(140,constants.Concrete("C25/30"))
    slab = CompositeSlab(100,constants.Concrete("C25/30"),50)
    
    p = CompositeIBeam(12000,IPE(450),slab)
    p.draw()
    
    print(p.slab.material.creep_coefficient(50,227,t=1e8,t0=28))
    
    print("Effective width of concrete flange {0:4.2f} mm".format(p.beff))    
    print("Centroid of concrete {0:4.2f} mm".format(p.yc))
    print("Centroid of steel part {0:4.2f} mm".format(p.ya))
    print("Centroid {0:4.2f} mm".format(p.centroid()))
    print("Bending stiffness of steel part: {0:4.2f} MNm".format(p.EIa*1e-12))
    print("Bending stiffness of concrete part: ",p.EIc*1e-12)
    print("Distance between centroids [mm]: ", p.distance_centroids())
    print("Liittovaikutuskerroin: ", p.composite_coefficient())
    print("Bending stiffness of composite section: ",p.EIcom()*1e-12)
    print("Bending resistance: MplRd = {0:4.2f}".format(p.MplRd()*1e-6))
    print("Bending resistanc of steel: Mpl,a,Rd = {0:4.2f}".format(p.steel.plastic_bending_resistance()*1e-6))
    """
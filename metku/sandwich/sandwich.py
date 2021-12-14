# -*- coding: utf-8 -*-
"""
Calculation methods for sandwich panels with steel faces

Created on Wed Aug  7 09:15:25 2019

@author: kmela
"""

import math
from materials.steel_data import Steel
from eurocodes.en1993.constants import gammaM2

fasteners = {'SXC14-S19-5.5':{'d':5.5,'dS':5.1,'d1':4.59}}

def ruukki_panel(name='SPA E', thickness=100, length=6000):
    
    if name == 'SPA E':
        TopFace = Face(tF=0.6,material="S280GD")
        BottomFace = Face(tF=0.5,material="S280GD")
        core = Core(tC=thickness-0.6-0.5,shear_modulus=3.7,density=110e-9)
        width = 1200
    
    
    panel = SandwichPanel([TopFace,BottomFace],core,width,length)
    
    return panel

class SandwichScrew:
    """ Class for sandwich panel fasteners """
    
    def __init__(self,d=5.5,dS=5.1,d1=4.59):
        """
        Constructor

        Parameters
        ----------
        d : float, optional
            Diameter [mm]. The default is 5.5.
        dS : TYPE, optional
            Diameter [mm]. The default is 5.1.
        d1 : TYPE, optional
            Diameter [mm]. The default is 4.59.

        Returns
        -------
        None.

        """
        
        self.d = d
        self.dS = dS
        self.d1 = d1
    

class Face:
    """ Face of sandwich panel """
    
    def __init__(self,tF,material="S280GD",tcoat=0.04):
        """ Constructor

            Parameters:
            -----------
                :param tF: thickness of face [mm]
                :param material: face material [mm]
                :param ctcoat: thickness of coating [mm]

                Variables:
                ----------
                :ivar t: thickness of face
                :ivar mat: material
                :ivar tcoat: coating thickness
                :ivar fy: yield strength [MPa]
                :ivar fu: ultimate strength [MPa]
        """
        
        self.t = tF
        self.mat = Steel(material)
        self.tcoat = tcoat        
        
    @property
    def fy(self):
        return self.mat.fy
    
    @property
    def fu(self):
        return self.mat.fu
    
    @property
    def E(self):
        return self.mat.E
    
    @property
    def tcor(self):
        """ Thickness of face core """
        return self.t-self.tcoat
    
    def A(self,b):
        """ Cross-sectional area of the face 
            :param b: width of panel        
        """
        
        return self.t*b
    
    def I(self,b):
        """ Second moment of area of the face with respect to its own
            neutral axis
            :param b: width of panel
        """
        return 1/12*b*self.t**3
       
    def EI(self,b):
        """ Bending stiffness of the face with respect ot its own neutral axis
        """
        return self.E*self.I(b)
    
    def EA(self,b):
        """ Axial stiffness of the face """
        return self.E*self.A(b)
            
class Core:
    """ Core material """

    def __init__(self,tC,shear_modulus,density):
        
        self.t = tC
        self.G = shear_modulus
        self.density = density
        
        

class SandwichPanel:
    """ Class for sandwich panels """
    
    def __init__(self,faces,core,width=1200,length=6000,c=[1000],
                 fastener=SandwichScrew(d=5.5,dS=5.1,d1=4.59)):
        """ Constructor

            Parameters:
            -----------
                :param faces: list of class Face objects
                :param core: class core object
                :param width, length: self-evident [mm]
                :param c: list of screw pair distances [mm], starting from
                        outermost screws

                Variables:
                ----------
                :ivar t: thickness of face
                :ivar mat: material
                :ivar tcoat: coating thickness
                :ivar fy: yield strength [MPa]
                :ivar fu: ultimate strength [MPa]
        """
        self.faces = faces
        self.core = core
        self.width = width
        self.length = length
        self.c = c
        self.fastener = fastener
        
    @property
    def thickness(self):
        """ Total thickness of panel """
        return self.faces[0].t + self.faces[1].t + self.core.t
        
    @property
    def d(self):
        """ distance between center lines of faces """
        return 0.5 * self.faces[0].t + 0.5 * self.faces[1].t + self.core.t

    @property
    def Atop(self):        
        return self.faces[0].A(self.width)

    @property
    def Abottom(self):        
        return self.faces[1].A(self.width)

    @property
    def S(self):
        """ Shear stiffness [N] """
        return self.core.G*self.width*self.d**2/self.core.t
    
    @property
    def EIo(self):
        """ Bending stiffness of faces """
        b = self.width
        return self.faces[0].EI(b) + self.faces[1].EI(b)
    
    @property
    def EIs(self):
        """ 'sandwich part' of the bending stiffness """
        EA = [face.EA(self.width) for face in self.faces]
        
        return EA[0]*EA[1]/sum(EA)*self.d**2
    
    @property
    def tF2_cor(self):
        """ Thickness of steel core of face 2 """
        
        return self.faces[1].tcor
    
    @property
    def fu_F2(self):
        """ Ultimate strength of face 2 """
        return self.faces[1].fu

    def fastener_transverse_stiffness(self,tsup,verb=False):
        """ Transverse stiffness of the fastener according to
            ECCS recommendations
            
            :param tsup: thickness of supporting structure [mm]
        """
        d1 = self.fastener.d1
        dS = self.fastener.dS
        
        # bending stiffness of fastener
        EI = 200000*math.pi*dS**4/64
        
        # stiffness of clamping in the supporting structure
        Csup = 2400*math.sqrt(tsup*d1**5)
        
        # Stiffness of internal face (hole elongation)
        if self.faces[1].tcor < 0.7:
            kF2 = 6.93*self.faces[1].fu*math.sqrt(self.faces[1].tcor**3*d1)/(0.26 + 0.8 * self.faces[1].t)
        else:
            kF2 = 4.2*self.faces[1].fu*math.sqrt(self.faces[1].tcor**3*d1)/0.373

        
        D = self.thickness
        xF = 1 - (1/kF2-0.5*D*tsup/Csup - D*tsup**2/8/EI)/(1/kF2+D**2/Csup+D**2*(2*D+3*tsup)/6/EI)
        
        kv = 1/(xF/kF2 + (tsup**2 + 2*(1-xF)*D*tsup)/4/Csup + (3*(1-xF)*D*tsup**2 + tsup**3)/24/EI)
        
        if verb:
            print("EI = {0:4.2f}".format(EI))
            print("Csup = {0:4.2f}".format(Csup))
            print("kF2 = {0:4.2f}".format(kF2))
            print("xF = {0:4.4f}".format(xF))
            
            print("1/k1 = {0:4.4g}".format(xF/kF2))
            print("1/k2 = {0:4.4g}".format((tsup**2 + 2*(1-xF*D*tsup))/4/Csup))
            print("1/k3 = {0:4.4g}".format((3*(1-xF)*D*tsup**2 + tsup**3)/24/EI))
            
        return kv
    
    def shear_stiffness_for_support(self,tsup,L,n):
        """
        Shear stiffness of the panel that it provides for stabilisation of an individual member

        Returns
        -------
        None.

        """
        kv = self.fastener_transverse_stiffness(tsup)
        
        Si = 0.5*kv*n/L*sum([ci**2 for ci in self.c])
        
        return Si
        
    def fastener_shear_resistance(self):
        """
        Shear resistance of fastener

        Returns
        -------
        None.

        """
        VRk = 4.2*math.sqrt(self.tF2_cor**3*self.fastener.d1)*self.fu_F2
        
        return VRk/gammaM2
        
    
def example():
    
    TopFace = Face(tF=0.6,material="S280GD")
    BottomFace = Face(tF=0.5,material="S280GD")
    core = Core(tC=200-0.6-0.5,shear_modulus=3.7,density=110e-9)
    
    panel = SandwichPanel([TopFace,BottomFace],core,
                           width=1200,length=6000)

    # Tuulikuorma
    qw = 0.6 # kN/m2
    
    qU = 1.5*qw*panel.width*1e-3
    MEd = qU*panel.length**2/8
    VEd = qU*panel.length/2
    d = panel.d
    
    sF1 = MEd/d/panel.Atop
    sF2 = MEd/d/panel.Abottom
    tauC = VEd/d/panel.width
    Ls = 100
    sC = VEd/Ls/panel.width
    
    print("Kuorma murtorajatilassa qU = {0:4.2f} kN/m".format(qU))
    print("Taivutusmomentti MEd = {0:4.2f} kNm".format(MEd*1e-6))
    print("Taivutusmomentti VEd = {0:4.2f} kN".format(VEd*1e-3))
    print("d = {0:4.2f} mm".format(d))
    print("Atop = {0:4.2f} mm2".format(panel.Atop))
    print("Abottom = {0:4.2f} mm2".format(panel.Abottom))
    print("sF1 = {0:4.2f} MPa".format(sF1))
    print("sF2 = {0:4.2f} MPa".format(sF2))
    print("tauC = {0:4.3f} MPa".format(tauC))
    print("SC = {0:4.3f} MPa".format(sC))

    return panel

if __name__ == '__main__':
    
    #f = Face(0.5,material="S250GD")
            
    p = example()
    
    #TopFace = Face(tF=0.5,material="S280GD")
    #BottomFace = Face(tF=0.5,material="S280GD")
    #Core = Core(tC=97,shear_modulus=3.7,density=110e-9)
    
    #SPA100 = SandwichPanel([TopFace,BottomFace],Core,
    #                       width=1200,length=6000)
    
    #kv = SPA100.fastener_transverse_stiffness(tsup=8)    
    #print("kv = {0:4.2f}".format(kv))
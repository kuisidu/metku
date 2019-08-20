# -*- coding: utf-8 -*-
"""
Calculation methods for sandwich panels with steel faces

Created on Wed Aug  7 09:15:25 2019

@author: kmela
"""

import math
from materials.steel_data import steel

fasteners = {'SXC14-S19-5.5':{'d1':5.5,'dS':5.1}}

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
        self.mat = material
        self.tcoat = tcoat
        try:
            self.fy = steel[material]['fy']
            self.fu = steel[material]['fu']
            self.E = steel[material]['E']
        except KeyError:
            print("Class 'Face' constructor: Material '{0:s}' not valid.".format(material))
            
    
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
    
    def __init__(self,faces,core,width=1200,length=6000):
        """ Constructor

            Parameters:
            -----------
                :param faces: list of class Face objects
                :param core: class core object
                :param width, length: self-evident [mm]

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
        
    @property
    def thickness(self):
        """ Total thickness of panel """
        return self.faces[0].t +self.faces[1].t + self.core.t
        
    @property
    def d(self):
        """ distance between center lines of faces """
        return 0.5*self.faces[0].t + 0.5*self.faces[1].t + self.core.t

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
        
        return EA[0]*EA[1]/sum(EA)*self.d
        

    def fastener_transverse_stiffness(self,tsup,fastener='SXC14-S19-5.5'):
        """ Transverse stiffness of the fastener according to
            ECCS recommendations
            
            :param tsup: thickness of supporting structure [mm]
            :param fastener: name of the fastener (see fasteners list)
        """
        d1 = fasteners[fastener]['d1']
        dS = fasteners[fastener]['dS']
        
        # bending stiffness of fastener
        EI = 200000*math.pi*dS**4/64
        
        # stiffness of clamping in the supporting structure
        Csup = 2400*math.sqrt(tsup*d1**5)
        
        # Stiffness of internal face (hole elongation)
        if self.faces[1].tcor < 0.7:
            kF2 = 6.93*self.faces[1].fu*math.sqrt(self.faces[1].tcor**3*d1)/(0.26+0.8*self.faces[1].t)
        else:
            kF2 = 4.2*self.faces[1].fu*math.sqrt(self.faces[1].tcor**3*d1)/0.373

        
        D = self.thickness
        xF = 1 - (1/kF2-0.5*D*tsup/Csup - D*tsup**2/8/EI)/(1/kF2+D**2/Csup+D**2*(2*D+3*tsup)/6/EI)
        
        kv = 1/(xF/kF2 + (tsup**2 + 2*(1-xF)*D*tsup)/4/Csup + (3*(1-xF)*D*tsup**2 + tsup**3)/24/EI)
        
        print("EI = {0:4.2f}".format(EI))
        print("Csup = {0:4.2f}".format(Csup))
        print("kF2 = {0:4.2f}".format(kF2))
        print("xF = {0:4.2f}".format(xF))
        
        print("1/k1 = {0:4.4g}".format(xF/kF2))
        print("1/k2 = {0:4.4g}".format((tsup**2 + 2*(1-xF*D*tsup))/4/Csup))
        print("1/k3 = {0:4.4g}".format((3*(1-xF)*D*tsup**2 + tsup**3)/24/EI))
        
        return kv

if __name__ == '__main__':
    
    #f = Face(0.5,material="S250GD")
            
    TopFace = Face(tF=0.5,material="S280GD")
    BottomFace = Face(tF=0.5,material="S280GD")
    Core = Core(tC=97,shear_modulus=3.7,density=110e-9)
    
    SPA100 = SandwichPanel([TopFace,BottomFace],Core,
                           width=1200,length=6000)
    
    kv = SPA100.fastener_transverse_stiffness(tsup=8)    
    print("kv = {0:4.2f}".format(kv))
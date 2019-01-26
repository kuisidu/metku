""" Class for steel members """

import sys
import math
import numpy as np

from sections.steel.steel_section import SteelSection
from eurocodes.en1993 import en1993_1_1
from eurocodes.en1993 import constants

class SteelMember:
    """ Steel members accroding to EN 1993
        
              
        Methods implemented in separate m files       
        % Computes the cost
        c = MemberCost(self)
        % Determine optimum profile
        %[C,sOpt] = OptimizeProfile(self,S)
        [C,sOpt] = OptimizeProfile(self,S,NLC,forceFun)


    """
    
    def __init__(self,profile,length,Lcr = [1.0,1.0],mtype="beam"):
        """ Constructor
            
            profile -- member profile (of cross_section class)
            length -- member length
            lcr -- array with length 2 including the buckling
                    length factors
            type -- member type (beam, column or truss)
            
            Other attributes:
            ned -- axial force in various sections (array)
            myed -- bending moment in various sections 
                    with respect to major axis (array)
            mzed -- bending moment in various sections
                    with respect to minor axis (array)
            vyed -- shear force in various sections
                    with respect to major axis (array)
            vzed -- shear force in various sections
                    with respect to minor axis (array)
            ted -- torsion moment in various sections (array)
            loc -- locations of sections along the members in
                    local coordinates (0 <= loc[i] <= 1)
        """
            
        self.profile = profile
        self.length = length
        self.lcr = np.array(Lcr)
        # self.lcr[0] = Lcr[0]*length
        # self.lcr[1] = Lcr[1]*length
        self.type = mtype
        self.ned = []
        self.myed = []
        self.mzed = []
        self.vyed = []
        self.vzed = []
        self.ted = []
        self.loc = []
    
    def nsect(self):
        """ Number of sections """
        return len(self.ned)
    
    def fy(self):
        """ Yield strength of the member """
        return self.profile.fy

    def ncrit(self):
        """ Buckling force according to Euler """
        #ncrit = [0.0,0.0]
        C = math.pi**2*self.profile.E
        ncrit = C*np.array(self.profile.I)/(self.length*self.lcr)**2
        #ncrit[0] = C*self.profile.I[0]/self.lcr[0]**2
        #ncrit[1] = C*self.profile.I[1]/self.lcr[1]**2
        return ncrit

    def slenderness(self):
        """ Non-dimensional slenderness according to EN 1993-1-1 """
        slend = [0,0]
        NRd = self.profile.A*self.fy()
        Ncr = self.ncrit()
        slend = np.sqrt(NRd/Ncr)
        #slend[0] = math.sqrt(NRd/Ncr[0])
        #slend[1] = math.sqrt(NRd/Ncr[1])
        return slend

    def LT_slenderness(self,Mcr):
        """ Non-dimensional slenderness for lateral-torsional buckling.
        """
        MRd = self.profile.bending_resistance()[0]
        #print("MRd: {}, Mcr: {}".format(MRd, Mcr))
        lambdaLT = math.sqrt(MRd/Mcr)
        return lambdaLT
    

    def buckling_strength(self):
        """ Member buckling strength according to EN 1993-1-1 """
        slend = self.slenderness()
        NbRd = []
        NRd = self.profile.A*self.fy()
        #self.profile.imp_factor()
        #print(slend)
        #p = 0.5*(1+np.array(self.profile.imp_factor)*(slend-0.2)+slend**2)
        #r = 1/(p+np.sqrt(p**2-slend**2))
        
        for i in range(len(slend)):
            p = 0.5*(1+self.profile.imp_factor[i]*(slend[i]-0.2)+slend[i]**2)
            r = min(1/(p+math.sqrt(p**2-slend[i]**2)),1.0)            
            NbRd.append(NRd*r)
        
        #NbRd = self.profile.A*self.fy()*r;
        return NbRd
    
    def LT_buckling_strength(self,Mcr):
        """ Member lateral-torsional bending strength """
        lambdaLT = self.LT_slenderness(Mcr)
        #self.profile.define_imp_factor_LT()
        #alpha = self.profile.ImperfectionLT
        alpha = self.profile.imp_factor
        # alpha[0] is the imperfection factor about y-axis
        p = 0.5*(1+alpha[0]*(lambdaLT-0.2)+lambdaLT**2)
        chiLT = min(1./(p+math.sqrt(p**2-lambdaLT**2)),1.0)
        MbRd = chiLT*self.profile.bending_resistance()[0]
        return MbRd
        
    def weight_per_length(self):
        """ Weight per unit length
            
            Default units: kg/mm
        """
        w = self.profile.A*constants.density
        return w
        
    def weight(self):
        """ Weight of the member (kg) """
        w = self.weight_per_length*self.length
        return w

    def mcrit(self,C=[1.12,0.45],za=0,k=[1,1]):
        """ Critical moment for lateral-torsional buckling
            
            For double-symmetric profiles
            C(1), C(2) -- coefficients
            za -- position of application of load
            k(1) -- kz
            k(2) -- kw
        """
                
        zs = 0.0
        zg = za-zs
        kz = k[0];
        kw = k[1];
            
        Iz = self.profile.I[1];
        Iw = self.profile.Iw;
        It = self.profile.It;
        E = eurocode3.young;
        G = eurocode3.G;
        L = self.length;
            
        Mcr = C[0]*math.pi**2*E*Iz/((kz*L)**2)*(math.sqrt((kz/kw)**2*Iw/Iz +
               (kz*L)**2*G*It/(math.pi**2*E*Iz) + (C[1]*zg)**2)-C[1]*zg)
        return Mcr

    def add_section(self,ned=0.0,myed=0.0,mzed=0.0,vyed=0.0,vzed=0.0,ted=0.0,loc=0.0):
        """ Adds new section with internal forces and location """

        self.ned.append(ned)
        self.myed.append(myed)
        self.mzed.append(mzed)
        self.vyed.append(vyed)
        self.vzed.append(vzed)
        self.ted.append(ted)
        self.loc.append(loc)

    def check_section(self,n=0,class1or2=False):
        """ Verify resistance of section 'n' """
        self.profile.Ned = self.ned[n]
        self.profile.Med = self.myed[n]
        self.profile.Ved = self.vzed[n]    
        
        r = self.profile.section_resistance(class1or2)
        return r
        
    def check_sections(self,class1or2=True):
        """ Verify resistance of all sections """
        r = []
        for n in range(self.nsect()):
            r.append(self.check_section(n,class1or2))
        return r
        
            
    def check_buckling(self):
        """ Verify buckling resistance
            output: r .. utilization ratio for buckling        
        """
        
        """ Find greatest compression force in the member
            clip replace all positive values with 0.0        
        """
        NEd = np.min(np.array(self.ned).clip(max=0.0))
        if NEd >= 0.0:
            """ if all axial forces are non-negative,
                buckling resistance need not be checked
            """
            r = [0.0,0.0]
        else:
            NbRd = self.buckling_strength()            
            r = abs(NEd)/NbRd
        
        return r
    
    def check_LT_buckling(self):
        """ Verify lateral torsional buckling resistance
            output: r .. utilization ratio for buckling        
        """
      
        
        MEd1 = np.max(np.array(self.myed))
        MEd2 = abs(np.min(np.array(self.myed)))
        MEd = max(MEd1, MEd2)
        MbRd = self.LT_buckling_strength(self.mcrit())            
        r = abs(MEd)/MbRd
        
        return r
        

    
    def check_beamcolumn(self,Cmy=0.9,class1or2=False):
        """ Verify stability for combined axial force and bending moment
        """
"""        
# For testing
if __name__ == "__main__":
    import sys
    
    from hollow_sections import SHS
    
    p = SHS(100.0,4.0)
    m = SteelMember(p,2e3,[0.7,1.2])
    
    #print(m.ncrit()*1e-3)
    #print(m.slenderness())
    
    m.add_section(ned=-250.0e3,myed=34.0e6,loc=0.0)
    m.add_section(ned=540.0e3,myed=-23.1e6,loc=0.5)
    
    #m.check_section(0)
    #m.check_section(1)
    
    #m.check_sections()
    m.check_buckling()
    """
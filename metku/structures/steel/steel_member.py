""" Class for steel members """

import math

import numpy as np

try:
    from metku.eurocodes.en1993 import constants, en1993_1_1
except:
    from eurocodes.en1993 import constants, en1993_1_1


class SteelMember:
    """ Steel members accroding to EN 1993
        
              
        Methods implemented in separate m files       
        % Computes the cost
        c = MemberCost(self)
        % Determine optimum profile
        %[C,sOpt] = OptimizeProfile(self,S)
        [C,sOpt] = OptimizeProfile(self,S,NLC,forceFun)


    """

    def __init__(self, profile, length, Lcr=[1.0, 1.0], mtype="beam",
                 LT_buckling=False, symmetry='dual'):
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
        self.LTB = LT_buckling
        self.symmetry = symmetry
        
        self.Cmy = 0.9
        self.Cmz = 0.9
        self.CmLT = 1.0
        
        
        self.utilization = []
        """ Utilization ratio for stability
            buck_y .. flexural buckling with respect to y axis
            buck_z .. flexural buckling with respect to z axis
            LTB .. lateral-torsional buckling
            bc_z .. beam-column stability, y axis
            bc_z .. beam-column stability, z axis
        """
        self.stability = {'buck_y':0.0,'buck_z':0.0,'LTB':0.0,'bc_y':0.0,'bc_z':0.0}

    @property
    def NbRd(self):
        return np.asarray(self.buckling_strength())

    @property
    def MbRd(self):
        Mcr = self.mcrit()
        return np.asarray(self.LT_buckling_strength(Mcr))
    
    @property
    def MEdY(self):
        return max(abs(np.array(self.myed)))
    
    @property
    def MEdZ(self):
        return max(abs(np.array(self.mzed)))

    @property
    def NEd(self):
        """ Maximum axial force (absolute value) """
        return max(abs(np.array(self.ned)))
    
    @property
    def NcEd(self):
        """ Maximum compressive axial force """
        return np.min(np.array(self.ned).clip(max=0.0))

    def nsect(self):
        """ Number of sections """
        return len(self.ned)

    def fy(self):
        """ Yield strength of the member """
        return self.profile.fy

    def ncrit(self, verb=False):
        """ Buckling force according to Euler """
        C = math.pi ** 2 * self.profile.E

        ncrit = C * np.array(self.profile.I) / (self.length * self.lcr) ** 2

        if verb:
            print("Ncr,y = {0:4.2f} kN".format(ncrit[0] * 1e-3))
            print("Ncr,z = {0:4.2f} kN".format(ncrit[1] * 1e-3))

        return ncrit

    def ncrit_T(self, verb=False):
        """ Critical laod for torsional buckling """
        GIt = self.profile.G*self.profile.It
        EIw = math.pi**2*self.profile.E*self.profile.Iw/self.length**2
        ncrit_T = 1/self.profile.i0()**2*(GIt + EIw)
        
        if verb:
            print("Ncr,T = {0:4.2f} kN".format(ncrit_T * 1e-3))            

        return ncrit_T

    def slenderness(self, verb=False):
        """ Non-dimensional slenderness according to EN 1993-1-1 """
        NRd = self.profile.A * self.fy()
        Ncr = self.ncrit(verb)
        # if NRd <= 1e-6:
        # print(self.profile.A)
        # print(self.profile.h)
        # print(NRd, Ncr)

        slend = np.sqrt(NRd / Ncr)

        if verb:
            print("lambda,y = {0:4.2f}".format(slend[0]))
            print("lambda,z = {0:4.2f}".format(slend[1]))

        return slend

    def LT_slenderness(self, Mcr, verb=False):
        """ Non-dimensional slenderness for lateral-torsional buckling.
        """
        MRd = self.profile.MRd
        lambdaLT = np.sqrt(MRd[0] / Mcr)
        if verb:
            # print("MRd: {}, Mcr: {}".format(MRd, Mcr))
            # print("MRd0 = {0:4.2f}".format(MRd[0] * 1e-6))
            if MRd[0] < 0:
                print(self.profile.h, self.profile.tw, self.profile.tf, self.profile.b)
            # print("MRd1 = {0:4.2f}".format(MRd[1] * 1e-6))
            # print("Mcr = {0:4.2f}".format(Mcr * 1e-6))
        return lambdaLT

    def buckling_strength(self, verb=False):
        """ Member buckling strength according to EN 1993-1-1 """

        if verb:
            print("** Buckling strength** ")

        slend = self.slenderness(verb)
        NbRd = []
        NRd = self.profile.A * self.fy()
        # self.profile.imp_factor()
        # print(slend)
        # p = 0.5*(1+np.array(self.profile.imp_factor)*(slend-0.2)+slend**2)
        # r = 1/(p+np.sqrt(p**2-slend**2))

        for i in range(len(slend)):
            p = 0.5 * (1 + self.profile.imp_factor[i] * (slend[i] - 0.2) +
                       slend[i] ** 2)
            r = min(1 / (p + math.sqrt(p ** 2 - slend[i] ** 2)), 1.0)
            if verb:
                print("Psi,{0:1} = {1:4.2f}".format(i, p))
                print("chi,{0:1} = {1:4.2f}".format(i, r))
                print("NbRd,{0:1} = {1:4.2f}".format(i, NRd * r * 1e-3))

            NbRd.append(NRd * r)

        # NbRd = self.profile.A*self.fy()*r;
        return NbRd

    def LT_buckling_strength(self, Mcr, axis='y', method='specific',
                             verb=False):
        """ Member lateral-torsional bending strength """

        MRd = self.profile.bending_resistance()[0]
        """
        if axis == 'y':
            idx = 0
        else:
            idx = 1
        """
        lambdaLT = self.LT_slenderness(Mcr) # [idx]
        
        # self.profile.define_imp_factor_LT()
        # alpha = self.profile.ImperfectionLT

        if method == 'general':
            alphaLT = self.profile.imp_factor_LT_gen
            lambdaLT0 = 0.2
            beta = 1.0
            chiLTmax = 1.0
        else:
            alphaLT = self.profile.imp_factor_LT
            lambdaLT0 = 0.4
            beta = 0.75
            chiLTmax = min(1.0, 1/lambdaLT**2)
        
        if lambdaLT <= lambdaLT0 or self.MEdY/Mcr <= lambdaLT0**2:
            chiLT = 1.0
        else:
            p = 0.5 * (1 + alphaLT * (lambdaLT - lambdaLT0) + beta*lambdaLT ** 2)
            chiLT = min(1. / (p + math.sqrt(p ** 2 - beta*lambdaLT ** 2)), chiLTmax)
        MbRd = chiLT * MRd
        
        if verb:
            print("lambdaLT = {0:5.4f}".format(lambdaLT))
            print("alphaLT = {0:5.4f}".format(alphaLT))
            # print("phi = {0:4.2f}".format(p))
            print("chi_LT = {0:5.4f}".format(chiLT))
            print("MbRd = {0:5.3f}".format(MbRd*1e-6))
            print("MRd = {0:5.3f}".format(MRd * 1e-6))

        return MbRd

    def weight_per_length(self):
        """ Weight per unit length
            
            Default units: kg/mm
        """
        w = self.profile.A * constants.density
        return w

    def weight(self):
        """ Weight of the member (kg) """
        w = self.weight_per_length() * self.length
        return w

    def mcrit(self, C=[1.12, 0.45, 0.525], k=[1, 1],za=0.0,verb=False):
        """ Critical moment for lateral-torsional buckling

            For mono- and double-symmetric profiles
            C(1), C(2), C(3) -- coefficients
            za -- position of application of load
            zs -- position of the shear centre S from centroid G
            zj -- ??
            k(1) -- kz
            k(2) -- kw
        """

        kz = k[0]
        kw = k[1]

        Iz = self.profile.I[1]
        Iw = self.profile.Iw
        It = self.profile.It

        E = constants.E
        G = constants.G
        L = self.length
        C1, C2, C3 = C

    

        if self.symmetry == 'dual':            
            zs = 0
            zg = za - zs
            part1 = C1 * (math.pi ** 2 * E * Iz) / (kz * L) ** 2
            part2 = np.sqrt((kz / kw) ** 2 * Iw / Iz + 
                            (kz * L) ** 2 * G * It /(math.pi ** 2 * E * Iz) + 
                            (C2 * zg) ** 2)
            part3 = C2 * zg
            
            #print("part 1 = {0:4.2f}".format(part1))
            #print("part 2 = {0:4.2f}".format(part2))
            #print("part 3 = {0:4.2f}".format(part3))

            Mcr = part1 * (part2 - part3)            

        elif self.symmetry == 'mono':
            za = 0
            zs = self.profile.zs
            zg = za - zs

            zj = self.profile.zj

            part1 = C1 * (math.pi ** 2 * E * Iz) / (kz * L) ** 2
            part2 = C3 * zj - C2 * zg
            part3 = np.sqrt((kz / kw) ** 2 * Iw / Iz + (kz * L) ** 2 * G * It /
                            (math.pi ** 2 * E * Iz) + (C2 * zg - C3 * zj) ** 2)

            Mcr = part1 * (part2 + part3)
            
        if verb:
            print("Critical moment for LT buckling")
            print("kz = {0:4.2f}; kw = {1:4.2f}".format(kz,kw))
            print("C1 = {0:4.2f}; C2 = {1:4.2f}".format(C1,C2))
            print("Mcr = {0:6.4f}".format(Mcr*1e-6))

        # print("Iz = {0:4.2f}".format(Iz*1e-4))
        # print("Iw = {0:4.2f}".format(Iw*1e-6))
        # print("It = {0:4.2f}".format(It*1e-4))
        # print("Mcr = {0:4.2f}".format(Mcr*1e-6))
        # print("kz =", kz, "kw =", kw, "L =", L, "part3 =", part3))

        return Mcr

    def add_section(self, ned=0.0, myed=0.0, mzed=0.0, vyed=0.0, vzed=0.0,
                    ted=0.0, loc=0.0):
        """ Adds new section with internal forces and location """

        self.ned.append(ned)
        self.myed.append(myed)
        self.mzed.append(mzed)
        self.vyed.append(vyed)
        self.vzed.append(vzed)
        self.ted.append(ted)
        self.loc.append(loc)
        self.utilization.append({})

    def clear_sections(self):
        self.ned.clear()
        self.myed.clear()
        self.mzed.clear()
        self.vyed.clear()
        self.vzed.clear()
        self.ted.clear()
        self.loc.clear()


    def check_section(self, n=0, verb=False):
        """ Verify resistance of section 'n' """
        self.profile.Ned = self.ned[n]
        self.profile.Med = self.myed[n]
        self.profile.Ved = self.vzed[n]
                

        r = self.profile.section_resistance(verb=verb)
        
        self.utilization[n] = {'n':r[0],'vz':r[1],'my':r[2],'mny':r[3]}
        return r

    def check_sections(self, class1or2=True, verb=False):
        """ Verify resistance of all sections """
        r = np.zeros_like(self.check_section())
        for n in range(self.nsect()):
            r = np.vstack((r, self.check_section(n)))
            # r.append(self.check_section(n))
        return np.max(r, axis=0)

    def check_buckling(self):
        """ Verify buckling resistance
            output: r .. utilization ratio for buckling        
        """

        """ Find greatest compression force in the member
            clip replace all positive values with 0.0        
        """
        NEd = self.NcEd
        if NEd >= 0.0:
            """ if all axial forces are non-negative,
                buckling resistance need not be checked
            """
            r = [0.0, 0.0]
        else:
            r = abs(NEd) / self.NbRd
        
        self.stability['buck_y'] = r[0]
        self.stability['buck_z'] = r[1]


        return r

    def check_LT_buckling(self):
        """ Verify lateral torsional buckling resistance
            output: r .. utilization ratio for buckling        
        """
        MEd = self.MEdY
        MbRd = self.LT_buckling_strength(self.mcrit())
        r = abs(MEd) / MbRd

        return r

    def check_beamcolumn(self, section_class=None):
        """ Verify stability for combined axial force and bending moment
        """
        # find the largest compression force in the member
        NEd = self.NcEd
        self.profile.Ned = NEd
        
        if section_class is None:
            cross_section_class = self.profile.section_class()
        else:
            cross_section_class = section_class
            
        #print("Check beam column:")
        #print("Section class:", cross_section_class)
        
        slend = self.slenderness()
    
        NbRd = self.buckling_strength()

        UNy = abs(NEd) / NbRd[0]
        UNz = abs(NEd) / NbRd[1]

        kyy = en1993_1_1.kyy(UNy, slend[0], self.Cmy,
                             section_class=cross_section_class)

        kzy = en1993_1_1.kzy(kyy, UNz, slend[1], CmLT=self.CmLT,
                             section_class=cross_section_class)

        """ Find maximum bending moment along the member """
        MyEd = self.MEdY

        """ For now, assume no bending with respect to z axis """
        # MzEd1 = np.max(np.array(self.mzed))
        # MzEd2 = abs(np.min(np.array(self.mzed)))
        # MzEd = max(MzEd1, MzEd2)

        #MRd = self.profile.bending_resistance(C=cross_section_class)
        MRd = self.profile.MRd
      
        com_comp_bend_y = -NEd / NbRd[0] + kyy * MyEd / MRd[0]
        com_comp_bend_z = -NEd / NbRd[1] + kzy * MyEd / MRd[0]

        com_comp_bend = [com_comp_bend_y, com_comp_bend_z]
        
        self.stability['bc_y'] = com_comp_bend_y
        self.stability['bc_z'] = com_comp_bend_z
        
        return com_comp_bend
    
    def required_shear_stiffness(self):
        """ EN 1993-1-1, Eq. (BB.2) """
        E = self.profile.E
        G = self.profile.G
        L = self.length
        h = self.profile.h
        
        Sreq = (E*self.profile.Iw*math.pi**2/L**2 + G*self.profile.It + E*self.profile.Iz*math.pi**2/L**2*0.25*h**2)*70/h**2
        
        return Sreq

    def design(self,verb=False):
        """ Carry out member design 
            
            Returns:
                rmax .. maximum utilization ratio
        
        """
        
        """ Design cross sections """
        rsec = self.check_sections(class1or2=False,verb=verb)
      
        rmax = np.max(rsec)  
      
        """ Check stability """
        if self.LTB:
            rLTB = self.check_LT_buckling()        
            rmax = max(rmax,rLTB)
        
        if self.NcEd < 0.0:
            rbuck = self.check_buckling()
            rbeam_column = self.check_beamcolumn()
            rmax = max(rmax,np.max(rbuck))
            rmax = max(rmax,np.max(rbeam_column))
            
        return rmax
        
        
        
        
        
        

if __name__ == "__main__":
    
    from metku.sections.steel.ISection import HEA
    from metku.sections.steel.RHS import RHS
    
    p = HEA(240)
    R = RHS(200,200,8)
    p = RHS(100,100,5)
    m = SteelMember(p,3000)
    mr = SteelMember(R,6000)    
    m.ncrit(True)
    print(m.weight())
    #m.ncrit_T(True)
    #mr.ncrit_T(True)
    #mr.ncrit(True)
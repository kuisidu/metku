# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
""" Class for steel members """

import math

import numpy as np

    
from metku.eurocodes.en1993 import constants, en1993_1_1
from metku.sections.steel.wi import WISection


class WebTaperedSteelMember:
    """ Web tapered steel members accroding to EN 1993                     
        TODO:
            - NbRd
            - MbRd            
            - Nurjahdus
            - Kiepahdus

    """

    def __init__(self, profile0, h1, length, Lcr=[1.0, 1.0], mtype="beam",
                 LT_buckling=False,taper_type="right"):
        """ Constructor
            
            profile0 -- smaller member profile
            h1 -- height of the profile at the other end
            length -- member length
            lcr -- array with length 2 including the buckling
                    length factors
            type -- member type (beam, column or truss)
            taper_type .. type of tapering
                        "right" right hand side of the profile is tapered
                        "left" left hand side of the profile is tapered
                        "dual" both sides are tapered
            
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
        
        # create 'top profile'
        profile1 = WISection(h1,profile0.tw,profile0.b,profile0.tf,profile0.fy,profile0.weld_throat)

        self.profiles = [profile0,profile1]
        self.length = length
        self.lcr = np.array(Lcr)
        self.h1 = h1
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
        self.LT_B = LT_buckling
        self.taper_type = taper_type

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
        return max(abs(np.array(self.ned)))

    @property
    def Iymax(self):
        return self.profiles[1].I[0]

    @property
    def Iymin(self):
        return self.profiles[0].I[0]

    @property
    def Ieq(self):
        """ Equivalent second moment of area with respect to major axis """
        
        Iymax = self.Iymax
        r = math.sqrt(self.Iymin/Iymax)
        C = 0.08+0.92*r
        
        return C*Iymax

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
        if axis == 'y':
            idx = 0
        else:
            idx = 1
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
            chiLT = min(1. / (p + math.sqrt(p ** 2 - beta*lambdaLT ** 2)), chiLTmax
                        )
        MbRd = chiLT * MRd
        
        if verb:
            print("lambdaLT = {0:4.2f}".format(lambdaLT))
            print("alphaLT = {0:4.2f}".format(alphaLT))
            # print("phi = {0:4.2f}".format(p))
            print("chi_LT = {0:4.2f}".format(chiLT))
            print("MbRd = {0:4.2f}".format(MbRd*1e-6))
            print("MRd = {0:4.2f}".format(MRd * 1e-6))

        return MbRd

    def weight_per_length(self):
        """ Weight per unit length
            
            Default units: kg/mm
        """
        w = self.profile.A * constants.density
        return w

    def weight(self):
        """ Weight of the member (kg) """
        w = self.weight_per_length * self.length
        return w

    def mcrit(self, C=[1.12, 0.45, 0.525], k=[1, 1]):
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
            za = 0
            zs = 0
            zg = za - zs
            part1 = C1 * (math.pi ** 2 * E * Iz) / (kz * L) ** 2
            part2 = np.sqrt((kz / kw) ** 2 * Iw / Iz + (kz * L) ** 2 * G * It /
                            (math.pi ** 2 * E * Iz) + (C2 * zg) ** 2)
            part3 = C2 * zg

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

    def check_section(self, n=0):
        """ Verify resistance of section 'n' """
        self.profile.Ned = self.ned[n]
        self.profile.Med = self.myed[n]
        self.profile.Ved = self.vzed[n]

        r = self.profile.section_resistance()
        return r

    def check_sections(self, class1or2=True):
        """ Verify resistance of all sections """
        r = np.zeros(len(self.check_section()))
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
        NEd = np.min(np.array(self.ned).clip(max=0.0))
        if NEd >= 0.0:
            """ if all axial forces are non-negative,
                buckling resistance need not be checked
            """
            r = [0.0, 0.0]
        else:
            r = abs(NEd) / self.NbRd

        return r

    def check_LT_buckling(self):
        """ Verify lateral torsional buckling resistance
            output: r .. utilization ratio for buckling        
        """

        MEd1 = np.max(np.array(self.myed))
        MEd2 = abs(np.min(np.array(self.myed)))
        MEd = max(MEd1, MEd2)
        MbRd = self.LT_buckling_strength(self.mcrit())
        r = abs(MEd) / MbRd

        return r

    def check_beamcolumn(self, Cmy=0.9, section_class=None):
        """ Verify stability for combined axial force and bending moment
            LT_buckling -- kiepahdus
                    True .. rajoitusehto huomioidaan
                    False .. rajoitusehtoa ei huomioida
        """
        # find the largest compression force in the member
        NEd = np.min(np.array(self.ned).clip(max=0.0))

        self.profile.Ned = NEd
        if section_class is None:
            cross_section_class = self.profile.section_class()
        else:
            cross_section_class = section_class
            
        #print("Check beam column:")
        #print("Section class:", cross_section_class)
        
        slend = self.slenderness()
        CmLT = 1  # Pitääkö käyttää jotakin kaavaa?

        # NRd = self.profile.A * self.fy()
        NbRd = self.buckling_strength()

        UNy = abs(NEd) / NbRd[0]
        UNz = abs(NEd) / NbRd[1]

        kyy = en1993_1_1.kyy(UNy, slend[0], Cmy,
                             section_class=cross_section_class)

        kzy = en1993_1_1.kzy(kyy, UNz, slend[1], CmLT=CmLT,
                             section_class=cross_section_class)

        MyEd1 = np.max(np.array(self.myed))
        MyEd2 = abs(np.min(np.array(self.myed)))
        MyEd = max(MyEd1, MyEd2)

        MRd = self.profile.bending_resistance(C=cross_section_class)    

        com_comp_bend_y = -NEd / NbRd[0] + kyy * MyEd / MRd[0]
        com_comp_bend_z = -NEd / NbRd[1] + kzy * MyEd / MRd[0]

        com_comp_bend = [com_comp_bend_y, com_comp_bend_z]

        return com_comp_bend

if __name__ == '__main__':
    
    p = WISection(200, 6.0, [140.0, 140.0], [8.0, 8.0],fy=355,weld_throat=4)
    
    m = WebTaperedSteelMember(p, 400, 8000, Lcr=[1.0, 1.0], mtype="column",LT_buckling=True,taper_type="right")
    
    print(m.Iymin,m.Ieq,m.Iymax)    
    
# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
""" Class for steel members """

import math
import numpy as np

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
        #self.ned = []
        #self.myed = []
        #self.mzed = []
        #self.vyed = []
        #self.vzed = []
        #self.ted = []
        #self.loc = []
        self.LTB = LT_buckling
        self.symmetry = symmetry
        
        # Dictionary of sections
        # key = local coordinate running from 0 to 1
        # value = dict with key, value pairs
        # {'ned':axial_force, 'myed':bending_moment_yaxis,
        #  'mzed': bending_moment_zaxis, 'vyed': shear force, y axis,
        #  'vzed': shear force z axis, 'ted': torsional moment,
        # 'utilization': utilization ratio in different loading cases}
        self.sections = {}
        
        self.__Cmy = 0.9
        self.__Cmz = 0.9
        self.__CmLT = 1.0
        
        
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
    def ned(self):
        """ Return a list of axial forces in the sections """
        return np.array([sec['N'] for sec in self.sections.values()])
    
    @ned.setter
    def ned(self,val):
        """ Set the same axial force for all sections """
        
        for sec in self.sections.values():
            sec['N'] = val        
        
    @property
    def myed(self):
        """ Return a list of bending moments MyEd in the sections """
        return np.array([sec['My'] for sec in self.sections.values()])
    
    @myed.setter
    def myed(self,val):
        for sec in self.sections.values():
            sec['My'] = val
    
    @property
    def mzed(self):
        """ Return a list of bending moments MzEd in the sections """
        return np.array([sec['Mz'] for sec in self.sections.values()])
    
    @mzed.setter
    def mzed(self,val):
        for sec in self.sections.values():
            sec['Mz'] = val
    
    @property
    def vyed(self):
        """ Return a list of shear forces VyEd in the sections """
        return np.array([sec['Vy'] for sec in self.sections.values()])
    
    @vyed.setter
    def vyed(self,val):
        for sec in self.sections.values():
            sec['Vy'] = val
    
    @property
    def vzed(self):
        """ Return a list of shear forces VzEd in the sections """
        return np.array([sec['Vz'] for sec in self.sections.values()])

    @vzed.setter
    def vzed(self,val):
        for sec in self.sections.values():
            sec['Vz'] = val

    @property
    def ted(self):
        """ Return a list of torques TEd in the sections """
        return np.array([sec['T'] for sec in self.sections.values()])
    
    @ted.setter
    def ted(self,val):
        for sec in self.sections.values():
            sec['T'] = val
    
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
        #return len(self.ned)
        return len(self.sections)

    def fy(self):
        """ Yield strength of the member """
        return self.profile.fy
    
    def Cmy(self,sway=True,load="uniform"):
        """ Equivalent moment factor """
        if sway:
            return 0.9
        else:
            try:
                # Get end moments
                Mend = [self.sections[0]['My'],self.sections[1]['My']]
            except:
                # If no end moments are given, use the value stored in Cmy variable
                return self.__Cmy
            
            # Get span moment. If it is not given, assume linear distribution
            try:                    
                Mspan = self.sections[0.5]['My']
            except:
                Mspan = None
            
            self.__Cmy = en1993_1_1.equivalent_moment_factor(Mend,Mspan,load)
                            
        return self.__Cmy

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
            print("lambda,y = {0:4.3f}".format(slend[0]))
            print("lambda,z = {0:4.3f}".format(slend[1]))

        return slend

    def slenderness_torsional(self, verb=False):
        """ Non-dimensional slenderness for torsional buckling
            according to EN 1993-1-1 
        """
        NRd = self.profile.A * self.fy()
        Ncr = self.ncrit_T(verb)
        # if NRd <= 1e-6:
        # print(self.profile.A)
        # print(self.profile.h)
        # print(NRd, Ncr)

        slend = np.sqrt(NRd / Ncr)

        if verb:
            print("lambda_T = {0:4.2f}".format(slend))
           
        return slend

    def LT_slenderness(self, Mcr, verb=False):
        """ Non-dimensional slenderness for lateral-torsional buckling.
        """
        MRd = self.profile.MRd
        lambdaLT = np.sqrt(MRd / Mcr)
        if verb:
            # print("MRd: {}, Mcr: {}".format(MRd, Mcr))
            # print("MRd0 = {0:4.2f}".format(MRd[0] * 1e-6))
            if MRd < 0:
                print(self.profile.h, self.profile.tw, self.profile.tf, self.profile.b)
            # print("MRd1 = {0:4.2f}".format(MRd[1] * 1e-6))
            # print("Mcr = {0:4.2f}".format(Mcr * 1e-6))
        return lambdaLT

    def buckling_strength(self, verb=False):
        """ Member buckling strength according to EN 1993-1-1 """

        if verb:
            print("** Buckling strength ** ")

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
                print("Psi,{0:1} = {1:4.3f}".format(i, p))
                print("chi,{0:1} = {1:4.3f}".format(i, r))
                print("NbRd,{0:1} = {1:4.2f}".format(i, NRd * r * 1e-3))

            NbRd.append(NRd * r)

        return NbRd
    
    def torsional_buckling_strength(self, verb=False):
        """ Member torsional buckling strength according to EN 1993-1-1 
            EN 1993-1-1, 6.3.1.4        
        """
        if verb:
            print("** Torsional buckling strength** ")
        
        slend = self.slenderness_torsional(verb)        
        NRd = self.profile.A * self.fy()
                
        p = 0.5*(1 + self.profile.imp_factor[1]*(slend-0.2) + slend**2)
        r = min(1/(p + math.sqrt(p**2-slend**2)), 1.0)
        NbRd = r*NRd
        if verb:      
            print(f"lambda_t = {slend:4.3f}")
            print(f"Psi = {p:4.2f}")
            print(f"chi = {r:4.2f}")
            print(f"NbRd = {NbRd*1e-3:4.2f}")
        
        return NbRd

    def LT_buckling_strength(self, Mcr, method='specific',verb=False):
        """ Member lateral-torsional buckling strength """

        MRd = self.profile.bending_resistance()[0]
        
        lambdaLT = self.LT_slenderness(Mcr)
        
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
        return self.profile.weight()

    def weight(self):
        """ Weight of the member (kg) """
        w = self.weight_per_length() * self.length
        return w
    
    def cost(self):
        """ Cost of member (€, $, etc.).
            By default, the cost includes only material cost.
            
            More specific cost function should be implemented for
            particular steel member using the class of that member
        """
        
        return self.profile.cost() * self.length
        

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
            print("Mcr = {0:6.4f} kNm".format(Mcr*1e-6))

        # print("Iz = {0:4.2f}".format(Iz*1e-4))
        # print("Iw = {0:4.2f}".format(Iw*1e-6))
        # print("It = {0:4.2f}".format(It*1e-4))
        # print("Mcr = {0:4.2f}".format(Mcr*1e-6))
        # print("kz =", kz, "kw =", kw, "L =", L, "part3 =", part3))

        return Mcr
    
    def mcrit_cantilever(self,load="uniform",loc="middle",warping_free=True,verb=False):
        """ Critical buckling moment for cantilevers with I-sections 
            Based on Andrade et al. (2007), J Constr Steel Res.
            
            load .. type of load "uniform" or "point"
            loc .. location of load "middle", "top" or "bottom"
            warping .. is warping of the fixed end free
        """
        C1 = 1.0
        C2 = 0.0
        
        # hS is the distance between center lines of the flanges
        hS = self.profile.h-self.profile.tf
        
        # K is a constant that is used in the equations.
        K = np.pi/self.length*np.sqrt(self.profile.E*self.profile.Iz*hS**2/4/self.profile.G/self.profile.It)
        
        if K < 0.1 or K > 2.5:
            print(f"Warning: constant K = {K:4.2f} is out of scope")
        
        denom = np.sqrt(1+K**2)
        
        if warping_free:
            if load == 'point':
                C1 = 2.437/denom + 0.613*K/denom - 0.105*K**2/denom
                
                if loc == "top":
                    C2 = 0.409 + 1.444*K + 0.070*K**2
                elif loc == "bottom":
                    C2 = 0.529 + 0.234*K + 0.149*K**2
            else:
                C1 = 3.840/denom + 1.496*K/denom - 0.247*K**2/denom
                
                if loc == "top":
                    C2 = 0.984 + 1.420*K + 0.165*K**2
                elif loc == "bottom":
                    C2 = 1.028 + 0.388*K + 0.150*K**2
        else:
            # Warping fully restrained
            if load == 'point':
                C1 = 2.462/denom + 2.383*K/denom
                
                if loc == "top":
                    C2 = 0.380 + 2.092*K - 0.318*K**2
                elif loc == "bottom":
                    C2 = 0.512 + 0.370*K -0.033*K**2
            else:
                C1 = 3.962/denom + 5.531*K/denom
                
                if loc == "top":
                    C2 = 1.130 + 1.539*K - 0.176*K**2
                elif loc == "bottom":
                    C2 = 1.049 + 0.234*K - 0.020*K**2
        
        if loc == "top":
            za = 0.5*self.profile.h
        elif loc == "bottom":
            za = -0.5*self.profile.h
        else:
            za = 0.0
        
        Mcr = self.mcrit([C1,C2,0],[2,1], za,verb)
        
        if verb:
            print(f'K = {K:4.3f}')
        
        return Mcr
                
    
    def mcrit_buckling(self,verb=False):
        """ Evaluate elastic critical moment by buckling of compressed
            flange
        """
        
        # Determine compressed flange and its
        # buckling load
        _, Ifz = self.profile.compressed_flange()
        Ncrz = np.pi**2*self.profile.E*Ifz/self.length**2
        
        Mcr = Ncrz * self.profile.h
        
        if verb:
            print("Critical moment for LT buckling using compressed flange")
            print("Ifz = {0:4.2f} mm4".format(Ifz))
            print("Ncr_z = {0:4.2f} kN".format(Ncrz*1e-3))
            print("Mcr = {0:6.4f} kNm".format(Mcr*1e-6))

        
    def add_section(self, ned=0.0, myed=0.0, mzed=0.0, vyed=0.0, vzed=0.0,
                    ted=0.0, loc=0.0):
        """ Adds new section with internal forces and location """

        self.sections[loc] = {'N': ned, 'My': myed, 'Mz': mzed,
                              'Vy': vyed, 'Vz': vzed, 'T': ted,
                              'utilization': {}}

        """
        self.ned.append(ned)
        self.myed.append(myed)
        self.mzed.append(mzed)
        self.vyed.append(vyed)
        self.vzed.append(vzed)
        self.ted.append(ted)
        self.loc.append(loc)
        self.utilization.append({})
        """
        
    def clear_sections(self):
        
        self.sections.clear()
        
        """
        self.ned.clear()
        self.myed.clear()
        self.mzed.clear()
        self.vyed.clear()
        self.vzed.clear()
        self.ted.clear()
        self.loc.clear()
        self.utilization.clear()
        """
    
    def check_section(self, loc=0, verb=False):
        """ Verify resistance of section at location 'loc' """
        self.profile.Ned = self.sections[loc]['N']
        self.profile.Med[0] = self.sections[loc]['My']
        self.profile.Ved[1] = self.sections[loc]['Vz']
                

        r = self.profile.section_resistance(verb=verb)
        
        self.sections[loc]['utilization'] = {'n':r[0],'vz':r[1],'my':r[2],'mny':r[3]}
        return r

    """
    def check_section(self, n=0, verb=False):
        # Verify resistance of section 'n'
        self.profile.Ned = self.ned[n]
        self.profile.Med = self.myed[n]
        self.profile.Ved = self.vzed[n]
                

        r = self.profile.section_resistance(verb=verb)
        
        self.utilization[n] = {'n':r[0],'vz':r[1],'my':r[2],'mny':r[3]}
        return r
    """
    def check_sections(self, class1or2=True, verb=False):
        """ Verify resistance of all sections """
        r = 0.0
        for loc in self.sections.keys():
            r = max(r,max(self.check_section(loc,verb)))
        
        return r
        
        #r = np.zeros_like(self.check_section())
        #for n in range(self.nsect()):
        #    r = np.vstack((r, self.check_section(n)))
            # r.append(self.check_section(n))
        #return np.max(r, axis=0)

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

    def check_beamcolumn(self, sway=False, load="uniform", section_class=None,verb=False):
        """ Verify stability for combined axial force and bending moment
        """
        
        """ Find maximum bending moment along the member """        
        MyEd = self.MEdY
        
        if abs(MyEd) < 1e-6:
            # If there is no bending moment, there is no need to
            # continue
            return [0,0]
            
        
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

        # Determine equivalent moment factor
        Cmy = self.Cmy(sway,load)        

        kyy = en1993_1_1.kyy(UNy, slend[0], Cmy,
                             section_class=cross_section_class)

        if self.LTB:
            # Member is susceptible to lateral-torsional buckling
            kzy = en1993_1_1.kzy(kyy, UNz, slend[1], CmLT=self.CmLT,
                                 section_class=cross_section_class)
        else:
            MzEd = self.MEdZ
            
            if MzEd == 0:
                # Uniaxial bending: according to current EN 1993-1-1,
                # kzy = 0 is possible
                # In the revised eurocode, this is no longer valid.
                kzy = 0
            else:
                if cross_section_class < 3:
                    kzy = 0.6*kyy
                else:
                    kzy = 0.8*kyy
        
        if verb:
            print("Beam-column stability:")
            print(f'Cmy = {Cmy:4.2f}')
            print(f'kyy = {kyy:4.3f}')
            print(f'kzy = {kzy:4.3f}')

        """ For now, assume no bending with respect to z axis """
        # MzEd1 = np.max(np.array(self.mzed))
        # MzEd2 = abs(np.min(np.array(self.mzed)))
        # MzEd = max(MzEd1, MzEd2)

        #MRd = self.profile.bending_resistance(C=cross_section_class)
        MRd = self.profile.MRd
      
        com_comp_bend_y = -NEd / NbRd[0] + kyy * MyEd / MRd
        com_comp_bend_z = -NEd / NbRd[1] + kzy * MyEd / MRd

        com_comp_bend = [com_comp_bend_y, com_comp_bend_z]
        
        self.stability['bc_y'] = com_comp_bend_y
        self.stability['bc_z'] = com_comp_bend_z
        
        if verb:
            if com_comp_bend_y <= 1.0:
                yOK = 'OK!'
            else:
                yOK = 'NOT OK!'
                
            if com_comp_bend_z <= 1.0:
                zOK = 'OK!'
            else:
                zOK = 'NOT OK!'
            print(f"y-axis: {com_comp_bend_y:4.3f} <= 1.0 ({yOK})")
            print(f"z-axis: {com_comp_bend_z:4.3f} <= 1.0 ({zOK})")
        
        return com_comp_bend
    
    def required_shear_stiffness(self):
        """ EN 1993-1-1, Eq. (BB.2) """
        E = self.profile.E
        G = self.profile.G
        L = self.length
        h = self.profile.h
        
        Sreq = (E*self.profile.Iw*math.pi**2/L**2 + G*self.profile.It + E*self.profile.Iz*math.pi**2/L**2*0.25*h**2)*70/h**2
        
        return Sreq

    def design(self,sway=False,verb=False):
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
            #print("Check stability")
            rbuck = self.check_buckling()            
            rbeam_column = self.check_beamcolumn(sway,verb=verb)
            rmax = max(rmax,np.max(rbuck))
            rmax = max(rmax,np.max(rbeam_column))

        return rmax
        
        
        
"""

if __name__ == "__main__":
    
    from metku.sections.steel.ISection import HEA
    from metku.sections.steel.RHS import RHS
    
    #p = HEA(240)
    #R = RHS(200,200,8)
    p = RHS(160,160,8)
    #m = SteelMember(p,3000)
    mr = SteelMember(p,3680,Lcr=[0.9,1.0])
    mr.add_section(loc=0.0,ned=-725.57e3,myed=-22.12e6)
    mr.add_section(loc=0.5,ned=-725.57e3,myed=17.07e6)
    mr.add_section(loc=1.0,ned=-725.57e3,myed=-24.95e6)
    
    #r = mr.check_section()
    #m.ncrit(True)
    #print(m.weight())
    #m.ncrit_T(True)
    #mr.ncrit_T(True)
    #mr.ncrit(True)
"""
# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 22:42:48 2018

Class for general steel cross sections 
    NOTE: This could be a package with the following structure
        cross_section
            hollow_sections
                RHS
                SHS
                CHS
            I_sections


@author: kmela
"""

from abc import abstractclassmethod, ABCMeta

import numpy as np

from metku.eurocodes.en1993 import constants
from metku.eurocodes.en1993 import en1993_1_1
from metku.materials.steel_data import Steel

INFEASIBLE = 999


class SteelSection(metaclass=ABCMeta):
    """ Base class for cross-sections
    
    Specific cross sections are subclasses of this class
    
    """

    # Constructor
    def __init__(self,
                 material,
                 A,
                 I,
                 Au,
                 Wpl,
                 Wel,
                 Ashear,
                 Ned=0.0,
                 Med=0.0,
                 Ved=0.0,
                 Ted=0.0,
                 design_code=en1993_1_1):
        """

        Parameters
        ----------
        material : int or str
            Either yields strength or name of the material.
        A : float
            Cross-sectional area.
        I : List
            second moment of area with respect to principal axes.
            I[0] .. major axis
            I[1] .. minor axis
        Au : float
            circumferential area per unit length.
        Wpl : list
            plastic section modulus with respect to principal axes.
        Wel : list
            elastic section modulus with respect to principal axes..
        Ashear : float
            shear area.
        Ned : float, optional
            design value of axial force. The default is 0.0.
        Med : list or float, optional
            bending moment. The default is 0.0.
        Ved : list or float, optional
            shear force. The default is 0.0.
            Ved[0] .. shear force is along the major axis
            Ved[1] .. shear force is along the minor axis
        Ted: float, optional
            torque, The default is 0.0
        design_code : TYPE, optional
            DESCRIPTION. The default is en1993_1_1.

        Returns
        -------
        None.

        """

        # print(material)
        # if isinstance(material,str):
        #    self.material = Steel(material)
        # elif isinstance(material,(float,int)):
        #    self.material = Steel("S" + str(int(material)))

        self.material = material

        self.E = constants.E
        self.A = A
        self.I = np.asarray(I)
        self.Au = Au
        self.Wpl = np.asarray(Wpl)
        self.Wel = np.asarray(Wel)

        if isinstance(Ashear, list):
            self.Ashear = Ashear
        else:
            self.Ashear = [0.0, Ashear]

        self.Ned = Ned
        if not isinstance(Med, list):
            self.Med = [Med, 0.0]
        else:
            self.Med = Med

        if not isinstance(Ved, list):
            self.Ved = [0.0, Ved]
        else:
            self.Ved = Ved

        self.Ted = Ted

        self.imp_factor = None
        self.code = design_code
        self.density = constants.density
        self.rotate = False

        # Cost of section, €/kg (or $/kg, or...)
        self.unit_cost = 1.0

    def __repr__(self):
        return type(self).__name__

    def cost(self):
        """ Material cost of section, units: €/m """
        return self.weight() * self.unit_cost

    @property
    def material(self):
        return self.__material

    @material.setter
    def material(self, val):
        """ Sets the material of section """
        if isinstance(val, str):
            self.__material = Steel(val)
        elif isinstance(val, (float, int)):
            self.__material = Steel("S" + str(int(val)))

    @property
    def G(self):
        """ Shear modulus """
        return self.material.G

    @property
    def fy(self):
        """ Yield strength """

        return self.material.fy

    @property
    def fu(self):
        """ Ultimate strength """

        return self.material.fu

    @property
    def Iy(self):
        if self.rotate:
            return self.I[1]
        return self.I[0]

    @Iy.setter
    def Iy(self, val):
        if self.rotate:
            self.I[1] = val
        self.I[0] = val

    @property
    def Iz(self):
        if self.rotate:
            return self.I[0]
        return self.I[1]

    @Iz.setter
    def Iz(self, val):
        if self.rotate:
            self.I[0] = val
        self.I[1] = val

    @property
    def Wply(self):
        if self.rotate:
            return self.Wpl[1]
        return self.Wpl[0]

    @Wply.setter
    def Wply(self, val):
        if self.rotate:
            self.Wpl[1] = val
        self.Wpl[0] = val

    @property
    def Wplz(self):
        if self.rotate:
            return self.Wpl[0]
        return self.Wpl[1]

    @Wplz.setter
    def Wplz(self, val):
        if self.rotate:
            self.Wpl[0] = val
        self.Wpl[1] = val

    @property
    def Wely(self):
        if self.rotate:
            return self.Wel[1]
        return self.Wel[0]

    @Wely.setter
    def Wely(self, val):
        if self.rotate:
            self.Wel[1] = val
        self.Wel[0] = val

    @property
    def Welz(self):
        if self.rotate:
            return self.Wel[0]
        return self.Wel[1]

    @Welz.setter
    def Welz(self, val):
        if self.rotate:
            self.Wel[0] = val
        self.Wel[1] = val

    @property
    def section_factor(self):
        """ Section factor for unprotected profile,
            exposed to fire on all sides
        """
        return self.Au / self.A * 1e3

    @property
    def box_section_factor(self):
        """ Section factor when considering the
            profile as a box
        """
        return self.section_factor

    def shadow_effect(self):
        """ Shadow effect for fire, EN 1993-1-2 """

        return 1.0

    def iy(self):
        """ Radius of gyration of the major axis """
        return np.sqrt(self.Iy / self.A)

    def iz(self):
        """ Radius of gyration of the minor axis """
        return np.sqrt(self.Iz / self.A)

    def i0(self):
        """ Polar radius """
        return np.sqrt(self.iy() ** 2 + self.iz() ** 2)

    @property
    def eps(self):
        return self.code.epsilon(self.fy)

    @property
    def NRd(self):
        if self.Ned > 0:
            return self.code.compression_resistance(self.A, self.fy)

        else:
            return self.code.compression_resistance(self.A, self.fy)

    @property
    def VRd(self):
        if abs(self.Ashear[1]) < 1e-6 or self.fy < 1e-3:
            print("STEEL SECTION, ASHEAR: ", self.Ashear[1])
        return self.code.shear_resistance(self.Ashear[1], self.fy)

    @property
    def VyRd(self):
        if abs(self.Ashear[0]) < 1e-6 or self.fy < 1e-3:
            print("STEEL SECTION, ASHEAR: ", self.Ashear[0])
        return self.code.shear_resistance(self.Ashear[0], self.fy)

    @property
    def VzRd(self):
        if abs(self.Ashear[1]) < 1e-6 or self.fy < 1e-3:
            print("STEEL SECTION, ASHEAR: ", self.Ashear[1])
        return self.code.shear_resistance(self.Ashear[1], self.fy)

    @property
    def MRd(self):
        M, C = self.bending_resistance()
        return M

        """
        if self.C < 3:
            return self.code.bending_resistance(self.Wpl[0], self.fy)
        elif self.C == 3:
            return self.code.bending_resistance(self.Wel[0], self.fy)
        else:

            #  print("USING ELASTIC PROPERTIES FOR CROSS-SECTION CLASS 4")
            return self.code.bending_resistance(self.Wel[0], self.fy)

            #raise NotImplemented("Calculation of cross-section class 4 is not implemented yet")
        """

    @property
    def MyRd(self):
        M, C = self.bending_resistance(axis="y")
        return M

        """
        if self.C < 3:
            return self.code.bending_resistance(self.Wpl[0], self.fy)
        elif self.C == 3:
            return self.code.bending_resistance(self.Wel[0], self.fy)
        else:

            #  print("USING ELASTIC PROPERTIES FOR CROSS-SECTION CLASS 4")
            return self.code.bending_resistance(self.Wel[0], self.fy)
        """

    @property
    def MzRd(self):
        M, C = self.bending_resistance(axis="z")
        return M

        """
        if self.C < 3:
            return self.code.bending_resistance(self.Wpl[1], self.fy)
        elif self.C == 3:
            return self.code.bending_resistance(self.Wel[1], self.fy)
        else:

            #  print("USING ELASTIC PROPERTIES FOR CROSS-SECTION CLASS 4")
            return self.code.bending_resistance(self.Wel[1], self.fy)
        """

    @property
    def C(self):
        C = self.section_class()
        return C

    @property
    def epsilon(self):
        e = self.code.epsilon(self.fy)
        return e

    @abstractclassmethod
    def flange_class(self):
        """ Determines flange's class in compression

            Returns:
            ---------
            :return C: flange's cross-section class in compression
            :rtype C: int
        """

    @abstractclassmethod
    def web_class_comp(self):
        """ Determines web's class in compression

            Returns:
            ---------
            :return C: web's cross-section class
            :rtype C: int

        """

    @abstractclassmethod
    def web_class_bend(self):
        """ Setermines web's class in bending

        Returns:
        ---------
        :return C: web's cross-section class in bending
        :rtype C: int
        """

    @abstractclassmethod
    def web_class_comp_bend(self, Ned):
        """ Determines web's class in combined bending and compression

            Parameters:
            -----------
            :param Ned: compressive force [kN]
            :type Ned: float

            Returns:
            --------
            :return C: web's cross-section class in combined bending and compression
            :rtype C: int

        """

    def robot(self):
        # returns the name of the section that Robot Structural Analysis can identify.
        pass

    def abaqus(self):
        # writes cross-section data to a file that abaqus can read
        pass

    def weight(self):
        """ Weight per unit length kg/mm """
        w = self.A * constants.density
        return w

    def self_weight(self):
        """ self-weight of the section (N/mm) """
        return self.weight() * 9.81

    def stress_ratio(self):
        """ computes the stress ratio Psi required
        in cross-section classification of parts in
        compression and bending
        See. EN 1993-1-1, Table 5.2
        NOTE: compression is positive, thus
        the '-' sign in front of Ned.
        """

        stmax = -self.Ned / self.A - self.Med / self.Wel[0]
        stmin = -self.Ned / self.A + self.Med / self.Wel[1]

        p = stmin / stmax

        return p

    def sigma_com(self):
        """ Elastic compressive stress 
            This works for double symmetric cross-sections
        """
        sN = -self.Ned / self.A
        sMy = self.Med[0] / self.Wel[0]

        return sN + sMy

    def sigmaN(self):
        """ Axial stress from normal force """

        return self.Ned / self.A

    def sigmaM(self):
        """ Axial stress from bending moment """
        return self.Med[0] / self.Wel[0]

    def section_class(self, verb=False):
        """ Determines cross-section class """

        C = 1
        Med = self.Med[0]
        
        if self.Ned >= 0.0 and abs(Med) < 1e-4:
            # No bending and zero axial force or tension
            if verb:
                print("No axial force or tension and no bending")
            C = 1
        
        # Pure bending        
        elif abs(self.Ned) < 1e-4 and abs(Med) > 1e-4:
            if verb:
                print("Pure bending")
            Cflange = self.flange_class(verb)
            Cweb = self.web_class_bend(verb)
            C = max(Cweb, Cflange)

        # Compression
        elif self.Ned < 0.0 or Med == 0:
            # Pure compression
            if abs(Med) < 1e-4:
                if verb:
                    print("Pure compression")
                Cflange = self.flange_class(verb)
                Cweb = self.web_class_comp(verb)
                C = max(Cweb, Cflange)
            # Bending and compression
            else:
                if verb:
                    print("Bending and compression")
                # Flange is in compression
                Cflange = self.flange_class(verb)
                # Classify web as a part in compression and bending
                Cweb = self.web_class_comp_bend(self.Ned, verb)
                C = max(Cweb, Cflange)

        # Axial tension and bending
        elif self.Ned > 0.0 and abs(Med) > 1e-4:
            C = 1

        if verb:
            print("Section class = {0}".format(C))

        return C

    def axial_force_resistance(self, verb=False):
        if self.Ned >= 0:
            NRd = self.code.tension_resistance(self.A, self.fy)
        else:
            NRd = self.code.compression_resistance(self.A, self.fy)

        if verb:
            print("NRd = {0:4.2f} kN".format(NRd * 1e-3))

        return NRd

    def bending_resistance(self, C=0, axis="y", verb=False):
        # Bending resistance, Nmm
        if C == 0:
            C = self.section_class()

        if axis == "y":
            n = 0
        else:
            n = 1

        if C < 3:
            WRd = self.Wpl[n]
        elif C == 3:
            WRd = self.Wel[n]
        else:
            raise NotImplementedError(f"Cross-section class 4 not implemented in: {self}")
            # WRd = self.Wel[n]

        MRd = self.code.bending_resistance(WRd, self.fy)

        if verb:
            print("MRd = {0:4.2f} kNm".format(MRd * 1e-6))

        return MRd, C

    def elastic_bending_resistance(self, axis="y"):
        """ Elastic bending resistance """
        if axis == "y":
            n = 0
        else:
            n = 1
        return self.code.bending_resistance(self.Wel[n], self.fy)

    def plastic_bending_resistance(self, axis="y"):
        """ Plastic bending resistance """
        if axis == "y":
            n = 0
        else:
            n = 1
        return self.code.bending_resistance(self.Wpl[n], self.fy)

    def shear_force_resistance(self, C=0, axis="z", verb=False):
        if C == 0:
            C = self.section_class()

        if axis == "y":
            n = 0
        else:
            n = 1

        if C < 3:
            VRd = self.code.shear_resistance(self.Ashear[n], self.fy)
        elif C == 3:
            VRd = self.code.shear_resistance(self.Ashear[n], self.fy)
        else:
            VRd = self.code.shear_resistance(self.Ashear[n], self.fy)

        if verb:
            print("Shear area = {0:4.2f} mm2".format(self.Ashear[n]))
            print("VRd = {0:4.2f} kN".format(VRd * 1e-3))

        return VRd

    def moment_axial_force_interact(obj, UN, MRd):
        return None

        """
        def a = increase_M_and_N(self):
            "" Increases bending moment and axial force proportionally,
            until cross-section resistance is obtained. Used for classification
            
            Needs to be modified for Python!
            ""
            
            a0 = 1.0
            
            #%a = fzero(@MNInc,a0)
            aIter = []
            amax = 5
            while isempty(aIter)
                amax = 2*amax
                [f,aIter] = secant(@MNInc,[1.0,amax])
                   
            a = aIter(end)                        
        
        
        def res = MNInc(k):
            selfN.Med = k*self.Med
            selfN.Ned = k*self.Ned
                
            MNRd = selfN.MomentAxialForceInteract
            if MNRd > 0:
                res = selfN.Med/MNRd-1
            else
                res = 1
        """

    def section_resistance(self, axis='y', return_list=True, verb=False):
        """ Calculates resistance of cross-section
            Checks the following:
                Axial force
                Shear force
                Bending moment
                Interaction of axial force and bending moment
            output: r .. maximum of different utilization ratios
        """
        if axis == 'y':
            idx = 0
            MRd = self.MyRd
            Med = self.Med[0]
            Ved = self.Ved[1]
            VRd = self.VzRd
        else:
            idx = 1
            Med = self.Med[1]
            Ved = self.Ved[0]
            MRd = self.MzRd
            VRd = self.VyRd

        # Normal force
        UN = abs(self.Ned) / self.NRd
        # Shear force
        UV = abs(Ved) / VRd
        # Bending moment
        UM = abs(Med) / MRd

        if self.C < 3:
            MNRd = self.moment_axial_force_interact(UN, MRd)
            if MNRd > 0.0:
                UMN = abs(Med) / MNRd
            else:
                UMN = UN + UM
                # UMN = INFEASIBLE
        else:
            UMN = UN + UM

        if return_list:
            r = [UN, UV, UM, UMN]
        else:
            r = max([UN, UV, UM, UMN])

        if verb:
            print("Cross-section design: " + self.__repr__() + " S" + str(self.fy))
            print("NEd = {0:4.2f} kN; NRd = {1:4.2f} kN => UN = {2:4.2f}".format(self.Ned * 1e-3, self.NRd * 1e-3, UN))
            print("VEd = {0:4.2f} kN; VRd = {1:4.2f} kN => UV = {2:4.2f}".format(Ved * 1e-3, VRd * 1e-3, UV))
            print("MEd = {0:4.2f} kNm; MRd = {1:4.2f} kNm => UM = {2:4.2f}".format(Med * 1e-6, MRd * 1e-6, UM))

        return r

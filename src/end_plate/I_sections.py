# Imported libraries
from tables_and_tuples import *
from math import sqrt, pi
import eurocode3 as ec3


class Beam:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, size="", material="S355", L=1000,
                 h=1.0e-6, b=1.0e-6, t_f=1.0e-6, t_w=1.0e-6, r=1.0e-6,
                 N=1.0e-6, Vy=1.0e-6, Vz=1.0e-6, Mt=1.0e-6, My=1.0e-6, Mz=1.0e-6):

        self.material = material
        self.L = L
        self.stress = []
        self.c_t_part = []
        self.PL = 0

        if size != "Welded" and size != "":
            # Standard size beam
            self.size = size

            self.h = float(profile[size]["h"])
            self.b = float(profile[size]["b"])
            self.t_f = float(profile[size]["t_f"])
            self.t_w = float(profile[size]["t_w"])
            self.r = float(profile[size]["r"])
        else:
            if size == "Welded":
                # Welded beam
                self.size = "Welded"
            else:
                # Size not defined
                self.size = ""

            self.h = float(h)
            self.b = float(b)
            self.t_f = float(t_f)
            self.t_w = float(t_w)
            self.r = float(r)

        if self.material not in mat:
            print("Defined material not in material table!\n"
                  "Material S355 used instead!")
            self.material = "S355"

        # Material properties
        self.f_y = float(mat[self.material]["f_y"])             # Yield strength
        self.f_u = float(mat[self.material]["f_u"])             # Ultimate strength
        self.E = float(mat[self.material]["E"])                 # Young's modulus
        self.v = float(mat[self.material]["v"])                 # Poisson's ratio
        self.rho = float(mat[self.material]["rho"])             # Density of material

        # Maximum web axial stress
        self.sigma_com_Ed = 0.0

        self.cross_section_properties()
        self.loading(N, Vy, Vz, Mt, My, Mz)
        self.stress_points()

    def cross_section_properties(self):
        # Cross-section properties
        if self.size == "Welded":
            print(NotImplemented + "Calculation of cross-section properties for welded I-beam.")
        else:
            # Calculation according to http://sections.arcelormittal.com/fileadmin/redaction/4-Library/
            # 1-Sales_programme_Brochures/Sales_programme/Sections_MB-ArcelorMittal_FR_EN_DE-V2017-3.pdf
            # page 193 ->

            # Area of section [mm^2]
            self.A = 2.0*self.t_f*self.b + (self.h - 2.0*self.t_f)*self.t_w + (4.0 - pi)*self.r**2.0

            # Length of straight portion of outstand part of flange [mm]
            self.d_f = 0.5*(self.b - self.t_w) - 0.8*self.r
            # Length of straight portion of web [mm]
            self.d_w = self.h - 2.0*(self.t_f + self.r)

            # Perimeter of cross-section
            self.P = 2.0*self.b + 4.0*self.d_f + 2.0*self.d_w + 2.0*pi*self.r

            # Shear area [mm^2]
            self.A_v = self.A - 2.0*self.b*self.t_f + (self.t_w + 2.0*self.r)*self.t_f

            # Second moment of area, I_y is stronger direction [mm^4]
            self.I_y = 1.0/12.0*(self.b*self.h**3.0 - (self.b - self.t_w)*(self.h - 2.0*self.t_f)**3.0) + \
                       0.03*self.r**4 + 0.2146*self.r**2.0*(self.h - 2.0*self.t_f - 0.4468*self.r)**2.0
            self.I_z = 1.0/12.0*(2.0*self.t_f*self.b**3.0 + (self.h - 2.0*self.t_f)*self.t_w**3.0) + \
                       0.03*self.r**4.0 + 0.2146*self.r**2.0*(self.t_w + 0.4468*self.r)**2.0

            # Torsional constant [mm^4]
            self.I_t = 2.0/3.0*(self.b - 0.63*self.t_f)*self.t_f**3.0 + \
                       1.0/3.0*(self.h - 2.0*self.t_f)*self.t_w**3.0 + \
                       2.0*(self.t_w/self.t_f)*(0.145 + 0.1*self.r/self.t_f)*\
                       (((self.r + self.t_w/2.0)**2.0 + (self.r + self.t_f)**2.0 - self.r**2.0)/(2.0*self.r + self.t_f))

            # Warping constant [mm^4]
            self.I_w = (self.t_f*self.b**3.0)/24.0*(self.h - self.t_f)**2.0

            # Plastic section modulus [mm^3]
            self.W_pl_y = (self.t_w*self.h**2.0)/4.0 + (self.b - self.t_w)*(self.h - self.t_f)*self.t_f + \
                          (4.0 - pi)/2.0*self.r**2.0*(self.h - 2.0*self.t_f) + \
                          (3.0*pi - 10.0)/3.0*self.r**3.0
            self.W_pl_z = (self.b**2.0*self.t_f)/2.0 + \
                          (self.h - 2.0*self.t_f)/4.0*self.t_w**2.0 + self.r**3.0*(10.0/3.0 - pi) + \
                          (2.0 - pi/2.0)*self.t_w*self.r**2.0

        # Mass per unit length [kg/mm]
        self.g = self.A*self.rho

        # Variable epsilon used in design
        self.eps = sqrt(235.0/self.f_y)

        # Elastic section modulus [mm^3]
        self.W_el_y = self.I_y/(0.5*self.h)
        self.W_el_z = self.I_z/(0.5*self.b)

    def loading(self, N=1.0e-6, Vy=1.0e-6, Vz=1.0e-6, Mt=1.0e-6, My=1.0e-6, Mz=1.0e-6):
        self.N = float(N)
        self.Vy = float(Vy)
        self.Vz = float(Vz)
        self.Mt = float(Mt)
        self.My = float(My)
        self.Mz = float(Mz)

        # Normal stress on stress points
        self.stress_points()
        # Cross-section class
        self.cross_section_class()

        if abs(Mt) > 1.0e-6:
            print("Effects of torsion to capacity not implemented!")

        # Calculating capacity of cross-section
        self.capacity()

    def cross_section_class(self):
        self.c_t_part = []
        if self.size == "Welded":
            print(self.NOT_IMPLEMENTED + "Cross-section classification of welded cross-sections")
            self.PL = 0
        elif len(self.size) != 0:
            # Web cross-section class, I-profile
            self.c_t_part.append(c_t_class(self, c=self.d_w, t=self.t_w, c_t_type="internal",
                                           c_t_stress1=self.stress[0].sigma_x, c_t_stress2=self.stress[1].sigma_x))
            # Upper flange cross-section classes
            self.c_t_part.append(c_t_class(self, c=self.d_f, t=self.t_f, c_t_type="outstand", point1='out', point2='in',
                                           c_t_stress1=self.stress[2].sigma_x, c_t_stress2=self.stress[3].sigma_x))
            self.c_t_part.append(c_t_class(self, c=self.d_f, t=self.t_f, c_t_type="outstand", point1='in', point2='out',
                                           c_t_stress1=self.stress[5].sigma_x, c_t_stress2=self.stress[6].sigma_x))
            # Lower flange cross-section classes
            self.c_t_part.append(c_t_class(self, c=self.d_f, t=self.t_f, c_t_type="outstand", point1='out', point2='in',
                                           c_t_stress1=self.stress[7].sigma_x, c_t_stress2=self.stress[8].sigma_x))
            self.c_t_part.append(c_t_class(self, c=self.d_f, t=self.t_f, c_t_type="outstand", point1='in', point2='out',
                                           c_t_stress1=self.stress[10].sigma_x, c_t_stress2=self.stress[11].sigma_x))
            try:
                self.PL = max(c_t.PL for c_t in self.c_t_part if c_t.PL != 0)
            except ValueError:
                self.PL = 0
        else:
            self.PL = 0

    def capacity(self):
        # Capacities
        self.N_Rd = self.A*self.f_y/gamma_M[0]              # Axial force capacity
        self.V_y_Rd = self.A_v*self.f_y/gamma_M[0]          # Shear capacity

        if self.PL == 1 or self.PL == 2:
            self.M_Rd_y = self.f_y*self.W_pl_y/gamma_M[0]           # Plastic moment capacity, stronger direction
            self.M_Rd_z = self.f_y*self.W_pl_z/gamma_M[0]           # Plastic moment capacity, weaker direction
        else:
            self.M_Rd_y = self.f_y*self.W_el_y/gamma_M[0]       # Elastic moment capacity, stronger direction
            self.M_Rd_z = self.f_y*self.W_el_z/gamma_M[0]       # Elastic moment capacity, weaker direction
            if self.PL == 4:
                print('Cross-section class PL = 4!')
                print('Calculation of effective cross-section values for PL 4 not implemented!!!')

        if abs(self.Vy) > 0.5*self.V_y_Rd and len(self.size) != 0:
            print("V_y_Rd = " + str(self.V_y_Rd))
            print("Vy = " + str(self.Vy))
            print("Shear force over 50% of plastic shear capacity! Moment capacity reduced!")

            rho_V = (2.0*self.Vy/self.V_y_Rd - 1.0)**2.0            # Reduction factor due to shear

            if self.PL == 1 or self.PL == 2:
                self.M_Rd_y = rho_V*self.f_y*self.W_pl_y/gamma_M[0]  # Plastic moment capacity, stronger direction
                self.M_Rd_z = rho_V*self.f_y*self.W_pl_z/gamma_M[0]  # Plastic moment capacity, weaker direction
            else:
                self.M_Rd_y = rho_V*self.f_y*self.W_el_y/gamma_M[0]  # Elastic moment capacity, stronger direction
                self.M_Rd_z = rho_V*self.f_y*self.W_el_z/gamma_M[0]  # Elastic moment capacity, weaker direction

    def stress_points(self):
        self.stress = []
        # Stress points calculation [MPa]
        # Web stress points
        self.stress.append(stress_point(self, z=-0.5*self.d_w, y=0.0))                    # Stress point 1
        self.stress.append(stress_point(self, z=0.5*self.d_w, y=0.0))                     # Stress point 2
        # Upper flange stress points
        self.stress.append(stress_point(self, z=-0.5*self.h, y=-0.5*self.b))              # Stress point 3
        self.stress.append(stress_point(self, z=-0.5*self.h, y=-0.5*self.b + self.d_f))   # Stress point 4
        self.stress.append(stress_point(self, z=-0.5*self.h, y=0.0))                      # Stress point 5
        self.stress.append(stress_point(self, z=-0.5*self.h, y=0.5*self.b - self.d_f))    # Stress point 6
        self.stress.append(stress_point(self, z=-0.5*self.h, y=0.5*self.b))               # Stress point 7
        # Lower flange stress points
        self.stress.append(stress_point(self, z=0.5*self.h, y=-0.5*self.b))               # Stress point 8
        self.stress.append(stress_point(self, z=0.5*self.h, y=-0.5*self.b + self.d_f))    # Stress point 9
        self.stress.append(stress_point(self, z=0.5*self.h, y=0.0))                       # Stress point 10
        self.stress.append(stress_point(self, z=0.5*self.h, y=0.5*self.b - self.d_f))     # Stress point 11
        self.stress.append(stress_point(self, z=0.5*self.h, y=0.5*self.b))                # Stress point 12

        # Maximum web axial stress [MPa]
        if abs(self.stress[0].sigma_x) > abs(self.stress[1].sigma_x):
            self.sigma_com_Ed = abs(self.stress[0].sigma_x)
        else:
            self.sigma_com_Ed = abs(self.stress[1].sigma_x)

    def info(self, fform):
        print(self.size + "\n---------------------------------------")
        print("h = {0} [mm]\nb = {1} [mm]\nt_f = {2} [mm]\nt_w = {3} [mm]\nr = {4} [mm]\n"
              .format(self.h,self.b,self.t_f,self.t_w,self.r))

        print("Material: " + self.material)
        print("f_y = {0} [MPa], f_u = {1} [MPa]\n".format(self.f_y,self.f_u))

        print("Cross-section properties")
        print("A = {0} [mm^2]\nA_v = {1} [mm^2]\nd_w = {2} [mm]".format(self.A, self.A_v, self.d_w))
        print("W_el_y = {0:{ff}} mm^3, W_el_z = {1:{ff}} mm^3".format(self.W_el_y,self.W_el_z, ff=fform))
        print("W_pl_y = {0:{ff}} mm^3, W_pl_z = {1:{ff}} mm^3".format(self.W_pl_y,self.W_pl_z, ff=fform))

        print('\nCross-section class: ' + str(self.PL))
        print('\nLoads:')
        print('   N = ' + str(round(self.N, 5)))
        print('   Vy = ' + str(round(self.Vy, 5)))
        print('   Vz = ' + str(round(self.Vz, 5)))
        print('   Mt = ' + str(round(self.Mt, 5)))
        print('   My = ' + str(round(self.My, 5)))
        print('   Mz = ' + str(round(self.Mz, 5)))

        n = 0
        print("")
        for i in self.stress:
            n = n + 1
            print("Stress point {0:2d}: sigma_x = {1:7{fm}} [MPa]   coords: (z, y) = ({2:8{fm}}, {3:8{fm}}) [mm]"
                  .format(n, i.sigma_x, i.z, i.y, fm='.16f'))#fform))


class stress_point:

    def __init__(self, beam, z=0.0, y=0.0):
        self.z = z              # [mm]
        self.y = y              # [mm]

        # Loads have been allocated using value 1.0e-6 to avoid dividing with zero,
        # Rounded up for stress point calculations
        N = round(beam.N, 5)              # [kN]
        Vy = round(beam.Vy, 5)            # [kN]
        Vz = round(beam.Vz, 5)            # [kN]
        Mt = round(beam.Mt, 5)            # [kNm]
        My = round(beam.My, 5)            # [kNm]
        Mz = round(beam.Mz, 5)            # [kNm]

        A = beam.A              # [mm^2]
        I_y = beam.I_y          # [mm^4]
        I_z = beam.I_z          # [mm^4]

        # Calculation of axial stress on arbitrary point in the cross-section [MPa]
        self.sigma_x = (N*1.0e3/A + (My*1.0e6/I_y)*self.z + (Mz*1.0e6/I_z)*self.y)


class c_t_class:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, beam, c, t, c_t_type, c_t_stress1, c_t_stress2, point1='', point2=''):
        self.c = float(c)
        self.t = float(t)
        c_t = c/t
        self.c_t = c/t
        self.c_t_type = c_t_type
        self.point1 = point1
        self.point2 = point2
        self.PL = 0

        if c_t_stress1 <= c_t_stress2:         # Compression defined to be negative in main program
            self.c_t_stress1 = c_t_stress1      # [MPa]
            self.c_t_stress2 = c_t_stress2      # [MPa]
        else:
            self.c_t_stress1 = c_t_stress2      # [MPa]
            self.c_t_stress2 = c_t_stress1      # [MPa]
            self.point1 = point2
            self.point2 = point1

        c_t_stress1 = -1.0*self.c_t_stress1          # Compression positive in SFS EN 1993-1-1 table 5.2
        c_t_stress2 = -1.0*self.c_t_stress2          # Compression positive in SFS EN 1993-1-1 table 5.2

        # Plastic stress ratio
        # alfa = NOT DEFINED
        # Elastic stress ratio
        # psi = NOT DEFINED

        if c_t_type == "internal":
            # Internal part, SFS EN 1993-1-1 table 5.2 (Part 1)
            if c_t_stress1 == -c_t_stress2:
                self.loading = "Pure bending"
                self.PL = ec3.internal_part_in_bending(c/t, beam.eps)
            elif c_t_stress1 == c_t_stress2:
                if c_t_stress1 > 0.0:                      # Compression positive in SFS EN 1993-1-1 table 5.2
                    self.loading = "Pure compression"
                    self.PL = ec3.internal_part_in_compression(c/t, beam.eps)
                else:
                    self.loading = "Pure tension"
                    self.PL = 1
            else:
                if c_t_stress1 <= 0.0 and c_t_stress2 <= 0.0:   # Compression positive in SFS EN 1993-1-1 table 5.2
                    self.loading = 'c/t part in tension'
                    self.PL = 1
                else:
                    self.loading = "Axial force and bending"
                    # Plastic and elastic ratios not defined! Used compression limits, which are conservative
                    # self.PL = ec3.internal_part_comp_bend(c/t, beam.eps, alfa, psi)
                    self.PL = ec3.internal_part_in_compression(c / t, beam.eps)

        elif c_t_type == "outstand":
            # External part, SFS EN 1993-1-1 table 5.2 (Part 2)
            if c_t_stress1 == c_t_stress2:
                if c_t_stress1 > 0.0:  # Compression positive in SFS EN 1993-1-1 table 5.2
                    self.loading = "Pure compression"
                    self.PL = ec3.outstand_part_in_compression(c/t, beam.eps)
                else:
                    self.loading = "Pure tension"
                    self.PL = 1
            else:
                if c_t_stress1 <= 0.0 and c_t_stress2 <= 0.0:  # Compression positive in SFS EN 1993-1-1 table 5.2
                    self.loading = 'c/t part in tension'
                    self.PL = 1
                else:
                    self.loading = "Axial force and bending"
                    # Plastic and elastic ratios not defined! Used compression limits, which are conservative
                    self.PL = ec3.outstand_part_in_compression(c/t, beam.eps)
        else:
            print("Type of c_t part must be (internal) or (outstand)")
            print("Given c_t part type is: " + c_t_type)


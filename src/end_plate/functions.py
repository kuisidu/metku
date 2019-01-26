# Contains functions needed to calculate end-plate connection capacity and stiffness

# Imported libraries
from math import sqrt, pi, ceil
import itertools
from tables_and_tuples import *

NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

# -----------------------------------------------------------------------------------------------
# Functions for calculation of basic joint components of end-plate connection

# --------------------------------------------------------------------------------------------
# Reduction factor, table 6.3
def min_coeff_w(beta, b_eff_c_wc, t_wc, A_vc):
    w = 0.0

    w_1 = 1.0/(sqrt(1.0 + 1.3*(b_eff_c_wc*t_wc/A_vc)**2.0))
    w_2 = 1.0/(sqrt(1.0 + 5.2*(b_eff_c_wc*t_wc/A_vc)**2.0))

    if 0.0 <= beta <= 0.5:
        w = 1.0
    elif 0.5 < beta < 1.0:
        w = w_1 + 2 * (1.0 - beta) * (1.0 - w_1)
    elif beta == 1.0:
        w = w_1
    elif 1.0 < beta < 2.0:
        w = w_1 + 2.0 * (beta - 1.0) * (1.0 - w_1)
    elif beta == 2.0:
        w = w_2
    else:
        print("Error: Reduction factor beta = {0} > 2".format(beta))

    return w


# Function defining alfa parameter figure (for l_eff of endplate)
def lambda11(lambda2, alfa):
    lambda1_lim = round(1.25/(alfa - 2.75), 12)
    lambda2_lim = round((alfa*lambda1_lim)/2.0, 12)
    if lambda2 >= lambda2_lim:
        return lambda1_lim
    else:
        return lambda1_lim + (1.0 - lambda1_lim)*((lambda2_lim - lambda2)/lambda2_lim)**(alfa/sqrt(2.0))


# Finding value for parameter alfa, SFS EN 1993-1-8, figure 6.11
def par_alfa(lambda1, lambda2):
    alfa = 6.5          # Initial guess for alfa
    diff = (lambda11(lambda2, alfa) - lambda1)

    max_iter = 1000
    iters = 0
    while abs(diff) > 1.0e-3:
        iters += 1

        diff = (lambda11(lambda2, alfa) - lambda1)
        alfa = alfa + diff

        if iters == max_iter:
            print("Determination of alfa (for end-plate bending) did not converge in " + str(iters) +
                  " iterations! Check results!")
            print("lambda1 = {0:{fm}}, lambda2 = {1:{fm}}".format(lambda1, lambda2, fm='.2f'))
            print("Used value of alfa = " + str(alfa))
            break

    return alfa


def define_welds(sd, beam, equ=True):
    plate_material = sd.plate_material
    if plate_material not in mat:
        print("Defined plate material not in material table!\n"
              "Material S355 used instead!")
        plate_material = "S355"

    if beam.f_y > mat[plate_material]["f_y"]:
        material = plate_material
        f_y = mat[plate_material]["f_y"]
        f_u = mat[plate_material]["f_u"]
    else:
        material = beam.material
        f_y = beam.f_y
        f_u = beam.f_u

    # Correlation factor
    if material == "S235":
        beta_w = 0.8
    elif material == "S275":
        beta_w = 0.85
    elif material == "S355":
        beta_w = 0.9
    elif material == "S420":
        beta_w = 1.0
    elif material == "S460":
        beta_w = 1.0
    else:
        print("Correlation factor beta_w not defined for material " + str(beam.material) + "!" +
              "\nCorrelation factor beta_w = 1.0 used")
        beta_w = 1.0

    t_f = beam.t_f
    t_w = beam.t_w

    # Effective lengths of welds
    if beam.size == "Welded":
        l_eff_f = 2.0*beam.b - (2.0*sqrt(2.0)*beam.r + t_w)   # Effective length of flange weld
    else:
        l_eff_f = 2.0*beam.b - (2.0*beam.r + t_w)               # Effective length of flange weld
    l_eff_w = 2.0*beam.d_w                                        # Effective length of web weld

    # Equivalent strength welds of flange and web according to Ruukki-Hitsatut-Profiilit-Kasikirja-2010, page 354
    # Calculation is based on simplified method, where all welds are taken to be sqrt(3) welds
    a_flange_equ = ((sqrt(3.0)*beta_w*gamma_M[2]*f_y)/(2.0*gamma_M[0]*f_u))*t_f
    a_web_equ = ((sqrt(3.0)*beta_w*gamma_M[2]*f_y)/(2.0*gamma_M[0]*f_u))*t_w

    sd.weld_f = float(ceil(a_flange_equ))
    sd.weld_w = float(ceil(a_web_equ))
    sd.welds = "End-plate welds defined to be equivalent strength with beam"

    if (not equ) and (abs(beam.N) > 1.0e-6 or abs(beam.Vy) > 1.0e-6 or abs(beam.My) > 1.0e-6):

        if abs(beam.Vz) > 1.0e-6 or abs(beam.Mt) > 1.0e-6 or abs(beam.Mz) > 1.0e-6:
            print('Effects of shear force and moment in weaker direction or torsion not taken into account!')

        # Flange and web weld loads from connection loads
        # My [kNm], N [kN], Vy [kN], h [mm], t_f [mm]
        F_weld_f = abs(beam.My)/((beam.h - t_f)*1.0e-3) + abs(beam.N)       # Flange weld load [kN]
        F_weld_w = abs(beam.Vy)                                             # Web weld load [kN]

        # Flange and web weld sizes according to connection loads
        a_flange_load = max((2.0*(F_weld_f*1.0e3)*beta_w*gamma_M[2])/(sqrt(2.0)*l_eff_f*f_u),
                            ((F_weld_f*1.0e3)*gamma_M[2])/(sqrt(2.0)*l_eff_f*0.9*f_u))
        a_web_load = (sqrt(3.0)*(F_weld_w*1.0e3)*beta_w*gamma_M[2])/(l_eff_w*f_u)

        n_extra = 1.2       # Coefficient to make sure welds do not define connection capacity
        sd.weld_f = float(ceil(min(max(n_extra*a_flange_load, 4.0), a_flange_equ)))
        sd.weld_w = float(ceil(min(max(n_extra*a_web_load, 4.0), a_web_equ)))
        sd.welds = "End-plate welds defined from beam loads"

# Design resistance F_t_Rd of a T-Stub flange, Table 6.2
def F_t_Rd_T_stub(column, sd, rg, component, A_s=0.0, sum_F_T_Rd=0.0):
    # Initializing values
    sum_l_eff_1 = 0.0
    sum_l_eff_2 = 0.0
    t_f = 0.0
    f_y = 0.0
    m = 0.0
    e_min = 0.0

    F_T_1_Rd = 0.0
    F_T_2_Rd = 0.0
    F_T_3_Rd = 0.0
    F_T_1_2_Rd = 0.0

    if component == "flange":
        sum_l_eff_1 = rg.l_eff_1f
        sum_l_eff_2 = rg.l_eff_2f
        t_f = column.t_f
        f_y = column.f_y
        m = rg.mc
        e_min = rg.e_min
    elif component == "plate":
        sum_l_eff_1 = rg.l_eff_1p
        sum_l_eff_2 = rg.l_eff_2p
        t_f = sd.plate_thickness
        f_y = mat[sd.plate_material]["f_y"]
        if type(rg).__name__ == "RowGroup":
            m = rg.m
            e_min = rg.e_min
        else:
            if rg.location_p == "Bolt-row outside tension flange of beam":
                m = rg.m_x
                e_min = rg.e_x
            else:
                m = rg.m
                e_min = rg.e_min

    d_w = rg.d_w

    t_bp = sd.t_bp
    f_y_bp = sd.f_y_bp

    e_w = d_w/4.0              # Used in alternative method
    n = min(e_min, 1.25*m)

    # Prying force check
    if sum_l_eff_1 == 0.0:
        L_b_star = 0.0
    else:
        L_b_star = (8.8e9*A_s)/(sum_l_eff_1*t_f**3.0)

    # Equivalent moment capacities of T_stub
    M_pl_1_Rd = 0.25*sum_l_eff_1*t_f**2.0*f_y/gamma_M[0]
    M_pl_2_Rd = 0.25*sum_l_eff_2*t_f**2.0*f_y/gamma_M[0]
    M_bp_Rd = 0.25*sum_l_eff_1*t_bp**2.0*f_y_bp/gamma_M[0]

    if rg.L_b <= L_b_star:
        rg.prying = "Prying forces may develop"

        # Mode 1
        if t_bp == 0.0:         # Without backing plate
            F_T_1_Rd = 4.0*M_pl_1_Rd/m
        else:                   # With backing plate
            F_T_1_Rd = (4.0*M_pl_1_Rd + 2.0*M_bp_Rd)/m

        # Mode 2
        F_T_2_Rd = (2.0*M_pl_2_Rd + n*sum_F_T_Rd)/(m + n)
    else:
        rg.prying = "No prying forces"
        F_T_1_2_Rd = 2.0*M_pl_1_Rd/m

    # Mode 3
    F_T_3_Rd = sum_F_T_Rd

    if component == "flange":
        rg.F_T_Rdf = [F_T_1_Rd, F_T_2_Rd, F_T_3_Rd, F_T_1_2_Rd]
    elif component == "plate":
        rg.F_T_Rdp = [F_T_1_Rd, F_T_2_Rd, F_T_3_Rd, F_T_1_2_Rd]


# Comment information
# Component, design resistance reference, stiffness coefficient reference

# Column web panel in shear, 6.2.6.1, 6.3.2
def cweb_shear(cn, column):
    """ Input:
        cn ... Connection object of class EndPlate
        column ... Profile of the column
    """
    d_wc = column.d_w                                     # Column web height
    t_wc = column.t_w                                     # Column web thickness
    f_y_wc = column.f_y                                   # Yield strength of column material
    A_vc = column.A_v                                     # Shear area of column web

    eps = sqrt(235.0/f_y_wc)                                     # Material reference value

    # Checking column web slenderness
    if (d_wc/t_wc) >= (69.0*eps):
        cn.warnigs.append("Column web too slender d/t_w <= 69*eps: " + str(d_wc/t_wc) + " <= " + str(69.0*eps))

    # Calculating column web panel shear capacity
    cn.V_wp_Rd = (0.9*f_y_wc*A_vc)/(sqrt(3.0)*gamma_M[0])


# Column web in transverse compression, 6.2.6.2, 6.3.2
def cweb_trv_comp(cn, column, sd, sigma_com_Ed=1.0):
    # Connection using screws

    # Column information
    r_c = column.r                                     # Column rounding, VALSSATTU
    t_fc = column.t_f                                  # Column flange thickness
    t_wc = column.t_w                                  # Column web thickness
    d_wc = column.d_w                                  # Column web height
    A_vc = column.A_v                                  # Shear area of column web
    f_y_wc = column.f_y                                # Yield strength of column material
    E = column.E                                       # Young's modulus of column material

    # Beam information
    t_fb = sd.beam.t_f                           # Beam flange thickness

    # Plate information
    t_p = sd.plate_thickness                     # Plate thickness
    lo_p = sd.lower_overhang                     # Lower overhang of plate
    a_p = sd.weld_f                              # Weld size of plate on flanges

    # Width of compressed zone under compressed flange using 45 deg stress distribution
    s_p = min(2.0*t_p, t_p+lo_p-sqrt(2.0)*a_p)           # lo_p is overhang of end-plate over lower flange of beam

    # Effective length of compressed web
    s = r_c
    b_eff_c_wc = t_fb+2.0*sqrt(2.0)*a_p+5.0*(t_fc+s)+s_p

    # Slenderness of column web
    lam_p = 0.932*sqrt((b_eff_c_wc*d_wc*f_y_wc)/(E*t_wc**2.0))

    # Reduction factor of plate buckling for web
    if lam_p <= 0.72:
        rho = 1.0
    else:
        rho = (lam_p-0.2)/(lam_p**2.0)

    # Reduction factor
    w = min_coeff_w(sd.beta, b_eff_c_wc, t_wc, A_vc)

    # Checking maximum longitudinal compressive stress duo to column axial force and bending moment, 6.2.6.2 (2)
    # Maximum axial stress on column web
    if sigma_com_Ed <= 0.7*f_y_wc:
        k_wc = 1.0
    else:
        k_wc = 1.7-sigma_com_Ed/f_y_wc

    # Calculating capacity of column web in transverse compression
    sd.F_c_wc_Rd = min((w*k_wc*b_eff_c_wc*t_wc*f_y_wc)/gamma_M[0],(w*k_wc*rho*b_eff_c_wc*t_wc*f_y_wc)/gamma_M[1])

    sd.b_eff_c_wc = b_eff_c_wc


# Column web in transverse tension, 6.2.6.3, 6.3.2
def cweb_trv_tens(column, sd, rg):# l_eff_1f, l_eff_2f):

    # Connection using screws
    # Column information
    t_wc = column.t_w  # Column web thickness
    A_vc = column.A_v  # Shear area of column web
    f_y_wc = column.f_y  # Yield strength of column material

    # Effective width of column web in tension
    b_eff_t_wc = min(rg.l_eff_1f, rg.l_eff_2f)

    # Reduction factor
    w = min_coeff_w(sd.beta, b_eff_t_wc, t_wc, A_vc)

    # Calculating design resistance of unstiffened column web subject to transverse tension
    rg.F_t_wc_Rd = (w*b_eff_t_wc*t_wc*f_y_wc)/gamma_M[0]


# Column flange in bending, 6.2.6.4, 6.3.2
def cflange_bending(column, sd, rg):
    if rg.thread_in_shear_plane:
        A_s = rg.A_s
    else:
        A_s = rg.A

    sum_F_T_Rd = rg.n_rows*rg.F_T_Rd

    # Capacity of T-stub
    F_t_Rd_T_stub(column, sd, rg, "flange", A_s, sum_F_T_Rd)

    try:
        rg.F_t_f_Rd = min(i for i in rg.F_T_Rdf if i != 0.0)
    except:
        rg.F_t_f_Rd = 0.0

    if rg.F_t_f_Rd == rg.F_T_Rdf[0]:
        rg.mode_f = "Mode 1, flange yielding"
    elif rg.F_t_f_Rd == rg.F_T_Rdf[1]:
        rg.mode_f = "Mode 2, flange yielding and bolt failure"
    elif rg.F_t_f_Rd == rg.F_T_Rdf[2]:
        rg.mode_f = "Mode 3, bolt failure"
    elif rg.F_t_f_Rd == rg.F_T_Rdf[3]:
        rg.mode_f = "Mode 1-2"
    else:
        rg.mode_f = "None"

# End-plate in bending, 6.2.6.5, 6.3.2
def end_plate_bending(column, sd, rg):
    if rg.thread_in_shear_plane:
        A_s = rg.A_s
    else:
        A_s = rg.A

    sum_F_T_Rd = rg.n_rows*rg.F_T_Rd

    # Capacity of T-stub
    F_t_Rd_T_stub(column, sd, rg, "plate", A_s, sum_F_T_Rd)

    try:
        rg.F_t_p_Rd = min(i for i in rg.F_T_Rdp if i != 0.0)
    except:
        rg.F_t_p_Rd = 0.0

    if rg.F_t_p_Rd == rg.F_T_Rdp[0]:
        rg.mode_p = "Mode 1, end-plate yielding"
    elif rg.F_t_p_Rd == rg.F_T_Rdp[1]:
        rg.mode_p = "Mode 2, end-plate yielding and bolt failure"
    elif rg.F_t_p_Rd == rg.F_T_Rdp[2]:
        rg.mode_p = "Mode 3, bolt failure"
    elif rg.F_t_p_Rd == rg.F_T_Rdp[3]:
        rg.mode_p = "Mode 1-2"
    else:
        rg.mode_p = "None"

# Beam flange and web in compression, 6.2.6.7, 6.3.2
# BEAM REINFORCEMENT WITH HAUNCHES NOT ADDED
def fc_comp(cn, sd):

    h = sd.beam.h            # Height of beam
    t_fb = sd.beam.t_f       # Thickness of beam flange

    if h > 600:
        cn.warnings.append(NOT_IMPLEMENTED + "Lisättävä ehto: jos h > 600 mm niin uuman vaikutus"
                                                     " puristuskestävyyteen saa olla enintään 20% (6.2.6.7 (1))")

    M_c_Rd = sd.beam.M_Rd_y      # Design moment resistance of beam cross-section, shear considered

    sd.F_c_fb_Rd = M_c_Rd/(h-t_fb)


# Beam web in tension, 6.2.6.8, 6.3.2
def bweb_tens(sd, rg):

    # Connection using screws

    # Beam information
    t_wb = sd.beam.t_w  # Beam flange thickness
    f_y_wb = sd.beam.f_y  # Yield strength of beam material

    # Effective width of column web in tension
    b_eff_t_wb = min(rg.l_eff_1p, rg.l_eff_2p)

    # Calculating design resistance of unstiffed column web subject to transverse tension
    rg.F_t_wb_Rd = (b_eff_t_wb * t_wb * f_y_wb) / gamma_M[0]


# Calculation of stiffness coefficient for each row, table 6.11
def stiffness_coeff(cn, sd, row):

    column = cn.column

    # Values
    t_fc = column.t_f                                  # Flange thickness column
    t_wc = column.t_w                                  # Web thickness of column
    d_c = column.d_w                                   # Height of straight portion of web in column
    t_p = sd.plate_thickness                              # End-plate thickness

    # Smallest effective length of row as part of row group
    # Flange
    try:
        l_eff_cpf_group = min(i for i in row.l_eff_cpf if i != 0.0)
    except ValueError:
        l_eff_cpf_group = 0.0
    try:
        l_eff_ncf_group = min(i for i in row.l_eff_ncf if i != 0.0)
    except ValueError:
        l_eff_ncf_group = 0.0
    # Plate
    try:
        l_eff_cpp_group = min(i for i in row.l_eff_cpp if i != 0.0)
    except ValueError:
        l_eff_cpp_group = 0.0
    try:
        l_eff_ncp_group = min(i for i in row.l_eff_ncp if i != 0.0)
    except ValueError:
        l_eff_ncp_group = 0.0

    if sd.stiffener == "unstiffened":
        # Column web in tension
        b_eff_t_wc = min(i for i in [row.l_eff_1f, row.l_eff_2f, l_eff_cpf_group, l_eff_ncf_group] if i != 0.0)
        row.k_3 = 0.7*b_eff_t_wc*t_wc/d_c
    else:
        row.k_3 = 1.0e12

    # Column flange in bending
    l_eff_f = min(i for i in [row.l_eff_1f, row.l_eff_2f, l_eff_cpf_group, l_eff_ncf_group] if i != 0.0)
    row.k_4 = (0.9*l_eff_f*t_fc**3.0)/(row.mc**3.0)

    # Endplate in bending
    l_eff_p = min(i for i in [row.l_eff_1p, row.l_eff_2p, l_eff_cpp_group, l_eff_ncp_group] if i != 0.0)
    if row.location_p == "Bolt-row outside tension flange of beam":
        row.k_5 = (0.9*l_eff_p*t_p**3.0)/(row.m_x**3.0)
    else:
        row.k_5 = (0.9*l_eff_p*t_p**3.0)/(row.m**3.0)

    # Flange cleat in bending
    # k_6 = 0.9*l_eff*t_a**3.0/m**3.0

    # Bolts in tension (for a single bolt-row)
    row.k_10 = 1.6*row.A_s/row.L_b

    # Bolts in shear
    # Bolts in bearing (for each component j on which the bolts bear)

    # Concrete in compression (including grout)
    # Plate in bending under compression
    # Base plate in bending under tension (for a single bolt row in tension)
    # Anchor bolts in tension

# --------------------------------------------------------------------------------

# Defining possible row groups
def groups(sd, rows):
    lst = list(itertools.product([0, 1], repeat=rows))

    i = 0
    while i < len(lst):
        # Taking into account shear rows
        j = 0
        while j < len(lst[i]):
            if lst[i][j] == 1 and sd.row[j].shear_row == True:
                lst[i] = list(lst[i])
                lst[i][j] = 0
                lst[i] = tuple(lst[i])
            j = j + 1

        # Removing groups with only 1 or less rows
        if sum(lst[i]) <= 1:
            lst.remove(lst[i])
            continue


        # Removing groups with blank rows between full rows
        n1 = lst[i].index(1)
        j = 0
        n0 = 0
        while j < len(lst[i]):
            if j > n1 and lst[i][j] == 0:
                n0 = 1
            elif j > n1 and n0 == 1 and lst[i][j] == 1:
                lst.remove(lst[i])
                if i > 0:
                    i = i - 1
                break
            j = j + 1

        # Removing duplicate row groups
        if lst.count(lst[i]) > 1:
            lst.remove(lst[i])
            continue

        i = i + 1

    return lst

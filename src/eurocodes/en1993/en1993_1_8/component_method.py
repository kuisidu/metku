# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 19:21:50 2019

Functions of the component method as described in EN 1993-1-8

@author: kmela
"""

from math import sqrt, pi

from eurocodes.en1993.constants import gammaM0, gammaM1

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

""" COMPONENTS """


def col_web_shear(column):
    """ Column web panel in shear, 6.2.6.1, 6.3.2
        Input: column .. steel_section class object
    """
    f_y_wc = column.fy                                   # Yield strength of column material
    A_vc = column.Ashear                                     # Shear area of column web

    # Calculating column web panel shear capacity
    return  0.9*f_y_wc*A_vc/sqrt(3.0)/gammaM0

def beam_web_compression(beam):
    """ Beam web in transverse compression:
        EN 1993-1-8, 6.2.6.7
        
        input:
            beam .. Beam profile as a SteelSection object
    """
    McRd = beam.MRd[0]
    return McRd/(beam.h-beam.tf)
    
def beam_web_tension(beam,beff):
    """ Beam web in tension:
        EN 1993-1-8, 6.2.6.8
        
        input:
            beam .. Beam profile as a SteelSection object
            beff .. Effective width of T-stub according to end plate in bending
                    (6.2.6.5)
    """
    return beff*beam.tw*beam.fy/gammaM0
    
def Tstub_compression_flange(fjd,beff,leff):
    """ Compression resistance of T-stub flange
        EN 1993-1-8, Eq. (6.4)
        input:
            fjd .. design bearing strength (6.2.5(7))
            beff .. effective width of T-stub flange
            leff .. effective length of T-stub flange
    """
    return fjd*beff*leff


def col_web_trv_comp(column, beam, t_p, ebottom, a_p, beta, sigma_com_Ed=0.0):
    """ Column web in transverse compression, 6.2.6.2, 6.3.2
        Input: 
            column .. SteelSection class object
            beam .. SteelSection class object
            t_p .. thickness of end plate
            ebottom .. lower overhang of the plate
            a_p .. Weld size between end plate and beam
            beta .. factor taking into account shear of the column web
            sigma_com_Ed .. compression stress in the column
    """
    # Connection using screws

    # Column information
    r_c = column.r                                     # Column rounding, VALSSATTU
    t_fc = column.tf                                  # Column flange thickness
    t_wc = column.tw                                  # Column web thickness
    d_wc = column.hw                                  # Column web height
    A_vc = column.Ashear                                  # Shear area of column web
    f_y_wc = column.fy                                # Yield strength of column material
    E = column.E                                       # Young's modulus of column material

    # Beam information
    t_fb = beam.tf                           # Beam flange thickness

    # Plate information
    lo_p = ebottom                  

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
    w = min_coeff_w(beta, b_eff_c_wc, t_wc, A_vc)

    # Checking maximum longitudinal compressive stress duo to column axial force and bending moment, 6.2.6.2 (2)
    # Maximum axial stress on column web
    if sigma_com_Ed <= 0.7*f_y_wc:
        k_wc = 1.0
    else:
        k_wc = 1.7-sigma_com_Ed/f_y_wc

    # Calculating capacity of column web in transverse compression
    F_c_wc_Rd = min((w*k_wc*b_eff_c_wc*t_wc*f_y_wc)/gammaM0,(w*k_wc*rho*b_eff_c_wc*t_wc*f_y_wc)/gammaM1)

    return F_c_wc_Rd, b_eff_c_wc


def column_web_tension(column, l_eff, beta):# l_eff_1f, l_eff_2f):
    """ Column web in transverse tension, 6.2.6.3, 6.3.2
        
    """
    # Connection using screws
    # Column information
    t_wc = column.tw  # Column web thickness
    A_vc = column.Ashear  # Shear area of column web
    f_y_wc = column.fy  # Yield strength of column material

    # Effective width of column web in tension
    b_eff_t_wc = min(l_eff)

    # Reduction factor
    w = min_coeff_w(beta, b_eff_t_wc, t_wc, A_vc)

    # Calculating design resistance of unstiffened column web subject to transverse tension
    F_t_wc_Rd = (w*b_eff_t_wc*t_wc*f_y_wc)/gammaM0
    
    return F_t_wc_Rd
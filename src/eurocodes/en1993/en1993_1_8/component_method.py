# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 19:21:50 2019

Functions of the component method as described in EN 1993-1-8

@author: kmela
"""

from math import sqrt, pi
from scipy.optimize import root_scalar

from eurocodes.en1993.constants import gammaM0, gammaM1
from eurocodes.en1993.en1993_1_1 import buckling_reduction_factor

# Reduction factor, table 6.3
def min_coeff_w(beta, b_eff_c_wc, t_wc, A_vc):
    w = 0.0

    #print("beff_c_wc = ",b_eff_c_wc)
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

def find_alpha(lambda1,lambda2):
    """ Determine alpha parameter of Fig. 6.11 by a root finding method """
    
    def a_zero(a):
        l1_lim = 1.25/(a-2.75)
        l2_lim = 0.5*a*l1_lim
        
        return l1_lim + (1-l1_lim)*((l2_lim-lambda2)/l2_lim)**(0.185*a**1.785) - lambda1
    
    r = root_scalar(a_zero,x0=6,bracket=[4.45,8])

    if r.converged:
        res = r.root
    else:
        res = 4.45
        print("Warning: alpha parameter iteration did not converge.")

    return res
    

def new_alpha(e,m,m2):
    """ Alpha parameter according to the final document of prEN 1993-1-8:2019 """
    
    return min(max(4+1.67*e/m*(m/m2)**0.67,4+1.25*e/m),8.0)

    
""" COMPONENTS """


def col_web_shear(column):
    """ Column web panel in shear, 6.2.6.1, 6.3.2
        Input: column .. steel_section class object
    """
    f_y_wc = column.fy    # Yield strength of column
    A_vc = column.Ashear  # Shear area of column web

    # Calculating column web panel shear capacity
    return  0.9*f_y_wc*A_vc/sqrt(3.0)/gammaM0

def col_web_shear_k(column,beta,z):
    """ Stiffness factor for Column web panel in shear
        Input: column .. steel_section class object
               beta .. transformation parameter
               z .. moment arm of the connection               
    """   
    return 0.38*column.Ashear/beta/z    

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

    """
    print("tbf = {0:4.2f}".format(t_fb))
    print("weld = {0:4.2f}".format(2*sqrt(2.0)*a_p))
    print("5.0*(t_fc+s) = {0:4.2f}".format(5.0*(t_fc+s)))
    print("s_p= {0:4.2f}".format(s_p))
    """
    # Slenderness of column web
    lam_p = 0.932*sqrt((b_eff_c_wc*d_wc*f_y_wc)/(E*t_wc**2.0))

    # Reduction factor of plate buckling for web
    if lam_p <= 0.72:
        rho = 1.0
    else:
        rho = (lam_p-0.2)/(lam_p**2.0)

    #print("rho = {0:4.3f}".format(rho))
    # Reduction factor
    w = min_coeff_w(beta, b_eff_c_wc, t_wc, A_vc)
        
    # Checking maximum longitudinal compressive stress duo to column axial force and bending moment, 6.2.6.2 (2)
    # Maximum axial stress on column web
    if sigma_com_Ed <= 0.7*f_y_wc:
        k_wc = 1.0
    else:
        k_wc = 1.7-sigma_com_Ed/f_y_wc

    """
    print("d_wc = ",d_wc)
    print("lam_p = ",lam_p)
    print("omega = ",w)
    print("rho = ",rho)
    """
    # Calculating capacity of column web in transverse compression
    F_c_wc_Rd = min((w*k_wc*b_eff_c_wc*t_wc*f_y_wc)/gammaM0,(w*k_wc*rho*b_eff_c_wc*t_wc*f_y_wc)/gammaM1)

    return F_c_wc_Rd, b_eff_c_wc

def col_web_trv_comp_k(column,beff):
    """ Stiffness of column web in compression """    
    
    return 0.7*beff*column.tw/column.hw

def col_web_trv_comp_stiffened(column, plate):
    """ Resistance of stiffened column web """
    e = column.eps
    tw = column.tw
    bsg = plate.b
    ts = plate.t
    
    fy = plate.material.fy
    As_eff = (30*e*tw+ts)*tw + 2*bsg*ts
    Is = (2*bsg + tw)**3*ts/12
    
    lcr = plate.h
    lambda_1 = 93.9*e
    i_s = sqrt(Is/As_eff)
    lambda_s = lcr/i_s/lambda_1
    
    if lambda_s > 0.2:
        chi = buckling_reduction_factor(lambda_s,0.49)
        gammaM = gammaM1
    else:
        chi = 1.0
        gammaM = gammaM0
    
    return chi*As_eff*fy/gammaM
        

def column_web_tension(column, l_eff, beta,verb=False):# l_eff_1f, l_eff_2f):
    """ Column web in transverse tension, 6.2.6.3, 6.3.2
        NOTE! l_eff is the effective length of column flange in bending component
    """
    # Connection using screws
    # Column information
    t_wc = column.tw  # Column web thickness
    A_vc = column.Ashear  # Shear area of column web
    f_y_wc = column.fy  # Yield strength of column material

    # Effective width of column web in tension
    # This is equal to the effective length of the T-stub
    # representing column flange.
    b_eff_t_wc = l_eff

    # Reduction factor
    w = min_coeff_w(beta, b_eff_t_wc, t_wc, A_vc)

    if verb:
        print("w = {0:5.5f}".format(w))
        print("t_wc = {0:4.3f}".format(t_wc))
        print("b_eff_tw = {0:4.3f}".format(b_eff_t_wc))
        print("f_y_wc = {0:4.3f}".format(f_y_wc))
        

    # Calculating design resistance of unstiffened column web subject to transverse tension
    F_t_wc_Rd = (w*b_eff_t_wc*t_wc*f_y_wc)/gammaM0
    
    return F_t_wc_Rd

def column_web_tension_k(column,beff):
    """ Stiffness of column web in compression """
    d_wc = column.hw
    t_wc = column.tw
    
    return 0.7*beff*t_wc/d_wc

def column_flange_bending(Tstub,verb=False):
    """ Column flange in bending, 6.2.6.4
        input:
            Tstub .. TStubColumnFlange object containing all the relevant information
            
    """
    
    return Tstub.FT_Rd(verb)

def column_flange_bending_k(column,leff,m):
    """ Column flange in bending
    
    """
    return 0.9*leff*column.tf**3/m**3

def end_plate_bending(Tstub,verb=False):
    """ End plate in bending, 6.2.6.5
        input:
            Tstub .. TStubEndPlate object containing all the relevant information
            
    """
    
    return Tstub.FT_Rd(verb)

def end_plate_bending_k(tp,leff,m):
    """ End plate in bending
    
    """
    return 0.9*leff*tp**3/m**3

def bolt_row_tension_k(As,Lb):
    """ Bolt row in tension (for two bolts in a row)
    """
    
    return 1.6*As/Lb
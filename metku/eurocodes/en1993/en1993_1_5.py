# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-5 Rules of plated structures

@author: kmela
"""

import math

try:
    from metku.eurocodes.en1993.constants import E, gammaM1
except:
    from eurocodes.en1993.constants import E, gammaM1


def shear_eta(fy):
    """ Computes the 'eta' according to EN 1993-1-5, clause 5.1(2) """

    eta = 1.0

    if fy <= 460.0:
        eta = 1.2

    return eta


def reduction_factor_internal(lp, psi):
    """ Reduction factor for internal compression elements
        4.4(2), and AC:2009
        input: lp .. slenderness (lambda_p)
              psi .. stress ratio
    """

    if lp <= 0.5 + math.sqrt(0.085 - 0.055 * psi):
        rho = 1.0
    else:
        rho = min((lp - 0.055 * (3 + psi)) / lp ** 2, 1.0)

    return rho


def reduction_factor_outstand(lp):
    """ Reduction factor for outstand compression elements
        4.4(2), and AC:2009
        
        input: lp .. slenderness (lambda_p)
    """

    if lp <= 0.748:
        rho = 1.0
    else:
        rho = min((lp - 0.188) / lp ** 2, 1.0)

    return rho


def lambda_p(b, t, epsilon, ksigma):
    """ Slenderness according to 4.4(2)
        input:
            b .. width of the element
            t .. thickness of the element
            epsilon = sqrt(235/fy)
            ksigma .. buckling factor
    """
    return b / t / 28.4 / epsilon / math.sqrt(ksigma)


def buckling_factor_internal(psi=1.0):
    """ Table 4.1 """

    if psi == 1.0:
        ksigma = 4.0
    elif psi > 0:
        ksigma = 8.2 / (1.05 + psi)
    elif psi == 0:
        ksigma = 7.81
    elif psi > -1:
        ksigma = 7.81 - 6.29 * psi + 9.78 * psi ** 2
    elif psi == -1.0:
        ksigma = 23.9
    elif psi > -3:
        ksigma = 5.98 * (1 - psi) ** 2

    return ksigma


def buckling_factor_outstand(psi=1.0, sigma1="tip"):
    """ Table 4.2 
        input:
            sigma1 .. location of the maximum compressive stress
                      "tip" = unsupported end
                      "support" = supported end
    """

    if sigma1 == "tip":
        if psi == 1.0:
            ksigma = 0.43
        elif psi == 0.0:
            ksigma = 0.57
        elif psi == -1:
            ksigma = 0.85
        else:
            ksigma = 0.57 - 0.21 * psi + 0.07 * psi ** 2
    else:
        if psi == 1.0:
            ksigma = 0.43
        if psi < 1.0 and psi > 0.0:
            ksigma = 0.578 / (psi + 0.34)
        elif psi == 0.0:
            ksigma = 1.7
        elif psi < 0 and psi > -1:
            ksigma = 1.7 - 5 * psi + 17.1 * psi ** 2
        else:
            ksigma = 23.8

    return ksigma


def effective_width_internal(b, psi, rho):
    """ EN 1993-1-5 Table 4.1
        Input:
            b .. width of the original plate
            psi .. stress ratio
            rho .. reduction factor for plate buckling
        
        Output:
            beff .. effective width
            be .. be[0] = b_e1, be[1] = b_e2
    """
    be = [0, 0]

    if abs(psi - 1.0) < 1e-6:
        # Psi = 1.0
        beff = rho * b
        be[0] = 0.5 * beff
        be[1] = 0.5 * beff
    elif psi > 0.0:
        beff = rho * b
        be[0] = 2 / (5 - psi) * beff
        be[1] = beff - be[0]
    else:
        # psi < 0
        beff = b * rho / (1 - psi)
        be[0] = 0.4 * beff
        be[1] = 0.6 * beff

    return beff, be


def effective_width_outstand(b, psi, rho):
    """ EN 1993-1-5 Table 4.2
        Input:
            b .. width of the original plate
            psi .. stress ratio
            rho .. reduction factor for plate buckling            
    """

    if psi < 0.0:
        beff = b * rho / (1 - psi)
    else:
        beff = rho * b

    return beff


def transverse_force_resistance(fyw, hw, tw, fyf, bf, tf, ss, a=0, type="a",
                                c=0):
    """ Clause 6 
        input:
            fyw .. yield strength of the web [MPa]
            hw .. height of the web (mm)
            tw .. web thickness (mm)
            fyf .. yield strength of the flange
            bf .. width of flange
            tf .. thicknes of flange
            a .. distance between transverse stiffeners (a=0 for no stiffeners)
            ss .. length of stiff bearing
            type .. type of load application
                    "a" .. through one flange at the middle of the member
                    "b" .. through two flanges at the middle of the member
                    "c" .. through one side of the flange at the end of the member
            c .. distance of the edge of the load application surface from the end
                of the member (for type "c")
    """
    if a == 0:
        a = 1000 * hw

    ss = min(ss, hw)

    if type == "a":
        kF = 6 + 2 * (hw / a) ** 2
    elif type == "b":
        kF = 3.5 + 2 * (hw / a) ** 2
    elif type == "c":
        kF = min(2 + 6 * (ss + c) / hw, 6.0)

    # Next, the reduction factor for effective length (Clause 6.4) and
    # Effective loaded length are determined. This may require iteration
    m1 = fyf * bf / fyw / tw
    m2 = 0

    Fcr = 0.9 * kF * E * tw ** 3 / hw

    if type == "a" or type == "b":
        ly = min(ss + 2 * tf * (1 + math.sqrt(m1 + m2)), a)
    elif type == "c":
        le = min(0.5 * kF * E * tw ** 2 / fyw / hw, ss + c)
        ly = min(le + tf * math.sqrt(0.5 * m1 + (le / tf) ** 2 + m2),
                 le + tf * math.sqrt(m1 + m2))

    lambdaF = math.sqrt(ly * tw * fyw / Fcr)

    if lambdaF > 0.5:
        # in this case, the initial assumption of m2=0 was wrong, and
        # ly needs to be corrected:
        m2 = 0.02 * (hw / tf) ** 2

    if type == "a" or type == "b":
        ly = min(ss + 2 * tf * (1 + math.sqrt(m1 + m2)), a)
    elif type == "c":
        le = min(0.5 * kF * E * tw ** 2 / fyw / hw, ss + c)
        ly = min(le + tf * math.sqrt(0.5 * m1 + (le / tf) ** 2 + m2),
                 le + tf * math.sqrt(m1 + m2))

    lambdaF = math.sqrt(ly * tw * fyw / Fcr)

    # Reduction factor (Eq. (6.3))
    chiF = min(0.5 / lambdaF, 1.0)

    Leff = chiF * ly

    FRd = fyw * Leff * tw / gammaM1

    return FRd


def shear_buckling_reduction_factor(slend_w, eta, end_post="non-rigid",
                                    verb=False):
    """ EN 1993-1-5, Table 5.1 """

    if slend_w < 0.83 / eta:
        chiW = eta
    elif slend_w < 1.08:
        chiW = 0.83 / slend_w
    else:
        if end_post == "non-rigid":
            chiW = 0.83 / slend_w
        elif end_post == "rigid":
            chiW = 1.37 / (0.7 + slend_w)
        else:
            raise ValueError(
                "en1993_1_5.shear_buckling_reduction_factor: bad end_post argument value.")

    if verb:
        print("Shear buckling reduction factor (EN 1993-1-5, Table 5.1)")
        print("End post:" + end_post)
        print("Slenderness: lambda_w = {0:4.3f}".format(slend_w))
        print("Eta: eta = {0:4.3f}".format(eta))
        print("Reduction factor: chi_w = {0:4.3f}".format(chiW))

    return chiW


def shear_buckling_coefficient(hw, a, verb=False):
    """ EN 1993-1-5, A.3
        input: 
            hw .. height of web
            a .. distance between transverse stiffeners
            
        It is assumed that there are no longitudinal stiffeners.
        This has to be modified, if longitudinal stiffeners are included!
    """
    rw = a / hw

    if rw >= 1.0:
        ktau = 5.34 + 4.0 / (rw ** 2)
    else:
        ktau = 4.0 + 5.34 / (rw ** 2)

    if verb:
        print("* Shear buckling coefficient (EN 1993-1-5, A.3): *")
        print("Distance between stiffeners a = {0:4.3f} mm".format(a))
        print("Web height hw = {0:4.3f} mm".format(hw))
        print("Buckling coefficient k_tau = {0:4.3f}".format(ktau))

    return ktau


def sigma_E(t, b):
    """ Critical stress of plate
        EN 1993-1-5, A.1
    """
    return 190_000 * (t / b) ** 2


def tau_crit(hw, a, t, b, verb):
    """ Critical shear stress for shear buckling
        EN 1993-1-5, Eq. (5.4)
    """

    if verb:
        print(
            "* Critical shear stress for shear buckling (EN 1993-1-5, Eq. (5.4)) *")
        print("Web height hw = {0:4.3f}".format(hw))
        print("Distance between stiffeners a = {0:4.3f}".format(a))

    ktau = shear_buckling_coefficient(hw, a, verb)
    sE = sigma_E(t, b)
    tau = ktau * sE

    if verb:
        print("Buckling stress sigma_E = {0:4.3f} MPa".format(sE))
        print("Shear buckling stress tau_crit = {0:4.3f} MPa".format(tau))

    return tau


def shear_buckling_slenderness(fyw, tau_cr, verb=False):
    """ Slenderness for shear buckling
        EN 1993-1-5, Eq. (5.3)
    """

    lambda_w = 0.76 * math.sqrt(fyw / tau_cr)

    if verb:
        print("* Shear buckling slenderness (EN 1993-1-5, Eq. (5.3)) *")
        print("Yield strength of the web: fyw = {0:4.3g}".format(fyw))
        print("Shear buckling stress: tau_crit = {0:4.3f}".format(tau_cr))
        print("Slenderness: lambda_w = {0:4.3f}".format(lambda_w))

    return lambda_w


def shear_buckling_web(chi_w, fyw, hw, t, verb=False):
    """ Shear buckling: contribution from the web
        EN 1993-1-5, Eq. (5.2)
    """

    VbwRd = chi_w * fyw * hw * t / (math.sqrt(3) * gammaM1)

    if verb:
        print("* Part of web: EN 1993-1-5, Eq. (5.2) *")
        print("Web height hw = {0:4.3f} mm".format(hw))
        print("Web thickness tw = {0:4.3f} mm".format(t))
        print("Web yield strength fyw = {0:4.3g} MPa".format(fyw))
        print("Resistance Vbw,Rd = {0:4.3f} kN".format(VbwRd * 1e-3))

    return VbwRd


def shear_buckling_flanges(bf, tf, fyf, a, hw, t, fyw, rM, verb=False):
    """ Shear buckling: contribution from the flanges
        EN 1993-1-5, Eq. (5.2)
    """

    c = a * (0.25 + 1.6 * bf * tf ** 2 * fyf / t / (hw ** 2) / fyw)

    VbfRd = bf * tf ** 2 * fyf / c / gammaM1 * (1 - rM ** 2)

    if verb:
        print("* Part of flanges: EN 1993-1-5, Eq. (5.2) *")
        print("Flange width bf = {0:4.3f} mm".format(bf))
        print("Flange thickness tf = {0:4.3f} mm".format(tf))
        print("Flange yield strength fyw = {0:4.3g} MPa".format(fyf))
        print("c = {0:4.3f} ".format(c))
        print("Resistance Vbf,Rd = {0:4.3f} kN".format(VbfRd * 1e-3))

    return VbfRd

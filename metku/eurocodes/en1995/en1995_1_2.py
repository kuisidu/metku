"""
@author Viktor Haimi
"""
import math
from eurocodes.en1995.fire_protection import FireProtection


def charring_speed(rhok: (int, float), type: str, hardness: str, h_p: (int, float)) -> list[float]:
    """
    Ensin tarkistaa RIL-205-2-2009 Taulukon 3.2 mukaisesti ja
    jos ei löydy tai paksuus on alle 20mm lasketaan beta_0_rho_t Kaavan 3.4 mukaan
    käyttää timber_datan rhok:ta eikä rhomeaniä

    @param rhok: ominaistiheys
    @param type: vaihtoehdot: solid_timber, glt, lvl
    @param hardness: vaihtoehdot: hardwood, softwood
    @param h_p: paksuus
    @return: list [beta_0, beta_n]
    """
    beta = None

    if type == 'solid_timber' or type == 'glt':
        if hardness == 'softwood':
            if type == 'solid_timber':
                if rhok >= 290:
                    beta = [0.65, 0.8]
            elif type == 'glt':
                if rhok >= 290:
                    beta = [0.65, 0.7]
        elif hardness == 'hardwood':
            if rhok == 290:
                beta = [0.65, 0.7]
            if rhok >= 450:
                beta = [0.5, 0.55]
    elif type == 'lvl':
        if rhok >= 480:
            beta = [0.65, 0.7]
        elif rhok >= 410:
            beta = [0.7, 0.75]
    else:
        beta = [1.0, 1.0]

    if h_p < 20:
        krho = math.sqrt(450 / rhok)
        k_h = math.sqrt(20 / h_p)
        beta = [krho * k_h * beta[0], beta[1]]

    return beta


def d_char_n(fp: FireProtection, t: (int, float), beta_n: float) -> float:
    """
    Returns the burned depth of single side
    @param fp: fire protection on the side
    @param t: how long the fire has burned
    @param beta_n: charring speed
    @return: depth of charring
    """

    k0 = 1.0
    d0 = 7
    # RIL-205-2-2009 Taulukko 4.1
    if fp:
        if fp.t_ch <= 20:
            if t < 20:
                k0 = t / 20
        else:
            if t <= fp.t_ch:
                k0 = t/fp.t_ch
    else:
        if t < 20:
            k0 = t / 20

    if not fp:
        return t * beta_n + k0 * d0
    # RIL-205-2-2009 kuvat 3.3 - 3.5
    if t < fp.t_ch:
        return 0
    if fp.t_ch < fp.t_f:
        if t < fp.t_f:
            return (t - fp.t_ch) * 0.73 * beta_n + k0 * d0
        elif fp.t_f < t < fp.t_a:
            return (fp.t_f - fp.t_ch) * 0.73 * beta_n + (t - fp.t_f) * 2 * beta_n + k0 * d0
        else:
            return (fp.t_f - fp.t_ch) * 0.73 * beta_n + (fp.t_a - fp.t_f) * 2 * beta_n + (t - fp.t_a) * beta_n + k0 * d0
    if t < fp.t_a:
        dt = t - fp.t_f
        return 2 * beta_n * dt + k0 * d0
    else:
        t_fa = fp.t_a - fp.t_f
        dt = t - fp.t_a
        return 2 * beta_n * t_fa + beta_n * dt + k0 * d0
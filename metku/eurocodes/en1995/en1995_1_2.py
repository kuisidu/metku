"""
@author Viktor Haimi
"""
import math

def charring_speed(rhok, type, hardness, h_p):
    """
    Ensin tarkistaa Taulukon 3.2 mukaisesti ja jos ei löydy sieltä lasketaan beta_0_rho_t
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

    # olettaa että tiheydet on oikein
    if h_p < 20:
        krho = math.sqrt(450 / rhok)
        k_h = math.sqrt(20 / h_p)
        beta = [krho * k_h * beta[0], beta[1]]

    if beta is not None:
        return beta
    return beta

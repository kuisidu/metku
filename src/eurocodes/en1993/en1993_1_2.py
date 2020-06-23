# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:46 2020

EN 1993-1-2 Supplementary rules for fire design

@author: kmela
"""

import math
import numpy as np

from eurocodes.en1993.constants import gammaM0, gammaM1, gammaM2

# EN 1993-1-2: Table 3.1
REDUCTION_FACTORS_CS = {20:{'ky':1.0,'kp':1.0,'kE':1.0},
                        100:{'ky':1.0,'kp':1.0,'kE':1.0},
                        200:{'ky':1.0,'kp':0.807,'kE':0.9},
                        300:{'ky':1.0,'kp':0.613,'kE':0.8},
                        400:{'ky':1.0,'kp':0.420,'kE':0.7},
                        500:{'ky':0.780,'kp':0.360,'kE':0.6},
                        600:{'ky':0.470,'kp':0.180,'kE':0.310},
                        700:{'ky':0.230,'kp':0.075,'kE':0.130},
                        800:{'ky':0.110,'kp':0.050,'kE':0.090},
                        900:{'ky':0.060,'kp':0.0375,'kE':0.0675},
                        1000:{'ky':0.040,'kp':0.0250,'kE':0.0450},
                        1100:{'ky':0.020,'kp':0.0125,'kE':0.0225},
                        1200:{'ky':0.000,'kp':0.0,'kE':0.0},
                        }

# EN 1993-1-2, Table C.1
REDUCTION_FACTORS_14301 = {20:{'ky':0.26,'k02':1.0,'ku':1.0,'kE':1.0,'kEct':0.11,'epsU': 0.40},
                           100:{'ky':0.24,'k02':0.82,'ku':0.87,'kE':0.96,'kEct':0.05,'epsU': 0.40},
                           200:{'ky':0.19,'k02':0.68,'ku':0.77,'kE':0.92,'kEct':0.02,'epsU': 0.40},
                           300:{'ky':0.19,'k02':0.64,'ku':0.73,'kE':0.88,'kEct':0.02,'epsU': 0.40},
                           400:{'ky':0.19,'k02':0.60,'ku':0.72,'kE':0.84,'kEct':0.02,'epsU': 0.40},
                           500:{'ky':0.19,'k02':0.54,'ku':0.67,'kE':0.80,'kEct':0.02,'epsU': 0.40},
                           600:{'ky':0.22,'k02':0.49,'ku':0.58,'kE':0.76,'kEct':0.02,'epsU': 0.35},
                           700:{'ky':0.26,'k02':0.40,'ku':0.43,'kE':0.71,'kEct':0.02,'epsU': 0.30},
                           800:{'ky':0.35,'k02':0.27,'ku':0.27,'kE':0.63,'kEct':0.02,'epsU': 0.20},
                           900:{'ky':0.38,'k02':0.14,'ku':0.15,'kE':0.45,'kEct':0.02,'epsU': 0.20},
                           1000:{'ky':0.040,'k02':0.06,'ku':0.07,'kE':0.2,'kEct':0.02,'epsU': 0.20},
                           1100:{'ky':0.040,'k02':0.03,'ku':0.03,'kE':0.1,'kEct':0.02,'epsU': 0.20},
                           1200:{'ky':0.040,'k02':0.0,'ku':0.00,'kE':0.0,'kEct':0.02,'epsU': 0.20},
                           }

SPECIAL_TEMPS = np.array([20,100,200,300,400,500,600,700,800,900,1000,1100,1200])




def standard_fire(t):
    """ Returns the temperature of standard fire ISO 834 
        input:
            t: time in minutes
        output:
            Temperature in Celsius
    """
    return 20 + 345*math.log(8*t+1,10)

def hydrocarbon_fire(t):
    """ Returns the temperature of hydrocarbon fire
        input:
            t: time in minutes
        output:
            Temperature in Celsius
    """
    return 20 + 1080*(1-0.325*math.exp(-0.167*t)-0.675*math.exp(-2.5*t))

def linear_interpolation(p1,p2,x):
    """ Performs linear interpolation.
        p1 = (x1,y1), p2 = (x2,y2)
        
        The function finds the value y corresponding to x1 <= x <= x2:
        Such that y(x) = k*x + b
    """
    
    y1 = p1[1]
    x1 = p1[0]
    y2 = p2[1]
    x2 = p2[0]
    
    k = (y2-y1)/(x2-x1)
    b = y1-k*x1
    
    return k*x+b

def reduce_property(t,prop,steel='carbon'):
    """ Function for determining the reduction factor at temperature 't'
        of property 'prop'
    """
    # Check if 't' matches any of the special temperatures
    if any(SPECIAL_TEMPS == t):
        if steel == 'carbon':
            rf = REDUCTION_FACTORS_CS[t][prop]
        elif steel =='1.4301':
            rf = REDUCTION_FACTORS_14301[t][prop]
    else:
        # t does not match any of the special temperatures;
        # Do linear interpolation
        
        """ Find the largest smaller and smallest larger
            temperatures and the corresponding reduction factors
        """
        t_smaller = np.max(SPECIAL_TEMPS[SPECIAL_TEMPS < t])
        t_larger = np.min(SPECIAL_TEMPS[SPECIAL_TEMPS > t])
        
        if steel == 'carbon':
            p_smaller = REDUCTION_FACTORS_CS[t_smaller][prop]
            p_larger = REDUCTION_FACTORS_CS[t_larger][prop]
        elif steel == '1.4301':
            p_smaller = REDUCTION_FACTORS_14301[t_smaller][prop]
            p_larger = REDUCTION_FACTORS_14301[t_larger][prop]
        
        rf = linear_interpolation((t_smaller,p_smaller), (t_larger,p_larger), t)
    
    return rf

""" Chapter 3: Materials """
def EaT(t,Ea=210000,steel='carbon'):
    """ Temperature-dependent Young's modulus
        for carbon steel or stainless steel
        input:
            t .. temperature [Celsius]
    """
    kE = reduce_property(t,'kE',steel)
    return kE*Ea

def fyT(t,fy=355,steel='carbon'):
    """ Temperature-dependent yield strength
        for carbon steel
        input:
            t .. temperature [Celsius]
    """
    ky = reduce_property(t,'ky',steel)
    return ky*fy

def fpT(t,fp=355,steel='carbon'):
    """ Temperature-dependent yield strength
        for carbon steel
        input:
            t .. temperature [Celsius]
    """
    kp = reduce_property(t,'kp',steel)
    return kp*fp

def f02pT(t,fy=210,steel="1.4301"):
    """ Temperature-dependent value for proof strength at 0.2 plastic strain """
    
    k = reduce_property(t,'k02',steel)
    
    return k*fy

def fuT(t,fy_T):
    """ EN 1993-1-2, Annex A(2) 
        Ultimate strength of carbon steel, when taking
        into account work hardening
    """
    
    if t < 300:
        fu_T = 1.25*fy_T
    elif t < 400:
        fu_T = fy_T*(2-0.0025*t)
    else:
        fu_T = fy_T
        
    return fu_T
    

def stress(strain,t,Ea,fya,hardening=False):
    """ Returns the stress of carbon steel for given
        strain and
        temperature 't'
        
        EN 1993-1-2, Fig. 3.1
    """
    eps_yt = 0.02
    eps_tt = 0.15
    eps_ut = 0.20
    
    # Proportional limit and Young's modulus at temperature 't'
    fp_t = fpT(t,fya)
    Ea_t = EaT(t,Ea)
    
    eps_pt = fp_t/Ea_t
    
    if strain <= eps_pt:
        st = Ea_t*strain
    elif strain < eps_yt:        
        fy_t = fyT(t,fya)
        c = (fy_t-fp_t)**2/((eps_yt-eps_pt)*Ea_t-2*(fy_t-fp_t))
        b = math.sqrt(c*(eps_yt-eps_pt)*Ea_t+c**2)
        a = math.sqrt((eps_yt-eps_pt)*(eps_yt-eps_pt+c/Ea_t))
        st = fp_t - c + (b/a)*math.sqrt(a**2-(eps_yt-strain)**2)
    elif strain <= eps_tt:
        fy_t = fyT(t,fya)
        if hardening:
            fu_T = fuT(t,fy_t)
            if strain < 0.04:
                st = 50*(fu_T-fy_t)*strain + 2*fy_t - fu_T
            else:
                st = fu_T
        else:
            st = fy_t            
    elif strain < eps_ut:
        fy_t = fyT(t,fya)
        if hardening:
            fu_T = fuT(t,fy_t)
            st = fu_T*(1-20*(strain-0.15))
        else:
            st = fy_t*(1-(strain-eps_tt)/(eps_ut-eps_tt))
    else:
        st = 0.0
        
    return st

def stress_stainless(strain,t,Ea,fya,steel='1.4031'):
    """ Stress-strain relationship of stainless steels at elevated temperatures """
    
    # Proportional limit and Young's modulus at temperature 't'
    f02p_t = f02pT(t,fya,steel)
    Ea_t = EaT(t,Ea,steel)
    
    eps_ct = f02p_t/Ea_t + 0.002
    kE_ct = reduce_property(t,'kEct',steel)
    Ee_ct = kE_ct*Ea
    
    eps_u = reduce_property(t,'epsU',steel)
     
    
    if strain <= eps_ct:
        b = (1-eps_ct)
        a = (Ea_t*eps_ct - f02p_t)/f02pt/pow(eps_ct,b)
        st = Ea_t*strain
    elif strain < eps_yt:        
        fy_t = fyT(t,fya)
        c = (fy_t-fp_t)**2/((eps_yt-eps_pt)*Ea_t-2*(fy_t-fp_t))
        b = math.sqrt(c*(eps_yt-eps_pt)*Ea_t+c**2)
        a = math.sqrt((eps_yt-eps_pt)*(eps_yt-eps_pt+c/Ea_t))
        st = fp_t - c + (b/a)*math.sqrt(a**2-(eps_yt-strain)**2)
    elif strain <= eps_tt:
        fy_t = fyT(t,fya)
        if hardening:
            fu_T = fuT(t,fy_t)
            if strain < 0.04:
                st = 50*(fu_T-fy_t)*strain + 2*fy_t - fu_T
            else:
                st = fu_T
        else:
            st = fy_t            
    elif strain < eps_ut:
        fy_t = fyT(t,fya)
        if hardening:
            fu_T = fuT(t,fy_t)
            st = fu_T*(1-20*(strain-0.15))
        else:
            st = fy_t*(1-(strain-eps_tt)/(eps_ut-eps_tt))
    else:
        st = 0.0
        
    return st
    

def therm_elongation(t):
    """ Relative thermal elongation of steel
        EN 1993-1-2, 3.4.1.1
    """
    
    if t < 750:
        dL = 1.2e-5*t + 0.4e-8*t**2-2.416e-4
    elif t <= 860:
        dL = 1.1e-2
    else:
        dL = 2e-5*t -6.2e-3
    
    return dL

def therm_conductivity(t):
    """ Thermal conductivity of steel [W/mK]
        EN 1993-1-2, 3.4.1.3
    """
    if t < 800:
        lambda_a = 54-3.33e-2*t
    else:
        lambda_a = 27.3
        
    return lambda_a
    

def spec_heat(t):
    """ Specific heat of steel [J/kgK] 
        EN 1993-1-2, 3.4.1.2
    """
    
    if t < 600:
        ca = 425+7.73e-1*t - 1.69e-3*t**2+2.22e-6*t**3
    elif t < 735:
        ca = 666 + 13002/(738-t)
    elif t < 900:
        ca = 545+17820/(t-731)
    else:
        ca = 650
    
    return ca

def therm_elongation_stainless(t):
    """ Relative thermal elongation of stainless steel
        EN 1993-1-2, C.3.1
    """
    
    return (16+4.79e-3*t-1.243e-6*t**2)*(t-20)*1e-6
    

def therm_conductivity_stainless(t):
    """ Thermal conductivity of stainless steel [W/mK]
        EN 1993-1-2, C.3.3
    """
    return 14.6+1.27e-2*t
    

def spec_heat_stainless(t):
    """ Specific heat of stainless steel [J/kgK] 
        EN 1993-1-2, C.3.2
    """
    return 450 + 0.280*t-2.91e-4*t**2 + 1.34e-7*t**3
    

def steel_temp(ksh,AmV,tmax,gas_temp=standard_fire,dt=5,Tsteel_ini=20):
    """ Determine temperature of unprotected steel at time 'T'
        starting from t=0 in time intervals 'dT'
    """
    rho_a = 7850 # kg/m3
    
    t = 0
    Tsteel = [Tsteel_ini]
    Time = [0]
    
    while t < tmax:
        """ Assume that gas_temp is a function that takes
            the time in minutes and returns the temperature
            in Celsius.
            
            Alternatively 'gas_temp' could be a table of
            pre-calculated temperatures. Then, dt should match
            the time interval used to produce the table.
            This feature is NOT yet implemented.
        """
        t_gas = gas_temp(t/60)
        dTsteel = ksh*AmV/spec_heat(Tsteel[-1])/rho_a*hnet(t_gas,Tsteel[-1])*dt
        Tsteel.append(Tsteel[-1]+dTsteel)
        t += dt
        Time.append(t/60)
    
    return Tsteel, Time

def hnet(t_gas,t_steel,phi=1.0,eps_f=1.0,eps_m=0.7,ac=25):
    """ Net heat flux per unit area [W/m2] """
    s = 5.67e-8 # Stefan-Boltzmann -constant [W/(m2K4)]

    """ heat flux from conductivity """
    hnet_c = ac*(t_gas-t_steel)
    
    """ heat flux from radiation """
    hnet_r = phi*eps_m*eps_f*s*((t_gas+273)**4-(t_steel+273)**4)
    
    return hnet_c+hnet_r


def stress_strain_curves(temp=[400,500,600,700],hardening=False,write=False):
    """ Plots several stree-strain curves """
    import matplotlib.pyplot as plt
    
    strain = np.linspace(0,0.2,200)
    st = []
    
    for idx, t in enumerate(temp):
        st.append([])
        for e in strain:
            st[idx].append(stress(e,t,210000,355,hardening)/355)
     
    if write:
        f = open(r'C:\Users\kmela\OneDrive - TUNI.fi\TRY\LaTeX\figures\luku6\jannitys_venyma_hiiliteras_lujittuminen.dat','w')
        top_row = 'Venymä'
        for t in temp:
            top_row += ' {0:4.0f}C'.format(t)
        f.write(top_row + '\n')
        
        for idx, e in enumerate(strain):
            row = "{0:5.3f} ".format(e)
            for s in st:
                row += ' {0:5.3f}'.format(s[idx])
            row += '\n'
            
            f.write(row)
        
    for stresses in st:
        plt.plot(strain,stresses)
    plt.grid(True)
    plt.xticks([0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2])
    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

    return strain, st

def unprotected_sections(write=False):
    import matplotlib.pyplot as plt
    from sections.steel.ISection import IPE, HEA, HEB, HEM
    from sections.steel.RHS import SHS
    
    """
    p1 = IPE(300)
    p2 = HEA(300)
    p3 = HEB(300)
    p4 = HEM(300)
    profiles = [p1,p2,p3,p4]
    """
    p1 = SHS(180,5)
    p2 = SHS(180,8)
    p3 = SHS(200,12.5)
    
    profiles = [p1,p2,p3]
    Tsteel = []
    
    for idx, p in enumerate(profiles):
        Tsteel.append([])
        ksh = p.shadow_effect()
        AmV = p.section_factor
        
        Tsteel[idx], times = steel_temp(ksh,AmV,30*60,gas_temp=standard_fire)
        Tgas = []
        for t in times:
            Tgas.append(standard_fire(t))
        
        plt.plot(times,Tsteel[idx])
        plt.plot(times,Tgas)
        plt.grid(True)
        #plt.xticks(np.arange(0,65,5))
        #plt.yticks(np.arange(0,1050,50))
        plt.xticks(np.arange(0,35,5))
        plt.yticks(np.arange(0,900,50))
    
    if write:
        f = open(r'C:\Users\kmela\OneDrive - TUNI.fi\TRY\LaTeX\figures\luku6\suojaamattomat_hiilivety.dat','w')
        """
        top_row = 'Venymä'
        for t in temp:
            top_row += ' {0:4.0f}C'.format(t)
        f.write(top_row + '\n')
        """
        for idx, t in enumerate(times):
            row = "{0:5.3f} ".format(t)
            for T in Tsteel:
                row += ' {0:5.3f}'.format(T[idx])
            row += '\n'
            
            f.write(row)
    
if __name__ == "__main__":
    
    #unprotected_sections(write=False)
    #strain,st = stress_strain_curves(write=True)
    strain,st = stress_strain_curves(temp=[300,350,400],hardening=True,write=True)
    

    
    
    
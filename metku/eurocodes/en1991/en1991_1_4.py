# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 15:16:22 2020

Calculations from EN 1991-1-4, Wind loads

@author: kmela
"""

import math

# Fundamental value of the basic wind velocity, m/s
vb0 = 21
# Directional factor
cdir = 1.0
# Seasonal factor
cseason = 1.0

# Density of air (kg/m3)
rho_air = 1.25

zmax = 200

terrain_category = {0: {"z0":0.003, "zmin":1.0},
                    1: {"z0":0.01, "zmin":1.0},
                    2: {"z0":0.05, "zmin":2.0},
                    3: {"z0":0.3, "zmin":5.0},
                    4: {"z0":1.0, "zmin":10.0}}

""" EN 1991-1-4, Table 7.1 """
ext_pressure_walls = {"A": {5: {"cpe10":-1.2, "cpe1":-1.4},
                            1: {"cpe10":-1.2, "cpe1":-1.4},
                            0.25: {"cpe10":-1.2, "cpe1":-1.4}},
                      "B": {5: {"cpe10":-0.8, "cpe1":-1.1},
                            1: {"cpe10":-0.8, "cpe1":-1.1},
                            0.25: {"cpe10":-0.8, "cpe1":-1.1}},
                      "C": {5: {"cpe10":-0.5, "cpe1":-0.5},
                            1: {"cpe10":-0.5, "cpe1":-0.5},
                            0.25: {"cpe10":-0.5, "cpe1":-0.5}},
                      "D": {5: {"cpe10":+0.8, "cpe1":1.0},
                            1: {"cpe10":+0.8, "cpe1":1.0},
                            0.25: {"cpe10":+0.7, "cpe1":1.0}},
                      "E": {5: {"cpe10":-0.7, "cpe1":-0.7},
                            1: {"cpe10":-0.5, "cpe1":-0.5},
                            0.25: {"cpe10":-0.3, "cpe1":-0.3}},
                      }

""" EN 1991-1-4, Table 7.2 """
ext_pressure_roof_sharp_eaves = {"F": {"cpe10":-1.8, "cpe1":-2.5},
                            "G": {"cpe10":-1.2, "cpe1":-2.0},
                            "H": {"cpe10":-0.7, "cpe1":-1.2},
                            "I": {"cpe10":[0.2,-0.2], "cpe1":[0.2,-0.2]}}

""" EN 1991-1-4, Table 7.3a """
ext_pressure_monopitch_roof = {"F": {"cpe10":-1.8, "cpe1":-2.5},
                            "G": {"cpe10":-1.2, "cpe1":-2.0},
                            "H": {"cpe10":-0.7, "cpe1":-1.2},
                            "I": {"cpe10":[0.2,-0.2], "cpe1":[0.2,-0.2]}}                      


def basic_wind_velocity(vb0=21.0,cdir=1.0,cseason=1.0):
    """ EN 1991-1-4, Eq. (4.1) """
    return vb0*cdir*cseason

def terrain_roughness(z,cat=2):
    """ EN 1991-1-4, Eq. (4.4) 
    
        :param z: height (m)
        :param cat: terrain category (0,..,4)
    """    
    z02 = 0.05
    z0  = terrain_category[cat]["z0"]
    zmin = terrain_category[cat]["zmin"]
    
    if cat == 0:
        """ Finnish Natioanl Annex """
        kr = 0.18
    else:
        kr = 0.19*pow(z0/z02,0.07)
    
    return kr*math.log(max(zmin,z)/z0)
    
    
def mean_wind_velocity(vb,z,cat=2):
    """ EN 1991-1-4, Eq. (4.3) """
    cr = terrain_roughness(z,cat)
    # Orography factor
    co = 1.0
    
    return vb*cr*co

def basic_velocity_pressure(vb,rho=1.25):
    """ EN 1991-1-4, Eq. (4.10) """
    return 0.5*rho*vb**2

def turbulence_intensity(z,cat=2):
    """ EN 1991-1-4, Eq. (4.7) """
    # Turbulence factor 
    kI = 1.0
    # Orography factor
    co = 1.0
    
    z0  = terrain_category[cat]["z0"]
    zmin = terrain_category[cat]["zmin"]
    
    return kI/co/(math.log(max(zmin,z)/z0))

def peak_velocity_pressure(z,cat,vb):
    """ EN 1991-1-4, Eq. (4.8) 
        :param z: height (m)
        :param cat: terrain category
        :param vb: basic wind velocity
        
        Output: peak velocity pressure (N/m2)
    """
    
    Iv = turbulence_intensity(z,cat)
    vm = mean_wind_velocity(vb,z,cat)
    
    return (1 + 7*Iv)*0.5*rho_air*vm**2
    
def wind_pressure_ext(qp,cpe):
    """ EN 1991-1-4, Eq. (5.1) """
    
    return qp*cpe

def wind_pressure_int(qp,cpi):
    """ EN 1991-1-4, Eq. (5.2) """
    
    return qp*cpi


if __name__ == "__main__":
    
    vb = basic_wind_velocity()
    
    qp = peak_velocity_pressure(8, 2, vb)

    print(qp)    
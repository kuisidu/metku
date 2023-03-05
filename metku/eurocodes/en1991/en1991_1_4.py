# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
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

def interpolate(X,Y,x0):
    """ Linear interpolation between two points """
    
    return Y[0] + (Y[1]-Y[0])/(X[1]-X[0])*(x0-X[0])

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
# TODO!
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

def rect_building_wall_pressure(length,width,h,qp,wind_dir="long"):
    """ Calculates external wind pressures in various zones
        for a rectangular building.
        
        EN 1991-1-4, Figure 7.5
    
        :param length: length of the building
        :param width: width of the building
        :param h: height of the building
        :param qp: peak velocity pressure
        :param wind_dir: wind direction.
                        "long" .. wind to the long side
                        "short" .. wind to the short side
    """
    
    if wind_dir == "long":
        b = length
        d = width
    elif wind_dir == "short":
        b = width
        d = length
    
    e = min(b,2*h)
    
    """ Define a dict for the wind pressures in the various
        zones of the wall
    """
    zones = {"A":{"length":0.0,'cpe':0.0,'we':0.0},
             "B":{"length":0.0,'cpe':0.0,'we':0.0},
             "C":{"length":0.0,'cpe':0.0,'we':0.0},
             "D":{"length":b,'cpe':0.0,'we':0.0},
             "E":{"length":b,'cpe':0.0,'we':0.0},
             }
    
    if e < d:
        zones["A"]['length'] = e/5
        zones["B"]['length'] = 0.8*e
        zones["C"]['length'] = d-e
    elif e >= 5*d:
        zones["A"]['length'] = d
    else:
        """ e >= d, but e < 5*d """
        zones["A"]['length'] = 0.2*e
        zones["B"]['length'] = d-0.2*e
    
    """ Area of the wind """
    A = b*h
    
    hd_ratio = h/d
    
    print("*** Wind load ***")
    print("Wind to " + wind_dir + " side of the building")
    print("b = {0:4.2f} m".format(b))
    print("d = {0:4.2f} m".format(d))
    print("h = {0:4.2f} m".format(h))
    print("e = {0:4.2f} m".format(e))
    print("h/d = {0:4.2f} ".format(hd_ratio))
    
    if A > 10:
        zones["A"]['cpe'] = ext_pressure_walls["A"][5]["cpe10"]
        zones["B"]['cpe'] = ext_pressure_walls["B"][5]["cpe10"]
        zones["C"]['cpe'] = ext_pressure_walls["C"][5]["cpe10"]

        if hd_ratio >= 5:
            zones["D"]['cpe'] = ext_pressure_walls["D"][5]["cpe10"]
            zones["E"]['cpe'] = ext_pressure_walls["E"][5]["cpe10"]
        elif hd_ratio >= 1:
            """ Linear interpolation for E """
            zones["D"]['cpe'] = ext_pressure_walls["D"][1]["cpe10"]
            Y = [ext_pressure_walls["E"][1]["cpe10"],
                 ext_pressure_walls["E"][5]["cpe10"]]
            X = [1,5]            
            zones["E"]['cpe'] = interpolate(X,Y,hd_ratio)
        elif hd_ratio >= 0.25:
            """ Linear interpolation for D and E """            
            X = [0.25,1]                        
            YD = [ext_pressure_walls["D"][0.25]["cpe10"],
                 ext_pressure_walls["D"][1]["cpe10"]]            
            YE = [ext_pressure_walls["E"][0.25]["cpe10"],
                 ext_pressure_walls["E"][1]["cpe10"]]            
            zones["D"]['cpe'] = interpolate(X,YD,hd_ratio)
            zones["E"]['cpe'] = interpolate(X,YE,hd_ratio)
            
            #zones["E"]['cpe'] = ext_pressure_walls["E"][1]["cpe10"]
        else:
            zones["D"]['cpe'] = ext_pressure_walls["D"][0.25]["cpe10"]
            zones["E"]['cpe'] = ext_pressure_walls["E"][0.25]["cpe10"]

    for key, value in zones.items():
        value['we'] = qp*value['cpe']
        
    return zones

def monopitch_roof_pressure(length,width,h,qp,wind_angle=180,alpha=5):
    """ Calculates external wind pressures in various zones
        for a rectangular building.
        
        EN 1991-1-4, Figure 7.5
    
        :param length: length of the building
        :param width: width of the building
        :param h: height of the building
        :param qp: peak velocity pressure
        :param wind_anlge: wind direction.
                        0 .. up the slope
                        180 .. against the slope
                        90 .. perpendicular to the slope
    """
    
    if wind_angle == 0 or wind_angle == 180:
        b = length
        d = width
    else:
        b = width
        d = length
    
    e = min(b,2*h)
    
    """ Define a dict for the wind pressures in the various
        zones of the wall
    """
    zones = {"F":{"length":0.25*e,'width':0.1*e,'cpe':0.0,'we':0.0},
             "G":{"length":0.5*e,'width':0.1*e,'cpe':0.0,'we':0.0},
             "H":{"length":b,'width':0.9*e,'cpe':0.0,'we':0.0},
             "I":{"length":b,'width':0.0,'cpe':0.0,'we':0.0},
             "Fup":{"length":0.25*e,'width':0.1*e,'cpe':0.0,'we':0.0},
             "Flow":{"length":0.25*e,'width':0.1*e,'cpe':0.0,'we':0.0},
             }
    
    if wind_angle == 90:
        zones['H']['width'] = 0.4*e
        zones['I']['width'] = length-0.5*e
            
    
    """ Area of the wind """
    A = b*h
    
    hd_ratio = h/d
    
    print("*** Wind load ***")
    print("Wind to " + wind_dir + " side of the building")
    print("b = {0:4.2f} m".format(b))
    print("d = {0:4.2f} m".format(d))
    print("h = {0:4.2f} m".format(h))
    print("e = {0:4.2f} m".format(e))
    print("h/d = {0:4.2f} m".format(hd_ratio))
    
    if A > 10:
        zones["A"]['cpe'] = ext_pressure_walls["A"][5]["cpe10"]
        zones["B"]['cpe'] = ext_pressure_walls["B"][5]["cpe10"]
        zones["C"]['cpe'] = ext_pressure_walls["C"][5]["cpe10"]

        if hd_ratio >= 5:
            zones["D"]['cpe'] = ext_pressure_walls["D"][5]["cpe10"]
            zones["E"]['cpe'] = ext_pressure_walls["E"][5]["cpe10"]
        elif hd_ratio >= 1:
            """ Linear interpolation for E """
            zones["D"]['cpe'] = ext_pressure_walls["D"][1]["cpe10"]
            Y = [ext_pressure_walls["E"][1]["cpe10"],
                 ext_pressure_walls["E"][5]["cpe10"]]
            X = [1,5]            
            zones["E"]['cpe'] = interpolate(X,Y,hd_ratio)
        elif hd_ratio >= 0.25:
            """ Linear interpolation for D and E """            
            X = [0.25,1]                        
            YD = [ext_pressure_walls["D"][0.25]["cpe10"],
                 ext_pressure_walls["D"][1]["cpe10"]]            
            YE = [ext_pressure_walls["E"][0.25]["cpe10"],
                 ext_pressure_walls["E"][1]["cpe10"]]            
            zones["D"]['cpe'] = interpolate(X,YD,hd_ratio)
            zones["E"]['cpe'] = interpolate(X,YE,hd_ratio)
            
            #zones["E"]['cpe'] = ext_pressure_walls["E"][1]["cpe10"]
        else:
            zones["D"]['cpe'] = ext_pressure_walls["D"][0.25]["cpe10"]
            zones["E"]['cpe'] = ext_pressure_walls["E"][0.25]["cpe10"]

    for key, value in zones.items():
        value['we'] = qp*value['cpe']
        
    return zones

def flat_roof_pressure(d,b,h,qp,hp=0):
    """ Calculates external wind pressures in various zones
        for a flat roof.
        
        EN 1991-1-4, Figure 7.6
    
        :param d: length of the building in wind direction
        :param b: length of the building in direction perpendicular to wind
        :param h: height of the building (to eaves)
        :param qp: peak velocity pressure
        :param hp: length of parapets (hp = 0 means no parapets)
    """
    
    e = min(b,2*h)
    
    """ Define a dict for the wind pressures in the various
        zones of the wall
    """
    zones = {"F":{"length":0.25*e,'width':0.1*e,'cpe':0.0,'we':0.0},
             "G":{"length":b-0.5*e,'width':0.1*e,'cpe':0.0,'we':0.0},
             "H":{"length":b,'width':0.4*e,'cpe':0.0,'we':0.0},
             "I":{"length":b,'width':d-0.5*e,'cpe':0.0,'we':0.0}
             }
    
    hd_ratio = h/d
    
    print("*** Wind load ***")
    print("Wind to flat roof")
    print("qp = {0:4.2f} ".format(qp))
    print("b = {0:4.2f} m".format(b))
    print("d = {0:4.2f} m".format(d))
    print("h = {0:4.2f} m".format(h))
    print("e = {0:4.2f} m".format(e))
    print("h/d = {0:4.2f} m".format(hd_ratio))
    
    if hp == 0:
        for key, zone in zones.items():
            zones[key]['cpe'] = ext_pressure_roof_sharp_eaves[key]["cpe10"]
            if isinstance(zone['cpe'],list):
                zones[key]['we'] = [qp*cpe for cpe in zone['cpe']]
            else:
                zones[key]['we'] = qp*zone['cpe']
        
    return zones

#if __name__ == "__main__":

    #wind = rect_building_wall_pressure(length=72,width=24,h=8,qp=qp,wind_dir="long")
    #wind = rect_building_wall_pressure(length=15,width=12,h=7,qp=qp,wind_dir="long")
    #wind = rect_building_wall_pressure(length=30,width=18,h=H,qp=qp,wind_dir="long")
    #wind = rect_building_wall_pressure(length=12*6,width=5*4.8,h=H,qp=qp,wind_dir="short")

    #print(qp)    
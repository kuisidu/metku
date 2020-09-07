# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 22:50:51 2018

Constants for EN 1992

@author: kmela
"""

from math import sqrt


# Partial safety factors
gammaC = 1.5    # Concrete
gammaS = 1.15   # Steel

# EN 1992-1-1 3.1.6
alpha_cc = 0.85 


""" Material data
fck .. lieriölujuus
fck_cube .. kuutiolujuus
fcm
Ecm
"""
concrete = {"C12/15": {"fck": 12.0, "fck_cube": 15.0, "fcm": 20, "Ecm": 27000.0},
            "C16/20": {"fck": 16.0, "fck_cube": 20.0, "fcm": 24, "Ecm": 29000.0},
            "C20/25": {"fck": 20.0, "fck_cube": 25.0, "fcm": 28, "Ecm": 30000.0},
            "C25/30": {"fck": 25.0, "fck_cube": 30.0, "fcm": 33, "Ecm": 31000.0},
            "C30/37": {"fck": 30.0, "fck_cube": 37.0, "fcm": 38, "Ecm": 33000.0},
            "C35/45": {"fck": 35.0, "fck_cube": 45.0, "fcm": 43, "Ecm": 34000.0},
            "C40/50": {"fck": 40.0, "fck_cube": 50.0, "fcm": 48, "Ecm": 35000.0},
            "C45/55": {"fck": 45.0, "fck_cube": 55.0, "fcm": 53, "Ecm": 36000.0},
            "C50/60": {"fck": 50.0, "fck_cube": 60.0, "fcm": 58, "Ecm": 37000.0},
            }

class Concrete:
    
    def __init__(self,material="C25/30"):
        """ Constructor for Concrete material class """

        self.name = material
        self.fck = concrete[material]["fck"]
        self.fck_cube = concrete[material]["fck_cube"]
        
        self.fcm = self.fck + 8.0
        self.Ecm = 22000*(self.fcm/10)**0.3
        
        # density [kn/m3], from EN 1991-1-1, Table A.1
        self.density = 24
        self.RH = 50
        
        #self.fcm = concrete[material]["fcm"]
        #self.Ecm = concrete[material]["Ecm"]
        
    def fcd(self,a_cc=alpha_cc,gC=gammaC):
        """ Design value for compressive strength """
        return a_cc*self.fck/gC
    
    def nu_cracked(self):
        """ Reduction factor for concrete cracked due to shear """
        return 0.6*(1-self.fck/250)
    
    def creep_coefficient(self,h0,t,t0=None):
        """ According to EN 1992-1-1, Annex B """
        
        RH = self.RH
        
        """ Eq. (B.8c) """
        alpha = [0,0,0]
        alpha[0] = (35/self.fcm)**0.7
        alpha[1] = (35/self.fcm)**0.2
        alpha[2] = (35/self.fcm)**0.5
        
        """ Eq. (B.3) """
        if self.fcm <= 35:
            phiRH = 1 + (1-RH/100)/0.1/h0**(1/3)
            betaH =  min(1.5*(1+(0.012*RH)**18)*h0+250,1500)
        else:
            phiRH = (1 + (1-RH/100)/0.1/h0**(1/3)*alpha[0])*alpha[1]
            betaH =  min(1.5*(1+(0.012*RH)**18)*h0+250*alpha[2],1500*alpha[2])
        
        """ Eq. (B.4) """
        beta_cm = 16.8/sqrt(self.fcm)
        
        """ Eq. (B.5) """
        beta_t0 = 1/(0.1+t0**0.2)
        
        """ Eq. (B.2) """
        phi0 = phiRH*beta_cm*beta_t0
                
        
        """ Eq. (B.7) """
        beta_c = ((t-t0)/(betaH + t-t0))**0.3
                
        phi = phi0 * beta_c
        
        return phi
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 09:52:00 2018

@author: mercader
"""



import numpy as np
from pyswarm import pso




 def costf_weight(self, X=[], *args):
        """ Cost function used for optimization
            
            Parameters
            ----------
            X : 1-D-array of float, optional
                List containing rotational stiffness values [Sj1, Sj2]
                and h and b values in order
                [alpha11, alpha12,..., alphan2, h1, b1, h2,...,bn]
            Return:
                weight : float
                    Weight of the members
        """ 
        # Total number of beams
        num_beams = self.bays * self.storeys
        # Number of symmetrical beams
        sym_beams = int(self.storeys*(self.bays -self.bays%2)/2)
        # Number of beams to be designed
        beams = num_beams - sym_beams
        # Number of connections
        joints = num_beams - 2*sym_beams
        
        
        
        # If list X is empty, create a list x with current frame's dimensions
        if len(X) == 0:
            X = self.initial_x0()
        # If joints aren't optimized, no need to change the list
        if not self.optimize_joints:
            x = X      
        else:           
            # remove joint rigidities from list
            x = X[joints:]
        
        
        rho = 7850e-9
        # IPE values
        K = [0.0155, 2.4921, 0.0256, 3.2831, 0.3486, 30.178, 0.0155, 2.4921]
       
        # HEA values
        Q = [0.014, 4.2949, 0.0294, 5.7651, 1.0323, 2.5368, 0.0781, 2,9206]
       
        
        weight = 0
        index = 0
      
        # weight of the beams
        for i in range(beams):
            h = x[i]
            t_w = K[0] * h + K[1]
            t_f = K[2] * h + K[3]
            b = K[4] * h + K[5]
            r = K[6] * h + K[7]
          
            # Area of I-profile
            A = 2.0*t_f*b + (h - 2.0*t_f)*t_w + (4.0 - math.pi)*r**2.0
            weight += A * rho * self.bay_length * 1000  # kg's
        x = x[beams+2:]
         
        # weight of the columns
        for i in range(int(len(x)/2)):
            h = x[i]
            t_w = Q[0] * h + Q[1]
            t_f = Q[2] * h + Q[3]
            
            if h<=270:
                b = Q[4] * h + Q[5]
                r = Q[6] * h + Q[7]
            else:
                b=300
                if h<=690:
                    r=27
                else:
                    r=30
                
            
            # Area of HEA-profile
            A = 2.0*t_f*b + (h - 2.0*t_f)*t_w + (4.0 - math.pi)*r**2.0
            weight += A * rho * self.storey_height * 1000 # kg's
        
        # Add penalty for every non-satisfactory constraint
        if self.optimizer == "DIF_EVOL":
            value = self.penalty_val
            a = np.asarray(self.constraints(X)) < 0
            a = sum(a)
            penalty = a * value
            return weight + penalty
        
        return weight





# Define the objective (to be minimize)
def weight(x, *args):
    r,h = x
    B, rho, E, P = args
    return rho*2*np.pi*d*t*np.sqrt((B/2)**2 + H**2)

# Setup the constraint functions
def yield_stress(x, *args):
    H, d, t = x
    B, rho, E, P = args
    return (P*np.sqrt((B/2)**2 + H**2))/(2*t*np.pi*d*H)

def buckling_stress(x, *args):
    H, d, t = x
    B, rho, E, P = args
    return (np.pi**2*E*(d**2 + t**2))/(8*((B/2)**2 + H**2))

def deflection(x, *args):
    H, d, t = x
    B, rho, E, P = args
    return (P*np.sqrt((B/2)**2 + H**2)**3)/(2*t*np.pi*d*H**2*E)

def constraints(x, *args):
    strs = yield_stress(x, *args)
    buck = buckling_stress(x, *args)
    defl = deflection(x, *args)
    return [100 - strs, buck - strs, 0.25 - defl]

# Define the other parameters
B = 60  # inches
rho = 0.3  # lb/in^3
E = 30000  # kpsi (1000-psi)
P = 66  # kip (1000-lbs, force)
args = (B, rho, E, P)

# Define the lower and upper bounds for H, d, t, respectively
lb = [10, 1, 0.01]
ub = [30, 3, 0.25]

xopt, fopt = pso(weight, lb, ub, f_ieqcons=constraints, args=args)

# The optimal input values are approximately
#     xopt = [29, 2.4, 0.06]
# with function values approximately
#     weight          = 12 lbs
#     yield stress    = 100 kpsi (binding constraint)
#     buckling stress = 150 kpsi
#     deflection      = 0.2 in
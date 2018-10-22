# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 09:18:15 2018

@author: mercader
"""

from frame2d import Frame2D, SteelBeam, SteelColumn, PointLoad,\
     FixedSupport, XHingedSupport, LineLoad, YHingedSupport, Hinge


import sys
# path needs to be added to use tables_and_tuples and I_sections
sys.path.append('S:\91202_METKU\Kristo\Python\src\End-plate')
sys.path.append('S:\91202_METKU\Kristo\Python\src\sections\steel')

from timeit import default_timer as timer

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fmin_slsqp, differential_evolution

from steel_section import SteelSection
from ISection import ISection, HEA, IPE


import frame_fem as fem
#import I_sections
#import anna_HEA_sections
from pso import pso
from steel_member import SteelMember

import tables_and_tuples







     
     
OPTIMIZERS = ["PSO", "SLSQP", "DIF_EVOL"]

#frame2d=Frame2D()

class optimization:
    def __init__(self,h):
        self.h=h
        self.frame2d=Frame2D()
        self.is_optimized = False
        self.self_weight = False
        self.penalty_val = 1000
        self.optimizer = "PSO"
        self.optimize_joints = True
        
    def costf_weight(self, X=[]):
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
            num_beams = frame2d.bays * frame2d.storeys
            # Number of symmetrical beams
            sym_beams = int(frame2d.storeys*(frame2d.bays -frame2d.bays%2)/2)
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
            
            #rho in kg/mm^3
            rho = 7850e-9
            
            """
            # IPE values
            K = [0.0256, 3.2831, 0.0155, 2.4921,  0.3486, 30.178, 0.0155, 2.4921]
           
            # HEA values
            Q = [0.0294, 5.7651, 0.014, 4.2949,  1.0323, 2.5368, 0.0781, 2,9206]
           
            
            weight = 0
            index = 0
            """
            
            # weight of the beams
            
            
            for key in frame2d.members.keys():
                member = frame2d.members[key]
                h=member.h
                b=member.b
                t_f = member.tf
                t_w = member.tw
                
           
                
                # Area profile
                A = 2.0*t_f*b + (h - 2.0*t_f)*t_w + (4.0 - math.pi)*r**2.0
                #weight += A * rho * self.storey_height * 1000 # kg's
                weight += A * rho * frame2d.storey_height # kg's
            
            # Add penalty for every non-satisfactory constraint
            if self.optimizer == "DIF_EVOL":
                value = self.penalty_val
                a = np.asarray(self.constraints(X)) < 0
                a = sum(a)
                penalty = a * value
                print ("pes penalty:",weight)
                return weight + penalty
           
            
            return weight
    
    
    
    def func2(self, x):
        print(len(x))
        for i, member in enumerate(self.members.values()):
            member.h = x[i*2]
            print(member.h)
        
        #self.generate_frame()
        self.calculate()
       
    def func(self, x):
        """ Function, that assigns x's-values to frame's members
            
            Parameters
            ----------
            x : 1-D-array of float
                Array containing rotational stifness values [alpha1, alpha2]
                and h and b values
        """ 
        num_beams = self.bays * self.storeys
        sym_beams = int(self.storeys*(self.bays -self.bays%2)/2)
        index = 0
        if self.optimize_joints:
            # Change beam's joint rigidities
            # Symmetrical beams
            i = 0
            for i in range(sym_beams):
                key = "Beam" + str(i)
                member = self.members[key]
                member.alpha1 = x[index]
                member.alpha2 = x[index+1]
                index += 2
            # Beams on symmetry axis  
            for j in range(num_beams - 2*sym_beams):
                key = "Beam" + str(i+j)
                member = self.members[key]
                member.alpha1 = x[index]
                member.alpha2 = x[index]
                index += 1            
        # Change beam's profile
        # Symmetrical beams
        value = 0
        for _ in range(math.floor(self.bays/2)):
            for i in range(self.storeys):
                key1 = "Beam" + str(i + value)
                key2 = "Beam" + str(i - value 
                                    + int(self.storeys*(self.bays-1)))              
                member1 = self.members[key1]
                member2 = self.members[key2]                
                #member1.profile = "Custom"
                #member2.profile = "Custom"  
                member1.h = x[index]
                member2.h = x[index]
                #member1.custom_section_properties_IPE([x[index], x[index+1]])
                #member2.custom_section_properties_IPE([x[index], x[index+1]])                
                # Symmetrical joint rigidities
                member2.alpha2 = member1.alpha1
                member2.alpha1 = member1.alpha2                
                member1.alpha_to_Sj()
                member2.alpha_to_Sj()
                index += 2
            value += self.storeys           
        # If there're beams on the symmetry axis
        if self.bays % 2 != 0:
            for i in range(self.storeys):
                key = "Beam" + str(i + value)
                member = self.members[key]
                #member.profile = "Custom"   
                member.h = x[index]
                #member.custom_section_properties_IPE([x[index], x[index+1]])                
                # Symmetrical joint rigidities
                member.alpha1 = member.alpha2
                member.alpha_to_Sj()
                index += 2                   
        # Change column's profile
        # Symmetrical columns
        value = 0
        for j in range(math.ceil(self.bays/2)):
            for i in range(self.storeys):
                key1 = "Column" + str(i + value)
                key2 = "Column" + str(i + self.storeys * (self.bays - j))                
                member1 = self.members[key1]
                member2 = self.members[key2]                
                #member1.profile = "Custom"
                #member2.profile = "Custom"  
                member1.h = x[index]
                member2.h = x[index]
                #member1.custom_section_properties_HEA([x[index], x[index+1]])
                #member2.custom_section_properties_HEA([x[index], x[index+1]])
                index += 2                
            value += self.storeys            
        # If there're columns on the symmetry axis
        # uses i from previous loop
        if self.bays % 2 == 0:
                for _ in range(self.storeys):
                    i += 1
                    key = "Column" + str(i)
                    member = self.members[key]
                    member.profile = "Custom"
                    member.h=x[index]
                    #member.custom_section_properties_HEA([x[index], x[index+1]])
                    index += 2            
        #self.generate_frame()
        self.calculate()


    def constraints(self, x=[]):
        """ Constraints for optimization
            Satisfactory when all values in the array >= 0
            
            Parameters
            ----------
            x : 1-D-array of float
                Array containing optimized members values
            x = [] in the function call allows to check these constraints
                without x-array
                
            Returns
            -------
            cons : 2-D-array of float
                Array of utilization ratio's
        """
        if len(x) != 0:
            self.func(x)
            #self.func2(x)
        # Sway constraint 
        top_node = self.num_elements * self.storeys
        height = self.storey_height * self.storeys
        max_sway = height / 400 # SWAY LIMIT
        sway = self.nodal_displacements[top_node][0]
        cons1 = max_sway - sway        
        # Beam deflection
        cons2 =[]
        num_beams = self.bays * self.storeys
        for i in range(num_beams):
            member = self.members["Beam" + str(i)]
            max_v = 0
            for node in member.nodal_displacements.keys():
                v_val = member.nodal_displacements[node][1]
                if abs(v_val) > max_v:
                    max_v = abs(v_val)
            cons2.append(max_v)
        cons2 = np.asarray(cons2)
        max_deflection = self.bay_length / 300 # DEFLECTION LIMIT
        cons2 = max_deflection - cons2             
        cons = [cons1] + list(cons2)
        # Cross-section, buckling and lateral-buckling strength
        for key in self.members.keys():
            member = self.members[key]
            for r in member.r:
                cons.append(1-r)                           
        return cons
    
    
    
    def log_optimization(self, parameters, result, time):
        """ Function for logging optimization results
        
            Parameters
            ----------
            parameters : dict
                Parameters used to get optimization result
            result : float
                Optimized frame's weight
            time : float
                Time elapsed for optimization
        """
        file_name = self.optimizer + "_optimization_log.txt"
        print("Optimization results logged to {}".format(file_name))
        with open(file_name, 'a') as f:          
            f.write("Frame properties: \n")
            f.write("    Bays: {} \n".format(self.bays))
            f.write("    Bay length: {} \n".format(self.bay_length))
            f.write("    Storeys: {} \n".format(self.storeys))
            f.write("    Storey height: {} \n".format(self.storey_height))
            f.write("    Elements: {} \n".format(self.num_elements))
            if self.line_loads != {}:
                f.write("Loads: \n")
                for loadid in self.line_loads.keys():
                    f.write("    LoadID: {} \n".format(loadid))
                    f.write("    Line loads: \n")
                    for loads in self.line_loads[loadid]:
                        member = loads[0]
                        q_value = loads[1]
                        direction = loads[2]
                        f.write("        Member: {}\n".format(member))
                        f.write("        q value [kN/m]: {}\n".format(q_value))
                        f.write("        Direction: {}\n".format(direction))
                    if self.point_loads[loadid]:
                        f.write("    Point loads: \n")
                        for point_loads in self.point_loads[loadid]:
                            node = point_loads[0]
                            value = point_loads[1]
                            factor = point_loads[2]
                            f.write("        Node: {}\n".format(node))
                            f.write("        Value [Fx,Fy,Mz] [kN]: {}\n"
                                    .format(value))
                            f.write("        Factor: {}\n".format(factor))
            f.write("Optimizer : {}\n".format(self.optimizer))
            if not self.optimize_joints:
                f.write("Joint rigidities NOT optimized")
            f.write("Parameters: \n")
            for key in parameters.keys():
                f.write("    {} : {} \n".format(key, parameters[key]))
            f.write("Optimization result: {0:.2f} kg \n".format(result))
            f.write("Optimization result is feasible: {}\n"
                    .format(((min(self.constraints()) > 0))))
            f.write("    MAX r: {:.3f}% \n"
                    .format((1-min(self.constraints()))*100))
            f.write("Optimized profiles: \n")
            for key in self.members.keys():
                member = self.members[key]
                f.write("    {} : {} r:{:.2f} %\n"
                        .format(key,member.profile, (max(member.r)*100)))
                if member.mtype == 'beam':
                    f.write("    Sj1 : {:.2f}  Sj2 : {:.2f} \n"
                            .format(member.Sj1, member.Sj2))
            f.write("Time elapsed: {0:.2f} s \n".format(time))
            f.write("{}\n".format('-'*81))
    
    
    def initial_x0(self):
        """ Creates initial guess x0 from current frame's member's properties
        """        
        num_beams = self.bays * self.storeys
        sym_beams = int(self.storeys*(self.bays -self.bays%2)/2)
        x0 = []
        index = 0
        # Append beam's joint rigidities
        # Symmetrical beams
        i = 0
        for i in range(sym_beams):
            key = "Beam" + str(i)
            member = self.members[key]
            x0.append(member.alpha1)
            x0.append(member.alpha2)
            index += 2
        # Beams on symmetry axis  
        for j in range(num_beams - 2*sym_beams):
            key = "Beam" + str(i+j)
            member = self.members[key]
            x0.append(member.alpha1)
            index += 1
        # Append beam's h and b
        # Symmetrical beams
        value = 0
        for _ in range(math.floor(self.bays/2)):
            for i in range(self.storeys):
                key = "Beam" + str(i + value)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2
            value += self.storeys            
        # If there're beams on the symmetry axis
        if self.bays % 2 != 0:
            for i in range(self.storeys):
                key = "Beam" + str(i + value)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2
        # Append column's h and b
        # Symmetrical columns
        value = 0
        for _ in range(math.ceil(self.bays/2)):
            for i in range(self.storeys):
                key = "Column" + str(i + value)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2               
            value += self.storeys           
        # If there're columns on the symmetry axis
        if self.bays % 2 == 0:
            for i in range(self.storeys):
                key = "Column" + str(i + self.bays)
                member = self.members[key]
                x0.append(member.h)
                x0.append(member.b)
                index += 2   
        return x0
    
    
    
    
    
    
        
         
    def optimize_frame(self, debug=True, print_data=False, optimizer="",
                       swarmsize=50, maxiter=20,omega=0.5, phip=0.5,
                       phig=0.5, penalty_val=1000, joints=False, log_data=False,
                       draw=True, lb=[80], ub=[500]):
        
        """ Optimizes frame's members and beam's joints rotational stifness
            
            Parameters
            ----------
            debug : bool (Default: True)
                boolean value, if True prints optimization data 
            print_data : bool (Default: False)
                boolean value, if True prints members' info
            optimizer : string
                Optimization algorithm used to optimize the frame
                    "PSO" -- Particle Swarm optimization
                    "SLSQP" -- Sequential least squares
                    "DIF-EVOL" -- Differential evolution
            swarmsize : int (Default: 50)
                PSO number of particles used for optimization
            maxiter : int (Default: 20)
                number of maximum iterations
            omega : float (Default: 0.5)
                PSO Particle velocity scaling factor 
            phip : float (Default: 0.5)
                PSO Scaling factor to search away from the particle's best
                known position 
            phig : float (Default: 0.5)
                PSO Scaling factor to search away from the swarm's
                best known position
            joints : bool (Default: True)
                boolean value, if True optimizes joint rigidities
            penalty_val : float (Default: 1000)
                DIF-EVOL Penalty value that's added to cost function's value 
                for every non-satisfactory constraint
            log_data : bool (Default: False)
                Logs optimization data
            draw : bool (Default: True)
                Draws frame after optimization
            lb : 1-D-array of float (Default: [80, 46])
                array of lower boundary values [h, b] [mm, mm]
            ub : 1-D-array of float (Default: [500, 246])
                array of upper boundary values [h, b] [mm, mm]               
        """
        self.optimize_joints = joints
        # If no optimizer is given in the fuction call, use last used optimizer
        # Also checks if given optimizer is in the optimizer list and raises
        # error if not, otherwise changes frame's optimizer to that algorithm
        if optimizer == "":
            optimizer = self.optimizer
         
        else:
            optimizer = optimizer.upper()
            if optimizer not in OPTIMIZERS:
                raise ValueError("There's no optimizer called {}"
                                 .format(optimizer))
            else:
                self.optimizer = optimizer
           
        # Set penalty value
        if self.penalty_val != penalty_val:
            self.penalty_val = penalty_val
        # Set boundary arrays
        #if len(ub) != 2 or len(lb) != 2:
        #   raise ValueError("Boundary arrays must have exactly two values.")           
        # Set true for printing optimized values
        self.is_optimized = True
        # Start timer
        start = timer()
        # Number of beams and symmetrical beams
        num_beams = self.bays * self.storeys
        sym_beams = self.storeys*(self.bays -self.bays%2)/2       
        # Objective function
        obj_func = self.costf_weight()
        # Constraint function
        cons = self.constraints        
        # Symmetrical beams
        # Joint rigidity boundaries according to SFS-EN 1993-1-8
        lb_j = [0.5]*2*int(sym_beams)
        ub_j = [25]*2*int(sym_beams)        
        # Beams on symmetry axis
        lb_j += [0.5]*1*(num_beams - 2*int(sym_beams))
        ub_j += [25]*1*(num_beams - 2*int(sym_beams))
        if not joints:
            lb_j = []
            ub_j = []
        # Symmetrical members
        # [h, b]
        num_sym_members = self.bays*self.storeys + self.storeys
        lb = lb_j + lb*num_sym_members
        ub = ub_j + ub*num_sym_members  
        # Particle Swarm Optimization
        if optimizer == "PSO": 
            """
            # Designs frame and uses those values as initial guess
            # also scales upper and lower boundary values according to
            # designed members
            x0 = np.asarray(self.initial_x0())
            lb_x0 = x0[len(lb_j):]*0.5
            ub_x0 = x0[len(ub_j):]*1.5
            lb = lb_j + list(lb_x0)
            ub = ub_j + list(ub_x0)
            """           
            xopt, fopt = pso(obj_func,
                             lb,
                             ub,
                             swarmsize=swarmsize,
                             maxiter=maxiter,
                             f_ieqcons=cons,
                             debug=debug,
                             omega=omega,
                             phip=phip,
                             phig=phig)
            if xopt == []:
                return
            self.func(xopt)
            # Parameters saved to dict for logger
            parameters = {"swarmsize": swarmsize,
                          "maxiter": maxiter,
                          "omega": omega,
                          "phip": phip,
                          "phig": phig,
                          "result": xopt}
        # Sequential Least SQuares 
        if optimizer == "SLSQP":     
            """
            # Designs frame and uses those values as initial guess
            # also scales upper and lower boundary values according to
            # designed members
            x0 = np.asarray(self.initial_x0())
            lb_x0 = x0[len(lb_j):]*0.5
            ub_x0 = x0[len(ub_j):]*1.5
            lb = lb_j + list(lb_x0)
            ub = ub_j + list(ub_x0)
            """
            bounds = []
            for i in range(len(lb)):
                bounds.append((lb[i], ub[i]))
            # initial guess
            x0 = np.random.rand(len(ub)) * np.asarray(ub)
            #x0 = self.initial_x0()           
            if debug:
                print("x0: ",x0)
                iprint = 2
            else:
                iprint = 0
            out, fx, its, imode, smode = fmin_slsqp(obj_func,
                                                    x0,
                                                    f_ieqcons=cons,
                                                    bounds=bounds,
                                                    iprint=iprint,
                                                    iter=maxiter,
                                                    full_output=True)
            
            self.constraints(out)
            parameters = {
                    "x0": x0,
                    "iterations": its,
                    "maxiter": maxiter,
                    "exit mode": smode,
                    "result": out
                    }
        # Differential evolution
        if optimizer == "DIF_EVOL":
            """
            x0 = np.asarray(self.initial_x0())
            lb_x0 = x0[len(lb_j):]-100
            ub_x0 = x0[len(ub_j):]+100
            lb = lb_j + list(lb_x0)
            ub = ub_j + list(ub_x0)
            """
            x0 = np.random.rand(len(ub)) * np.asarray(ub)
            bounds = []
            
            if debug:
                print("x0: ",x0)
                disp = True
            else:
                disp = False
            
            for i in range(len(lb)):
                bounds.append((lb[i], ub[i]))
            result = differential_evolution(obj_func,
                                            bounds,
                                            disp=disp,
                                            maxiter=maxiter)
            self.constraints(result.x)
            parameters = {
                    "x0": x0,
                    "maxiter": maxiter,
                    "penalty": self.penalty_val
                    }            
        # Change optimized continuous members to discrete
        self.continuous_to_discrete() 
        # End timer
        end = timer()
        print("Time elapsed: ", end - start, "s")
        if log_data:
            self.log_optimization(parameters, self.weight, end - start)       
        if print_data:
            self.print_member_info()
        if draw:
            self.draw_frame()
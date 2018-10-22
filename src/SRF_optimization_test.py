# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 10:41:23 2018

@author: huuskoj
"""

def alpha_func(self, x):
    """ Cost function for optimization
        Input:
        x -- list containing rotational stifness values [Sj1, Sj2]
        Return:
        self.weight -- weight of the frame
    """        
    num_beams = self.bays * self.storeys
    index = 0
    for i in range(num_beams):
        key = "Beam" + str(i)
        member = self.members[key]
        member.alpha_to_Sj(x)
        index += 2
    
    for key in self.members.keys():
        member = self.members[key]
        profile_idx = math.floor(x[index])
        member.profile = profiles_list[profile_idx]
        member.profile_idx = profile_idx
        index += 1
        
    self.generate_frame()
    self.calculate()   
    return self.weight

def constraints(self, x, *args):
    
    # Sway constraint 
    top_node = self.num_elements * self.storeys
    height = self.storey_height * self.storeys
    max_sway = height / 400
    sway = self.nodal_displacements[top_node][0]
    cons1 = 1 - sway / max_sway
    
    # Cross-section strength check
    if self.is_strong_enough:
        cons2 = 1
    else:
        cons2 = -1
    return [cons1, cons2]

def optimize_frame(self, debug=False, print_data=False):
    """ Optimizes frame's members and beam's joints rotational stifness
        kwargs:
            debug -- boolean value, if True prints PSO debug data 
            print_data -- boolean value, if True prints members' info
    """
    self.is_optimized = True
    start = timer()
    num_beams = self.bays * self.storeys
    # boudaries for joint rigidness values
    # alpha boundaries
    lb = [0]*2*num_beams
    ub = [1]*2*num_beams
    # add profile indices to boundaries
    lb += [0]*self.num_members
    ub += [len(profiles_list)]*self.num_members

    pso(self.func, lb, ub, swarmsize=50, maxiter=50, f_ieqcons=self.constraints, debug=debug)
    end = timer()
    print("Time elapsed: ", end - start, "s")
    if print_data:
        self.print_member_info()
    self.draw_frame()
    
    
    # ANOTHER IDEA FOR OPTIMIZATION
    def optimize_frame(self, debug=False, print_data=False):
    """ Optimizes frame's members and beam's joints rotational stifness
        kwargs:
            debug -- boolean value, if True prints PSO debug data 
            print_data -- boolean value, if True prints members' info
    """
    self.is_optimized = True
    start = timer()
    num_beams = self.bays * self.storeys
    # boudaries for joint rigidness values
    # alpha boundaries
    lb = [0, 0]*num_beams #  [Sj1, Sj2]
    ub = [1, 1]*num_beams #  [Sj1, Sj2]
    # add profile indices to boundaries
    lb += [0, 0, 0, 0, 0, 0]*self.num_members # [A,Iy,Iz,Au,Wply, Wplz ,Wely, Welz, Ashear]
    ub += [1e20, 1e20, 1e20, 1e20, 1e20, 1e20]*self.num_members # [A,I,Au,Wpl,Wel,Ashear]
    pso(self.func, lb, ub, swarmsize=50, maxiter=50, f_ieqcons=self.constraints, debug=debug)
    end = timer()
    print("Time elapsed: ", end - start, "s")
    if print_data:
        self.print_member_info()
    self.draw_frame()
    
    

    # SYMMETRICAL COLUMNS
###############################################################################
    def func(self, x):
        """ Cost function for optimization
            Input:
            x -- list containing rotational stifness values [Sj1, Sj2]
            Return:
            self.weight -- weight of the frame
        """        
        num_beams = self.bays * self.storeys
        index = 0
        # Change beam's joint rigidities
        for i in range(num_beams):
            key = "Beam" + str(i)
            member = self.members[key]
            # SYMMETRICAL JOINT RIGIDITIES!
            member.alpha_to_Sj([x[index], x[index]])       
            index += 1
        # Change beam's profile
        for j in range(num_beams):
            member = "Beam" + str(j)
            self.members[member].profile = profiles_list[math.floor(x[index])]
            index += 1
        # Symmetrical columns
        value = 0
        while self.bays - value >= 0:     
            for i in range(self.storeys):
                member1 = "Column" + str(i + value)
                member2 = "Column" + str(i + value + self.storeys * (self.bays - value))        
                self.members[member1].profile = profiles_list[math.floor(x[index])]
                self.members[member2].profile = profiles_list[math.floor(x[index])] 
                index += 1
            value += 2  
        self.generate_frame()
        self.calculate()
        
        return self.weight
    
    def constraints(self, x, *args):
        
        # Sway constraint 
        top_node = self.num_elements * self.storeys
        height = self.storey_height * self.storeys
        max_sway = height / 400 # SWAY LIMIT
        sway = self.nodal_displacements[top_node][0]
        cons1 = 1 - sway / max_sway
        
        # Cross-section strength check
        if self.is_strong_enough:
            cons2 = 1
        else:
            cons2 = -1
            
        # Beam deflection
        cons3 =[]
        num_beams = self.bays * self.storeys
        for i in range(num_beams):
            member = self.members["Beam" + str(i)]
            max_v = 0
            for node in member.nodal_displacements.keys():
                v_val = member.nodal_displacements[node][1]
                if abs(v_val) > max_v:
                    max_v = abs(v_val)
            cons3.append(max_v)
        cons3 = np.asarray(cons3)
        max_deflection = self.bay_length / 300 # DEFLECTION LIMIT
        cons3 = 1 - max(cons3 / max_deflection)
        
        
        return [cons1, cons2, cons3]
    
    def optimize_frame(self, debug=False, print_data=False):
        """ Optimizes frame's members and beam's joints rotational stifness
            kwargs:
                debug -- boolean value, if True prints PSO debug data 
                print_data -- boolean value, if True prints members' info
        """
        self.is_optimized = True
        start = timer()
        num_beams = self.bays * self.storeys
        # boudaries for joint rigidness values
        # alpha boundaries
        lb = [0]*num_beams
        ub = [1]*num_beams
        # add profile indices to boundaries
        lb += [0]*(self.num_members -self.bays*self.storeys)
        ub += [len(profiles_list)-1]*(self.num_members -self.bays*self.storeys)
    
        pso(self.func, lb, ub, swarmsize=50, maxiter=20,
            f_ieqcons=self.constraints, debug=debug)
        end = timer()
        print("Time elapsed: ", end - start, "s")
        if print_data:
            self.print_member_info()
        self.draw_frame()

###############################################################################
    
    
""" OLD OPTIMIZATION FUNCTIONS
    def func(self, x):
        """ Cost function for optimization
            Input:
            x -- list containing rotational stifness values [Sj1, Sj2]
            Return:
            self.weight -- weight of the frame

        # Sj values between kHinge and kRigid
        num_beams = self.bays * self.storeys
        index = 0
        for i in range(num_beams):
            key = "Beam" + str(i)
            member = self.members[key]
            member.Sj1 = x[index]
            member.Sj2 = x[index+1]
            index += 2
        """          
        # fixity factor alpha -values
        num_beams = self.bays * self.storeys
        index = 0
        for i in range(num_beams):
            key = "Beam" + str(i)
            member = self.members[key]
            member.alpha_to_Sj(x)
            index += 2

        for key in self.members.keys():
            member = self.members[key]
            profile_idx = math.floor(x[index])
            member.profile = profiles_list[profile_idx]
            member.profile_idx = profile_idx
            index += 1
        """   

        self.generate_frame()
        self.calculate()
        self.design_frame()    
        return self.weight
    
    
    def optimize_frame(self, debug=False, print_data=False):
        """ Optimizes frame's members and beam's joints rotational stifness
            kwargs:
                debug -- boolean value, if True prints PSO debug data 
                print_data -- boolean value, if True prints members' info
        """
        self.is_optimized = True
        start = timer()
        num_beams = self.bays * self.storeys
        # boundaries for joint rigidness values
        lb = [fem.kHinge]*2*num_beams  
        ub = [fem.kRigid]*2*num_beams
        """
        # alpha boundaries
        lb = [0]*2*num_beams
        ub = [1]*2*num_beams
        
        # add profile indices to boundaries
        lb += [0]*self.num_members
        ub += [len(profiles_list)-2]*self.num_members
        
        cons = [self.constraints]
        """
        pso(self.func, lb, ub, swarmsize=50, maxiter=50, debug=debug)
        end = timer()
        print("Time elapsed: ", end - start, "s")
        if print_data:
            self.print_member_info()
        self.draw_frame()
"""
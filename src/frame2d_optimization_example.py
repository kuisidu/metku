from frame2d import Frame2D, SteelBeam, SteelColumn, PointLoad,\
     FixedSupport, XHingedSupport, LineLoad, YHingedSupport, Hinge

import math
import numpy as np

from scipy.optimize import fmin_slsqp, differential_evolution

from timeit import default_timer as timer


def create_frame():
    """
    Creates frame
    :return: Frame2D object
    """
    # Create your frame here
    # You can add parameters to the function call so your frame is created
    # according to those parameters

    # Coordinates of members, loads and supports
    coord1 = [[0.0,0], [0.0, 1]]
    coord2 = [[0,1], [1,1]]
    coord3 = [[1,0], [1,1]]
    
    supp_coord1 = [0.0,0]
    supp_coord2 = [1,0]

    # Loads and members
    col1 = SteelColumn(coord1)
    col2 = SteelColumn(coord3)
    beam1 = SteelBeam(coord2)
    #load1 = PointLoad(coord1[1], [50, 0,0])
    load2 = LineLoad(beam1, [-50, -50], 'y')

    # Create empty frame 'envelope'
    frame = Frame2D(num_elements=3)

    # Add members
    frame.add(col1)
    frame.add(col2)
    frame.add(beam1)

    # Add loads
    #frame.add(load1)
    frame.add(load2)

    # Add supports
    frame.add(FixedSupport(supp_coord1, supp_id=1))
    frame.add(FixedSupport(supp_coord2, supp_id=2))




    # Generate frame
    frame.generate()

    # Return generated frame
    

    #change column profiles to HEA
    """    
    for key in frame.members.keys():
            member = frame.members[key]
            if member.mtype=="column":
                member.profile = "he 100 a"
    """



   



    frame.plot()

    # Calculate result
    frame.calculate()
    #print(frame.nodal_forces[2])
    
    
    #frame.f.draw()
    frame.bmd(20)
    frame.plot_deflection()

    print ("weight =", frame.weight,"kg")
   
    # Print nodal forces
    """
    for member in frame.members.values():
            for element in member.elements.values():
                print(element.bending_moment)
    """
    
    
    """
    for key in frame.members.keys():
        member = frame.members[key]
           
        print ("\n",member.profile, member.mtype,"h =",member.cross_section.h, "b=",member.cross_section.b, "tf=",
               member.cross_section.tf,"tw=",member.cross_section.tw,"r=",member.cross_section.r )
    print(frame.nodes)
               
    for i in range (len(frame.nodes)):
        print("\n","node",i,",",
              "nodal forces",frame.nodal_forces[i],
              "nodal displacements",frame.nodal_displacements[i])  
    """ 
    
    
    
    
    
    
    
        #h= frame.members[key].profile.split()[1]
        #print ("h=",h)
    
    #print ("nodal displacements", frame.nodal_displacements)
    
   
    #print (self.storey_height())
    #print(frame.bays())
        
    
    return frame
    
    """
        h=member.__h
        b=member.b
        t_f = member.tf
        t_w = member.tw
        
        print ("h=",h,"b=",b,"tf=",t_f,"tw",t_w)
        """       
    

def cost_function(X):
    """
    Cost function used in optimization
    :param X: list of float values, e.g. member's dimensions
    :return cost_value: value of cost_function
    """
    # Implement cost function
    # Assign new optimized X values to the frame here
    # e.g. frame.members[1].h = X[1]
    
    
    #rho in kg/mm^3
   
    for i, member in enumerate(frame.members.values()):
        member.h = X[i]
  
        #print ("bucle optimitzacio", member.profile)
     
    return frame.weight



def constraint_function(X):
    """
    Constraints used in the optimization
    :param X: list of float values, e.g. member's dimensions
    :return constraints: list of stress ratios
    """
    # Implement constraints for the optimization
    # Calculate new stress ratios with the new optimized dimensions X
    # e.g. stress_ratios = np.asarray(frame.r)
    # e.g. constraints = 1 - stress_ratios
    
    
    
    
    # Iterate through every member and change their height's
    
    """
    for i, member in enumerate(frame.members.values()):
        member.h = X[i]
    frame.calculate()
    """
    # Constraint values need to be <= 0
    # Stress indices must be less than 1
    stress_indices = np.asarray(frame.r)
   
    #constf_values = 1 - stress_indices
    constf_values = stress_indices - 1
    
    #print("constf_values", constf_values)
    
    return constf_values.flatten()
    

    
    
    """
    if len(X) != 0:
        self.func(X)
            #self.func2(x)
        # Sway constraint 
        top_node = frame.num_elements * frame.storeys
        height = self.storey_height * self.storeys
        max_sway = height / 400 # SWAY LIMIT
        sway = self.nodal_displacements[top_node][0]
        cons1 = max_sway - sway        
        # Beam deflection
        cons2 =[]
        num_beams = frame.bays * self.storeys
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
    """
    
    
    #pass

def optimize_frame():
    """
    Optimizes frame
    :param frame: Frame2D object
    :return:
    """
    # Implement function that optimizes the frame
    # Optimization functions usually call constraint_function and cost_function
    # in their own code
    
    
    #self.optimize_joints = joints
        # If no optimizer is given in the fuction call, use last used optimizer
        # Also checks if given optimizer is in the optimizer list and raises
        # error if not, otherwise changes frame's optimizer to that algorithm


    
    #"SLSQP":     
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
    
    start = timer()
    
    lb=[80]
    ub=[500]
    bounds = [(80,500)]
    

    
    
    #x0 = self.initial_x0()           
    
    debug=True
    maxiter=20
    
    
    lb=lb*len(frame.members)
    ub=ub*len(frame.members)
    #x0=x0*len(frame.members)
    bounds=bounds*len(frame.members)
    
    # initial guess
    x0 = np.random.rand(len(ub)) * np.asarray(ub)
    
    
    if debug:
        print("x0: ",x0)
        iprint = 2
    else:
        iprint = 0
   
    
        
    out, fx, its, imode, smode = fmin_slsqp(cost_function,
                                            x0,
                                            f_ieqcons=constraint_function,
                                            bounds=bounds,
                                            iprint=iprint,
                                            iter=maxiter,
                                            full_output=True)
        
    end = timer()
    print("Time elapsed: ", end - start, "s")     
        
 #no se si cal aquest tros
    

    constraint_function(out)
    parameters = {
            "x0": x0,
            "iterations": its,
            "maxiter": maxiter,
            "exit mode": smode,
            "result": out
            }
    
    

def continuous_to_discrete(frame):
    """
    Changes frame's member's continuous dimensions to discrete profiles
    :param frame: optimized Frame2D object
    """
    # Implement function that changes continuously optimized members to discrete profiles
    
    
        
    pass


if __name__ == '__main__':
    # Implement here function calls to optimize the frame
    """
    e.g.
    frame = create_frame()
    optimize_frame(frame)
    continuous_to_discrete(frame)
    frame.plot()
    """
    
    frame=create_frame()
    
    #print(cost_function(X=[100,80,100]))  
    #frame.plot_deflection(10)
    
    
    
    
    
    
    """
    optimize_frame()
   
    
    frame.plot()
    frame.bmd(20)
    frame.plot_deflection(10)
    
    print ("weight optimized frame =",frame.weight)
    for key in frame.members.keys():
        member = frame.members[key]
        print(member.profile, member.mtype)       
        
        #h= frame.members[key].profile.split()[1]
        #print ("h=",h)
        
        print ("h =",member.cross_section.h, "b=",member.cross_section.b, "tf=",
               member.cross_section.tf,"tw=",member.cross_section.tw,"r=",member.cross_section.r )

        #print ("vinclament", frame.r)
   """
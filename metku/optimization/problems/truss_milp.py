# -*- coding: utf-8 -*-
"""
Created on Sat Apr  3 20:51:03 2021

Discrete truss optimization by mixed-integer linear programming.

TODO:
    1) Modify equilibrium equations for multiple loading conditions
    2) Kinematic stability constraints for topology optimization (auxiliary load case)
    3) Chains for topology optimization
    4) Member grouping (sizing and topology)
    5) Profile counting
290
@author: kmela
"""

import numpy as np
import re

try:
    from metku.frame2d.frame2d import SteelBeam
    from metku.frame2d.simple_truss import SimpleTruss, three_bar_truss, ten_bar_truss
    from metku.sections.steel import SteelSection
    import metku.optimization.structopt as sop
    from metku.sections.steel.catalogue import ipe_profiles, shs_profiles, hea_profiles
    from metku.sections.steel import ISection, WISection, RHS, SHS
    from metku.optimization.solvers.lp import LP
    from metku.optimization.solvers.milp import MILP
except:
    from frame2d.frame2d import SteelBeam
    from sections.steel import SteelSection
    import optimization.structopt as sop
    from sections.steel.catalogue import ipe_profiles, shs_profiles, hea_profiles
    from sections.steel import ISection, WISection, RHS
    

class TrussMILP:
    """ Class for mixed-integer linear programming formulation for truss structures
    """
    
    def __init__(self,name,
                 truss,profiles,ulb=1e-3,uub=1e3,
                 scaling={'N':1.0, 'u':1.0},
                 groups=None,
                 member_profiles=None):
        """ Constructor
            Parameters:
            ----------
            :param name: name of the problem
            :param truss: truss structure to be optimized (ground structure for topology optimization)
            :param profiles: list of profiles from which the optimum is sought
            :param ulb: lower bound for nodal displacements
            :param uub: upper bound for nodal displacements
            :param scaling: scaling factors for axial force and displacements (lengths)
            :param groups: member groups
            :param member_profiles: for each truss member, a subset of the list 'profiles' can be given to limit
                                    the available profiles for that member.
        """
        self.name = name
        self.structure = truss
        self.profiles = profiles
        
        self.ulb = ulb*scaling['u']
        self.uub = uub*scaling['u']
        
        # Scaling
        self.scaling = scaling
        self.scaling['A'] = self.scaling['u']**2
        self.scaling['E'] = self.scaling['N']/(self.scaling['u']**2)
        
        # Design variables
        self.vars = {'profile':[[] for _ in range(self.nmem)],
                     'force':[[[] for _ in range(self.nlc)] for _ in range(self.nmem)],
                     'disp':[[] for _ in range(self.nlc)],'member':[],'node':[],'group':[]}
        
        self.__nvars = 0        
        
        self.milp = None
        
        # Resistance of cross-sections
        self.NRd = [np.floor(p.NRd)*self.scaling['A'] for p in profiles]
    
        # Member data
        self.members = [{'L':0.0,'k0':0.0,'profiles':[],
                         'NRd':[],'NbRd':[],'Nlb':[],'Nub':[],'w':[]} for _ in range(self.nmem)]
        
        # Node data:
        # keep track of supported nodes and connectivity. This data is needed
        # in topology optimization.
        self.nodes = [{'supported':0,'members':[]} for _ in range(self.nnodes)]
        
        # Allocate available member profiles
        # Assume initially that all members can have all profiles
        self.member_profiles = [list(range(self.nprof))] * self.nmem
        
        """ If member groups are given, allow for each member of a given group
            the indicated profiles
        """
        self.groups = groups
        
        if groups is not None:
            """ groups is a dict with keys: 'members','profiles' """
            for group in groups:
                for mem in group['members']:
                    self.member_profiles[mem] = group['profiles']
        
        if member_profiles is not None:
            """ member_profiles is list of dicts with keys: 'member', 'profiles' """
            for mem in member_profiles:
                self.member_profiles[mem['member']] = mem['profiles']
        
        # If no selection is given for member profiles, initialize all profiles to all members
        #if member_profiles is None:
        #    self.member_profiles = [list(range(self.nprof))] * self.nmem
        #else:
            # If member profiles is given, it must be a list of lists,
            # with length self.nmem. For each member, a list of available profile
            # indices is given.
        #    self.member_profiles = member_profiles
        
        # Groups here
        # group is a dict with keys 'members', 'profiles', which are lists
        # these lists contain indices of members that belong to the same group
        # and the corresponding profile alternatives.
        
        # Initialize member profiles
        for mem, profs in zip(self.members,self.member_profiles):
            mem['profiles'] = [self.profiles[i] for i in profs]
            mem['Nlb'] = np.zeros(len(profs))
            mem['Nub'] = np.zeros(len(profs))
            mem['NRd'] = [self.NRd[i] for i in profs]
            mem['NbRd'] = np.zeros(len(profs))
        
        #self.NbRd = np.zeros([self.nprof,self.nmem])
        
        B = self.statics_matrix
        # Initialize member data and variables that are common in all problems
        for i, mem in self.structure.members.items():
            # Calculate length and "axial stiffness"
            self.members[i]['L'] = mem.length*self.scaling['u']
            self.members[i]['k0'] = mem.E/self.members[i]['L']*self.scaling['E']
            # For each profile of the member:
            # i) calculate buckling resistance
            # ii) create profile selection variable
            # iii) calculate force variable bounds
            nlb = np.sum(np.maximum(B[:,i],0)*ulb) + np.sum(np.minimum(B[:,i],0)*uub)
            nub = np.sum(np.minimum(B[:,i],0)*ulb) + np.sum(np.maximum(B[:,i],0)*uub)
            for j, p in enumerate(self.members[i]['profiles']):
                mem.profile = p
                # Calculate buckling resistance for section 'p'
                self.members[i]['NbRd'][j] = np.floor(mem.NbRd[0])*self.scaling['N']
                
                # Add profile selection variable
                self.vars['profile'][i].append(sop.BinaryVariable("y({0:.0f},{1:.0f})".format(i,j),
                                                section=p,
                                                target={"objects":[mem],"property":"profile"}))
                self.__nvars += 1
                # Force variable bounds
                self.members[i]['Nlb'][j] = self.members[i]['k0']*p.A*nlb*self.scaling['A']
                self.members[i]['Nub'][j] = self.members[i]['k0']*p.A*nub*self.scaling['A']
        
                # Weight of member for profile 'p'
                self.members[i]['w'].append(p.weight()/self.scaling['u']*self.members[i]['L'])
                
        """ Design variables: member forces 
        Force variables are self.vars['force'][i][k][j], where i is member, j is profile and k is load case
    
        """
        for k, lc in enumerate(self.structure.load_cases):
            for i, mem in self.structure.members.items():                            
                for j, p in enumerate(self.members[i]['profiles']):
                    self.vars['force'][i][k].append(
                        sop.Variable(
                            name="N({0:.0f},{1:.0f},{2:.0f})".format(i,j,lc),
                            lb=self.members[i]['Nlb'][j],
                            ub=self.members[i]['Nub'][j],
                            target=None,value=1))
                    self.__nvars += 1
                    #c.append(0)
            """ Design variables: nodal displacements 
                u[j,k]
            """
            for j in range(self.ndof):
                self.vars['disp'][k].append(
                    sop.Variable(name="u({0:.0f},{1:.0f})".format(j,k),
                                 lb=ulb,ub=uub,target=None,value=0.0))
                self.__nvars += 1
                #c.append(0)
        
        
    @property
    def nmem(self):
        """ Number of truss members """
        
        return len(self.structure.members)

    @property
    def nnodes(self):
        """ Number of truss nodes """
        return len(self.structure.nodes)

    @property
    def nprof(self):
        """ Number of profiles """
        
        return len(self.profiles)
    
    @property
    def ndof(self):
        """ number of degrees of freedom """
        return self.structure.f.nfree_dofs
    
    @property
    def nlc(self):
        """ Number of load cases """
        return len(self.structure.load_cases)
    
    @property
    def nvars(self):
        """ Number of variables """
        return self.__nvars
    
    @property
    def nprof_vars(self):
        """ Total number of profile selection variables """    
        return sum([len(y) for y in self.vars['profile']])
    
    @property
    def nstate_vars(self):
        """ Total number of state variables, i.e. the sum of force and
            displacement variables
        """
        nu = sum([len(u) for u in self.vars['disp']])
        nf = sum([len(Nlc) for Nmem in self.vars['force'] for Nlc in Nmem])
        """
        nf = 0
        for Nmem in self.vars['force']:
            for Nlc in Nmem:
                nf += len(Nlc)
        """
        return nu + nf
    
    @property
    def statics_matrix(self):
        """ Statics matrix of the truss """
        return self.structure.f.statics_matrix()
    
    def start_ndx(self,var_type):
        """ Returns the start index of the first variable of 'var_type' """
        ndx = 0
        
        if var_type == 'profiles':
            ndx = 0
        elif var_type == 'member':
            ndx = self.nprof_vars + self.nstate_vars
        elif var_type == 'node':
            ndx = self.nprof_vars + self.nstate_vars + self.nmem
        elif var_type == 'group':
            ndx = self.nprof_vars + self.nstate_vars + self.nmem + self.nnodes

        return ndx
    
    def weight_obj(self):
        """ Returns the vector for weight of the truss """
        c = np.zeros(self.nvars)
        start = 0
        # It is assumed that the profile variables are first
        # in the vector of design variables
        for mem in self.members:
            # Find the end position for the current member profile variables.
            end = start + len(mem['w'])
            # The weights of the member for each profile are stored in 'w'.
            c[start:end] = mem['w']
            start = end
        
        return c
    
    def compatibility_constraints(self,B):
        """ Constraints: compatibility conditions """
        ny_vars = self.nprof_vars    
        ndof = self.ndof
        nx = self.nvars
        Nndx = self.nprof_vars
        for k, lc in enumerate(self.structure.load_cases):
            yndx = 0
            # Starting variable index for displacement variables
            #          y            N               u 
            ustart = ny_vars + (k+1) * ny_vars + k*ndof            
            uend = ustart +ndof
            for i, mem in enumerate(self.members):                            
                for j, p in enumerate(mem['profiles']):
                    a = np.zeros(nx)
                    """ Upper bound """
                    a[yndx] = mem['Nub'][j]              # y variable
                    a[Nndx] = -1  # N variable
                    
                    a[ustart:uend] = mem['k0']*p.A*B[:,i]*self.scaling['A']                  # u variable
                    self.milp.add(sop.LinearConstraint(a,mem['Nub'][j],con_type="<", name="COMP-UB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
        
                    """ Lower bound """
                    a[Nndx] = 1  # N variable
                    a[yndx] = -mem['Nlb'][j]             # y variable
                    a[ustart:uend] *= -1    # u variable
                    self.milp.add(sop.LinearConstraint(a,-mem['Nlb'][j],con_type="<", name="COMP-LB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
                    
                    yndx += 1                
                    Nndx +=1
                    
    def strength_constraints(self):
        """ Constraints: member strength (stress constraints)
            These are expressed in terms of member forces, i.e.
            -NRd[j]*y[i,j] <= N[i,j,k] <= NRd[j]*y[i,j]
        """
        nx = self.nvars
        # Starting index for force variables
        Nndx = self.nprof_vars
        for k, lc in enumerate(self.structure.load_cases):
            yndx = 0
            for i, mem in enumerate(self.members):
                for j, Nrd in enumerate(mem['NRd']):
                    a = np.zeros(nx)
                    # Tension resistance 
                    a[yndx] = -Nrd # -Nub[j,i] #               # y variable                
                    a[Nndx] = 1  # N variable
                    self.milp.add(sop.LinearConstraint(a,0,con_type="<", name="STR-UB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
                
                    #a[i] = -slb[i]
                    # Compression resistance, including buckling 
                    a[Nndx] = -1
                    a[yndx] = -min(Nrd,mem['NbRd'][j])# Nlb[j,i] # 
                    self.milp.add(sop.LinearConstraint(a,0,con_type="<", name="STR-LB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
        
                    Nndx += 1
                    yndx += 1
            # Skip the displacement variables in variable indexing
            Nndx += self.ndof
            
    def add_top_opt_vars(self):
        """ Adds binary variables for topology optimization """
        # Add member and node existence variables
        for i, mem in self.structure.members.items():
            self.vars['member'].append(sop.BinaryVariable("Y({0:.0f})".format(i),
                                                          target={"objects":[mem],"property":"member"}))
            self.__nvars += 1
        
        for i, node in enumerate(self.structure.nodes):
            self.vars['node'].append(sop.BinaryVariable("z({0:.0f})".format(i),
                                                          target={"objects":[node],"property":"node"}))
            self.__nvars += 1
    
        if self.groups is not None:
            """ In case of member grouping, add corresponding binary variables """
            for g, group in enumerate(self.groups):
                for p in group['profiles']:
                    self.vars['group'].append(sop.BinaryVariable("w({0:.0f},{1:0f})".format(g,p)))
                    self.__nvars += 1
    
    
    def gather_node_data(self):
        """ Collects needed topological data for nodes.
            i) which nodes are supported
            ii) which members are connected to each node
        """
        for i, mem in self.structure.members.items():        
            for n in mem.nodes.keys():
                self.nodes[n]['members'].append(i)
        
        for i, supp in self.structure.supports.items():            
            self.nodes[supp.node_id]['supported'] = len(supp.dofs)
    
    def add_topology_constraints(self):
        """ Adds the following constraints:
            1) members attached to vanishing node vanish
            2) minimum number of members meeting at existing node            
            3) necessary condition for kinematic stability
            4) minimum number of supported nodes
            5) constraints for groups (optional)
        """
        nx = self.nvars
        
        Ystart = self.start_ndx('member')
        Zstart = self.start_ndx('node')
        # 1) Members attached to vanishing node must vanish
        for j, node in enumerate(self.nodes):            
            for m in node['members']:
                a = np.zeros(nx)
                a[Zstart+j] = 1
                a[Ystart+m] = -1
                self.milp.add(sop.LinearConstraint(a,0,con_type="<",
                                                   name="Y[{0:.0f}] <= Z[{1:.0f}]".format(m,j)))
        
        # 2) Minimum number of members meeting at existing node
        for j, node in enumerate(self.nodes):            
            a = np.zeros(nx)
            mndx = [m + Ystart for m in node['members']]
            a[mndx] = -1
            if node['supported']:
                a[Zstart+j] = 1
            else:
                a[Zstart+j] = 2
            
            self.milp.add(sop.LinearConstraint(a,0,con_type="<",
                                                   name="Sum Y[i] => C*Z[{0:.0f}]".format(j)))
        
        # 3) Necessary condition for kinematic stability
        a = np.zeros(nx)
        a[[Ystart + i for i in range(self.nmem)]] = -1
        ND = 2
        a[[Zstart + i for i in range(self.nnodes)]] = ND
        
        for j, node in enumerate(self.nodes):
            if node['supported']:
                a[Zstart + j] -= node['supported']
        
        self.milp.add(sop.LinearConstraint(a,0,con_type="<",
                                           name="Kinematic stability (necessary cond.)"))
        
        # 4) Minimum number of supported nodes is 2
        a = np.zeros(nx)
        for j, node in enumerate(self.nodes):
            if node['supported']:
                a[Zstart + j] = -1
        
        self.milp.add(sop.LinearConstraint(a,-2,con_type="<",
                                           name="At least 2 supported nodes"))
        
        # 5) Constraints for groups (if any)
        if self.groups is not None:
            Wstart = self.start_ndx('group')
            dW = 0
            for g, group in enumerate(self.groups):
                # First constraints: sum of group profile variables can be at most 1
                aw = np.zeros(nx)
                aw[[Wstart + dW + i for i in range(len(group['profiles']))]] = 1
                
                self.milp.add(sop.LinearConstraint(aw,1,con_type="<",
                                                   name="SUM W[{0:.0f},j] <= 1".format(g)))
                
                # Second constraints:
                # For all members of the group:
                # y(i,j) <= w(j) for all j
                for mem in group['members']:
                    # Find index of the first member profile variable
                    yndx = self.milp.vars.index(self.vars['profile'][mem][0])
                    
                    for j, p in enumerate(group['profiles']):
                        a = np.zeros(nx)
                        a[Wstart + dW + j] = -1
                        a[yndx+j] = 1
                        
                        self.milp.add(sop.LinearConstraint(a,0,con_type="<",
                                                   name="y({0:.0f},{1:.0f}) <= W[{0:.0f},j]".format(mem,p)))
                
                # Third constraints:
                # SUM y[i,j] >= w[j] for all j
                for j, p in enumerate(group['profiles']):
                    a = np.zeros(nx)
                    a[Wstart + dW + j] = 1
                    
                    for mem in group['members']:
                        yndx = self.milp.vars.index(self.vars['profile'][mem][j])
                        a[yndx] = -1
                    
                    self.milp.add(sop.LinearConstraint(a,0,con_type="<",
                                                   name="SUM y(i,{0:.0f}) >= W[{0:.0f},{1:.0f}]".format(j,g)))
                
                
                dW += len(group['profiles'])
                        
                
            #print(Wstart)
        
                          
    def create_milp(self,problem_type='sizing',constraints={'strength':True,'buckling':True}):
        """ Create OptimizationProblem class instance that contains
            the mixed-integer linear program of the truss.
            
             Parameters:
            ------------
            :param problem_type: 'sizing' or 'topology'
            :param constraints: different kinds of constriaints to be included in the problem
    
            :type problem_type: string
            :type constraints: dict
        """
        dvars = []
        
        # Add profile selection variables
        for y in self.vars['profile']:
            dvars.extend(y)
        
        # For each load case, append the force and displacement variables
        # Design variable vector is then
        # x = [y, N[0] u[0], N[1] u[1], ... N[nlc] N[nlc]]
        for k, disp_vars in enumerate(self.vars['disp']):
            for mem_force_var in self.vars['force']:
                dvars.extend(mem_force_var[k])
            dvars.extend(disp_vars)
        
        if problem_type == 'topology':
            self.add_top_opt_vars()
            dvars.extend(self.vars['member'])
            dvars.extend(self.vars['node'])
            
            if self.groups is not None:
                dvars.extend(self.vars['group'])
            
    
        """ Initialize optimization problem with the created variables """
        self.milp = sop.OptimizationProblem(name=self.name,variables=dvars,structure=self.structure)
    
        """ Objective function: truss weight.
            This is the same for all problem types.
        """
        obj = sop.LinearObjective("Weight", self.weight_obj())
        self.milp.add(obj)
        
        ndof = self.ndof
        nvars = self.nvars
        
        """ Constraints: unique profile selection """
        Aeq = np.zeros([self.nmem,nvars])
        for i, var in enumerate(self.milp.vars):
            """ Find the profile selection variables based on their name """
            s = var.name
            if s[0] != 'y':
                break
            else:
                """ The first number in the name of the variable corresponds to the member index """
                m = re.findall(r'\d+',s)
                Aeq[int(m[0]),i] = 1
        
        if problem_type == 'topology':
            lhs = 0.0
            Yndx = self.start_ndx('member')
        else:
            lhs = 1.0            
        for i, a in enumerate(Aeq):                    
            if problem_type == 'topology':
                a[Yndx + i] = -1
                con_name = "Sum y({0:.0f},j)=Y({0:.0f})".format(i)                
            else:
                con_name = "Sum y({0:.0f},j)=1".format(i)
                
            self.milp.add(sop.LinearConstraint(a,lhs,con_type="=", name=con_name))
        
        if problem_type == 'topology':
            self.gather_node_data()
            self.add_topology_constraints()
        
        """ Constraints: equilibrium equations """
        # Create augmented statics matrix, whose columns correspond to
        # the columns of the original statics matrix copied as many times
        # as there are profile variables for each member.        
        B0 = np.zeros([ndof,nvars])
        B = self.statics_matrix
        
        # TODO! This has to be accommodated for several load cases!
        start = self.nprof_vars
        for i, yvars in enumerate(self.vars['profile']):
            end = start + len(yvars)
            B0[:,start:end] += np.reshape(B[:,i],(ndof,1))
            start = end
                            
        for j, lc in enumerate(self.structure.load_cases):
            # Load vector for the load case, including only unsupported dofs """
            p = np.delete(self.structure.f.global_load_vector(lc),self.structure.f.supp_dofs)*self.scaling['N']
            # a is a placeholder for the coefficient vector of a linear constraint """
            a = np.zeros(nvars)
            for i, (b,pi) in enumerate(zip(B0,p)):
                # The first 'nmem' variables are the cross-sectional areas.
                #    Insert the current row of the statics matrix as the coefficient
                #    vector.            
                self.milp.add(sop.LinearConstraint(b,pi,con_type="=", name="EQ-{0:.0g}-LC{1:.0g}".format(i,lc)))
        
        # Constraints: compatibility and strength
        self.compatibility_constraints(B)
        self.strength_constraints()
                
        
        return lp
    
    def substitute_design(self,x):
        """ Substitutes the design in 'x' to the truss structure """
        
        self.milp.substitute_variables(x)
    
    def plot_design(self,x):
        """ Plots the design represented by the variable vector 'x' """
        self.substitute_design(x)
        self.milp.structure.plot()

    def plot_ground_structure(self):
        """ plots the ground structure with member numbers """
        
        for mem in self.structure.members.values():
            mem.active = True
        
        self.structure.plot(print_text=False)
        
if __name__ == '__main__':
    #t = three_bar_truss(L=3000,F1=-200e3,F2=-250e3,nlc=1)
    t = ten_bar_truss(L=3000,F=200e3)
    
    #lp = truss_lp(t,0,1e8)
    #LPsolver = LP()
    #sol = LPsolver.solve(lp,verb=True)
    
    P = []
    
    hmax = 120
    hmin = 50
    tmin = 4.0
    
    for key, val in shs_profiles.items():
        if val['h'] <= hmax and val['h'] >= hmin and val['t'] >= tmin:
            P.append(SHS(val['h'],val['t']))
    
    """
    member_profiles = [[],[],[]]
    
    member_profiles[0] = [0, 2, 4, 6, 8]
    member_profiles[1] = [3, 5, 7]
    member_profiles[2] = [0, 1, 2, 9, 10, 11]
    """
    groups = []
    groups.append({'members':[6,3],'profiles':[7,8,9,10,11]})
    groups.append({'members':[1,7],'profiles':[1,2,3,4,5,6,7]})
    
    member_profiles = None
    lp = TrussMILP('TrussMILP',t,P,scaling={'N':1e-3,'u':1e-1},ulb=-200,uub=200,member_profiles=member_profiles,groups=groups)
    
    lp.create_milp(problem_type='topology')
    MILPsolver = MILP()
    MILPsolver.solve(lp.milp,verb=True)
    x = MILPsolver.X
    
    lp.plot_design(x)
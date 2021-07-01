# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 17:03:17 2020

Truss optimization

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


class TrussOptimization(sop.OptimizationProblem):
    """
    Class for general structural problem
    
    13.7.2020
    
    Tehtävätyypit:
        'sizing'
        'geometry'
        'topology'
        
    Tehtävät, joita on tarkoitus ratkaista:
        'sizing':
            1) jatkuvina muuttujina vain poikkipinta-alat
            2) diskreetteinä muuttujina vain poikkipinta-alat
            3) profiilit katalogista: MILP -formulointi
            4) profiilit katalogista: NAND -formulointi
            5) SAND-formulointi
    
    Rakenne on lähtökohtaisesti nivelpäisten sauvojen joukko, jossa sauvat
    kantavan vain normaalivoimaa ja joissa siirtyminä ovat vain solmusiirtymät.
    
    Tätä voitaisiin vielä laajentaa ristikoille, jotka analysoidaan palkkielementtien
    avulla.
    
    Pitänee tehdä erilliset tiedostot erilaisille rakenteille:
        'PortalTruss'
        'ArchedTruss'
        'BeamColumnFrame'
    
    jne.

    """

    def __init__(self,
                 name="Truss Optimization Problem",
                 problem_type = "sizing",
                 variables=[], constraints=[], objective=None, \
                 gradient=None, hess=None, structure=None, profiles=None):
        
        """         problem_type="sizing",
                 formulation='NAND',
                 structure=None,
                 var_groups=None,
                 con_groups=None,
                 constraints=None,
                 profiles=None,
                 obj=None):
        """
        """ Constructor 
            Parameters:
            ------------
            :param name: Name of the problem (string)
            :param structure: Structure to be optimized (Frame2D)
            :param var_groups: grouping of design variables
            :param con_groups: constraint groups
            :param constraints: constraints of the problem
            :param profiles: available profiles
            :param obj: objective function
    
            :type name: string
            :type structure: Frame2D member?
            :type var_groups: list or dict?
            :type con_groups: list or dict?
            :type constraints: list or dict?
            :type profiles: list or dict?
            :type obj: callable or string?
    
    
            Variables:
            ----------
            :ivar vars: list of variables
            :ivar var_groups:
            :ivar con_groups:
            :ivar profiles:
            :ivar constraints:
            :ivar
        
        """
        
        """ 'constraints' can be a dict or list of constraint functions """
        if isinstance(constraints,dict):
            constraint_funs = []
        else:
            constraint_funs = constraints
        super().__init__(name,variables, constraint_funs, objective, gradient, hess, structure, profiles)
        
        self._ptype = problem_type
        #self.vars = []
        var_groups = []
        con_groups = []
        self.var_groups = var_groups
        self.con_groups = con_groups
        
        """ If no profile list is given, use the default SHS profiles """
        if profiles is None:
            profiles = list(shs_profiles.keys())
            
        self.profiles = profiles
        
        """ Generate variables.
            By default, index variables are generated. This can be overridden
            by var_groups, which gives the variable data.
        """
        #self.create_variables()
        
        """ Generate constraints 
            By default, cross-section resistance constraints are generated.
            This can be overridden by 'constraints' dict.
        """
        if constraints is None or len(constraints) == 0:
            constraints = {
                'stress': True,
            }
        
        self.constraints = constraints
    
        self.create_constraints(**self.constraints)
        #self.group_constraints()
        if objective is not None:
            self.obj = objective
        else:
            self.create_objective()

    @property
    def available_constraints(self):
        """
        Returns all available constraints and their corresponding keyword
        :return:
        """
        return dict(stress=False,
                    strength=False,
                    buckling=False)
    
    @property
    def available_formulations(self):
        """
        Returns all available formulations and their corresponding keyword
        :return:
        """
        return dict(lp=False,
                    nand=False,
                    sand=False,
                    discrete=False,
                    sizing=False,
                    geometry=False,
                    topology=False)
    
    def create_variables(self,var_type="area",var_cat="continuous",**kwargs):
        """
        Creates design variables
            :param var_type: variable type. Possible values: 'area', 'profile', 'index'
            :param var_cat: variable category. Posstible values: 'continuous', 'discrete', 'integer'
            :param **kwargs: additional parameters to complete the definition, depending on the variable type and category
        """
        self.vars.clear()

        if var_type == "area":
            var_prop = 'A'
            if var_cat == "continuous":
                for key, value in kwargs.items():
                    if key == 'lb':
                        lb = value
                    elif key == 'ub':
                        ub = value
                
                x0 = 0.5*(lb+ub)

        # If no grouping is given, create a group for each member
        if self.var_groups is None:
            self.var_groups = []
            for mem in self.structure.members.values():
                group = {
                    "name": f"{mem.mtype.capitalize()} {mem.mem_id}",
                    "objects": [mem],
                    "value": x0,
                    "profiles": self.profiles,
                    "property": var_prop,
                    "var_type": 'index',
                }
                self.var_groups.append(group)

        # Main loop
        for i, group in enumerate(self.var_groups):
            if group["var_type"] == 'discrete':
                var = sop.DiscreteVariable(
                    name=group['name'],
                    values=group['values'],
                    value=group['value'],
                    target={"property": group["property"],
                            "objects": group["objects"]}
                )
                self.add(var)

            elif group["var_type"] == 'continuous':
                # Multiple properties
                if "properties" in group.keys():
                    bounds = group["bounds"]
                    if len(bounds) != len(group["properties"]):
                        raise ValueError(
                            "There must be same number of bounds as properties!"
                            f"{group['bounds']} != {group['properties']}")
                    for bounds, prop, value in zip(group["bounds"],
                                                   group["properties"],
                                                   group["values"]):
                        lb, ub = bounds
                        var = sop.Variable(
                            name=group['name'],
                            lb=lb,
                            ub=ub,
                            value=value,
                            target={"property": prop,
                                    "objects": group["objects"]}
                        )
                        self.add(var)

                # Single property
                else:
                    var = sop.Variable(
                        name=group['name'],
                        lb=group['lb'],
                        ub=group['ub'],
                        value=group['value'],
                        target={"property": group["property"],
                                "objects": group["objects"]}
                    )
                    self.add(var)

            elif group["var_type"] == 'index':
                var = sop.IndexVariable(
                    name=group['name'],
                    value=group['value'],
                    values=group['values'],
                    target={"property": group["property"],
                            "objects": group["objects"]}
                )
                self.add(var)

            elif group["var_type"] == 'binary':
                # TODO
                pass
            else:
                raise ValueError("var_type must be either 'discrete',"
                                 " 'continuous', 'index' or 'binary")
    
    def create_objective(self,scaling=10.0):
        """ By default, weight of the structure is taken as the
            objective function.
        """
        def obj_fun(x):
            """ Variable substitution is not shown """
            self.substitute_variables(x)
            return self.structure.weight/scaling

        obj = sop.ObjectiveFunction(name="Weight",
                                obj_fun=obj_fun,
                                obj_type='MIN')
        self.add(obj)
        
    def create_constraints(self,
                           stress=False,
                           strength=False,
                           buckling=False,
                           members=None
                           ):
        """
        Creates constraints

        NOTE! The name of these kwargs must match the key in the cons -dict
        created in constraint function creator method
        e.g.
        CORRECT:
        cons = {'compression': comp_fun, 'tension': ten_fun}
        FALSE:
        cons = {'comp_func': comp_fun, 'tension_con': ten_fun}

        :param stress:
        :param strength:
        :param buckling:
        :param members:
        :return:
        """

        if members is None:
            members = self.structure.members.values()
            # Clear constraint list
            self.cons.clear()

        for mem in members:
            if stress:
                stress_cons = self.stress_constraints(mem)
                for key, val in stress_cons.items():                    
                    con = sop.NonLinearConstraint(
                        name=f"{key}: {mem.mem_id}",
                        con_fun=val,
                        con_type='<',
                        fea_required=True,
                        vars=mem)
                    self.add(con)
            
            if strength:
                new_cons = self.strength_constraints(mem)
                for key, val in new_cons.items():                    
                    con = sop.NonLinearConstraint(
                        name=f"{key}: {mem.mem_id}",
                        con_fun=val,
                        con_type='<',
                        fea_required=True,
                        vars=mem)
                    self.add(con)
            
            if buckling:
                new_cons = self.buckling_constraints(mem,buckling)
                for key, val in new_cons.items():
                    con = sop.NonLinearConstraint(
                    name=f"{key}: {mem.mem_id}",
                    con_fun=val,
                    con_type='<',
                    fea_required=True,
                    vars=mem)
                
                self.add(con)
            

    
    def stress_constraints(self, mem, max_stress=None):
        """
        Axial stress constraints

        :param mem: FrameMember object
        :param max_stress: maximum allowable stress

        :return: dict of created functions
        """
        
        if max_stress == None:
            max_stress = mem.fy

        def compression(x):
            return -mem.ned / (max_stress*mem.A) - 1

        def tension(x):            
            return mem.ned / (max_stress*mem.A) - 1

        cons = {
            'compression': compression,
            'tension': tension,
        }

        return cons
    
    def strength_constraints(self, mem):
        """
        Constraints for member strength against axial force

        :param mem: FrameMember object

        :return: dict of created functions
        """
                
        def compression(x):
            self.substitute_variables(x)
            return -mem.ned / mem.NRd - 1

        def tension(x):            
            self.substitute_variables(x)
            return mem.ned / mem.NRd - 1

        cons = {
            'compression': compression,
            'tension': tension,
        }

        return cons

    def buckling_constraints(self, mem, buckling_type="EN 1993"):
        """
        Member buckling constraints
        
        :param mem: FrameMeber object
        :param buckling_type: 'ec3', 'euler', or something else
        
        :return: dict of constraints
        """
        
        def buckling(x):
            self.substitute_variables(x)
            if buckling_type == 'EN 1993':
                NRd = mem.NbRd[0]
            elif buckling_type == 'euler':
                NRd = mem.member.ncrit()[0]
            
            return -mem.ned/NRd - 1
        
        cons = {f"buckling ({buckling_type})": buckling}
        
        return cons
            
            
def truss_lp(truss,Alb=0,Aub=1e9,slb=None,sub=None):
    """ Creates minimum weight LP problem for the truss 
        Design variables:
            1) Cross-sectional areas of the members
            2) Axial forces of the members, for each load combination
        
        Constraints:
            1) Equilibrium equations, for each load combination
            2) Force constraints, tension and compression
            3) Box constraints for cross-sectional areas
    """
        
    """ Create variables """
    dvars = []
    c = []
    
    nmem = len(truss.members)
    
    """ Initialize cross-sectional area and stress bounds """
    if not (isinstance(Alb,np.ndarray) or isinstance(Alb,list)):        
        Alb = Alb*np.ones(nmem)
    
    if not (isinstance(Aub,np.ndarray) or isinstance(Aub,list)):
        Aub = Aub*np.ones(nmem)
        
    if slb is None:
        slb = [mem.fy for i, mem in t.members.items()]
    
    if sub is None:
        sub = [mem.fy for i, mem in t.members.items()]
    
    if not (isinstance(slb,np.ndarray) or isinstance(slb,list)):        
        slb = slb*np.ones(nmem)
    
    if not (isinstance(sub,np.ndarray) or isinstance(sub,list)):
        sub = sub*np.ones(nmem)
    
    
    
    """ Design variables: member areas 
        The assumption here is that each member has its own variable.
        If member grouping is used, this needs to be changed according to the
        grouping.
        
        Also the stress constraints must be changed to take into account groups.
    """
    for i, mem in truss.members.items():
        dvars.append(sop.Variable(name="A{0:.0g}".format(i),lb=Alb[i],ub=Aub[i],
                                  target={"objects":[mem],"property":["A"]},
                                  value=mem.A))
        c.append(mem.length)
    
    """ Design variables: member forces """
    for lc in truss.load_cases:
        for i, mem in truss.members.items():
            dvars.append(sop.Variable(name="N{0:.0g}{1:.0g}".format(lc,i),lb=-1e9,ub=1e9,
                                  target=None,
                                  value=1))
            c.append(0)
            
    
    """ Initialize optimization problem with the created variables """
    lp = sop.OptimizationProblem(name="Truss LP",variables=dvars,structure=truss)
    
    """ Objective function: truss weight """
    obj = sop.LinearObjective("Weight", c)
    lp.add(obj)
    
    """ Constraints: equilibrium equations """
    B = truss.f.statics_matrix()
    ndof = truss.f.nfree_dofs
        
    for j, lc in enumerate(truss.load_cases):
        """ Load vector for the load case, including only unsupported dofs """
        p = np.delete(truss.f.global_load_vector(lc),truss.f.supp_dofs)
        """ a is a placeholder for the coefficient vector of a linear constraint """
        a = np.zeros(lp.nvars())
        for i, (b,pi) in enumerate(zip(B,p)):
            """ The first 'nmem' variables are the cross-sectional areas.
                Insert the current row of the statics matrix as the coefficient
                vector. 
            """
            a[(j+1)*nmem:(j+2)*nmem] = b            
            lp.add(sop.LinearConstraint(a,pi,con_type="=", name="EQ-{0:.0g}-LC{1:.0g}".format(i,lc)))
    
    """ Constraints: stress constraints 
        These are expressed in terms of member forces, i.e.
        -sigma_lb*A[i] <= N[i] <= sigma_ub*A[i]
    """
    for j, lc in enumerate(truss.load_cases):
        for i, mem in truss.members.items():
            a = np.zeros(lp.nvars())
            a[i] = -sub[i]
            a[(j+1)*nmem+i] = 1
            lp.add(sop.LinearConstraint(a,0,con_type="<", name="STR-UB-{0:.0g}-LC{1:.0g}".format(i,lc)))
            
            a[i] = -slb[i]
            a[(j+1)*nmem+i] = -1
            lp.add(sop.LinearConstraint(a,0,con_type="<", name="STR-LB-{0:.0g}-LC{1:.0g}".format(i,lc)))
    
    return lp

def truss_milp(truss,profiles,ulb=-1e3,uub=1e3):
    """ Discrete sizing optimization of trusses using MILP formulation 
        
        It is assumed that all members have the same available set of profiles and
        each member has its own profile, i.e. no profile grouping.
    
        input:
            profiles .. list of available profiles. Profiles are RHS, SHS, IPE, etc. objects.
            ulb, uub .. bounds for nodal displacements
    """
    
    """ Do unit scaling. By default, the system is (N,mm). Perhaps
        (kN,cm) would be better    
    """
    Nscale = 1e-3
    uscale = 1e-1
    Escale = Nscale/(uscale**2)
    Ascale = uscale**2
    
    ulb *= uscale
    uub *= uscale
    
    """ Create variables """
    dvars = []
    c = []
    
    nmem = len(truss.members)   # number of members
    nprof = len(profiles)       # number of profiles
    ndof = truss.f.nfree_dofs   # number of degrees of freedom
    nlc = len(truss.load_cases) # number of load cases
    
    """ Statics matrix """
    B = truss.f.statics_matrix()  
    
    """ Calculate cross-sectional resistance of members and round down.
        By rounding, we try to avoid numerical issues in the computations.
    """
    NRd = [np.floor(p.NRd)*Nscale for p in profiles]
    
    """ Buckling resistance of members. This must be calculated for all profiles """
    NbRd = np.zeros([nprof,nmem])
    
    """ Bounds for force variables """
    Nlb = np.zeros([nprof,nmem])
    Nub = np.zeros([nprof,nmem])

    """ Design variables: member profile selection
        The assumption here is that all members can take all profiles from the list.
        If member grouping is used or if different profile alternatives are used for
        different members, this needs to be changed.
    """
    
    for i, mem in truss.members.items():
        L = mem.length*uscale
        k0 = mem.E/L*Escale
        nlb = np.sum(np.maximum(B[:,i],0)*ulb) + np.sum(np.minimum(B[:,i],0)*uub)
        nub = np.sum(np.minimum(B[:,i],0)*ulb) + np.sum(np.maximum(B[:,i],0)*uub)
        for j, p in enumerate(profiles):                        
            dvars.append(sop.BinaryVariable("y({0:.0f},{1:.0f})".format(i,j),
                                            section=p,
                                            target={"objects":[mem],"property":"profile"}))
            """ Calculate buckling resistance of the member for profile 'p' """
            mem.profile = p
            NbRd[j,i] = np.floor(mem.NbRd[0])*Nscale
            Nlb[j,i] = k0*p.A*nlb*Ascale
            Nub[j,i] = k0*p.A*nub*Ascale
            c.append(p.weight()/uscale*L)
    
  
    """ Design variables: member forces 
        Force variables are N[i,j,k], where i is member, j is profile and k is load case
    
    """
    for lc in truss.load_cases:
        for i, mem in truss.members.items():            
            for j in range(nprof):
                dvars.append(sop.Variable(name="N({0:.0f},{1:.0f},{2:.0f})".format(i,j,lc),lb=Nlb[j,i],ub=Nub[j,i],
                                          target=None,value=1))
                c.append(0)
    
    """ Design variables: nodal displacements 
        u[j,k]
    """
    for k in range(nlc):
        for j in range(ndof):
            dvars.append(sop.Variable(name="u({0:.0f},{1:.0f})".format(j,k),
                                      lb=ulb,ub=uub,target=None,value=0.0))
            c.append(0)
            
    """ Initialize optimization problem with the created variables """
    lp = sop.OptimizationProblem(name="Truss LP",variables=dvars,structure=truss)
    
    """ Objective function: truss weight """
    obj = sop.LinearObjective("Weight", c)
    lp.add(obj)
    
    """ Constraints: unique profile selection """
    Aeq = np.zeros([nmem,lp.nvars()])
    for i, var in enumerate(lp.vars):
        """ Find the profile selection variables based on their name """
        s = var.name
        if s[0] != 'y':
            break
        else:
            """ The first number in the name of the variable corresponds to the member index """
            m = re.findall(r'\d+',s)
            Aeq[int(m[0]),i] = 1
            
    for i, a in enumerate(Aeq):
        lp.add(sop.LinearConstraint(a,1,con_type="=", name="Sum y({0:.0f},j)=1".format(i)))
    
    """ Constraints: equilibrium equations """
    B0 = np.zeros([ndof,lp.nvars()])
    
    for i in range(nmem):
        B0[:,(nmem*nprof + i*nprof):(nmem*nprof + (i+1)*nprof)] += np.reshape(B[:,i],(ndof,1))
        
    
    for j, lc in enumerate(truss.load_cases):
        """ Load vector for the load case, including only unsupported dofs """
        p = np.delete(truss.f.global_load_vector(lc),truss.f.supp_dofs)*Nscale
        """ a is a placeholder for the coefficient vector of a linear constraint """
        a = np.zeros(lp.nvars())
        for i, (b,pi) in enumerate(zip(B0,p)):
            """ The first 'nmem' variables are the cross-sectional areas.
                Insert the current row of the statics matrix as the coefficient
                vector. 
            """
            lp.add(sop.LinearConstraint(b,pi,con_type="=", name="EQ-{0:.0g}-LC{1:.0g}".format(i,lc)))
    
    """ Constraints: compatibility conditions """
    nx = lp.nvars()
    Nndx = nmem*nprof
    for k, lc in enumerate(truss.load_cases):
        for i, mem in truss.members.items():            
            k0 = mem.E/mem.length*Escale/uscale
            for j, p in enumerate(profiles):
                a = np.zeros(nx)
                """ Upper bound """
                a[i*nprof+j] = Nub[j,i]              # y variable
                a[Nndx] = -1  # N variable
                ustart = nmem*nprof + nlc*nmem*nprof + k*ndof
                uend = ustart +ndof
                a[ustart:uend] = k0*p.A*B[:,i]*Ascale                  # u variable
                lp.add(sop.LinearConstraint(a,Nub[j,i],con_type="<", name="COMP-UB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
    
                """ Lower bound """
                a[Nndx] = 1  # N variable
                a[i*nprof+j] = -Nlb[j,i]             # y variable
                a[ustart:uend] *= -1    # u variable
                lp.add(sop.LinearConstraint(a,-Nlb[j,i],con_type="<", name="COMP-LB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
                
                Nndx +=1
    """ Constraints: stress constraints 
        These are expressed in terms of member forces, i.e.
        -NRd[j]*y[i,j] <= N[i,j,k] <= NRd[j]*y[i,j]
    """
    Nndx = nmem*nprof
    for k, lc in enumerate(truss.load_cases):
        for i, mem in truss.members.items():            
            for j, Nrd in enumerate(NRd):
                a = np.zeros(lp.nvars())
                # Tension resistance 
                a[i*nprof+j] = -Nrd # -Nub[j,i] #               # y variable                
                a[Nndx] = 1  # N variable
                lp.add(sop.LinearConstraint(a,0,con_type="<", name="STR-UB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
            
                #a[i] = -slb[i]
                # Compression resistance, including buckling 
                a[Nndx] = -1
                a[i*nprof+j] = -min(Nrd,NbRd[j,i])# Nlb[j,i] # 
                lp.add(sop.LinearConstraint(a,0,con_type="<", name="STR-LB({0:.0f},{1:.0f},{2:.0f})".format(i,j,k)))
    
                Nndx += 1
    
    return lp

if __name__ == '__main__':
    #t = three_bar_truss(L=3000,F1=-200e3,F2=-250e3)
    t = ten_bar_truss(L=3000,F=200e3)
    
    #lp = truss_lp(t,0,1e8)
    #LPsolver = LP()
    #sol = LPsolver.solve(lp,verb=True)
    
    P = []
    
    hmax = 110
    hmin = 50
    tmin = 3.0
    
    for key, val in shs_profiles.items():
        if val['h'] <= hmax and val['h'] >= hmin and val['t'] >= tmin:
            P.append(SHS(val['h'],val['t']))
    
    lp = truss_milp(t,P,ulb=-200,uub=200)
    
    MILPsolver = MILP()
    sol = MILPsolver.solve(lp,verb=True)
    
    
    
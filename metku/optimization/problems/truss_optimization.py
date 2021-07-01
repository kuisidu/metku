# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 17:03:17 2020

Truss optimization

@author: kmela
"""

try:
    from metku.frame2d.frame2d import SteelBeam
    from metku.sections.steel import SteelSection
    import metku.optimization.structopt as sop
    from metku.sections.steel.catalogue import ipe_profiles, shs_profiles, hea_profiles
    from metku.sections.steel import ISection, WISection, RHS
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
            
            
        
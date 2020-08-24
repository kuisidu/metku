# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 17:03:17 2020

Truss optimization

@author: kmela
"""

try:
    from src.frame2d.frame2d import SteelBeam
    from src.sections.steel import SteelSection
    from src.optimization.structopt import *
    from src.sections.steel.catalogue import ipe_profiles, shs_profiles, hea_profiles
    from src.sections.steel import ISection, WISection, RHS
except:
    from frame2d.frame2d import SteelBeam
    from sections.steel import SteelSection
    from optimization.structopt import *
    from sections.steel.catalogue import ipe_profiles, shs_profiles, hea_profiles
    from sections.steel import ISection, WISection, RHS


class TrussOptimization(OptimizationProblem):
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
                 problem_type="sizing",
                 formulation='NAND',
                 structure=None,
                 var_groups=None,
                 con_groups=None,
                 constraints=None,
                 profiles=None,
                 obj=None):
        """ Constructor 
            Parameters:
            ------------
            :param name: Name of the problem (tuple)
            :param structure: Structure to be optimized (Frame2D member?)
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
        super().__init__(name=name,
                         structure=structure)
        
        self._ptype = problem_type
        self.vars = []
        self.var_groups = var_groups
        self.con_groups = con_groups
        if profiles is None:
            profiles = list(ipe_profiles.keys())
        self.profiles = profiles
        
        """ Generate variables.
            By default, index variables are generated. This can be overridden
            by var_groups, which gives the variable data.
        """
        self.create_variables()
        
        """ Generate constraints 
            By default, cross-section resistance constraints are generated.
            This can be overridden by 'constraints' dict.
        """
        if constraints is None:
            constraints = {
                'compression': True,
                'tension': True,
                'shear': True,
                'bending_moment': True
            }
        self.constraints = constraints
        self.create_constraints(**self.constraints)
        self.group_constraints()
        if obj is not None:
            self.obj = obj
        else:
            self.create_objective()

    @property
    def available_constraints(self):
        """
        Returns all available constraints and their corresponding keyword
        :return:
        """
        return dict(joint_geometry_constraints=False,
                    joint_strength_constraints=False,
                    compression=False,
                    tension=False,
                    shear=False,
                    bending_moment=False,
                    buckling_y=False,
                    buckling_z=False,
                    LT_buckling=False,
                    compression_bending_y=False,
                    compression_bending_z=False,
                    deflection_y=None,
                    deflection_x=None,
                    alpha_cr=None)
    
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
    
    def create_variables(self):
        """
        Creates design variables
        """
        self.vars.clear()

        # If no grouping is given, create a group for each member
        if self.var_groups is None:
            self.var_groups = []
            for mem in self.structure.members.values():
                group = {
                    "name": f"{mem.mtype.capitalize()} {mem.mem_id}",
                    "objects": [mem],
                    "value": 0,
                    "profiles": self.profiles,
                    "property": "profile",
                    "var_type": 'index',
                }
                self.var_groups.append(group)

        # Main loop
        for i, group in enumerate(self.var_groups):
            if group["var_type"] == 'discrete':
                var = DiscreteVariable(
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
                        var = Variable(
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
                    var = Variable(
                        name=group['name'],
                        lb=group['lb'],
                        ub=group['ub'],
                        value=group['value'],
                        target={"property": group["property"],
                                "objects": group["objects"]}
                    )
                    self.add(var)

            elif group["var_type"] == 'index':
                var = IndexVariable(
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
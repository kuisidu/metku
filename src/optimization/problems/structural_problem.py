try:
    from src.optimization.structopt import *
    from src.sections.steel.catalogue import ipe_profiles
except:
    from optimization.structopt import *
    from sections.steel.catalogue import ipe_profiles


class StructuralProblem(OptimizationProblem):
    """
    Class for general structural problem

    """
    def __init__(self,
                 name="Structural Optimization Problem",
                 structure=None,
                 var_groups=None,
                 con_groups=None,
                 constraints=None,
                 profiles=None,
                 obj=None):
        super().__init__(name=name,
                         structure=structure)
        self.var_groups = var_groups
        self.con_groups = con_groups
        if profiles is None:
            profiles = list(ipe_profiles.keys())
        self.profiles = profiles
        self.create_variables()
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
        Returns all available constraints anf their corresponding keyword
        :return:
        """
        return dict(compression=False,
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

    def cross_section_constraints(self, mem, elem):
        """
        Creates cross-section constraints

        :param mem: FrameMember object
        :param elem: Element object

        :return: dict of created functions
        """

        def compression(x):
            N = elem.axial_force[0]
            return -N / mem.NRd - 1

        def tension(x):
            N = elem.axial_force[0]
            return N / mem.NRd - 1

        def shear(x):
            V = elem.shear_force[0]
            return abs(V) / mem.VRd - 1

        def bending_moment(x):
            # Moment about y
            # TODO: Moment about z
            M = elem.bending_moment[0]
            return abs(M) / mem.MRd[0] - 1

        cons = {
            'compression': compression,
            'tension': tension,
            'shear': shear,
            'bending_moment': bending_moment
        }

        return cons

    def stability_constraints(self, mem):
        """
        Creates stability constraint functions

        :param mem: FrameMember object
        :return:
        """

        def buckling_y(x):
            return -mem.NEd / mem.NbRd[0] - 1

        def buckling_z(x):
            return -mem.NEd / mem.NbRd[1] - 1

        def com_compression_bending_y(x):
            return mem.check_beamcolumn()[0] - 1

        def com_compression_bending_z(x):
            return mem.check_beamcolumn()[1] - 1

        cons = {
            'buckling_y': buckling_y,
            'buckling_z': buckling_z,
            'compression_bending_y': com_compression_bending_y,
            'compression_bending_z': com_compression_bending_z
        }

        return cons

    def deflection_constraints(self, mem, limit_x=None, limit_y=None):
        """
        Creates deflection constraint functions
        :param mem: FrameMember object
        :return:
        """

        cons = {}

        def deflection_y(x):
            displacements = mem.nodal_displacements.values()
            vals = [abs(l[1]) for l in displacements]
            abs_max = max(vals)
            return abs_max / limit_y - 1

        def deflection_x(x):
            displacements = mem.nodal_displacements.values()
            vals = [abs(l[0]) for l in displacements]
            abs_max = max(vals)
            return abs_max / limit_x - 1

        if limit_x:
            cons['deflection_x'] = deflection_x
        if limit_y:
            cons['deflection_y'] = deflection_y

        return cons

    def cross_section_class_constraints(self, mem):
        """
        Cross-section class constraints
        :param mem: FrameMember -object
        :return:
        """

        def section_class(x):
            # TODO
            return None

        return section_class

    def create_constraints(self,
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
                           section_class=None,
                           alpha_cr=None,
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

        :param compression:
        :param tension:
        :param shear:
        :param bending_moment:
        :param buckling_y:
        :param buckling_z:
        :param LT_buckling:
        :param compression_bending_y:
        :param compression_bending_z:
        :param deflection_y:
        :param deflection_x:
        :param alpha_cr:
        :return:
        """

        if members is None:
            members = self.structure.members.values()
            # Clear constraint list
            self.cons.clear()

        for mem in members:
            # STABILITY CONSTRAINTS
            for i, smem in enumerate(mem.steel_members):
                stability_cons = self.stability_constraints(smem)
                for key, val in stability_cons.items():
                    if eval(key.lower()):
                        con = NonLinearConstraint(
                            name=f"{key}: {mem.mem_id}|{i}",
                            con_fun=val,
                            con_type='<',
                            fea_required=True,
                            vars=mem)
                        self.add(con)

                # CROSS-SECTION CONSTRAINTS
                for eid, elem in mem.elements.items():
                    section_cons = self.cross_section_constraints(mem, elem)
                    for key, val in section_cons.items():
                        if eval(key.lower()):
                            con = NonLinearConstraint(
                                name=f"{key}: {mem.mem_id}|{eid}",
                                con_fun=val,
                                con_type='<',
                                fea_required=True,
                                vars=mem)
                            self.add(con)

            # CROSS-SECTION CLASS CONSTRAINTS
            # TODO

            # DEFLECTION CONSTRAINTS
            def_cons = self.deflection_constraints(mem, deflection_x,
                                                   deflection_y)
            for key, val in def_cons.items():
                if eval(key.lower()):
                    con = NonLinearConstraint(
                        name=f"{key}: {mem.mem_id}",
                        con_fun=val,
                        con_type='<',
                        fea_required=True,
                        vars=mem)
                    self.add(con)

    def create_variables(self):
        """
        Creates design variables
        """
        self.vars = []

        # If no grouping is given, create IndexVariables
        # with all IPE and HEA profiles
        # TODO
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

    def group_constraints(self):

        con_groups = [group for group in self.con_groups if "constraints" in group]

        for group in self.con_groups:
            cons = group['constraints']
            for obj in group['objects']:
                for con in self.cons:
                    if obj in con.vars:
                        # Remove old constraints
                        self.cons.remove(con)
                        del con
                init_cons = self.constraints.copy()
                init_cons.update(cons)
                self.create_constraints(**init_cons, members=group['objects'])



    def create_objective(self):

        def obj_fun(x):
            return self.structure.weight

        obj = ObjectiveFunction(name="Weight",
                                obj_fun=obj_fun,
                                obj_type='MIN')
        self.add(obj)


if __name__ == '__main__':

    try:
        from src.frame2d.frame2d import *
        from src.truss2d import Truss2D
        from src.optimization.solvers import *
    except:
        from frame2d.frame2d import *
        from truss2d import Truss2D
        from optimization.solvers import *

    n_bays = 1
    n_storeys = 1
    L_bay = 24000
    H_storey = 8000

    simple_frame = [n_storeys, n_bays, H_storey, L_bay]

    frame = Frame2D(simple=simple_frame, create_beams=False, supports='fixed')

    simple_truss = dict(
        H0=H_storey,
        H1=1800,
        H2=2400,
        L1=L_bay / 2,
        dx=0,
        n=14
    )

    truss = Truss2D(simple=simple_truss, fem_model=frame.f)
    frame.add(truss)

    for tc in truss.top_chords:
        frame.add(LineLoad(tc, [-30, -30], 'y'))
    for bc in truss.bottom_chords:
        frame.add(LineLoad(bc, [-1, -1], 'y'))
    col1, col2 = frame.columns
    frame.add(LineLoad(col1, [3.5, 3.5], 'x'))
    frame.add(LineLoad(col2, [0.2, 0.2], 'x'))
    truss.generate()
    # frame.calculate()
    truss.f.draw()

    TC_group = {
        'name': 'TopChords',
        'var_type': 'index',
        'value': 10,
        'values': list(rhs_profiles.keys()),
        'property': 'profile',
        'objects': truss.top_chords,
        'constraints': {
            'buckling_y': True,
            'buckling_z': False,
            'deflection_y': frame.L / 300,
            'deflection_x': frame.H / 300
        }
    }

    BC_group = {
        'name': 'BottomChords',
        'var_type': 'index',
        'value': 10,
        'values': list(rhs_profiles.keys()),
        'property': 'profile',
        'objects': truss.bottom_chords
    }

    WEB_group = {
        'name': 'Webs',
        'var_type': 'index',
        'value': 10,
        'values': list(rhs_profiles.keys()),
        'property': 'profile',
        'objects': list(truss.webs.values())
    }

    COL_group = {
        'name': 'Columns',
        'var_type': 'index',
        'value': 10,
        'values': list(ipe_profiles.keys()),
        # 'properties': ['h', 'b', 'tf', 'tw'],
        'objects': frame.columns
    }


    problem = StructuralProblem(name="TEST",
                                structure=frame,
                                var_groups=[TC_group, BC_group, WEB_group,
                                            COL_group],
                                con_groups=[TC_group],
                                constraints={
                                    'buckling_y': True,
                                    'buckling_z': True,
                                    'compression_bending_y': True,
                                    'compression_bending_z': True,
                                    'compression': True,
                                    'tension': True,
                                    'shear': True,
                                    'deflection_y': frame.L / 200,
                                    'deflection_x': frame.H / 300
                                })
    solver = GA(pop_size=50, mut_rate=0.15)
    x0 = [1 for var in problem.vars]
    fopt, xopt = solver.solve(problem, x0=x0, maxiter=10, plot=True)
    problem(xopt)
    frame.plot_deflection(100)

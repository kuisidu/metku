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

    def joint_geometry_constraints(self, joint):
        """
        Creates joint geometry constraint functions
        :param joint: TrussJoint
        :return: constraint functions in a dict
        """

        con_funs = {}

        if joint.joint_type == 'Y':
            w1 = list(joint.webs.values())[0]

            def b(x):
                """
                Web's breadth constraint
                """
                return w1.cross_section.B / joint.chord.cross_section.B - 1

            def b_min(x):
                """
                Minimum breadth for web 1
                """
                val = 0.35 * joint.chord.cross_section.B
                val /= w1.cross_section.B
                val -= 1
                return val

            def theta_min(x):

                val = 30 / np.degrees(joint.theta1) - 1
                # val = 30 - np.degrees(w1.angle - joint.chord.angle)
                return val

            con_funs = {
                f'Joint {joint.jid} b': b,
                f'Joint {joint.jid} b1_min': b_min,
                f'Joint {joint.jid} theta1_min': theta_min,
            }

        elif joint.joint_type == 'K':
            w1, w2 = joint.webs.values()

            def beta(x):
                """
                Beta constraint
                """
                b1 = w1.cross_section.B
                b2 = w2.cross_section.B
                h1 = w1.cross_section.H
                h2 = w2.cross_section.H
                b0 = joint.chord.cross_section.B

                return (b1 + b2 + h1 + h2) / (4 * b0) - 1
                # return (b1 + b2 + h1 + h2) - (4 * b0)

            def b1b0(x):
                b1 = w1.cross_section.B
                b0 = joint.chord.cross_section.B
                t0 = joint.chord.cross_section.T

                return (0.1 + 0.01 * b0 / t0) / (b1 / b0) - 1

            def b2b0(x):
                b2 = w2.cross_section.B
                b0 = joint.chord.cross_section.B
                t0 = joint.chord.cross_section.T

                return (0.1 + 0.01 * b0 / t0) / (b2 / b0) - 1

            def g_min(x):
                """
                Minimum value for gap
                """
                gap = w1.cross_section.T + w2.cross_section.T
                val = gap - joint.g1
                val /= abs(joint.g1)
                return val

            def g_min_beta(x):
                """
                Minimum value for gap
                """
                beta = joint.rhs_joint.beta()
                val = 0.5 * (1 - beta) * joint.chord.cross_section.B
                val /= abs(joint.g1)
                val -= 1
                return val

            def g_max_beta(x):
                """
                Maximum value for gap
                """
                beta = joint.rhs_joint.beta()
                val = joint.g1
                val /= (1.5 * (1 - beta) * joint.chord.cross_section.B)
                val -= 1
                return val

            def b1_min(x):
                """
                Minimum breadth for web 1
                """
                val = 0.35 * joint.chord.cross_section.B
                val /= w1.cross_section.B
                val -= 1
                return val

            def b1_max(x):
                val = -0.85 * joint.chord.cross_section.B
                val /= w1.cross_section.B
                val += 1
                return val

            def b2_min(x):
                val = 0.35 * joint.chord.cross_section.B
                val /= w2.cross_section.B
                val -= 1
                return val

            def b2_max(x):
                val = -0.85 * joint.chord.cross_section.B
                val /= w2.cross_section.B
                val += 1
                return val

            def theta1_min(x):

                val = 30 / np.degrees(joint.theta1) - 1
                # val = 30 - np.degrees(w1.angle - joint.chord.angle)
                return val

            def theta1_max(x):
                val = np.degrees(joint.theta1) / 90 - 1
                # val = 30 - np.degrees(w1.angle - joint.chord.angle)

                return val

            def theta2_min(x):
                val = 30 / np.degrees(joint.theta2) - 1
                # val = 30 - np.degrees(w2.angle - joint.chord.angle)
                return val

            def theta2_max(x):
                val = np.degrees(joint.theta2) / 90 - 1
                # val = 30 - np.degrees(w2.angle - joint.chord.angle)
                return val

            def e_pos(x):
                return joint.e / ((0.25 * joint.chord.cross_section.H)) - 1

                # return joint.e - ((0.25 * joint.chord.cross_section.H))

            def e_neg(x):
                return -1 - joint.e / (0.55 * joint.chord.cross_section.H)

                # return  - joint.e - (0.55 * joint.chord.cross_section.H)

            con_funs = {
                f'Joint {joint.jid} beta': beta,
                f'Joint {joint.jid} b1b0': b1b0,
                f'Joint {joint.jid} b2b0': b2b0,
                f'Joint {joint.jid} g_min': g_min,
                f'Joint {joint.jid} g_min_beta': g_min_beta,
                f'Joint {joint.jid} g_max_beta': g_max_beta,
                f'Joint {joint.jid} b1_min': b1_min,
                f'Joint {joint.jid} b1_max': b1_max,
                f'Joint {joint.jid} b2_min': b2_min,
                f'Joint {joint.jid} b2_max': b2_max,
                f'Joint {joint.jid} theta1_min': theta1_min,
                f'Joint {joint.jid} theta1_max': theta1_max,
                f'Joint {joint.jid} theta2_min': theta2_min,
                f'Joint {joint.jid} theta2_max': theta2_max,
                f'Joint {joint.jid} e_pos': e_pos,
                f'Joint {joint.jid} e_neg': e_neg}

        elif joint.joint_type == 'KT':
            # TODO
            pass

        return con_funs

    def joint_strength_constraints(self, joint):
        """
        Creates joint strength costraint functions

        :param joint: TrussJoint to be calculated
        :return: constraint functions in a dict
        """
        con_funs = {}

        if joint.joint_type == 'Y':
            w1 = list(joint.webs.values())[0]

            def chord_face_failure(x):
                NRd1 = joint.rhs_joint.chord_face_failure()
                return w1.ned / NRd1 - 1

            def chord_web_buckling(x):
                NRd1 = joint.rhs_joint.chord_web_buckling()
                return -w1.ned / NRd1 - 1

            def brace_failure(x):
                N1Rd = joint.rhs_joint.brace_failure()
                return abs(w1.ned) / N1Rd - 1

            con_funs = {
                f'Joint {joint.jid} chord face failure': chord_face_failure,
                f'Joint {joint.jid} chord web buckling': chord_web_buckling,
                f'Joint {joint.jid} brace failure': brace_failure,
            }


        elif joint.joint_type == 'K':
            w1, w2 = joint.webs.values()

            def chord_face_failure_1(x):

                NRd1, NRd2 = joint.rhs_joint.chord_face_failure()
                return abs(w1.ned) / NRd1 - 1
                # return abs(w1.ned) - NRd1

            def chord_face_failure_2(x):
                NRd1, NRd2 = joint.rhs_joint.chord_face_failure()
                return abs(w2.ned) / NRd2 - 1
                # return abs(w2.ned) - NRd2

            def punching_shear_1(x):
                NRd1 = joint.rhs_joint.punching_shear()[0]
                return abs(w1.ned) / NRd1 - 1
                # return abs(w1.ned) - NRd1

            def punching_shear_2(x):
                NRd2 = joint.rhs_joint.punching_shear()[1]
                return abs(w2.ned) / NRd2 - 1
                # return abs(w2.ned) - NRd2

            def chord_shear_0(x):
                joint.rhs_joint.V0 = joint.V0
                (N1Rd, N2Rd), N0Rd = joint.rhs_joint.chord_shear()

                return joint.N0 / N0Rd - 1
                # return joint.N0 - N0Rd

            def chord_shear_1(x):
                joint.rhs_joint.V0 = joint.V0
                (N1Rd, N2Rd), N0Rd = joint.rhs_joint.chord_shear()

                return abs(w1.ned) / N1Rd - 1

            def chord_shear_2(x):
                joint.rhs_joint.V0 = joint.V0
                (N1Rd, N2Rd), N0Rd = joint.rhs_joint.chord_shear()

                return abs(w2.ned) / N2Rd - 1

            def brace_failure_1(x):
                N1Rd = joint.rhs_joint.brace_failure()[0]
                return abs(w1.ned) / N1Rd - 1
                # return abs(w1.ned) - N1Rd

            def brace_failure_2(x):
                N2Rd = joint.rhs_joint.brace_failure()[1]
                return abs(w2.ned) / N2Rd - 1
                # return abs(w2.ned) - N2Rd

            con_funs = {
                f'Joint {joint.jid} chord face failure 1': chord_face_failure_1,
                f'Joint {joint.jid} chord face failure 2': chord_face_failure_2,
                f'Joint {joint.jid} punching shear 1': punching_shear_1,
                f'Joint {joint.jid} punching shear 2': punching_shear_2,
                f'Joint {joint.jid} chord shear 0': chord_shear_0,
                # f'Joint {joint.jid} chord shear 1': chord_shear_1,
                # f'Joint {joint.jid} chord shear 2': chord_shear_2,
                f'Joint {joint.jid} brace failure 1': brace_failure_1,
                f'Joint {joint.jid} brace failure 2': brace_failure_2,
            }

        elif joint.joint_type == 'KT':
            # TODO
            pass

        return con_funs

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
            return -min(mem.ned) / mem.NbRd[0] - 1

        def buckling_z(x):
            return -min(mem.ned) / mem.NbRd[1] - 1

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

        # JOINTS
        for truss in self.structure.truss:
            for joint in truss.joints.values():
                # JOINT GEOMETRY
                joint_geom_cons = self.joint_geometry_constraints(joint)
                for key, val in joint_geom_cons.items():
                    con = NonLinearConstraint(
                        name=f"{key}:",
                        con_fun=val,
                        con_type='<',
                        fea_required=False)
                    self.add(con)
                # JOINT STRENGTH
                joint_str_cons = self.joint_strength_constraints(joint)
                for key, val in joint_str_cons.items():
                    con = NonLinearConstraint(
                        name=f"{key}:",
                        con_fun=val,
                        con_type='<',
                        fea_required=True)
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

        if self.con_groups is not None:
            con_groups = [group for group in self.con_groups if
                          "constraints" in group]

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
                    self.create_constraints(**init_cons,
                                            members=group['objects'])

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
        dx=400,
        n=16
    )

    truss = Truss2D(simple=simple_truss, fem_model=frame.f)
    frame.add(truss)

    for tc in truss.top_chords:
        frame.add(LineLoad(tc, [-30, -30], 'y'))
    for bc in truss.bottom_chords:
        frame.add(LineLoad(bc, [-1, -1], 'y'))
    col1, col2 = frame.columns
    # frame.add(LineLoad(col1, [3.5, 3.5], 'x'))
    # frame.add(LineLoad(col2, [0.2, 0.2], 'x'))
    frame.generate()
    # truss.generate()
    frame.calculate()
    # truss.f.draw()

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
        'property': 'profile',
        # 'properties': ['h', 'b', 'tf', 'tw'],
        'objects': frame.columns
    }

    # # truss.f.draw()
    frame.to_robot("PROTAL_TEST")
    # frame.bmd(100)

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

    # solver = GA(pop_size=50, mut_rate=0.15)
    # x0 = [1 for var in problem.vars]
    # fopt, xopt = solver.solve(problem, x0=x0, maxiter=10, plot=True)
    # problem(xopt)
    # frame.plot_deflection(100)

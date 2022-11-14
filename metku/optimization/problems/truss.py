import metku.frame2d as f2d
import metku.truss2d as t2d
from metku.optimization.objective_enum import ObjectiveEnum
from metku.optimization.structopt import *
from metku.structures.steel.steel_member import SteelMember


class TrussProblem(OptimizationProblem):
    """
    Class for truss optimization problem
    """

    def __init__(self,
                 name='TrussProblem',
                 displacement_limit: float = 30,
                 variable_groups=None,
                 objective: callable or ObjectiveEnum = ObjectiveEnum.WEIGHT,
                 structure: f2d.Frame2D = None,
                 include_joint_strength: bool = False,
                 include_joint_geometry: bool = False,
                 include_section_strength: bool = True,
                 include_stability: bool = True,
                 include_section_class: bool = True,
                 include_displacement: bool = True,
                 load_id_SLS=LoadIDs.SLS_Charasteristic,
                 load_id_ULS=LoadIDs.ULS,
                 penalty_fun: callable = None):

        super().__init__(name, structure=structure)
        self.variable_groups = variable_groups
        self.objective = objective
        self.include_joint_strength = include_joint_strength
        self.include_joint_geometry = include_joint_geometry
        self.include_section_strength = include_section_strength
        self.include_stability = include_stability
        self.include_section_class = include_section_class
        self.include_displacement = include_displacement
        self.displacement_limit = displacement_limit
        self.load_id_SLS = load_id_SLS
        self.load_id_ULS = load_id_ULS
        # Penalty function
        if penalty_fun is None:
            penalty_fun = lambda *x: 0
        self.penalty_fun = penalty_fun
        # Create variables
        self.create_variables()
        # Create constraints
        self.create_constraints()

        # Assign objective
        match self.objective:
            case ObjectiveEnum.WEIGHT:
                self.obj = self.weight_fun
            case ObjectiveEnum.COST:
                raise NotImplementedError("Cost Objective onot implemented")
            case callable(self.objective):
                self.obj = self.objective
            case _:
                raise ValueError("Objective must be either "
                                 "ObjectiveTypeEnum or callable")

    @property
    def weight_fun(self):
        """
        Function for calculating structure's weight
        :return: ObjectiveFunction
        """
        if self.structure is not None:
            def obj_fun(x):
                return self.structure.weight + self.penalty_fun(x)

            obj = ObjectiveFunction(name="Objective ",
                                    obj_fun=obj_fun)
            obj.problem = self

            return obj

        else:
            raise ValueError("Structure is not defined!")

    def create_variables(self):
        """
        Creates design variables
        """
        self.vars = []
        for vg in self.variable_groups:
            for var in vg.variables:
                self.add(var)

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
                return w1.NEd[self.load_id_ULS] / NRd1 - 1

            def chord_web_buckling(x):
                NRd1 = joint.rhs_joint.chord_web_buckling()
                return -w1.NEd[self.load_id_ULS] / NRd1 - 1

            def brace_failure(x):
                N1Rd = joint.rhs_joint.brace_failure()
                return abs(w1.NEd[self.load_id_ULS]) / N1Rd - 1

            con_funs = {
                f'Joint {joint.jid} chord face failure': chord_face_failure,
                f'Joint {joint.jid} chord web buckling': chord_web_buckling,
                f'Joint {joint.jid} brace failure': brace_failure,
            }


        elif joint.joint_type == 'K':
            w1, w2 = joint.webs.values()

            def chord_face_failure_1(x):

                NRd1, NRd2 = joint.rhs_joint.chord_face_failure()
                return abs(w1.NEd[self.load_id_ULS]) / NRd1 - 1
                # return abs(w1.ned) - NRd1

            def chord_face_failure_2(x):
                NRd1, NRd2 = joint.rhs_joint.chord_face_failure()
                return abs(w2.NEd[self.load_id_ULS]) / NRd2 - 1
                # return abs(w2.ned) - NRd2

            def punching_shear_1(x):
                NRd1 = joint.rhs_joint.punching_shear()[0]
                return abs(w1.NEd[self.load_id_ULS]) / NRd1 - 1
                # return abs(w1.ned) - NRd1

            def punching_shear_2(x):
                NRd2 = joint.rhs_joint.punching_shear()[1]
                return abs(w2.NEd[self.load_id_ULS]) / NRd2 - 1
                # return abs(w2.ned) - NRd2

            def chord_shear_0(x):
                joint.rhs_joint.V0 = joint.V0[self.load_id_ULS]
                (N1Rd, N2Rd), N0Rd = joint.rhs_joint.chord_shear()

                return joint.N0[self.load_id_ULS] / N0Rd - 1
                # return joint.N0 - N0Rd

            def chord_shear_1(x):
                joint.rhs_joint.V0 = joint.V0[self.load_id_ULS]
                (N1Rd, N2Rd), N0Rd = joint.rhs_joint.chord_shear()

                return abs(w1.NEd[self.load_id_ULS]) / N1Rd - 1

            def chord_shear_2(x):
                joint.rhs_joint.V0 = joint.V0[self.load_id_ULS]
                (N1Rd, N2Rd), N0Rd = joint.rhs_joint.chord_shear()

                return abs(w2.NEd[self.load_id_ULS]) / N2Rd - 1

            def brace_failure_1(x):
                N1Rd = joint.rhs_joint.brace_failure()[0]
                return abs(w1.NEd[self.load_id_ULS]) / N1Rd - 1
                # return abs(w1.ned) - N1Rd

            def brace_failure_2(x):
                N2Rd = joint.rhs_joint.brace_failure()[1]
                return abs(w2.NEd[self.load_id_ULS]) / N2Rd - 1
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

    def cross_section_constraints(self, mem, elem, load_id=LoadIDs.ULS):
        """
        Creates cross-section constraints
        :return:
        """

        def absmax(vals):
            min_val = min(vals)
            max_val = max(vals)
            abs_max = max(abs(min_val), abs(max_val))
            if abs_max == abs(min_val):
                return min_val
            else:
                return max_val

        def compression(x):
            NEd = elem.fint[load_id]['fx']
            return -min(NEd) / mem.NRd - 1

        def tension(x):
            NEd = elem.fint[load_id]['fx']
            return max(NEd) / mem.NRd - 1

        def shear(x):
            VEd = elem.fint[load_id]['fz']
            V = absmax(VEd)
            return abs(V) / mem.VRd - 1

        def bending_moment(x):
            # Moment about y
            # TODO: Moment about z
            MEd = elem.fint[load_id]['my']
            M = absmax(MEd)
            return abs(M) / mem.cross_section.plastic_bending_resistance() - 1

        def shear_moment(x):

            gammaM0 = 1
            VEd = abs(absmax(elem.fint[load_id]['fz']))
            VRd = mem.cross_section.VRd
            if VEd <= 0.5 * VRd:
                rho = 0
            else:
                rho = (2 * VEd / VRd - 1) ** 2
            MRd = mem.MRd
            t = mem.cross_section.T
            Wply = mem.cross_section.Wply
            Av = mem.cross_section.Ashear[0]
            MVRd = min((Wply - rho * Av ** 2 / (8 * t) / gammaM0) * mem.fy, MRd)
            M = abs(absmax(elem.fint[load_id]['my']))
            return abs(M) / MVRd - 1

        def NVM(x):
            """
            Shear + Normal + Moment
            """
            gammaM0 = 1

            NEd = abs(absmax(elem.fint[load_id]['fx']))
            VEd = abs(absmax(elem.fint[load_id]['fz']))
            MEd = abs(absmax(elem.fint[load_id]['my']))
            VRd = mem.cross_section.VRd
            MRd = mem.cross_section.MRd
            A = mem.cross_section.A
            t = mem.cross_section.T
            b = mem.cross_section.B
            if VEd <= 0.5 * VRd:
                rho = 0
            else:
                rho = (2 * VEd / VRd - 1) ** 2
            Wply = mem.cross_section.Wply
            Av = mem.cross_section.Ashear[0]
            Awred = (1 - rho) * (A - 2 * b * t)
            Atotred = A - rho * (A - 2 * b * t)
            NVRd = Atotred * mem.fy / gammaM0
            MVRd = min((Wply - rho * Av ** 2 / (8 * t) / gammaM0) * mem.fy,
                       MRd)
            nV = NEd / NVRd
            aV = min(Awred / Atotred, 0.5)

            if NEd <= 0.25 * NVRd and NEd <= (0.5 * Awred * mem.fy / gammaM0):
                MNVRd = MVRd
            else:
                MNVRd = MVRd * (1 - nV) / (1 - 0.5 * aV)

            return abs(MEd) / MNVRd - 1

        # Cross-Section constraints
        return {
            'Ncd': compression,
            'Ntd': tension,
            'VEd': shear,
            'MEd Y': bending_moment,
            'VEd + MEd': shear_moment,
            'NEd + VEd + MEd': NVM
        }

    def stability_constraints(self, mem: f2d.FrameMember, smem: SteelMember,
                              load_id=LoadIDs.ULS) -> dict[str: callable]:
        """

        :param mem: FrameMember -object
        :param smem: SteelMember -object
        :param section_class:
        :param load_id:
        :return:
        """
        def buckling_y(x):
            return - mem.NEd[load_id] / mem.NbRd[0] - 1

        def buckling_z(x):
            return - mem.NEd[load_id] / mem.NbRd[1] - 1

        def com_compression_bending_y(x):
            section_class = mem.cross_section.section_class()
            return smem.check_beamcolumn(section_class=section_class)[0] - 1

        def com_compression_bending_z(x):
            section_class = mem.cross_section.section_class()
            return smem.check_beamcolumn(section_class=section_class)[1] - 1

        # ADD ANY MISSING CONSTRAINTS HERE
        def my_own_stability_con_fun(x):
            ...

        return {
            'Buckling Y': buckling_y,
            'Buckling Z': buckling_z,
            'Stab. NEd + MEd Y': com_compression_bending_y,
            'Stab. NEd + MEd Z': com_compression_bending_z,
            # 'My stability test': my_own_stability_con_fun
        }

    def deflection_constraints(self, truss: t2d.Truss2D):
        """
        Creates deflection constraint functions
        :param mem:
        :return:
        """

        def disp_fun(x):
            truss.calculate(load_id=self.load_id_SLS)
            max_deflection = truss.get_deflection(LoadIDs.SLS_Charasteristic)
            # To make sure correct values are used in stability design
            truss.calculate(load_id=self.load_id_ULS)
            return max_deflection - self.displacement_limit

        return disp_fun

    def cross_section_class_constraints(self, mem):
        """
        Cross-section class constraints
        :param mem:
        :return:
        """
        # SECTION IN PURE COMPRESSION
        def section_class(x):
            b_ = mem.cross_section.H - 2 * mem.cross_section.T - 2 * mem.cross_section.R
            c_t = b_ / mem.cross_section.T
            return c_t / (38 * mem.cross_section.epsilon) - 1

        return section_class

    def create_constraints(self):
        """
        Cretes constraints and saves them to self.cons
        """
        # Initialize cons as an empty list
        self.cons = []
        # MEMBER AND CROSS-SECTION CONSTRAINTS
        for mem in self.structure.members.values():
            # STABILITY CONSTRAINTS
            if self.include_stability:
                for j, smem in enumerate(mem.members):
                    stability_con_dict = self.stability_constraints(mem, smem)
                    for stab_con_name, stab_con_fun in stability_con_dict.items():
                        if "MEd" in stab_con_name and "chord" in mem.mtype:
                            self.non_linear_constraint(name=f"{stab_con_name} {mem.mtype} {mem.mem_id}|{j}",
                                                       con_fun=stab_con_fun,
                                                       fea_required=True)
                        elif "Buckling" in stab_con_name:
                            self.non_linear_constraint(name=f"{stab_con_name} {mem.mtype} {mem.mem_id}|{j}",
                                                       con_fun=stab_con_fun,
                                                       fea_required=True)
            # CROSS-SECTION CLASS CONSTRAINT
            if self.include_section_class:
                sect_fun = self.cross_section_class_constraints(mem)
                sect_con = self.non_linear_constraint(sect_fun,
                                                      name=f"Section Class {mem.mtype} {mem.mem_id} ")
            # CROSS-SECTION STRENGTH
            if self.include_section_strength:
                for i, elem in enumerate(mem.elements.values()):
                    con_fun_dict = self.cross_section_constraints(mem, elem)

                    for con_name, con_fun in con_fun_dict.items():
                        self.non_linear_constraint(con_fun, name=f"{con_name}: {mem.mtype} {mem.mem_id}|{i}", fea_required=True)
                    # Last element's end node
                    if i == len(mem.elements) - 1 and mem.mtype != 'web':
                        con_fun_dict = self.cross_section_constraints(mem, elem)

                        for con_name, con_fun in con_fun_dict.items():
                            self.non_linear_constraint(con_fun, name=f"{con_name}: {mem.mtype} {mem.mem_id}|{i + 1}")
        # DISPLACEMENT
        if self.include_displacement:
            disp_fun = self.deflection_constraints(self.structure)
            disp_con = self.non_linear_constraint(
                con_fun=disp_fun, name=" Truss Deflection")
            self.add(disp_con)
        # JOINT CONS
        if self.include_joint_geometry or self.include_joint_strength:
            for joint in self.structure.joints.values():
                # Geometry Constraints
                if self.include_joint_geometry:
                    con_funs = self.joint_geometry_constraints(joint)
                    for name, con_fun in con_funs.items():
                        self.non_linear_constraint(con_fun=con_fun,
                                                   name=name)
                # Strength Constraints
                if self.include_joint_strength:
                    strength_con_funs = self.joint_strength_constraints(joint)
                    for name, con_fun in strength_con_funs.items():
                        self.non_linear_constraint(con_fun=con_fun,
                                                   name=name)


def create_truss(h1: float,
                 h2: float,
                 L: float,
                 n: int,
                 dx: float,
                 qd: float = 0,
                 qk: float = 0) -> t2d.Truss2D:
    simple = dict(
        H0=0,
        H1=h2,
        H2=h1,
        L1=L / 2,
        n=n,
        dx=dx
    )
    truss = t2d.Truss2D(simple=simple)
    for mem in truss.members.values():
        if mem.mtype == 'top_chord':
            truss.add(f2d.LineLoad(mem, [qd, qd], 'y', load_id=LoadIDs.ULS))
            truss.add(f2d.LineLoad(mem, [qk, qk], 'y', load_id=LoadIDs.SLS_Charasteristic))
    truss.add(f2d.XYHingedSupport([0, truss.H1]))
    truss.add(f2d.YHingedSupport([truss.L, truss.H1]))
    truss.generate()
    truss.calculate('ALL')
    return truss


if __name__ == '__main__':
    from metku.optimization.variable_group import VariableGroup, VariableTypeEnum
    from metku.sections.steel.catalogue import shs_profiles
    from metku.optimization.solvers import SLP

    # THESE ARE NOT SORTED VALUES
    SHS_PROFILE_NAMES = list(shs_profiles.keys())

    truss = create_truss(h1=2000, h2=1200, L=20_000, dx=1200, n=16, qd=-25, qk=-10)

    top_chords = truss.top_chords
    bottom_chords = truss.bottom_chords
    webs = list(truss.webs.values())

    tc_sections = [mem.cross_section for mem in top_chords]
    bc_sections = [mem.cross_section for mem in bottom_chords]
    web_sections = [mem.cross_section for mem in webs]

    tc_group = VariableGroup(name="Top Chord",
                             var_type=VariableTypeEnum.INDEX,
                             attribute="profile",
                             objects=top_chords,
                             values=SHS_PROFILE_NAMES)

    bc_group = VariableGroup(name="Bottom Chord",
                             var_type=VariableTypeEnum.INDEX,
                             attribute="profile",
                             objects=bottom_chords,
                             values=SHS_PROFILE_NAMES)

    web_group = VariableGroup(name="Webs",
                              var_type=VariableTypeEnum.INDEX,
                              attribute="profile",
                              objects=webs,
                              values=SHS_PROFILE_NAMES)

    tc_HBT = VariableGroup(name="Top Chords",
                           var_type=VariableTypeEnum.CONTINUOUS,
                           attributes=[["H", "B"], "T"],
                           objects=tc_sections,
                           lower_bounds=[100, 3],
                           upper_bounds=[300, 15])
    bc_HBT = VariableGroup(name="Bottom Chords",
                           var_type=VariableTypeEnum.CONTINUOUS,
                           attributes=[["H", "B"], "T"],
                           objects=bc_sections,
                           lower_bounds=[100, 3],
                           upper_bounds=[300, 15])
    web_HBT = VariableGroup(name="Webs",
                            var_type=VariableTypeEnum.CONTINUOUS,
                            attributes=[["H", "B"], "T"],
                            objects=web_sections,
                            lower_bounds=[100, 3],
                            upper_bounds=[300, 15])

    groups = [tc_HBT, bc_HBT]

    # Number of webs
    nw = int(len(webs) / 2)
    for i in range(nw):
        w1 = webs[2 * i]
        w2 = webs[2 * i + 1]
        WEB_group = VariableGroup(name="Webs",
                                  var_type=VariableTypeEnum.CONTINUOUS,
                                  attributes=[["H", "B"], "T"],
                                  objects=[w1.cross_section, w2.cross_section],
                                  lower_bounds=[100, 3],
                                  upper_bounds=[300, 15])
        groups.append(WEB_group)

    problem = TrussProblem(name="Truss Optimization",
                           structure=truss,
                           displacement_limit=50,  # mm
                           include_joint_strength=False,
                           include_joint_geometry=False,
                           include_stability=True,
                           include_displacement=True,
                           include_section_class=True,
                           include_section_strength=True,
                           variable_groups=groups)

    x0 = [var.ub for var in problem.vars]

    solver = SLP([0.05, 0.05])
    fopt, xopt = solver.solve(problem, x0=x0, maxiter=50, verb=True)
    problem(xopt)
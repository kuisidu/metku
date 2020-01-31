"""
    Cost optimization of end plate connections
    
    @author: Kristo Mela
    @date: 31.1.2020
"""

import math

try:
    from src.frame2d.frame2d import SteelBeam
    from src.sections.steel import SteelSection
    from src.optimization.structopt import *
    from src.sections.steel.catalogue import ipe_profiles, shs_profiles, hea_profiles
    from src.sections.steel import ISection, WISection, RHS
except:
    #from optimization.structopt import 
    import optimization.structopt as sopt
    from structures.steel.plates import THICKNESSES


class EndPlateOpt(sopt.OptimizationProblem):
    """
    Class cost optimization of end plate connections

    """

    def __init__(self,
                 name="End-plate Optimization Problem",
                 structure=None,                 
                 constraints={'geometry':True},
                 plate_thicknesses=THICKNESSES,
                 obj=None):
        """ Constructor
            input:
                structure .. EndPlateJoint class object
                constraints .. list of constraints
                plate_thicknesses .. list of available plate thicknesses
                
        """
        super().__init__(name=name,
                         structure=structure)        
        self.plate_thicknesses = plate_thicknesses
        self.create_variables()
        if constraints is None:
            constraints = {
                'geometry': True,
                'bending_resistance': True,
                'shear_resistance': True,
                'stiffness': True
            }
                        
        self.constraints = constraints
        self.create_constraints(**self.constraints)

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
        return dict(geometry=False,
                    bending_resistance=False,
                    shear_resistance=False,
                    stiffness=False)

    def vertical_row_distance_constraints(self, i):
        """
        Creates a constraint that states that two adjacent
        bolt rows need to be at least 2.2*d0 apart from each other
        vertically.
        :param i: index of the row
        :return: constraint functions in a dict
        """

        con_funs = {}
        d0 = self.structure.bolt_rows[0].bolt.d0
        """ Condition: z[i+1]-z[i] >= 2.2*d0 """

        def p(x):
            return self.structure.bolt_rows[i+1].z_top - self.structure.bolt_rows[i].z_top - 2.2*d0

        con_funs = {f'Row {i+2} - Row{i+1} >= 2.2*d0': p}
              
        con_name = f'Row {i+2} - Row{i+1} >= 2.2*d0'
        con_fun = p
        
        return con_name, con_fun
    

    def create_constraints(self,
                           geometry=False,
                           bending_resistance=False,
                           shear_resistance=False,
                           stiffness=False
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

        :param geometry:        
        :return:
        """
        
        if geometry:
            for i in range(self.structure.nrows-1):
                p_name, p_fun = self.vertical_row_distance_constraints(i)
                con = sopt.NonLinearConstraint(name=p_name,
                                               con_fun=p_fun,
                                               con_type='>',
                                               fea_required=False)                
                            #vars=mem)
                self.add(con)

    def create_variables(self):
        """
        Creates design variables
        
        By default, the design variables are:
            tp .. end plate thickness
            x .. horizontal distance of all bolts from center axis of the beam
            zi .. vertical distance of bolt row i from top edge of the plate
            bp .. width of the end plate
            
        Alternatively, the zi could be replaced by:
            z1 .. distance of the top row from the top of the plate
            pi .. distance of adjacent rows
            -> this could be set to all rows such that the distance between
               selected rows would be equal. If zi variables are used,
               equality constraints need to be added to control the distance
               between rows.
        
        """
        self.vars = []

        """ End plate thickness """
        tp_var = sopt.Variable("tp",
                               lb=min(THICKNESSES),
                               ub=max(THICKNESSES),
                               target={"property":"t","objects":[self.structure.end_plate]},
                               value=self.structure.tp)
        
        self.add(tp_var)
        
        """ Distance from bolt rows from axis of the beam """
        d0 = self.structure.bolt_rows[0].bolt.d0
        dw = self.structure.bolt_rows[0].bolt.dw
        hb = self.structure.beam.h
        tb = self.structure.beam.tw
        bc = self.structure.col.b
        tc = self.structure.col.tw
        rc = self.structure.col.r
        af = self.structure.weld_f
        aw = self.structure.weld_w
        x_min = max(0.5*tb + math.sqrt(2)*aw + 0.5*dw,
                    0.5*tc + rc + 0.5*dw,
                    1.2*d0)
        x_max = min(0.5*bc-1.2*d0,
                    0.5*bc-0.5*dw)
        
        axis_x_var = sopt.Variable("x",
                                   lb=x_min,
                                   ub=x_max,
                                   target={"property":"xbolts","objects":[self.structure]},
                                   value=0.5*self.structure.bp-self.structure.ebolts)
        
        self.add(axis_x_var)
        
        """ Vertical distance of bolt rows """
        zlb = max(1.2*d0,0.5*dw)
        zub = hb
        for i, row in enumerate(self.structure.bolt_rows):
            self.add(sopt.Variable("z"+f"{i+1}",lb=zlb, ub=zub, target={"property":'z_top','objects':[row]},
                              value=0.5*hb+self.structure.etop-row.z))   
            
    def create_objective(self):

        def obj_fun(x):
            return self.structure.cost()
        
        obj = sopt.ObjectiveFunction(name="Cost",
                                     obj_fun=obj_fun,
                                     obj_type='MIN')
        self.add(obj)


if __name__ == '__main__':

    import structures.steel.end_plate_joint as ep    
    #from structures.steel.end_plate_joint import example_1
    
    conn = ep.example_1()

    problem = EndPlateOpt(name="TEST",
                          structure=conn,                                
                          constraints={
                                    'geometry': True,                                    
                                })

    # solver = GA(pop_size=50, mut_rate=0.15)
    # x0 = [1 for var in problem.vars]
    # fopt, xopt = solver.solve(problem, x0=x0, maxiter=10, plot=True)
    # problem(xopt)
    # frame.plot_deflection(100)

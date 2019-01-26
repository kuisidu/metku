# Imported libraries
from tables_and_tuples import *
from cost_data import *

# Cost calculations according to thesis:
# Feature-Based Costing Method for Skeletal Steel
# Structures based on the Process Approach, Jaakko Haapio, 2012

# Additional cost calculation functions according to Diaz et al.
# Optimum design of semi-rigid connections using metamodels, 2012


# Cost of material (steel parts, bolts, nuts and washers)
class MaterialCost:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self):
        self.plate_unit_cost = 0.0
        self.steel_grade_add_on = 0.0
        self.thickness_add_on = 0.0
        self.quantity_add_on = 0.0
        self.quantity_per_thickness_add_on = 1.0
        self.inspection_add_on = 1.0

        self.plate = 0.0                                # Material cost of plate
        self.bolts = 0.0                                # Material cost of bolts, nuts and washers'

        self.total = 0.0

    # Chapter 4.5.2 [Haapio]
    def cost_of_plate(self, cn, sd):

        b_p = sd.plate_width                    # [mm]
        h_p = sd.plate_height                   # [mm]
        t_p = sd.plate_thickness                # [mm]
        volume = b_p*h_p*t_p                    # [mm^3]
        material = sd.plate_material            # [S***]
        cost_function = cn.cost_function

        if cost_function == 'Diaz':
            k_s = Diaz['k_s']                   # [e/kg]
            rho = Diaz['rho']*1.0e-9            # [kg/mm^3]
            self.plate = k_s*volume*rho
        else:
            rho = mat[material]["rho"]                                      # Material density [kg/m^3]
            self.plate_unit_cost = c_SMBPl*(1.0/1000.0)                     # [e/kg]
            self.steel_grade_add_on = c_SMG(material)*(1.0/1000.0)          # [e/kg]
            self.thickness_add_on = c_SMT(t_p)*(1.0/1000.0)                 # [e/kg]

            # For whole structure, used add-on 0.0 for connection cost calculations
            self.quantity_add_on = 0.0                                      # [e/kg]
            self.quantity_per_thickness_add_on = 0.0                        # [e/kg]
            self.inspection_add_on = 0.0                                    # [e/kg]

            # Material cost of plate [e/kg]
            self.plate = volume*rho*(self.plate_unit_cost + self.steel_grade_add_on + self.thickness_add_on +
                                     self.quantity_add_on + self.quantity_per_thickness_add_on + self.inspection_add_on)


def bolt_cost(row, cost_function='Haapio'):

    n_bolts = row.n_row             # Number of bolts in row [pcs]

    d = 0
    if row.d in bolt_sizes:
        d = row.d                       # Nominal bolt size [mm]
    else:
        for i in range(0, len(bolt_sizes)):
            if bolt_sizes[i] > row.d and i != 1:
                d = bolt_sizes[i]
            elif bolt_sizes[i] > row.d and i > 1:
                d = bolt_sizes[i - 1]
            else:
                d = bolt_sizes[-1]
    L = row.L                       # Bolt length [mm]

    if cost_function == 'Diaz':
        k_b_mt = 3.076*(d*0.1)**2.0 - 7.373*(d*0.1) + 4.62      # Bolt unit cost + nut and two washers
        row.material_cost = n_bolts*k_b_mt
    else:
        if row.thread_in_shear_plane:
            row.bolt_unit_cost = bolt_unit_cost_full_thread[L][d]               # [e/pcs]
        else:
            row.bolt_unit_cost = bolt_unit_cost_part_thread[L][d]               # [e/pcs]
        if row.bolt_unit_cost == "N/A":
            row.bolt_unit_cost = bolt_unit_cost_part_thread['max']
            print("Size-length (M{0}-{1}) combination not in cost data tables, max bolt unit cost {2} used!"
                  .format(d, L, row.bolt_unit_cost))
        row.nut_unit_cost = nut_unit_cost[d]                                    # [e/pcs]
        if row.d_washer != 0.0:
            row.washer_unit_cost = washer_unit_cost[d]                          # [e/pcs]
        else:
            row.washer_unit_cost = 0.0

        # Total cost of bolts, nuts and washers
        row.material_cost = n_bolts*(row.bolt_unit_cost + row.nut_unit_cost + row.washer_unit_cost)


# Annual cost, Haapio, eq (2)
# Uniform series end-of-period installment [e/a]
def annual_cost(P, I, n):
    # P = a sum of money invested in the initial year less resale value [e]
    # I = interest rate (cost of capital) [%]
    # n = time, the number of units (years) over which interest accumulates

    return P*((I*(1.0 + I)**n)/((1.0 + I)**n - 1.0))       # [e/a]


# General form for cost of cost center
class CostCenter:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, T_Nk=0.0, T_Pk=0.0, c_Lk=0.0, c_Eqk=0.0, c_Mk=0.0, c_REk=0.0, c_Sek=0.0, c_Ck=0.0,
                 c_Enk=0.0, C_Ck=0.0, u_k=1.0):
        self.T_Nk = T_Nk
        self.T_Pk = T_Pk
        self.c_Lk = c_Lk
        self.c_Eqk = c_Eqk
        self.c_Mk = c_Mk
        self.c_REk = c_REk
        self.c_Sek = c_Sek
        self.c_Ck = c_Ck
        self.c_Enk = c_Enk
        self.C_Ck = C_Ck
        self.u_k = u_k

        # Total cost of cost center [e]
        self.cost = (T_Nk + T_Pk)*(c_Lk + c_Eqk + c_Mk + c_REk + c_Sek)/u_k + T_Pk*(c_Ck + c_Enk) + C_Ck

    def info(self):
        print(self.T_Nk)
        print(self.T_Pk)
        print(self.c_Lk)
        print(self.c_Eqk)
        print(self.c_Mk)
        print(self.c_REk)
        print(self.c_Sek)
        print(self.c_Ck)
        print(self.c_Enk)
        print(self.C_Ck)
        print(self.u_k)


# Cost center for blasting
class BlastingCost(CostCenter):

    def __init__(self, cn, sd):

        cost_function = cn.cost_function
        L_b = max(sd.plate_width, sd.plate_height)      # [mm] (max dimension of part to be blasted)

        if cost_function == 'Diaz':
            self.cost = 0.0
        else:
            L = 40.0                                # Blasting workspace length [m]
            B = 10.0                                # Blasting workspace width [m]
            Area = L*B                              # Workspace area [m^2]

            # Non-productive time
            T_Nk = 0.0

            # Productive time
            v_c = 3000.0                            # Conveyor speer [mm/m]
            T_Pk = L_b/v_c                          # Total processing time [min]

            # Labour
            c_hour = 27.68                          # Hourly cost [e/h]
            n_worker = 1.0                          # Number of workers [-]
            c_min = c_hour/60.0                     # Minute rate of worker [e/min]
            c_Lk = n_worker*c_min                   # Unit labor cost [e/min]

            # Equipment installment cost [e/min]
            P_invest_Eq = 200000.0                  # Investment of equipment [e]
            P_resale_Eq = 0.0                       # Resale value [e]
            I_Eq = 0.05                             # Interest rate (cost of capital) [%]
            n_Eq = 20.0                             # n = time [years]
            c_Eqk = annual_cost(P_invest_Eq - P_resale_Eq, I_Eq, n_Eq)/work_mins_in_a

            # Maintenance
            c_maintenance_Eq = 1000.0                   # Maintenance cost of equipment [e/a]
            c_Mk = c_maintenance_Eq/work_mins_in_a      # Maintenance cost of equipment [e/min]

            # Real estate installment cost [e/min]
            P_invest_RealS = Area*c_construction    # Investment cost of workspace area [e]
            P_resale_RealS = Area*c_land            # Resale value of workspace [e]
            I_RealS = 0.05                          # Interest rate (cost of capital) [%]
            n_RealS = 50.0                          # n = time (years)
            c_REk = annual_cost(P_invest_RealS - P_resale_RealS, I_RealS, n_RealS)/work_mins_in_a

            # Real estate maintenance cost
            c_maintenance_RealS = Area*RE_maintenance             # Real estate maintenance cost[e]
            c_Sek = c_maintenance_RealS/work_mins_in_a            # Real estate maintenance cost[e/min]

            # Consumables: grains
            grain_consump = 0.042                   # Grain consumption on blasting [kg/min], equipment utilization ratio 1
            grain_unit_cost = 0.50                  # Grain weight price [e/kg]
            c_Ck = grain_consump*grain_unit_cost    # Cost of consumables [e/min]

            # Energy
            P_total = 40                            # Total power consumption using four blasting turbines [kW]
            c_En_unit = 0.1                         # Energy unit cost [e/kWh]
            c_Enk = P_total*(c_En_unit/60.0)        # Cost of energy [e/min]

            C_Ck = 0.0                              # Total cost of non-time-related consumables

            u_k = 1.0                               # Workspace utilization ratio, set to 1

            CostCenter.__init__(self, T_Nk, T_Pk, c_Lk, c_Eqk, c_Mk, c_REk, c_Sek, c_Ck, c_Enk, C_Ck, u_k)


# Cost center for cutting
class CuttingCost(CostCenter):

    def __init__(self, cn, sd, n_holes=0, L_holes=0.0):

        t_p = sd.plate_thickness                            # [mm]
        t_f = cn.column.t_f                                 # [mm]
        d_0 = sd.d_0                                        # [mm]
        n_holes = n_holes                                   # [-]
        L_holes = L_holes                                   # [mm]
        L_plate = 2.0*(sd.plate_width + sd.plate_height)    # [mm]
        cost_function = cn.cost_function

        if cost_function == 'Diaz':
            # End-plate thickness, flange thickness and bolt nominal size in cm in calculations

            # Cost of cutting the end-plate to size
            k_c = Diaz['k_c']                               # Unit cost for cutting [e/min]
            f_c = Diaz['f_c']                               # Cutting factor that increases labor [-]
            L_c = L_plate*0.001                             # Cutting length [m]
            T_c_ex = Diaz['T_c_ex']                                 # Additional cutting time
            T_c = -0.015*(t_p*0.1)**2.0 + 0.421*(t_p*0.1) + 1.43    # Cutting time [min/m]
            plate_cutting = k_c*(f_c*T_c*L_c + T_c_ex)              # Total cost of cutting [e]

            # Cost of oxygen propane mixture used for cutting the end-plate
            M_c_o = 1.645*(t_p*0.1)**2.0 + 56.644*(t_p*0.1) - 6.73      # Material consumption of oxygen [l/m]
            M_c_p = 2.171*(t_p*0.1) + 7.87                              # Material consumption of propane [l/m]
            k_c_m_o = Diaz['k_c_m_o']                           # Unit cost of oxygen [e/l]
            k_c_m_p = Diaz['k_c_m_p']                           # Unit cost of propane [e/l]
            gas_cost = L_c*(k_c_m_o*M_c_o + k_c_m_p*M_c_p)      # Total cost of cutting gas [e]

            # Hole forming cost (drilling)
            V_d = 0.763*(d_0*0.1)**2.0 - 5.720*(d_0*0.1) + 20.96    # Drillig speed in progression [cm/min]
            t = (t_p + t_f)*0.1                                     # Total drilling thickness [cm]
            k_h = Diaz['k_h']                                       # Unit cost of hole forming [e]
            n_h = n_holes                                           # Number of holes [-]
            t_ex = Diaz['t_ex']                                     # Additional drilling path [cm]
            T_h_ex = Diaz['T_h_ex']                                 # Additional hole forming time [min]
            hole_forming = k_h*(n_h*(t + t_ex)/V_d + T_h_ex)        # Total cost of hole forming [e]

            self.cost = plate_cutting + gas_cost + hole_forming     # Total cost of end-plate shaping
        else:
            L = 20.0            # Workspace length [m]
            B = 10.0            # Workspace width [m]
            Area = L*B          # Workspace area [m^2]

            # Hole cutting on column flange taken into account with multiplier 2
            L_cut = 2.0*L_holes + L_plate

            # Non-productive time [min]
            T_Nk = 3.0

            # Productive time [min]
            # Cutting speed according to thickness [mm] of the plate
            if t_p <= 30.0:   # Plasma cutting
                v_cut = (8.9212*t_p**2.0 - 486.87*t_p + 8155.8)     # [mm/min]
            else:   # Flame cutting
                v_cut = -4.1939*t_p + 658.67                        # [mm/min]
            T_Pk = L_cut/v_cut                                      # [min]

            # Labour
            c_hour = 27.68                  # Hourly cost [e/h]
            n_worker = 2.0                  # Number of workers [-]
            c_min = c_hour / 60.0           # Minute rate of worker [e/min]
            c_Lk = n_worker * c_min         # Unit labor cost [e/min]

            # Equipment installment cost [e/min]
            P_invest_Eq = 280000.0          # Investment of equipment [e]
            P_resale_Eq = 0.0               # Resale value [e]
            I_Eq = 0.05                     # Interest rate (cost of capital) [%]
            n_Eq = 20.0                     # n = time [years]
            c_Eqk = annual_cost(P_invest_Eq - P_resale_Eq, I_Eq, n_Eq)/work_mins_in_a

            # Maintenance
            c_maintenance_Eq = 1000.0                       # Maintenance cost of equipment [e/a]
            c_Mk = c_maintenance_Eq/work_mins_in_a        # Maintenance cost of equipment [e/min]

            # Real estate installment cost [e/min]
            P_invest_RealS = Area * c_construction          # Investment cost of workspace area [e]
            P_resale_RealS = Area * c_land                  # Resale value of workspace [e]
            I_RealS = 0.05                                  # Interest rate (cost of capital) [%]
            n_RealS = 50.0                                  # n = time (years)
            c_REk = annual_cost(P_invest_RealS - P_resale_RealS, I_RealS, n_RealS)/work_mins_in_a

            # Real estate maintenance cost
            c_maintenance_RealS = Area * RE_maintenance         # Real estate maintenance cost[e]
            c_Sek = c_maintenance_RealS/work_mins_in_a        # Real estate maintenance cost[e/min]

            # Consumables: Gases, plasma electrodes and nozzles
            if t_p <= 30.0:   # Plasma cutting
                # Total consumables cost [e/min]
                cutting_oxygen = 0.071                                       # [m^3/min]
                oxygen_cost = 0.3                                            # [e/min]
                nozzle_unit_cost = 40                                        # [e]
                nozzle_life = 480                                            # [min]
                nozzle_cost = nozzle_unit_cost/nozzle_life                   # [e/min]
                c_Ck = cutting_oxygen*oxygen_cost + nozzle_cost

                # Plasma cutting, Total power consumption [kW]
                P_total = 72.0
            else:   # Flame cutting
                # Total consumables cost [e/min]
                propane_consump = 0.0063                                        # [m^3/min]
                propane_cost = 18.40                                            # [e/m^3]
                preheat_oxygen_consump = 0.025                                  # [m^3/min]
                cutting_oxygen_consump = 1.0e-5*t_p**2.0 + 0.001*t_p + 0.0224   # [m^3/min]
                oxygen_cost = 4.18                                              # [e/m^3]
                c_Ck = propane_consump*propane_cost + (preheat_oxygen_consump + cutting_oxygen_consump)*oxygen_cost

                # Flame cutting, Total power consumption [kW]
                P_total = 0.0

            c_Enk = P_total*(c_En_unit/60.0)        # Cost of energy [e/min]

            C_Ck = 0.0      # Total cost of non-time-related consumables

            u_k = 1.0       # Workspace utilization ratio, set to 1

            CostCenter.__init__(self, T_Nk, T_Pk, c_Lk, c_Eqk, c_Mk, c_REk, c_Sek, c_Ck, c_Enk, C_Ck, u_k)


class DrillingCost(CostCenter):
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, cn, sd):

        cost_function = cn.cost_function

        if cost_function == 'Diaz':
            self.cost = 0.0
        else:
            L = 35.0            # Workspace length [m]
            B = 15.0            # Workspace width [m]
            Area = L * B        # Workspace area [m^2]

            # Non-productive time of cost center [min], setup and moving
            T_Nk = 0.0
            T_Pk = 0.0
            c_Lk = 0.0
            c_Eqk = 0.0
            c_Mk = 0.0
            c_REk = 0.0
            c_Sek = 0.0
            c_Ck = 0.0
            c_Enk = 0.0
            C_Ck = 0.0
            u_k = 1.0

            CostCenter.__init__(self, T_Nk, T_Pk, c_Lk, c_Eqk, c_Mk, c_REk, c_Sek, c_Ck, c_Enk, C_Ck, u_k)


class PartAssemblyCost(CostCenter):

    def __init__(self, cn, sd):

        d_w = sd.beam.d_w                       # Length of staigth portion of the beam web [mm]
        P = sd.beam.P                           # Perimeter of beam cross-section [mm]
        a_web = sd.weld_w                       # [mm]
        a_flange = sd.weld_f                    # [mm]
        cost_function = cn.cost_function

        if cost_function == 'Diaz':
            # Manufacturing cost for welding the end-plate
            L_f = (2.0*sd.beam.b + 2.0*(sd.beam.b - 2.0*sd.beam.r - sd.beam.t_w))*1.0e-3      # Flange weld length [m]
            L_w = (2.0*(sd.beam.h - 2.0*sd.beam.r - 2.0*sd.beam.t_f))*1.0e-3                  # Web weld length [m]
            T_w_af = 9.03*(a_flange*0.1)**2.0 + 4.68*(a_flange*0.1) - 0.82    # Flange welding time [min/m]
            T_w_aw = 9.03*(a_web*0.1)**2.0 + 4.68*(a_web*0.1) - 0.82          # Web welding time [min/m]
            k_w = Diaz['k_w']                                                 # Unit cost for welding [e/min]
            f_w = Diaz['f_w']                                                 # Welding factor that increases labor [-]
            T_w_ex = Diaz['T_w_ex']                                           # Additional welding time [min]
            welding_cost = k_w*(f_w*(T_w_af*L_f + T_w_aw*L_w) + T_w_ex)       # Total cost of welding operation

            # Cost of the welding consumables
            M_w_mc_af = 0.66*(a_flange*0.1)**2.0 + 0.18*(a_flange*0.1) - 0.10   # Material cost of welding flange
            M_w_mc_aw = 0.66*(a_web*0.1)**2.0 + 0.18*(a_web*0.1) - 0.10         # Material cost of welding web
            k_w_mt = Diaz['k_w_mt']                                             # Unit cost of weld [e/kg]
            welding_materials = k_w_mt*(M_w_mc_af*L_f + M_w_mc_aw*L_w)          # Total cost of welding materials

            self.cost = welding_cost + welding_materials        # Total cost of welding
        else:
            L = 15.0            # Workspace length [m]
            B = 5.0             # Workspace width [m]
            Area = L * B        # Workspace area [m^2]

            # Non-productive time
            T_Nk = 0.0

            # Productive time
            T_PTa = 1.59                    # Tackwelding time [min]
            L_fw_web = 2.0*d_w              # Flange weld length [mm]
            L_fw_flange = P - 2.0*d_w       # Web weld length [mm]
            # Welding time of web weld [min]
            T_Pfw_web = (L_fw_web/1000.0)*(0.4988*a_web**2.0 - 0.0005*a_web + 0.0021)
            # Welding time of flange weld [min]
            T_Pfw_flange = (L_fw_flange/1000.0)*(0.4988*a_flange**2.0 - 0.0005*a_flange + 0.0021)
            T_Pk = T_PTa + T_Pfw_web + T_Pfw_flange         # Total processing time [min]

            # Labour
            c_hour = 27.42              # Hourly cost [e/h]
            n_worker = 1.0              # Number of workers [-]
            c_min = c_hour / 60.0       # Minute rate of worker [e/min]
            c_Lk = n_worker * c_min     # Unit labor cost [e/min]

            # Equipment installment cost [e/min]
            P_invest_Eq = 5000.0        # Investment of equipment [e]
            P_resale_Eq = 0.0           # Resale value [e]
            I_Eq = 0.05                 # Interest rate (cost of capital) [%]
            n_Eq = 10.0                 # n = time [years]
            c_Eqk = annual_cost(P_invest_Eq - P_resale_Eq, I_Eq, n_Eq) / work_mins_in_a

            # Maintenance
            c_maintenance_Eq = 1000.0                       # Maintenance cost of equipment [e/a]
            c_Mk = c_maintenance_Eq / work_mins_in_a        # Maintenance cost of equipment [e/min]

            # Real estate installment cost [e/min]
            P_invest_RealS = Area * c_construction          # Investment cost of workspace area [e]
            P_resale_RealS = Area * c_land                  # Resale value of workspace [e]
            I_RealS = 0.05                                  # Interest rate (cost of capital) [%]
            n_RealS = 50.0                                  # n = time (years)
            c_REk = annual_cost(P_invest_RealS - P_resale_RealS, I_RealS, n_RealS) / work_mins_in_a

            # Real estate maintenance cost
            c_maintenance_RealS = Area * RE_maintenance     # Real estate maintenance cost[e]
            c_Sek = c_maintenance_RealS / work_mins_in_a    # Real estate maintenance cost[e/min]

            # Total cost of non-time-related consumables [e]
            # Consumables: welding wire, shield gas
            c_Ck = 0.0                                                          # Total cost of time related consumables
            C_Ck_web = L_fw_web*a_web**2.0*0.00000785*(1.91 + 4.44)             # Consumables cost for web fillet weld
            C_Ck_flange = L_fw_flange*a_flange**2.0*0.00000785*(1.91 + 4.44)    # Consumables cost for web fillet weld
            C_Ck = C_Ck_web + C_Ck_flange

            # Energy
            P_total = 6.0                               # Power of welding equipment
            c_Enk = P_total * (c_En_unit / 60.0)        # Cost of energy [e/min]

            u_k = 1.0       # Workspace utilization ratio, set to 1

            CostCenter.__init__(self, T_Nk, T_Pk, c_Lk, c_Eqk, c_Mk, c_REk, c_Sek, c_Ck, c_Enk, C_Ck, u_k)


class PostTreatmentCost(CostCenter):

    def __init__(self, cn, sd):

        weld_test_length = sd.beam.P
        test_method = cn.test_method
        P = weld_test_length                        # Perimeter of beam cross-section [mm]
        cost_function = cn.cost_function

        if cost_function == 'Diaz':
            self.cost = 0.0
        else:
            L = 35.0            # Workspace length [m]
            B = 10.0            # Workspace width [m]
            Area = L * B        # Workspace area [m^2]

            # Non-productive time [min]
            T_Nk = 0.0

            # Productive time, P=testing length of welds [mm]
            if test_method == "VT":       # Visual testing
                T_Pk = 0.0*P    # [min]
            elif test_method == "UT":     # Ultrasonic testing
                T_Pk = 0.06*P   # [min]
            elif test_method == "MT":     # Magnetic particle testing
                T_Pk = 0.015*P  # [min]
            else:
                print("Defined testing method \"{}\" not implemented. Used testing time 0.0 min".format(test_method))
                T_Pk = 0.0

            # Labour
            c_hour = 27.32              # Hourly cost [e/h]
            n_worker = 1.0              # Number of workers [-]
            c_min = c_hour / 60.0       # Minute rate of worker [e/min]
            c_Lk = n_worker * c_min     # Unit labor cost [e/min]

            c_Eqk = 0.0         # Equipment installment cost [e/min]
            c_Mk = 0.0          # Maintenance cost of equipment [e/min]

            # Real estate installment cost [e/min]
            P_invest_RealS = Area * c_construction          # Investment cost of workspace area [e]
            P_resale_RealS = Area * c_land                  # Resale value of workspace [e]
            I_RealS = 0.05                                  # Interest rate (cost of capital) [%]
            n_RealS = 50.0                                  # n = time (years)
            c_REk = annual_cost(P_invest_RealS - P_resale_RealS, I_RealS, n_RealS) / work_mins_in_a

            # Real estate maintenance cost
            c_maintenance_RealS = Area * RE_maintenance     # Real estate maintenance cost[e]
            c_Sek = c_maintenance_RealS / work_mins_in_a    # Real estate maintenance cost[e/min]

            c_Ck = 0.0          # Total cost of time related consumables
            C_Ck = 0.0          # Total cost of non-time-related consumables
            c_Enk = 0.0         # Total energy cost

            u_k = 1.0           # Workspace utilization ratio, set to 1

            CostCenter.__init__(self, T_Nk, T_Pk, c_Lk, c_Eqk, c_Mk, c_REk, c_Sek, c_Ck, c_Enk, C_Ck, u_k)


class PaintingCost(CostCenter):
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, cn, sd):

        paint = cn.paint
        cost_function = cn.cost_function

        if cost_function == 'Diaz':
            A_p = 2.0*(sd.plate_width*sd.plate_height)      # Painting area [mm^2]

            k_p = Diaz['k_p']           # Unit cost for painting [e/min]
            T_p = Diaz['T_p']           # Time consumption of painting [min/m^2]
            k_p_mt = Diaz['k_p_mt']     # Unit cost of paint [e/l]
            M_p = Diaz['M_p']           # Paint consumption [l/m^2]
            self.cost = A_p*1.0e-6*(k_p*T_p + k_p_mt*M_p)      # Total cost of painting [e]
        else:
            L = 15.0        # Workspace length [m]
            B = 5.0         # Workspace width [m]
            Area = L * B    # Workspace area [m^2]

            # Non-productive time [min], drying
            T_Nk = 0.0

            # Painting area [mm], plate sides and beam cross-section taken into account
            A_p = 2.0*sd.plate_width*sd.plate_height + 2.0*sd.plate_thickness*(sd.plate_width + sd.plate_height) - \
                  sd.beam.A

            # Productive time [min], painting
            if paint == "alkyd":
                T_Pk = (0.513/900000)*A_p
            elif paint == "epoxy":
                T_Pk = (0.481/900000)*A_p
            elif paint == "polyurethane":
                T_Pk = (0.481/900000)*A_p
            elif paint == "acryl":
                T_Pk = (0.481/900000)*A_p
            else:
                T_Pk = 0.0
                print("Values for given paint \"{0}\" not implemented!".format(paint))

            # Labour
            c_hour = 27.59              # Hourly cost [e/h]
            n_worker = 1.0              # Number of workers [-]
            c_min = c_hour / 60.0       # Minute rate of worker [e/min]
            c_Lk = n_worker * c_min     # Unit labor cost [e/min]

            c_Eqk = 0.0     # Equipment installment cost [e/min]
            c_Mk = 0.0      # Maintenance cost of equipment [e/min]

            # Real estate installment cost [e/min]
            P_invest_RealS = Area * c_construction      # Investment cost of workspace area [e]
            P_resale_RealS = Area * c_land              # Resale value of workspace [e]
            I_RealS = 0.05                              # Interest rate (cost of capital) [%]
            n_RealS = 50.0                              # n = time (years)
            c_REk = annual_cost(P_invest_RealS - P_resale_RealS, I_RealS, n_RealS) / work_mins_in_a

            # Real estate maintenance cost
            c_maintenance_RealS = Area * RE_maintenance     # Real estate maintenance cost[e]
            c_Sek = c_maintenance_RealS / work_mins_in_a    # Real estate maintenance cost[e/min]

            # Consumables: Paint
            c_Ck = 0.0      # Total cost of time related consumables

            # Consumables: Paint
            # Total cost of non-time-related consumables [e]
            if paint == "alkyd":
                C_Ck = 3.87e-6*A_p
            elif paint == "epoxy":
                C_Ck = 4.11e-6*A_p
            elif paint == "polyurethane":
                C_Ck = 4.58e-6*A_p
            elif paint == "acryl":
                C_Ck = 5.22e-6*A_p
            else:
                C_Ck = 0.0
                print("Values for given paint \"{0}\" not implemented!".format(paint))

            c_Enk = 0.0     # Total energy cost

            u_k = 1.0       # Workspace utilization ratio, set to 1

            CostCenter.__init__(self, T_Nk, T_Pk, c_Lk, c_Eqk, c_Mk, c_REk, c_Sek, c_Ck, c_Enk, C_Ck, u_k)

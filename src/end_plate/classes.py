# Contains classes and methods for building end-plate object

# Imported libraries
from math import sqrt, pi, ceil

# Libraries for calculating strength and stiffness of connection
from functions import cweb_shear
from functions import cweb_trv_comp
from functions import cweb_trv_tens
from functions import bweb_tens
from functions import fc_comp
from functions import cflange_bending
from functions import end_plate_bending
from functions import groups
from functions import define_welds
from functions import par_alfa
from functions import lambda11
from functions import stiffness_coeff
from tables_and_tuples import *
from I_sections import Beam

# Libraries for cost calculations
from cost_data import *
from cost_calculation import MaterialCost
from cost_calculation import bolt_cost
from cost_calculation import BlastingCost
from cost_calculation import CuttingCost
from cost_calculation import DrillingCost
from cost_calculation import PartAssemblyCost
from cost_calculation import PostTreatmentCost
from cost_calculation import PaintingCost

# Libraries for plotting
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as patches


class EndPlate:

    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    warnings = []
    errors = []

    def __init__(self, column_config, column_end_length=0.0, sway="sway", cost_function='Haapio',
                 test_method="VT", paint="alkyd"):
        self.column_config = column_config
        if column_config == "continuous column" and column_end_length == 0.0:
            self.column_end_length = 4500.0
        elif column_config == "column end" and column_end_length > 1000.0:
            self.warnings.append("Column config = 'column end' and column_end_length > 1000.0 [mm]")
            self.column_end_length = column_end_length
        else:
            self.column_end_length = column_end_length
        self.sway = sway
        self.test_method = test_method
        self.paint = paint
        self.joint_config = ""

        # Sway, horizontal stiffening of structure
        if self.sway == "sway":
            self.k_b = 25.0  # Connections stiffness reduces horizontal displacements considerable
        else:
            self.k_b = 8.0  # Bracing reduces horizontal displacement by at least 80%

        self.analysis_method = 'elastic'        # Analysis method of global structure analysis, can be:
                                                # elastic, elastic-plastic or plastic

        self.column = Beam()                    # Initialize column
        self.side = [Side(), Side()]            # Initialize side

        self.V_wp_Rd = 0.0                      # Shear resistance of column web
        self.V_wp_Ed = 0.0                      # Resulting shear force to column web from connection
        self.n_V_wp = 0.0                       # Column shear utilization rate
        # Loading of column over and under the connection
        self.column_loads_over = [1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6]
        self.column_loads_under = [1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6]

        self.shear_block = None                 # Shear block, not implemented yet

        self.cost_function = cost_function
        self.total_cost = 0.0

        self.plot_nro = 0
        self.evaluations = 0

    def add_column(self, size="", material="S355", L=1000,
                   h=1.0e-6, b=1.0e-6, t_f=1.0e-6, t_w=1.0e-6, r=1.0e-6,
                   N=1.0e-6, Vy=1.0e-6, Vz=1.0e-6, Mt=1.0e-6, My=1.0e-6, Mz=1.0e-6):

        if len(self.column.size) != 0:
            print("Column already defined! \n"
                  "Beam before: {0}, Beam now: {1}!".format(self.column.size, size))

        self.column = Beam(size, material, L, h, b, t_f, t_w, r, N, Vy, Vz, Mt, My, Mz)

    def add_column_load(self, side, N=1.0e-6, Vy=1.0e-6, Vz=1.0e-6, Mt=1.0e-6, My=1.0e-6, Mz=1.0e-6):
        if side == 'over':
            self.column_loads_over = [N, Vy, Vz, Mt, My, Mz]
        elif side == 'under':
            self.column_loads_under = [N, Vy, Vz, Mt, My, Mz]
            self.column.loading(N=N, Vy=Vy, Vz=Vz, Mt=Mt, My=My, Mz=Mz)

    def add_beam(self, beam_id, size="", material="S355", L=1000,
                 h=1.0e-6, b=1.0e-6, t_f=1.0e-6, t_w=1.0e-6, r=1.0e-6,
                 N=1.0e-6, Vy=1.0e-6, Vz=1.0e-6, Mt=1.0e-6, My=1.0e-6, Mz=1.0e-6):

        if beam_id != 0 and beam_id != 1:
            print("Beam ID {0} must be 0 or 1!".format(beam_id))
            return
        if len(self.side[beam_id].beam.size) != 0:
            print("Beam {0} already defined! \n"
                  "Beam before: {1}, Beam now: {2}!".format(beam_id, self.side[beam_id].beam, size))

        self.side[beam_id].beam = Beam(size, material, L, h, b, t_f, t_w, r, N, Vy, Vz, Mt, My, Mz)

    def add_end_plate(self, beam_id, b, t, upper_overhang=0.0, lower_overhang=0.0, plate_height=0.0,
                      material="S355", equ=True):

        if beam_id != 0 and beam_id != 1:
            print("Beam ID must be 0 or 1! Given value: {0}".format(beam_id))
            return
        if self.side[beam_id].plate_width != 0.0:
            print("End-plate for beam {0} was already defined! \nEnd-plate redefined!".format(beam_id))

        sd = self.side[beam_id]

        sd.plate_width = float(b)
        sd.plate_thickness = float(t)
        sd.plate_material = material

        # Define end-plate welds
        # Using parameter value equ=True forces welds to be designed equivalent strength with beam
        define_welds(sd, sd.beam, equ=equ)

        # Upper overhang
        if upper_overhang > 0.0:
            sd.upper_overhang = max(float(upper_overhang), float(ceil(sd.weld_f*sqrt(2.0))))
        else:
            sd.upper_overhang = float(ceil(sd.weld_f*sqrt(2.0) + sd.plate_thickness))

        # Lower overhang
        if lower_overhang > 0.0:
            sd.lower_overhang = max(float(lower_overhang), float(ceil(sd.weld_f*sqrt(2.0))))
        else:
            sd.lower_overhang = float(ceil(sd.weld_f*sqrt(2.0) + sd.plate_thickness))

        # Plate height
        if plate_height > 0.0:
            sd.plate_height = plate_height
            sd.lower_overhang = max(sd.plate_height - sd.upper_overhang - sd.beam.h, float(lower_overhang),
                                    float(ceil(sd.weld_f*sqrt(2.0))))

        sd.plate_height = sd.upper_overhang + sd.beam.h + sd.lower_overhang


    def add_row(self, beam_id, z, e=0, d=0, bolt_material = 8.8, d_washer=0.0, t_washer=0.0,
                thread_in_shear_plane=True, shear_row=False):

        sd = self.side[beam_id]

        if 0.0 < z < sd.plate_height:
            temp_var1 = []
            # Checking duplicate rows (according to z coord of row), original gets overwritten!
            for i in sd.row:
                temp_var1.append(i.z)
            temp_var2 = [i for i, x in enumerate(temp_var1) if x == z]
            if temp_var2 != []:
                print("Row with same z-coordinate already defined! Original row overwritten!")
                for i in temp_var2:
                    del sd.row[i]
                    sd.row.append(Row(self, sd, z, e, d, bolt_material, d_washer=d_washer, t_washer=t_washer,
                                                  thread_in_shear_plane=thread_in_shear_plane, shear_row=shear_row))
            else:
                sd.row.append(Row(self, sd, z, e, d, bolt_material, d_washer=d_washer, t_washer=t_washer,
                                  thread_in_shear_plane=thread_in_shear_plane, shear_row=shear_row))

            # Sorting bolt rows from top to bottom according to z-value
            sd.row = sorted(sd.row, key=lambda rows: rows.z)

            # Adding id number to each row
            for i in range(0,len(sd.row)):
                sd.row[i].id = i

    def add_row_groups(self, beam_id):

        sd = self.side[beam_id]

        # Define row locations
        self.row_location(beam_id)

        if len(sd.row) < 1:
            if "No rows defined for side " + str(beam_id) not in self.errors:
                self.errors.append("No rows defined for side " + str(beam_id))
            return
        elif len(sd.row) == 1:
            if "No possible row groups, too few bolt rows for side " + str(beam_id) not in self.warnings:
                self.warnings.append("No possible row groups, too few bolt rows for side " + str(beam_id))
            return

        # Create row groups
        sd.groups = groups(sd, len(sd.row))

        # Adding id number to each row group
        for i in range(0, len(sd.groups)):
            sd.row_group.append(RowGroup(sd, i))

        # Allocate effective lengths for rows due to row group effect
        for row in sd.row:
            row.location_fg = ['']*len(sd.groups)
            row.location_pg = ['']*len(sd.groups)
            row.l_eff_cpf = [0.0]*len(sd.groups)
            row.l_eff_ncf = [0.0]*len(sd.groups)
            row.l_eff_cpp = [0.0]*len(sd.groups)
            row.l_eff_ncp = [0.0]*len(sd.groups)

    # Determining row locations
    def row_location(self, beam_id):

        sd = self.side[beam_id]

        # Defining row locations on column flange and end-plate side
        for i in range(0, len(sd.row)):
            # Checking if rows above are shear rows
            if sd.row[i].id == 0 and sd.row[i].shear_row:
                self.warnings.append("First row cannot be shear row! First row changed to tension row!")
                sd.row[i].shear_row = False
            elif sd.row[i].id > 0 and sd.row[i - 1].shear_row:
                sd.row[i].shear_row = True
        for i in range(0, len(sd.row)):
            # ----------------Row location in end-plate--------------------------- #
            z1 = sd.upper_overhang
            h_p = sd.plate_height
            h_b = sd.beam.h

            if 0.0 < sd.row[i].z < z1:          # Rows on upper overhang
                sd.row[i].location_p = "Bolt-row outside tension flange of beam"
                sd.row[i].e_x = sd.row[i].z
                sd.row[i].m_x = z1 - sd.row[i].e_x - 0.8*sqrt(2.0)*sd.weld_f
                if sd.row[i].id != 0 and sd.row[i - 1].location_p == "Bolt-row outside tension flange of beam":
                    self.errors.append('More than one bolt-row outside tension flange of beam!')
            elif z1 < sd.row[i].z < z1 + h_b:         # Rows between flanges
                if i == 0:
                    sd.row[i].location_p = "First bolt-row below tension flange of beam"
                    sd.row[i].m_2 = sd.row[i].z - z1 - sd.beam.t_f - 0.8*sqrt(2.0)*sd.weld_f
                    if sd.row[i].shear_row:
                        sd.row[i].location_p = "Shear row"
                elif i >= 1 and sd.row[i - 1].location_p == "Bolt-row outside tension flange of beam":
                    sd.row[i].location_p = "First bolt-row below tension flange of beam"
                    sd.row[i].m_2 = sd.row[i].z - z1 - sd.beam.t_f - 0.8*sqrt(2.0)*sd.weld_f
                    if sd.row[i].shear_row:
                        sd.row[i].location_p = "Shear row"
                elif sd.row[i].id == sd.row[-1].id:
                    sd.row[i].location_p = "Other end bolt-row"
                    if sd.row[i].shear_row:
                        sd.row[i].location_p = "Shear row"
                        if z1 < sd.row[i - 1].z < z1 + h_b and \
                           sd.row[i - 1].location_p != "First bolt-row below tension flange of beam" \
                           and not sd.row[i - 1].shear_row:
                            sd.row[i - 1].location_p = "Other end bolt-row"
                elif sd.row[i].shear_row:
                    sd.row[i].location_p = "Shear row"
                    if z1 < sd.row[i - 1].z < z1 + h_b and \
                       sd.row[i - 1].location_p != "First bolt-row below tension flange of beam" and \
                       not sd.row[i - 1].shear_row:
                        sd.row[i - 1].location_p = "Other end bolt-row"
                else:
                    sd.row[i].location_p = "Other inner bolt-row"
            elif z1 + h_b < sd.row[i].z < h_p:        # Bolt rows on lower overhang
                sd.row[i].location_p = "Bolt row in compression"
                sd.row[i].shear_row = True
                if sd.row[i].id == 0: pass
                elif z1 < sd.row[i - 1].z < z1 + h_b and \
                   sd.row[i - 1].location_p != "First bolt-row below tension flange of beam":
                    sd.row[i - 1].location_p = "Other end bolt-row"
                    if sd.row[i - 1].shear_row:
                        sd.row[i - 1].location_p = "Shear row"
            else:                   # Bolt rows outside end-plate
                sd.row[i].location_p = "Bolt row outside end-plate"
                sd.row[i].shear_row = True
                if sd.row[i].id == 0: pass
                elif z1 < sd.row[i - 1].z < z1 + h_b and sd.row[i - 1].location_p not in \
                        ["Bolt-row outside tension flange of beam",
                         "First bolt-row below tension flange of beam",
                         "Shear row"]:
                    sd.row[i - 1].location_p = "Other end bolt-row"
                    if sd.row[i - 1].shear_row:
                        sd.row[i - 1].location_p = "Shear row"

            # ----------------Row location in column flange--------------------------- #
            if sd.stiffener == "unstiffened":
                # Row in end-plate tension side
                if 0.0 <= sd.row[i].z <= z1 + h_b and (sd.row[i].id == 0 or sd.row[i].id == sd.row[-1].id):
                    if not sd.row[i].shear_row:
                        sd.row[i].location_f = "End bolt-row"
                elif 0.0 <= sd.row[i].z <= z1 + h_b and 0 < sd.row[i].id < sd.row[-1].id:
                    if not sd.row[i].shear_row:
                        sd.row[i].location_f = "Inner bolt-row"
                        if sd.row[i - 1].shear_row:
                            sd.row[i].location_f = "End bolt-row"
                # Row in end-plate compression side
                if z1 + h_b < sd.row[i].z <= h_p:
                    sd.row[i].location_f = "Bolt row in compression"
                    sd.row[i].shear_row = True
                    if sd.row[i].id == 0: pass
                    elif 0.0 <= sd.row[i-1].z <= z1 + h_b and not sd.row[i - 1].shear_row:
                        sd.row[i - 1].location_f = "End bolt-row"
                elif not 0.0 < sd.row[i].z < h_p:
                    sd.row[i].location_f = "Bolt row outside end-plate"
                    sd.row[i].shear_row = True
                    if sd.row[i].id == 0: pass
                    elif 0.0 <= sd.row[i-1].z <= z1 + h_b and not sd.row[i - 1].shear_row:
                        sd.row[i - 1].location_f = "End bolt-row"
                # Checking shear rows
                if sd.row[i].shear_row and sd.row[i].location_f not in ["Bolt row in compression",
                                                                        "Bolt row outside end-plate"]:
                    sd.row[i].location_f = "Shear row"
                    sd.row[i].shear_row = True
                    if sd.row[i].id == 0: pass
                    elif not sd.row[i - 1].shear_row and 0.0 <= sd.row[i - 1].z <= z1 + h_b:
                        sd.row[i - 1].location_f = "End bolt-row"
            else:
                print(self.NOT_IMPLEMENTED + "Column web stiffeners")

    # Determining effective lengths of rows and row groups
    def effective_lengths(self, beam_id):

        sd = self.side[beam_id]

        # Effective lengths of rows
        for i in range(0, len(sd.row)):
            T_stub(sd.row[i], component="flange", rg="row", stiffener=sd.stiffener, b_p=sd.plate_width)
            T_stub(sd.row[i], component="plate", rg="row", stiffener=sd.stiffener, b_p=sd.plate_width)

        x1 = "First bolt-row below tension flange of beam"
        x2 = "Other inner bolt-row"
        x3 = "Other end bolt-row"

        x4 = "Bolt-row outside tension flange of beam"
        x5 = "Bolt row outside end-plate"
        x6 = "Shear row"
        x7 = "Bolt row in compression"

        # Effective lengths of row groups
        for j in range(0, len(sd.groups)):
            sd.row_group[j].p = []

            # Center distances of rows in row group
            for i in range(0, len(sd.row_group[j].row)):
                if i == 0:
                    sd.row_group[j].p.append(sd.row_group[j].row[1].z - sd.row_group[j].row[0].z)
                elif i == len(sd.row_group[j].row) - 1:
                    sd.row_group[j].p.append(sd.row_group[j].row[-1].z - sd.row_group[j].row[-2].z)
                else:
                    sd.row_group[j].p.append((sd.row_group[j].row[i].z - sd.row_group[j].row[i - 1].z) / 2.0 +
                                             (sd.row_group[j].row[i + 1].z - sd.row_group[j].row[i].z) / 2.0)

            # Effective lengths of row groups on flange side
            l_eff_ncf_sum = 0.0
            l_eff_cpf_sum = 0.0
            temp_var_name = []
            for i in range(0, len(sd.row_group[j].row)):
                if i == 0:
                    sd.row_group[j].row[i].location_fg[sd.row_group[j].id] = "End bolt-row"
                elif i == len(sd.row_group[j].row) - 1:
                    sd.row_group[j].row[i].location_fg[sd.row_group[j].id] = "End bolt-row"
                else:
                    sd.row_group[j].row[i].location_fg[sd.row_group[j].id] = "Inner bolt-row"

                T_stub(sd.row_group[j].row[i], component="flange", rg="row group", stiffener=sd.stiffener,
                       b_p=sd.plate_width, p=sd.row_group[j].p[i], group_id=sd.row_group[j].id)
                l_eff_ncf_sum += sd.row_group[j].row[i].l_eff_ncf[sd.row_group[j].id]
                l_eff_cpf_sum += sd.row_group[j].row[i].l_eff_cpf[sd.row_group[j].id]
                temp_var_name.append(sd.row_group[j].row[i].location_f)

            # Removing redundant row groups from flange side
            if (x5 in temp_var_name) or (x6 in temp_var_name) or (x7 in temp_var_name):
                sd.row_group[j].l_eff_1f = 0.0
                sd.row_group[j].l_eff_2f = 0.0
                for row in sd.row_group[j].row:
                    row.l_eff_cpf[sd.row_group[j].id] = 0.0
                    row.l_eff_ncf[sd.row_group[j].id] = 0.0
            else:
                sd.row_group[j].l_eff_1f = min(l_eff_ncf_sum, l_eff_cpf_sum)
                sd.row_group[j].l_eff_2f = l_eff_ncf_sum

            # Effective lengths of row groups on plate side
            l_eff_ncp_sum = 0.0
            l_eff_cpp_sum = 0.0
            temp_var_name = []
            for i in range(0, len(sd.row_group[j].row)):
                if i == 0:
                    sd.row_group[j].row[i].location_pg[sd.row_group[j].id] = "Other end bolt-row"
                    if sd.row_group[j].row[i].location_p == "First bolt-row below tension flange of beam":
                        sd.row_group[j].row[i].location_pg[sd.row_group[j].id] = "First bolt-row below tension flange of beam"
                elif i == len(sd.row_group[j].row) - 1:
                    sd.row_group[j].row[i].location_pg[sd.row_group[j].id] = "Other end bolt-row"
                    if sd.row_group[j].row[i].location_p == "First bolt-row below tension flange of beam":
                        sd.row_group[j].row[i].location_pg[sd.row_group[j].id] = "First bolt-row below tension flange of beam"
                else:
                    sd.row_group[j].row[i].location_pg[sd.row_group[j].id] = "Other inner bolt-row"
                    if sd.row_group[j].row[i].location_p == "First bolt-row below tension flange of beam":
                        sd.row_group[j].row[i].location_pg[sd.row_group[j].id] = "First bolt-row below tension flange of beam"

                T_stub(sd.row_group[j].row[i], component="plate", rg="row group", stiffener=sd.stiffener,
                       b_p=sd.plate_width, p=sd.row_group[j].p[i], group_id=sd.row_group[j].id)
                l_eff_ncp_sum += sd.row_group[j].row[i].l_eff_ncp[sd.row_group[j].id]
                l_eff_cpp_sum += sd.row_group[j].row[i].l_eff_cpp[sd.row_group[j].id]
                temp_var_name.append(sd.row_group[j].row[i].location_p)

            # Removing redundant row groups from plate side
            if (x4 in temp_var_name) or (x5 in temp_var_name) or (x6 in temp_var_name) or (x7 in temp_var_name):
                sd.row_group[j].l_eff_1p = 0.0
                sd.row_group[j].l_eff_2p = 0.0
                for row in sd.row_group[j].row:
                    row.l_eff_cpp[sd.row_group[j].id] = 0.0
                    row.l_eff_ncp[sd.row_group[j].id] = 0.0
            else:
                sd.row_group[j].l_eff_1p = min(l_eff_ncp_sum, l_eff_cpp_sum)
                sd.row_group[j].l_eff_2p = l_eff_ncp_sum

    def calculate_capacity(self, beam_id):
        # Defining connection configuration
        if len(self.side[0].beam.size) != 0 and len(self.side[1].beam.size) != 0:
            self.joint_config = "Double-sided"
            if len(self.column.size) == 0:
                self.joint_config = "Beam splice"
        else:
            self.joint_config = "Single-sided"

        sd = self.side[beam_id]

        # Defining transformation parameter beta
        if sd.beam.My == 1.0e-6 and sd.beam.My == 1.0e-6:
            if "No bending loads in stronger axis defined for side " + \
                    str(beam_id) + "! Used value beta = 1.0!" not in self.warnings:
                self.warnings.append("No bending loads in stronger axis defined for side " + str(beam_id) +
                                     "! Used value beta = 1.0!")
            sd.beta = 1.0
        else:
            if beam_id == 0:
                sd.beta = min(abs(1.0 - self.side[1].beam.My/self.side[0].beam.My), 2.0 - 1.0e-12) + 1.0e-12
            elif beam_id == 1:
                sd.beta = min(abs(1.0 - self.side[0].beam.My/self.side[1].beam.My), 2.0 - 1.0e-12) + 1.0e-12

        # Capacity of compression side
        cweb_shear(self, self.column)                                               # Column web in shear
        cweb_trv_comp(self, self.column, sd, sigma_com_Ed=self.column.sigma_com_Ed) # Column web transverse compression
        fc_comp(self, sd)                                                           # Beam web and flange in compression
        sd.F_comp_Rd = min(self.V_wp_Rd/sd.beta, sd.F_c_wc_Rd, sd.F_c_fb_Rd)        # Total capacity

        # Defining effective lengths of rows
        self.effective_lengths(beam_id)

        # Distance of compression center from top of the end-plate
        sd.compression_center = sd.upper_overhang + sd.beam.h - 0.5*sd.beam.t_f

        # Capacities for each bolt row
        sum_F_t = 0.0
        for i in range(0, len(sd.row)):
            # Row distance from compression center
            if sd.row[i].z > sd.compression_center:
                sd.row[i].h_r = 0.0
            else:
                sd.row[i].h_r = sd.compression_center - sd.row[i].z

            if sd.row[i].location_p == "Bolt row outside end-plate":
                self.warnings.append("Row {0} on side {1} outside endplate!".format(i, beam_id))
                sd.row[i].mode = "None"
            elif sd.row[i].location_p == "Bolt row in compression":
                sd.row[i].mode = "None"
            else:
                cweb_trv_tens(self.column, sd, sd.row[i])
                bweb_tens(sd, sd.row[i])
                cflange_bending(self.column, sd, sd.row[i])
                end_plate_bending(self.column, sd, sd.row[i])

                # Defining total capacity of row from basic components
                sd.row[i].F_t_Rd = 0.0
                sd.row[i].F_t_Rd = min(i for i in [sd.row[i].F_t_f_Rd, sd.row[i].F_t_p_Rd,
                                                   sd.row[i].F_t_wc_Rd, sd.row[i].F_t_wb_Rd] if i > 0.0)

                # Defining row break mode
                if sd.row[i].F_t_Rd == sd.row[i].F_t_f_Rd:
                    sd.row[i].F_t_Rd = sd.row[i].F_t_f_Rd
                    sd.row[i].mode = "Flange: " + sd.row[i].mode_f + ", " + \
                                     sd.row[i].yield_mode_f + " yield lines"
                elif sd.row[i].F_t_Rd == sd.row[i].F_t_p_Rd:
                    sd.row[i].F_t_Rd = sd.row[i].F_t_p_Rd
                    sd.row[i].mode = "Plate: " + sd.row[i].mode_p + ", " + \
                                     sd.row[i].yield_mode_p + " yield lines"
                elif sd.row[i].F_t_Rd == sd.row[i].F_t_wc_Rd:
                    sd.row[i].F_t_Rd = sd.row[i].F_t_wc_Rd
                    sd.row[i].mode = "Column web in transverse tension"
                elif sd.row[i].F_t_Rd == sd.row[i].F_t_wb_Rd:
                    sd.row[i].F_t_Rd = sd.row[i].F_t_wb_Rd
                    sd.row[i].mode = "Beam web in tension"

                # Row capacity limited by capacity of compression side, 6.2.7.2 (7)
                if sd.row[i].F_t_Rd > sd.F_comp_Rd - sum_F_t:
                    sd.row[i].F_t_Rd = sd.F_comp_Rd - sum_F_t
                    sd.row[i].mode = "capacity of compression side"
                    if sd.row[i].F_t_Rd < 0.0:
                        sd.row[i].F_t_Rd = 0.0
                sum_F_t = sum_F_t + sd.row[i].F_t_Rd

        # Capacities for each bolt row group
        for i in range(0, len(sd.row_group)):
            if sd.row_group[i].l_eff_1f != 0.0 or sd.row_group[i].l_eff_2f != 0.0:
                cweb_trv_tens(self.column, sd, sd.row_group[i])
                cflange_bending(self.column, sd, sd.row_group[i])
            if sd.row_group[i].l_eff_1p != 0.0 or sd.row_group[i].l_eff_2p != 0.0:
                bweb_tens(sd, sd.row_group[i])
                end_plate_bending(self.column, sd, sd.row_group[i])

            # Defining total capacity of row group from basic components
            sd.row_group[i].F_t_Rd = 0.0
            sd.row_group[i].F_t_Rd = min(i for i in [sd.row_group[i].F_t_f_Rd, sd.row_group[i].F_t_p_Rd,
                                                     sd.row_group[i].F_t_wc_Rd, sd.row_group[i].F_t_wb_Rd] if i > 0.0)

            # Defining row group break mode
            if sd.row_group[i].F_t_Rd == sd.row_group[i].F_t_f_Rd:
                sd.row_group[i].F_t_Rd = sd.row_group[i].F_t_f_Rd
                sd.row_group[i].mode = "Flange: " + sd.row_group[i].mode_f
            elif sd.row_group[i].F_t_Rd == sd.row_group[i].F_t_p_Rd:
                sd.row_group[i].F_t_Rd = sd.row_group[i].F_t_p_Rd
                sd.row_group[i].mode = "Plate: " + sd.row_group[i].mode_p
            elif sd.row_group[i].F_t_Rd == sd.row_group[i].F_t_wc_Rd:
                sd.row_group[i].F_t_Rd = sd.row_group[i].F_t_wc_Rd
                sd.row_group[i].mode = "Column web in transverse tension"
            elif sd.row_group[i].F_t_Rd == sd.row_group[i].F_t_wb_Rd:
                sd.row_group[i].F_t_Rd = sd.row_group[i].F_t_wb_Rd
                sd.row_group[i].mode = "Beam web in tension"

            # Row capacity reduced by row group capacity
            sum_group_F_t = 0.0
            for row in sd.row_group[i].row:
                if row.F_t_Rd > sd.row_group[i].F_t_Rd - sum_group_F_t:
                    row.F_t_Rd = sd.row_group[i].F_t_Rd - sum_group_F_t
                    row.mode = str("capacity of row group " + str(sd.row_group[i].id + 1) +
                                   ": " + sd.row_group[i].mode)
                    if row.F_t_Rd < 0.0:
                        row.F_t_Rd = 0.0
                sum_group_F_t += row.F_t_Rd

        # Row capacity limited by effective design tension resistance of previous bolt rows
        for i in range(0, len(sd.row)):
            if sd.row[i].F_t_Rd > 1.9*(sd.row[i].F_T_Rd/sd.row[i].n_row):
                for j in range(i, len(sd.row)):
                    if sd.row[i].F_t_Rd*(sd.row[j].h_r/sd.row[i].h_r) < sd.row[j].F_t_Rd:
                        sd.row[j].F_t_Rd = sd.row[i].F_t_Rd*(sd.row[j].h_r/sd.row[i].h_r)
                        sd.row[j].mode = str("effective tension resistance of previous bolt-row "
                                             + str(sd.row[i].id + 1))
                break

        # Total bending capacity of connection in along stronger axis
        sd.M_j_Rd = 0.0
        for row in sd.row:
            if row.location_f != 'Shear row' and row.location_p != 'Shear row':
                sd.M_j_Rd += row.h_r*row.F_t_Rd         # [Nmm]

        # Connection strength limits
        if self.column_config == "column end":
            sd.M_rigid = max(sd.beam.W_pl_y*sd.beam.f_y, self.column.W_pl_y*self.column.f_y)
        elif self.column_config == "continuous column":
            sd.M_rigid = max(sd.beam.W_pl_y*sd.beam.f_y, 2.0*self.column.W_pl_y*self.column.f_y)
        sd.M_pinned = 0.25*sd.M_rigid

        # Connection strength classification
        if sd.M_j_Rd >= sd.M_rigid:
            sd.strength_class = "full-strength"
        elif sd.M_j_Rd <= sd.M_pinned:
            sd.strength_class = "nominally pinned"
        else:
            sd.strength_class = "partial-strength"

    def calculate_stiffness(self, beam_id):

        sd = self.side[beam_id]

        # Values
        A_vc = self.column.A_v                              # Shear area of column
        beta = sd.beta                                      # Beta factor of connection on beam_id side
        t_fb = sd.beam.t_f                                  # Flange thickness of beam
        t_fc = self.column.t_f                              # Flange thickness column
        s = self.column.r                                   # Rounding on column
        t_p = sd.plate_thickness                            # End-plate thickness
        lo_p = sd.lower_overhang                            # Length of lower overhang
        a_p = sd.weld_f                                     # Weld size of plate on flanges
        s_p = min(2.0*t_p, t_p + lo_p - sqrt(2.0)*a_p)      # Width of compressed zone under compressed flange
                                                            # using 45 deg stress distribution
        t_wc = self.column.t_w                              # Web thickness of column
        d_c = self.column.d_w                               # Height of straight portion of web in column
        E_b = sd.beam.E                                     # Youngs modulus of beam
        I_b = sd.beam.I_y                                   # second moment of inertia of beam
        L_b = sd.beam.L                                     # Length of beam
        I_c = self.column.I_y                               # second moment of inertia of column
        L_c = self.column.L                                 # Length of column

        # Calculating total spring stiffness of tension side
        sum1 = 0.0
        sum2 = 0.0
        for i in range(0, len(sd.row)):
            if sd.row[i].h_r > 0.0 and not sd.row[i].shear_row:
                stiffness_coeff(self, sd, sd.row[i])      # Defining spring stiffness of components for the row
                if self.joint_config == "Beam splice":
                    sd.row[i].k_eff = 1.0/(1.0/sd.row[i].k_5 + 1.0/sd.row[i].k_10)
                else:
                    sd.row[i].k_eff = 1.0/(1.0/sd.row[i].k_3 + 1.0/sd.row[i].k_4 +
                                           1.0/sd.row[i].k_5 + 1.0/sd.row[i].k_10)
                sum1 = sum1 + sd.row[i].k_eff*sd.row[i].h_r
                sum2 = sum2 + sd.row[i].k_eff*sd.row[i].h_r**2.0
        sd.z_eq = sum2/sum1
        sd.k_eq = sum1/sd.z_eq

        # Column web panel in shear
        if len(sd.row) > 1:
            z = sd.z_eq                               # 'Exact value'
            # z = (sd.row[0].h_r + sd.row[1].h_r)/2.0     # Approximation according to EN 1993-1-8 figure 6.15
            sd.k_1 = (0.38 * A_vc)/(beta * z)
        elif len(sd.row) == 1:
            sd.k_1 = (0.38*A_vc)/(beta*sd.row[0].h_r)
        else:
            sd.k_1 = 1.0e-12

        # Column web in compression
        b_eff_c_wc = t_fb + 2.0*sqrt(2.0)*a_p + 5.0*(t_fc + s) + s_p
        sd.k_2 = 0.7*b_eff_c_wc*t_wc/d_c

        # Calculating total spring stiffness of compression side
        if self.joint_config == "Double-sided":
            # Calculating initial rotational stiffness of connection
            sd.S_j_ini = (E_b*sd.z_eq**2.0)/(1.0/sd.k_2 + 1.0/sd.k_eq)  # [Nmm/rad]
        else:
            # Calculating initial rotational stiffness of connection
            sd.S_j_ini = (E_b*sd.z_eq**2.0)/(1.0/sd.k_1 + 1.0/sd.k_2 + 1.0/sd.k_eq)  # [Nmm/rad]

        # Limits for connection sfiffness
        sd.S_j_rigid = self.k_b*E_b*I_b/L_b              # [Nmm/rad]
        sd.S_j_pinned = 0.5*E_b*I_b/L_b                  # [Nmm/rad]

        # Rotational stiffness classification of connection
        if sd.S_j_ini >= sd.S_j_rigid:
            sd.stiffness_class = "rigid"
        elif sd.S_j_ini <= sd.S_j_pinned:
            sd.stiffness_class = "pinned"
        else:
            sd.stiffness_class = "semi-rigid"

        if (I_b/L_b)/(I_c/L_c) < 0.1:
            self.warnings.append("(I_b/L_b)/(I_c/L_c) < 0.1 => Semi-rigid joint!")
            sd.stiffness_class = "semi-rigid"

        # Stiffness ratio
        psi = 2.7  # Coefficient obtained from table 6.8, bolted end-plate
        eta = 2.0  # Stiffness modification coefficient for beam-to-column joints
        if abs(sd.beam.My) <= 2.0/3.0*sd.M_j_Rd:
            if self.analysis_method == 'plastic':
                sd.my = 1.0e-6
            else:
                sd.my = 1.0
        elif 2.0/3.0*sd.M_j_Rd < abs(sd.beam.My) <= sd.M_j_Rd:
            if self.analysis_method == 'elastic':
                sd.my = eta
            elif self.analysis_method == 'elastic-plastic':
                sd.my = (1.5*sd.beam.My/sd.M_j_Rd)**psi
            elif self.analysis_method == 'plastic':
                sd.my = 1.0e-6

        # Stiffness of connection, takes loading and global analysis method into account
        sd.S_j = sd.S_j_ini/sd.my          # [Nmm/rad]

    def calculate_utilization_rate(self):
        # Resulting shear force on column web [N]
        self.V_wp_Ed = abs(self.side[0].beam.My*1.0e6/self.side[0].z_eq -
                           self.side[1].beam.My*1.0e6/self.side[1].z_eq) + \
                       abs(self.column_loads_under[2] - self.column_loads_over[2])

        # Utilization rate of column web in shear
        self.n_V_wp = self.V_wp_Ed/self.V_wp_Rd

        for beam_id in range(0, 2):
            sd = self.side[beam_id]

            if sd.beam.size == "":
                continue

            sd.V_Rd = 1.0e-6
            for row in sd.row:
                # Tension load for row due to moment
                # NOT EXACT: GIVES ESTIMATION TO BE USED WHEN DEFINING SHEAR CAPACITY OF ROW!
                if row.location_f != 'Shear row' and row.location_p != 'Shear row':
                    row.F_t_Ed = (abs(sd.beam.My*1.0e6)/sd.M_j_Rd)*row.F_t_Rd        # [N]

                # NORMAL FORCE OF BEAM NOT TAKEN INTO ACCOUNT:
                # ACCORDING TO SFS EN 1993-1-8 NORMAL FORCE DOES NOT NEED TO BE TAKEN INTO ACCOUNT IF N < N_PL_BEAM
                # Tension load for row due to axial force, positive axial force is tension
                # if sd.beam.N > 1.0e-6:
                #     row.F_t_Ed += sd.beam.N/(len(sd.row)*row.n_row)     # [N]

                # Bearing resistance of row, SFS EN 1993-1-8, table 3.4
                # ONLY SHEAR IN STRONGER DIRECTION OF BEAM TAKEN INTO ACCOINT!!!
                # ALL BOLTS TAKEN TO BE EDGE BOLTS!!!
                f_ub = row.f_ub                         # Ultimate strength of bolt material
                f_u = mat[sd.plate_material]["f_u"]     # Ultimate strength of plate material

                # In shear load direction
                if row.location_p in ["Bolt-row outside tension flange of beam",
                                      "First bolt-row below tension flange of beam"]:       # End bolt-row
                    alfa_df = row.e_1/(3.0*row.d_0)
                    alfa_dp = row.z/(3.0*row.d_0)
                elif row.location_p not in ["Bolt row outside end-plate"]:                  # Inner bolt-row
                    alfa_df = (row.z - sd.row[row.id-1].z)/(3.0*row.d_0) - 1.0/4.0
                    alfa_dp = alfa_df
                else:                                                                       # Bolt-row outside plate
                    alfa_df = 0.0
                    alfa_dp = 0.0
                # Perpendicular to shear load direction
                k_1f = min(2.8*row.ec/row.d_0 - 1.7, 2.5)
                k_1p = min(2.8*row.e/row.d_0 - 1.7, 2.5)

                alfa_bf = min(alfa_df, f_ub/f_u, 1.0)
                alfa_bp = min(alfa_dp, f_ub/f_u, 1.0)

                row.F_b_Rdf = k_1f*alfa_bf*f_u*row.d*sd.plate_thickness/gamma_M[2]      # Bearing on flange side
                row.F_b_Rdp = k_1p*alfa_bp*f_u*row.d*sd.plate_thickness/gamma_M[2]      # Bearing on plate side

                # Total bearing resistance of row
                row.F_b_Rd = row.n_row*min(row.F_b_Rdf, row.F_b_Rdp)

                print('row.F_b_Rd = ' + str(row.F_b_Rd))
                print('row.id = ' + str(row.id))

                # Total shear resistance of row, taking into account bolt shear capacity and bolt bearing
                row.V_Rd = min(row.V_Rd, row.F_b_Rd)

                # Shear capacity of rows, taking into account tension on row
                # According to SFS EN 1993-1-8 chapter 6.2.2 (2)
                if row.F_t_Ed > 0.0:
                    row.V_Rd = (0.4/1.4)*row.V_Rd       # [N]
                if row.location_f == "Bolt row outside end-plate" or row.location_p == "Bolt row outside end-plate":
                    row.V_Rd = 0.0
                sd.V_Rd += row.V_Rd        # [N]
                if self.shear_block is not None:
                    print(self.NOT_IMPLEMENTED + ' Shear block!')

            # Utilization rate of connection under moment
            sd.n_M = abs(sd.beam.My*1.0e6)/sd.M_j_Rd

            # Utilization rate of side under shear
            sd.n_V = sd.beam.Vz*1.0e3/sd.V_Rd

    def cost_calculation(self):
        self.total_cost = 0.0
        # Cost calculation for each side
        for beam_id in range(0,2):

            sd = self.side[beam_id]

            if sd.beam.size != "":
                # Cost of plate material
                sd.material_cost.cost_of_plate(self, sd)

                # Calculating cost of bolts, nuts and washers
                temp_cost = 0.0
                L_holes = 0.0
                n_holes = 0
                for j in range(0, len(sd.row)):
                    # Disregarding rows outside end-plate
                    if sd.row[j].location_f == "Bolt row outside end-plate" or \
                                    sd.row[j].location_p == "Bolt row outside end-plate":
                        continue
                    # Adding bolt material cost for row
                    if sd.row[j].bolt_material != 8.8:
                        self.warnings.append("Price values defined only for bolt material 8.8! Given bolt material " +
                                             str(sd.row[j].bolt_material))
                    bolt_cost(sd.row[j], cost_function=self.cost_function)
                    temp_cost = temp_cost + sd.row[j].material_cost
                    L_holes = L_holes + pi*sd.row[j].d_0*sd.row[j].n_row
                    n_holes = n_holes + sd.row[j].n_row

                # Total cost of bolts, nuts and washers
                sd.material_cost.bolts = temp_cost

                # Total cost of material for side
                sd.material_cost.total = sd.material_cost.plate + sd.material_cost.bolts

                # Cost of cost centers
                sd.blasting = BlastingCost(self, sd)
                sd.cutting = CuttingCost(self, sd, n_holes=n_holes, L_holes=L_holes)
                sd.drilling = DrillingCost(self, sd)
                sd.part_assembly = PartAssemblyCost(self, sd)
                sd.post_treatment = PostTreatmentCost(self, sd)
                sd.painting = PaintingCost(self, sd)

                # Total cost of side [e]
                sd.cost = sd.material_cost.total + sd.blasting.cost + sd.cutting.cost + sd.drilling.cost + \
                          sd.part_assembly.cost + sd.post_treatment.cost + sd.painting.cost

                # Total cost of connection [e]
                self.total_cost = self.total_cost + sd.cost

    def update(self):
        for beam_id in range(0, 2):
            sd = self.side[beam_id]
            if sd.beam.size != '':
                # Update plate height
                sd.plate_height = sd.upper_overhang + sd.beam.h + sd.lower_overhang

                # Update row information
                for i in range(0, len(sd.row)):
                    sd.row[i].update(self, sd)

                # Sorting bolt rows from top to bottom according to z-value
                sd.row = sorted(sd.row, key=lambda rows: rows.z)

                # Adding id number to each row
                for i in range(0, len(sd.row)):
                    sd.row[i].id = i

                # Redefine row groups
                sd.groups = []
                sd.row_group = []
                self.add_row_groups(beam_id)

                if not (str("No rows defined for side " + str(beam_id)) in self.errors):
                    # Redefine row location and effective lengths
                    self.row_location(beam_id)
                    self.effective_lengths(beam_id)

                    # Calculate connection information again
                    self.calculate_capacity(beam_id)
                    self.calculate_stiffness(beam_id)
        self.calculate_utilization_rate()
        self.cost_calculation()

        self.evaluations += 1

    def print_geom(self, beam_id, side=True, front=True, savepath='', save_only=False):

        self.plot_nro = 1
        sd = self.side[beam_id]

        grey = '0.8'

        t_fc = self.column.t_f
        r_c = self.column.r
        h_c = self.column.h
        b_c = self.column.b
        h_b = sd.beam.h
        b_b = sd.beam.b
        t_fb = sd.beam.t_f
        t_wb = sd.beam.t_w
        d_wb = sd.beam.d_w
        a_f = sd.weld_f
        r_b = sd.beam.r
        h_p = sd.plate_height
        t_p = sd.plate_thickness
        b_p = sd.plate_width
        up_p = sd.upper_overhang

        if side:
            # Plotting side view
            plt.figure("side " + str(self.plot_nro))
            plt.clf()
            plt.xlabel("[mm]")
            plt.ylabel("[mm]")

            axis_lim = max(h_c + 500.0, h_b + 500.0, h_p + 200.0)

            # Plot column in figure
            column_y = [-100.0, axis_lim - 100.0]
            plt.plot([0.0, 0.0], column_y, 'k', zorder=2)
            plt.plot([-t_fc, -t_fc], column_y, 'k', zorder=2)
            plt.plot([-(t_fc + r_c), -(t_fc + r_c)], column_y, grey, zorder=1)
            plt.plot([-(h_c - (t_fc + r_c)), -(h_c - (t_fc + r_c))], column_y, grey, zorder=1)
            plt.plot([-(h_c - t_fc), -(h_c - t_fc)], column_y, 'k', zorder=2)
            plt.plot([-h_c, -h_c], column_y, 'k', zorder=2)

            # Plot plate in figure
            plate_x = [0.0, t_p, t_p, 0.0]
            plate_y = [0.0, 0.0, h_p, h_p]
            plt.plot(plate_x, plate_y, 'k', zorder=2)

            # Plot beam in figure
            beam_x = [t_p, axis_lim - (100.0 + h_c + t_p)]
            plt.plot(beam_x, [up_p, up_p], 'k', zorder=2)
            plt.plot(beam_x, [up_p + t_fb, up_p + t_fb], 'k', zorder=2)
            plt.plot(beam_x, [up_p + t_fb + r_b, up_p + t_fb + r_b], grey, zorder=1)
            plt.plot(beam_x, [up_p + h_b - (t_fb + r_b), up_p + h_b - (t_fb + r_b)], grey, zorder=1)
            plt.plot(beam_x, [up_p + h_b - t_fb, up_p + h_b - t_fb], 'k', zorder=2)
            plt.plot(beam_x, [up_p + h_b, up_p + h_b], 'k', zorder=2)

            # Plot welds
            plt.plot([t_p, t_p + sqrt(2.0)*a_f], [up_p - sqrt(2.0)*a_f, up_p], 'k', zorder=1)
            plt.plot([t_p, t_p + sqrt(2.0)*a_f], [up_p + t_fb + sqrt(2.0)*a_f, up_p + t_fb], 'k', zorder=1)
            plt.plot([t_p, t_p + sqrt(2.0)*a_f], [up_p + h_b - t_fb - sqrt(2.0)*a_f, up_p + h_b - t_fb], 'k', zorder=1)
            plt.plot([t_p, t_p + sqrt(2.0)*a_f], [up_p + h_b + sqrt(2.0)*a_f, up_p + h_b], 'k', zorder=1)

            # Plot rows on figure
            for i in range(0, len(sd.row)):
                z = sd.row[i].z
                d = sd.row[i].d
                h_nut = sd.row[i].h_nut
                h_bolt = sd.row[i].h_bolt
                d_m = sd.row[i].d_m
                L_b = sd.row[i].L_b
                t_washer = sd.row[i].t_washer
                d_washer = sd.row[i].d_washer

                if sd.row[i].location_p != "Bolt row outside end-plate":
                    # Plot bolt hole
                    plt.plot([-t_fc, t_p], [z + 0.5*d, z + 0.5*d], 'k--', linewidth=0.6)
                    plt.plot([-t_fc, t_p], [z - 0.5*d, z - 0.5*d], 'k--', linewidth=0.6)
                    # Plot washers
                    plt.gca().add_artist(plt.Rectangle((t_p, z - 0.5*d_washer), width=t_washer, height=d_washer,
                                                       edgecolor='k', facecolor="w", zorder=2))
                    plt.gca().add_artist(plt.Rectangle((-t_fc - t_washer, z - 0.5*d_washer), width=t_washer,
                                                       height=d_washer, edgecolor='k', facecolor="w", zorder=2))
                    # Plot nut and bolt head
                    plt.gca().add_artist(plt.Rectangle((t_p + t_washer, z - 0.5*d_m), width=h_bolt, height=d_m,
                                                       edgecolor='k', facecolor="w", zorder=2))
                    plt.gca().add_artist(plt.Rectangle((-t_fc - t_washer - h_nut, z - 0.5*d_m), width=h_nut, height=d_m,
                                                       edgecolor='k', facecolor="w", zorder=2))

                # Annotate row location on flange side
                plt.annotate(sd.row[i].location_f,
                             xy=(-25.0, z),
                             xycoords='data', xytext=(-160.0, 0.0),
                             textcoords='offset points', bbox=dict(boxstyle="square", fc="0.9"),
                             arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10,
                                             connectionstyle="angle,angleA=0,angleB=90,rad=10"))
                # Annotate row location on plate side
                plt.annotate(sd.row[i].location_p,
                             xy=(25.0, z),
                             xycoords='data', xytext=(25.0, 0.0),
                             textcoords='offset points', bbox=dict(boxstyle="square", fc="0.9"),
                             arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10,
                                             connectionstyle="angle,angleA=0,angleB=90,rad=10"))

                # Plot center line
                plt.plot([-0.5 * L_b, 0.5 * L_b], [z, z], 'b-.', linewidth=0.6, zorder=3)

            axes = plt.gca()
            axes.invert_yaxis()
            axes.set_aspect('equal', 'datalim')

            if savepath != '':
                plt.savefig(savepath + 'side.png', dpi=150)

        if front:
            # Plotting front view
            plt.figure("front " + str(self.plot_nro))
            plt.clf()
            plt.xlabel("[mm]")
            plt.ylabel("[mm]")

            axis_lim = max(b_c + 200.0, h_b + 200.0, h_p + 200.0)

            # Show column in figure
            plt.plot([-0.5*b_c, -0.5*b_c], [-100.0, axis_lim - 100.0], 'k', linewidth=1.0)
            plt.plot([0.5*b_c, 0.5*b_c], [-100.0, axis_lim - 100.0], 'k', linewidth=1.0)

            # Show plate in figure
            plt.gca().add_artist(plt.Rectangle((-0.5*b_p, 0.0), width=b_p, height=h_p, linewidth=1.0,
                                               edgecolor='k', facecolor="w", zorder=2))

            # Show beam in figure
            plt.plot([-0.5*b_b , -0.5*b_b, 0.5*b_b, 0.5*b_b], [up_p + t_fb, up_p, up_p, up_p + t_fb], 'k', zorder=3)
            plt.plot([-0.5*b_b, -0.5*t_wb - r_b], [up_p + t_fb, up_p + t_fb], 'k', zorder=3)
            plt.plot([0.5*b_b, 0.5*t_wb + r_b], [up_p + t_fb, up_p + t_fb], 'k', zorder=3)
            plt.plot([-0.5*b_b, -0.5*t_wb - r_b], [up_p + h_b - t_fb, up_p + h_b - t_fb], 'k', zorder=3)
            plt.plot([0.5*b_b, 0.5*t_wb + r_b], [up_p + h_b - t_fb, up_p + h_b - t_fb], 'k', zorder=3)
            plt.plot([-0.5*b_b , -0.5*b_b, 0.5*b_b, 0.5*b_b],
                     [up_p + h_b - t_fb, up_p + h_b, up_p + h_b, up_p + h_b - t_fb], 'k', zorder=3)
            plt.plot([-0.5*t_wb, -0.5*t_wb], [up_p + t_fb + r_b, up_p + t_fb + r_b + d_wb], 'k', zorder=3)
            plt.plot([0.5*t_wb, 0.5*t_wb], [up_p + t_fb + r_b, up_p + t_fb + r_b + d_wb], 'k', zorder=3)
            # Plot radius of root fillet
            plt.gca().add_patch(patches.Arc([-0.5 * t_wb - r_b, up_p + t_fb + r_b], 2.0 * r_b, 2.0 * r_b,
                                            theta1=270.0, theta2=360.0, color='k', linewidth=1.3, zorder=3))
            plt.gca().add_patch(patches.Arc([0.5 * t_wb + r_b, up_p + t_fb + r_b], 2.0 * r_b, 2.0 * r_b,
                                            theta1=180.0, theta2=270.0, color='k', linewidth=1.3, zorder=3))
            plt.gca().add_patch(patches.Arc([-0.5 * t_wb - r_b, up_p + h_b - t_fb - r_b], 2.0 * r_b, 2.0 * r_b,
                                            theta1=0.0, theta2=90.0, color='k', linewidth=1.3, zorder=3))
            plt.gca().add_patch(patches.Arc([0.5 * t_wb + r_b, up_p + h_b - t_fb - r_b], 2.0 * r_b, 2.0 * r_b,
                                            theta1=90.0, theta2=180.0, color='k', linewidth=1.3, zorder=3))

            # Show rows on figure
            for i in range(0, len(sd.row)):
                z = sd.row[i].z
                d = sd.row[i].d
                w = sd.row[i].w
                d_m = sd.row[i].d_m
                d_washer = sd.row[i].d_washer

                if sd.row[i].location_p != "Bolt row outside end-plate":
                    # Plot washers
                    plt.gca().add_artist(plt.Circle((-0.5*w, z), 0.5*d_washer, edgecolor='k', facecolor="w", zorder=4))
                    plt.gca().add_artist(plt.Circle((0.5*w, z), 0.5*d_washer, edgecolor='k', facecolor="w", zorder=4))
                    # Plot bolts
                    plt.gca().add_artist(plt.Circle((-0.5*w, z), 0.5*d_m, edgecolor='k', facecolor="w", zorder=4))
                    plt.gca().add_artist(plt.Circle((0.5*w, z), 0.5*d_m, edgecolor='k', facecolor="w", zorder=4))
                    plt.gca().add_artist(plt.Circle((-0.5*w, z), 0.5*d, linestyle='--', linewidth=0.6, edgecolor='k',
                                                    facecolor="w", zorder=4))
                    plt.gca().add_artist(plt.Circle((0.5*w, z), 0.5*d, linestyle='--', linewidth=0.6, edgecolor='k',
                                                    facecolor="w", zorder=4))

                    # Annotate row capacity
                    if sd.row[i].shear_row:
                        plt.annotate("Shear row, {0} = {1:6.2f} kN".format("$V_{Rd}$", sd.row[i].V_Rd * 1.0e-3),
                                     xy=(0.6 * b_p, z),
                                     xycoords='data', xytext=(25.0, 0.0),
                                     textcoords='offset points', bbox=dict(boxstyle="square", fc="0.9"),
                                     arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10,
                                                     connectionstyle="angle,angleA=0,angleB=90,rad=10"))
                    else:
                        plt.annotate("{0} = {1:6.2f} kN".format("$F_{t.Rd}$", sd.row[i].F_t_Rd * 1.0e-3),
                                     xy=(0.6 * b_p, z), xycoords='data', xytext=(25.0, 0.0),
                                     textcoords='offset points', bbox=dict(boxstyle="square", fc="0.9"),
                                     arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=10,
                                                     connectionstyle="angle,angleA=0,angleB=90,rad=10"))

                # Plot cross
                plt.plot([-0.5 * w, -0.5 * w], [z - 0.6 * d_m, z + 0.6 * d_m], 'b-.', linewidth=0.6, zorder=5)
                plt.plot([-0.5 * w - 0.6 * d_m, -0.5 * w + 0.6 * d_m], [z, z], 'b-.', linewidth=0.6, zorder=5)
                plt.plot([0.5 * w, 0.5 * w], [z - 0.6 * d_m, z + 0.6 * d_m], 'b-.', linewidth=0.6, zorder=5)
                plt.plot([0.5 * w - 0.6 * d_m, 0.5 * w + 0.6 * d_m], [z, z], 'b-.', linewidth=0.6, zorder=5)

            axes = plt.gca()
            axes.invert_yaxis()
            axes.set_aspect('equal', 'datalim')

            if savepath != '':
                plt.savefig(savepath + 'front.png', dpi=150)

        if save_only:
            plt.close('all')
        else:
            plt.draw()
            plt.pause(0.00000001)

    def info(self, values=True, cost_info=True, warnings_and_errors=True, figures=False, side=True, front=True):
        # Float formatting style
        fform = ".2f"

        if values:
            print(self.joint_config + ", " + self.column_config + ", " + self.sway)
            print('Global structure analysis method = ' + self.analysis_method)

            if self.column.size != "":
                print("\nColumn information:")
                self.column.info(fform)

            for beam_id in range(0, 2):

                sd = self.side[beam_id]

                if len(sd.beam.size) != 0:
                    print("\nBeam information:")
                    sd.beam.info(fform)

                    print("\nh_p = " + str(sd.plate_height) + "[mm]")
                    print(sd.welds)
                    print("a_f = {0} mm, a_w = {1} mm".format(sd.weld_f, sd.weld_w))
                    print("Compression center " + str(sd.compression_center) + " [mm] from top of the end-plate")

                    for j in range(0, len(sd.row)):
                        sd.row[j].info(fform, figures=figures)

                    print()
                    print("  V_wp_Rd = {0:{fm}} [kN] (Column web panel in shear)".format(self.V_wp_Rd*1.0e-3, fm='7.2f'))
                    print("F_c_wc_Rd = {0:{fm}} [kN] (Column web in transverse compression)"
                          .format(sd.F_c_wc_Rd*1.0e-3, fm='7.2f'))
                    print("F_c_fb_Rd = {0:{fm}} [kN] (Beam flange and web in compression)"
                          .format(sd.F_c_fb_Rd*1.0e-3, fm='7.2f'))
                    print('F_comp_Rd = {0:{fm}} [kN] (Design resistance of compression side, beta = {1:.2f})'
                          .format(sd.F_comp_Rd*1.0e-3, sd.beta, fm='7.2f'))

                    print("\nPossible row groups: " + str(sd.groups))

                    for j in range(0, len(sd.row_group)):
                        sd.row_group[j].info(fform)

                    print("\nEffective spring stiffness of components on compression side:")
                    print("          k_1 = {0:{fm}} [mm] (Column web panel in shear)\n"
                          "          k_2 = {1:{fm}} [mm] (Column web in compression)\n"
                          "   k_eff_comp = {2:{fm}} [mm] (Total effective spring stiffness of compression side)"
                          .format(sd.k_1, sd.k_2, sd.k_eff_comp, fm="8.4f"))
                    print("Effective spring stiffness of components on tension side:")
                    print("   k_eq = {0:{fm}} [mm] (Total effective spring stiffness of tension components)\n"
                          "   z_eq = {1:{fm}} [mm] (Total effective lever arm of tension components)"
                          .format(sd.k_eq, sd.z_eq, fm="8.4f"))
                    print()

                    print("Connection strength classification: " + sd.strength_class)
                    print("   Moment capacity of endplate connection M_j_Rd = {0:6.2f} [kNm]\n"
                          "   (M_rigid = {1:.2f} [kNm], M_pinned = {2:.2f} [kNm])"
                          .format(sd.M_j_Rd * 1.0e-6, sd.M_rigid * 1.0e-6, sd.M_pinned * 1.0e-6))
                    print("Connection stiffness classification: " + sd.stiffness_class)
                    print("   Rotational stiffness of endplate connection S_j_ini = {0:.2f} [kNm/rad]\n"
                          "   Rotational stiffness of endplate connection S_j = {1:.2f} [kNm/rad]"
                          .format(sd.S_j_ini * 1.0e-6, sd.S_j * 1.0e-6))
                    print("   (S_j_rigid = {0:.2f} [kNm/rad], S_j_pinned = {1:.2f} [kNm/rad])"
                          .format(sd.S_j_rigid * 1.0e-6, sd.S_j_pinned * 1.0e-6))

                    if cost_info:
                        print('\nUsed cost function: ' + self.cost_function)
                        print("Total costs on side {0}: {1:{fm}} [e]".format(beam_id, sd.cost, fm="6.2f"))
                        print("  Material costs:")
                        print("   Plate:          {0:{fm}} [e]".format(sd.material_cost.plate, fm="6.2f"))
                        print("   Bolts:          {0:{fm}} [e]".format(sd.material_cost.bolts, fm="6.2f"))
                        print("  Work step costs:")
                        print("   Blasting:       {0:{fm}} [e]".format(sd.blasting.cost, fm="6.2f"))
                        print("   Cutting:        {0:{fm}} [e]".format(sd.cutting.cost, fm="6.2f"))
                        # print("   Drilling:       {0:{fm}} [e]".format(sd.drilling.cost, fm="6.2f"))
                        print("   Part assembly:  {0:{fm}} [e]".format(sd.part_assembly.cost, fm="6.2f"))
                        print("   Post treatment: {0:{fm}} [e]".format(sd.post_treatment.cost, fm="6.2f"))
                        print("   Painting:       {0:{fm}} [e]".format(sd.painting.cost, fm="6.2f"))

                    print('\nConnection utilization ratios:')
                    print('Side ' + str(beam_id))
                    print('   M_Ed/M_Rd = {0:{fm}} (Connection moment utilization rate)'.format(sd.n_M, fm=fform))
                    print('   V_Ed/V_Rd = {0:{fm}} (Connection shear utilization rate)'.format(sd.n_V, fm=fform))
            print('V_wp_Ed/V_wp_Rd = {0:{fm}} (Column web shear utilization rate)'
                  .format(self.n_V_wp, fm=fform))
            print()

        if warnings_and_errors:
            if len(self.warnings) == 0:
                print("\nNo warnings encountered during calculations")
            else:
                print("\nWarnings encountered during calculations:")
                for i in self.warnings:
                    print("   " + i)
            if len(self.errors) == 0:
                print("\nNo errors encountered during calculations")
            else:
                print("\nErrors encountered during calculations:")
                for i in self.errors:
                    print("   " + i)

        if figures:
            if len(self.side[0].beam.size) != 0:
                self.print_geom(0, side, front)
            if len(self.side[1].beam.size) != 0:
                self.print_geom(1, side, front)

class Side:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self):
        self.beam = Beam()

        self.beta = 1.0

        self.welds = ""
        self.weld_f = 0.0
        self.weld_w = 0.0

        self.plate_width = 0.0
        self.plate_height = 0.0
        self.plate_thickness = 0.0
        self.upper_overhang = 0.0
        self.plate_material = ""
        self.lower_overhang = 0.0

        self.d = 0
        self.d_0 = 0.0
        self.d_w = 0.0
        self.e = 0.0
        self.e_min = 0.0

        self.compression_center = 0.0
        self.k_eff_comp = 1.0e-12
        self.k_1 = 1.0e-12
        self.k_2 = 1.0e-12
        self.k_eq = 1.0e-12
        self.z_eq = 1.0e12

        self.row = []
        self.groups = []
        self.row_group = []

        self.F_c_wc_Rd = 0.0                # Column web in transverse compression
        self.F_c_fb_Rd = 0.0                # Beam web/flange in compression
        self.F_comp_Rd = 0.0                # Design resistance of compression side

        self.M_j_Rd = 0.0                   # Bending capacity of connection, stronger axis
        self.M_rigid = 0.0                  # Full-strength connection limit
        self.M_pinned = 0.0                 # Pinned connection limit
        self.strength_class = ""            # Connection strength class, stronger axis

        self.S_j_ini = 0.0                  # Initial rotational stiffness of connection
        self.S_j_rigid = 0.0                # Rotational stiffness limit for rigid connection
        self.S_j_pinned = 0.0               # Rotational stiffness limit for pinned connection
        self.my = 0.0                       # Stiffness ratio, takes loading into account
        self.S_j = 0.0                      # Rotational stiffness of connection, takes loading into account
        self.stiffness_class = ""           # Rotational stiffness classification of connection

        self.n_M = 0.0                      # Utilization rate of side under moment
        self.n_V = 0.0                      # Utilization rate of side under shear

        self.stiffener = "unstiffened"

        self.t_bp = 0.0                     # Backing plate thickness
        self.f_y_bp = 0.0                   # Yield strength of backing plate material

        self.material_cost = MaterialCost()
        self.blasting = None
        self.cutting = None
        self.drilling = None
        self.part_assembly = None
        self.post_treatment = None
        self.painting = None

        self.cost = 0.0


class Row:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, cn, sd, z=0, e=0, d=0, bolt_material=8.8,
                 d_washer=0.0, t_washer=0.0, thread_in_shear_plane=True, shear_row=False):

        self.z = float(z)
        self.d = float(d)
        self.temp_d = 0.0
        self.e = float(e)
        self.bolt = ""

        self.bolt_material = bolt_material
        self.thread_in_shear_plane = thread_in_shear_plane
        self.shear_row = shear_row

        # Washer dimendions
        self.d_washer = d_washer
        self.t_washer = t_washer

        self.n_row = 2.0  # Number of bolts in row
        self.n_rows = 1.0  # Number of bolt rows, helping variable

        self.location_p = ''            # Row location individually
        self.location_f = ''            # Row location individually
        self.location_pg = []           # Row location in row group
        self.location_fg = []           # Row location in row group
        self.h_r = 0.0
        self.id = 0

        self.m_x = 0.0
        self.e_x = 0.0

        # Initializing effective length on flange side
        self.l_eff_1f = 0.0        # Individual row
        self.l_eff_2f = 0.0        # Individual row
        # Initializing effective length on plate side
        self.l_eff_1p = 0.0        # Individual row
        self.l_eff_2p = 0.0        # Individual row
        # Initializing lambda for calculation of alfa for endplate bending
        self.lambda1 = 0.0
        self.lambda2 = 0.0
        self.alfa = 0.0

        # Initializing effective lengths for yield line modes
        self.l_eff_cpf = []
        self.l_eff_ncf = []
        self.l_eff_cpp = []
        self.l_eff_ncp = []

        # Initializing capacities
        self.F_t_Rd = 0.0                           # Total capacity of row
        self.F_t_wc_Rd = 0.0                        # Capacity of column web tension
        self.F_t_wb_Rd = 0.0                        # Capacity of beam web/flange in tension
        self.F_t_f_Rd = 0.0                         # Capacity of flange in bending
        self.F_t_p_Rd = 0.0                         # Capacity of plate in bending
        self.F_T_Rdf = [0.0, 0.0, 0.0, 0.0]         # Mode capacities in flange bending
        self.F_T_Rdp = [0.0, 0.0, 0.0, 0.0]         # Mode capacities in plate bending
        self.mode_f = "None"                        # Failure mode on flange side
        self.mode_p = "None"                        # Failure mode on plate side
        self.mode = ""                              # Total failure mode
        self.F_t_Ed = 0.0                           # Design tension force of row
        self.n_F_t = 0.0                            # Tension utilization ratio of row

        # Initialize spring stiffnesses of components
        self.k_3 = 0.0                  # Spring stiffness of column web tension
        self.k_4 = 0.0                  # Spring stiffness of bending of column flange
        self.k_5 = 0.0                  # Spring stiffness of bending of end-plate
        self.k_10 = 0.0                 # Spring stiffness of bolts
        self.k_eff = 0.0                # Effective spring stiffness of tension row

        # Initialize cost values for the row
        self.bolt_unit_cost = 0.0
        self.nut_unit_cost = 0.0
        self.washer_unit_cost = 0.0
        self.material_cost = 0.0

        self.update(cn, sd)

    def update(self, cn, sd):

        column = cn.column

        # Bolt nominal size and edge distance according to side information,
        # overrides given value for row if information given for side != 0
        if sd.d != 0:
            self.d = float(sd.d)
        if sd.e != 0:
            self.e = float(sd.e)
        self.bolt = str("M" + str(int(self.d)))

        # Bolt material properties
        bolt_material = self.bolt_material
        self.f_yb = mat_bolt[bolt_material]["f_yb"]
        self.f_ub = mat_bolt[bolt_material]["f_ub"]

        #if self.d not in bolt_sizes:
        #    print("Given bolt size not standard size!\n"
        #          "Used bolt size is M{0} mm".format(self.d))

        if self.d <= 12.0:
            self.d_0 = self.d + 1.0
        elif 12.0 < self.d < 30.0:
            self.d_0 = self.d + 2.0
        else:
            self.d_0 = self.d + 3.0
        sd.d_0 = self.d_0

        # ISO metric hexagon nuts and bolts, structural, plain thread
        self.S = 1.8*self.d             # Max flats dimension of bolt/nut head (Avainväli)
        self.E = 2.0*self.d             # Max points dimension of bolt/nut head
        self.h_nut = 0.925*self.d         # Height of nut head
        self.h_bolt = 0.625*self.d        # Height of bolt head

        self.d_m = (self.S + self.E)/2.0            # Average bolt/nut diameter

        # Effective diameter for T-stub calculations
        self.d_w = max(self.d_m, self.d_washer)
        if sd.d_w == 0.0:
            sd.d_w = self.d_w
        else:
            sd.d_w = min(sd.d_w, self.d_w)

        # Bolt length
        self.L = column.t_f + sd.plate_thickness + 2.0 * self.t_washer + self.h_nut + self.h_bolt

        # Bolt elongation length
        self.L_b = column.t_f + sd.plate_thickness + 2.0 * self.t_washer + 0.5 * (self.h_nut + self.h_bolt)

        # Bolt length and cost, non standard sizes also valid
        if self.d in bolt_size:
            if self.thread_in_shear_plane:
                for i in range(0, len(bolt_unit_cost_full_thread) - 1):
                    if bolt_lengths[i] > self.L:
                        self.L = bolt_lengths[i]
                        if bolt_unit_cost_full_thread[self.L][self.d] != "N/A":
                            break
            else:
                for i in range(0, len(bolt_unit_cost_full_thread) - 1):
                    if bolt_lengths[i] > self.L:
                        self.L = bolt_lengths[i]
                        if bolt_unit_cost_part_thread[i][self.d] != "N/A":
                            break

            # Tensile stress area of bolt
            self.A_s = bolt_size[self.d]["A_s"]
        else:
            for i in range(0, len(bolt_sizes)):
                if bolt_sizes[i] > self.d:
                    if self.thread_in_shear_plane:
                        for j in range(0, len(bolt_unit_cost_full_thread) - 1):
                            if bolt_lengths[j] > self.L:
                                self.L = bolt_lengths[j]
                                if bolt_unit_cost_full_thread[self.L][bolt_sizes[i]] != "N/A":
                                    break
                    else:
                        for j in range(0, len(bolt_unit_cost_full_thread) - 1):
                            if bolt_lengths[j] > self.L:
                                self.L = bolt_lengths[j]
                                if bolt_unit_cost_part_thread[self.L][bolt_sizes[i]] != "N/A":
                                    break

            self.A_s = 0.8*pi*(0.5*self.d)**2.0

        # Nominal area of bolt
        self.A = pi*(self.d/2.0)**2.0

        # Bolt tension capacity, SFS 1993-1-8 table 3.4
        k_2 = 0.9      # Normaalikantainen ruuvi, uppokantaisille k_2=0.63
        if self.thread_in_shear_plane:
            if bolt_material in [4.6, 5.6, 8.8]:
                self.alfa_v = 0.6
            else:
                self.alfa_v = 0.5
            self.F_T_Rd = self.n_row*k_2*self.f_ub*self.A_s/gamma_M[2]     # Tension capacity of row
            self.V_Rd = self.n_row*self.alfa_v*self.f_ub*self.A_s/gamma_M[2]    # Shear capacity of row
        elif not self.thread_in_shear_plane:
            self.alfa_v = 0.6
            self.F_T_Rd = self.n_row*k_2*self.f_ub*self.A/gamma_M[2]       # Tension capacity of row
            self.V_Rd = self.n_row*self.alfa_v*self.f_ub*self.A/gamma_M[2]      # Shear capacity of row

        # Minimum and maximum spacing, end and edge distances, Table 3.3
        # Not slotted holes and steel not exposed to weather or other corrosive influences
        self.e_1_min = 1.2*self.d_0
        self.e_2_min = 1.2*self.d_0
        self.p_1_min = 2.2*self.d_0
        self.p_2_min = 2.4*self.d_0

        self.p_1_max = min(14.0*column.t_f, 14.0*sd.plate_thickness, 200.0)
        self.p_2_max = min(14.0*column.t_f, 14.0*sd.plate_thickness, 200.0)

        # Horizontal center distance of bolts in plate
        self.w = sd.plate_width - 2.0*self.e
        # Edge distance of bolts to web
        self.mc = 0.5*self.w - 0.5*column.t_w - 0.8*column.r                       # column
        self.m = 0.5*self.w - 0.5*sd.beam.t_w - 0.8*sqrt(2.0)*sd.weld_w            # beam
        # Horizontal edge distance of bolts in column flange
        self.ec = 0.5*(column.b - self.w)
        # Smallest edge distance (flange or plate)
        self.e_min = min(self.e, self.ec)

        if sd.e_min == 0.0:
            sd.e_min = self.e_min
        else:
            sd.e_min = min(sd.e_min, self.e_min)

        # Smallest edge distance in vertical direction in column flange
        self.e_1 = self.z + cn.column_end_length

    def geom_check(self):
        print(self.NOT_IMPLEMENTED + "Geometry check of row: edge distances")

    def info(self, fform=".2f", figures=False):
        print("\nRow {0}: z = {1:{fm}}, bolt = {2}"
              "\n-------------------------------------------\n"
              "   location on plate side = {3}\n"
              "   location on flange side = {4}"
              .format(self.id+1, self.z, self.bolt, self.location_p, self.location_f, fm=fform))
        print("   w = {0:{fm}}, h_r = {1:{fm}}".format(self.w, self.h_r, fm=fform))
        print("   e = {0}, m = {1:{fm}}".format(self.e, self.m, fm=fform))
        print("   ec = {0}, mc = {1:{fm}}".format(self.ec, self.mc, fm=fform))
        print("   e_min = {0:{fm}}".format(self.e_min, fm=fform))
        print("   e_x = {0:{fm}}, m_x = {1:{fm}}".format(self.e_x, self.m_x, fm=fform))
        print("   e_1 = {0:{fm}}".format(self.e_1, fm=fform))
        print("   L = {0:{fm}}, L_b = {1:{fm}}".format(self.L, self.L_b, fm=fform))
        print("   l_eff_1f = {0:{fm}} mm, l_eff_2f = {1:{fm}} mm".format(self.l_eff_1f,
                                                                         self.l_eff_2f, fm=fform))
        print("   l_eff_1p = {0:{fm}} mm, l_eff_2p = {1:{fm}} mm".format(self.l_eff_1p,
                                                                         self.l_eff_2p, fm=fform))
        print("   l_eff_cpf = {0} mm, l_eff_ncf = {1} mm".format(self.l_eff_cpf, self.l_eff_ncf, fm=fform))
        print("   l_eff_cpp = {0} mm, l_eff_ncp = {1} mm".format(self.l_eff_cpp, self.l_eff_ncp, fm=fform))

        if self.location_p == "First bolt-row below tension flange of beam":
            print("      lambda1 = {0:{fm}}, lambda2 = {1:{fm}}".format(self.lambda1, self.lambda2, fm="5.4f"))
            print("      alpha = {0:{fm}}".format(self.alfa, fm="5.4f"))

            if figures:
                plt.figure("alfa")
                plt.xlabel('$\\lambda$1')
                plt.ylabel('$\\lambda$2')

                for i in range(0,280):
                    plt.scatter(lambda11(i/200.0, self.alfa), i/200.0, s=1.0, c="r")
                plt.scatter(self.lambda1, self.lambda2, marker="o", s=5.0, c="b")

                # For image with alfa values shown
                img = mpimg.imread('plot_alfa_2.png')
                plt.imshow(img, extent=[0.0, 1.053, 0.0, 1.59])
                # For image with only contour lines of alfa shown
                #img = mpimg.imread('plot_alfa.png')
                #plt.imshow(img, extent=[0.0, 0.9, 0.0, 1.4])
                plt.draw()
                plt.pause(0.00000001)

        print("Total effective spring stiffness of row is k_eff = {0:{fm}} [mm]".format(self.k_eff, fm="8.4f"))
        print("   k_3 = {0:{fm}} [mm] (Column web tension)".format(self.k_3, fm="8.4f"))
        print("   k_4 = {0:{fm}} [mm] (Column flange bending)".format(self.k_4, fm="8.4f"))
        print("   k_5 = {0:{fm}} [mm] (End-plate bending)".format(self.k_5, fm="8.4f"))
        print("   k_10 = {0:{fm}} [mm] (Bolts)".format(self.k_10, fm="8.4f"))

        print("Shear capacity of row is {0:{fm}} kN".format(self.V_Rd*1.0e-3, fm=fform))
        print("Tension capacity of row is defined by " + self.mode)
        print("Total tension capacity is {0:{fm}} kN".format(self.F_t_Rd*1.0e-3, fm=fform))

        print("   F_t_wc_Rd = {0:{fm}} kN (Column web in transverse tension)".format(self.F_t_wc_Rd * 1.0e-3,
                                                                                  fm=fform))
        print("   F_t_f_Rd = {0:{fm}} kN (Column flange in bending)".format(self.F_t_f_Rd * 1.0e-3, fm=fform))
        print("      F_T_1_Rd = {0:{fm}} kN\n      F_T_2_Rd = {1:{fm}} kN\n      F_T_3_Rd = {2:{fm}} kN"
              "\n      F_T_1_2_Rd = {3:{fm}} kN".format(
                self.F_T_Rdf[0] * 1.0e-3, self.F_T_Rdf[1] * 1.0e-3,
                self.F_T_Rdf[2] * 1.0e-3, self.F_T_Rdf[3] * 1.0e-3, fm=fform))
        print("      Breaking mode of flange in bending: {0}".format(self.mode_f))

        print("   F_t_wb_Rd = {0:{fm}} kN (Beam web in tension)".format(self.F_t_wb_Rd * 1.0e-3, fm=fform))
        print("   F_t_p_Rd = {0:{fm}} kN (End-plate in bending)".format(self.F_t_p_Rd * 1.0e-3, fm=fform))
        print("      F_T_1_Rd = {0:{fm}} kN\n      F_T_2_Rd = {1:{fm}} kN\n      F_T_3_Rd = {2:{fm}} kN"
              "\n      F_T_1_2_Rd = {3:{fm}} kN".format(
                self.F_T_Rdp[0] * 1.0e-3, self.F_T_Rdp[1] * 1.0e-3,
                self.F_T_Rdp[2] * 1.0e-3, self.F_T_Rdp[3] * 1.0e-3, fm=fform))
        print("      Breaking mode of end-plate in bending: {0}\n".format(self.mode_p))


class T_stub:

    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature"
    ERROR = "Error in defining T-stub for bolt row bolt or row group!"

    def __init__(self, row, component, rg, group_id=0, stiffener="unstiffened", b_p=0.0, p=0.0):
        e = row.e
        ec = row.ec
        mc = row.mc
        e_1 = row.e_1
        w = row.w
        m_x = row.m_x
        e_x = row.e_x

        l_eff_nc = 0.0
        l_eff_cp = 0.0

        if component == "flange":
            m = row.mc
            # Effective length of column flange in bending
            if stiffener == "unstiffened" and rg == "row":
                location = row.location_f
                if location == "Bolt row outside end-plate":
                    l_eff_nc = 0.0
                    l_eff_cp = 0.0
                elif location == "Inner bolt-row":
                    l_eff_cp = 2.0*pi*mc
                    l_eff_nc = 4.0*mc + 1.25*ec
                elif location == "End bolt-row":
                    l_eff_cp = min(2.0*pi*mc, pi*mc+2.0*e_1)
                    l_eff_nc = min(4.0*mc+1.25*ec, 2.0*mc+0.625*ec+e_1)

                if l_eff_nc > l_eff_cp:
                    row.l_eff_1f = l_eff_cp
                    row.yield_mode_f = "Non-circular"
                else:
                    row.l_eff_1f = l_eff_nc
                    row.yield_mode_f = "Circular"
                row.l_eff_2f = l_eff_nc

            elif stiffener == "unstiffened" and rg == "row group":
                location = row.location_fg[group_id]
                if location == "Bolt row outside end-plate":
                    l_eff_nc = 0.0
                    l_eff_cp = 0.0
                elif location == "Inner bolt-row":
                    l_eff_cp = 2.0*p
                    l_eff_nc = p
                elif location == "End bolt-row":
                    l_eff_cp = min(pi*mc+p, 2.0*e_1+p)
                    l_eff_nc = min(2.0*mc+0.625*ec+0.5*p, e_1+0.5*p)

                row.l_eff_ncf[group_id] = l_eff_nc
                row.l_eff_cpf[group_id] = l_eff_cp

            elif stiffener == "stiffened":
                print(self.NOT_IMPLEMENTED + "Effective length of flange side with stiffeners.")

        elif component == "plate":
            m = row.m
            # Effective length of beam end-plate in bending
            if rg == "row":
                location = row.location_p
                if location == "Bolt row outside end-plate":
                    l_eff_nc = 0.0
                    l_eff_cp = 0.0
                elif location == "Bolt-row outside tension flange of beam":
                    l_eff_cp = min(2.0*pi*m_x, pi*m_x + w, pi*m_x + 2.0*e)
                    l_eff_nc = min(4.0*m_x + 1.25*e_x, e*2.0*m_x + 0.625*e_x, 0.5*b_p, 0.5*w + 2.0*m_x + 0.625*e_x)
                elif location == "First bolt-row below tension flange of beam":
                    row.lambda1 = m/(m+e)
                    row.lambda2 = row.m_2/(m+e)
                    row.alfa = par_alfa(row.lambda1, row.lambda2)
                    l_eff_cp = 2.0*pi*m
                    l_eff_nc = row.alfa*m
                elif location == "Other inner bolt-row":
                    l_eff_cp = 2.0*pi*m
                    l_eff_nc = 4.0*m + 1.25*e
                elif location == "Other end bolt-row":
                    l_eff_cp = 2.0*pi*m
                    l_eff_nc = 4.0*m + 1.25*e

                if l_eff_nc > l_eff_cp:
                    row.l_eff_1p = l_eff_cp
                    row.yield_mode_p = "Non-circular"
                else:
                    row.l_eff_1p = l_eff_nc
                    row.yield_mode_p = "Circular"
                row.l_eff_2p = l_eff_nc

            elif rg == "row group":
                location = row.location_pg[group_id]
                if location == "Bolt row outside end-plate":
                    l_eff_nc = 0.0
                    l_eff_cp = 0.0
                elif location == "Bolt-row outside tension flange of beam":
                    l_eff_cp = 0.0
                    l_eff_nc = 0.0
                elif location == "First bolt-row below tension flange of beam":
                    # Value for parameter alfa already defined in individual row information
                    l_eff_cp = pi*m + p
                    l_eff_nc = 0.5*p + row.alfa*m - (2*m + 0.625*e)
                elif location == "Other inner bolt-row":
                    l_eff_cp = 2.0 * p
                    l_eff_nc = p
                elif location == "Other end bolt-row":
                    l_eff_cp = pi * m + p
                    l_eff_nc = 2.0*m + 0.625*e + 0.5*p

                row.l_eff_ncp[group_id] = l_eff_nc
                row.l_eff_cpp[group_id] = l_eff_cp


class RowGroup:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature"

    def __init__(self, endplate, i):
        # Initializing effective lengths on flange side
        self.l_eff_1f = 0.0         # Row group
        self.l_eff_2f = 0.0         # Row group
        # Initializing effective lengths on plate side
        self.l_eff_1p = 0.0         # Row group
        self.l_eff_2p = 0.0         # Row group

        # Initializing capacities
        self.F_t_wc_Rd = 0.0
        self.F_t_wb_Rd = 0.0
        self.F_t_f_Rd = 0.0
        self.F_t_p_Rd = 0.0
        self.F_T_Rdf = [0.0, 0.0, 0.0, 0.0]
        self.F_T_Rdp = [0.0, 0.0, 0.0, 0.0]
        self.F_t_Rd = 0.0

        self.mode_f = "None"
        self.mode_p = "None"
        self.mode = "None"

        self.id = i
        self.row = []
        self.p = []

        tem_j = 0
        self.n_rows = 0
        self.group = endplate.groups[i]
        for j in range(0, len(self.group)):
            if self.group[j] == 1:
                self.n_rows = self.n_rows + 1.0
                tem_j = j
                self.row.append(endplate.row[j])

        # Row group info according to last bolt in group
        self.bolt = endplate.row[tem_j].bolt
        self.A_s = endplate.row[tem_j].A_s
        self.A = endplate.row[tem_j].A
        self.L_b = endplate.row[tem_j].L_b
        self.mc = endplate.row[tem_j].mc
        self.m = endplate.row[tem_j].m
        self.e_min = endplate.row[tem_j].e_min
        self.d_w = endplate.row[tem_j].d_w
        self.F_T_Rd = endplate.row[tem_j].F_T_Rd
        self.thread_in_shear_plane = endplate.row[tem_j].thread_in_shear_plane
        self.n_row = endplate.row[tem_j].n_row

    def info(self, fform=".2f"):
        print("\nRow group {0}, p = {1} mm, {2}".format(self.id + 1, self.p, self.group))
        print("-----------------------------------------------------------")
        for j in self.row:
            print("Row id = {0}, z = {1} mm, location_fg = {2}, location_pg = {3}"
                  .format(j.id, j.z, j.location_fg[self.id], j.location_pg[self.id]))
        print("   l_eff_1f = {0:{fm}} mm, l_eff_2f = {1:{fm}} mm".format(self.l_eff_1f,
                                                                         self.l_eff_2f, fm=fform))
        print("   l_eff_1p = {0:{fm}} mm, l_eff_2p = {1:{fm}} mm".format(self.l_eff_1p,
                                                                         self.l_eff_2p, fm=fform))

        print("Tension capacity of row group is defined by {0}\nTotal tension capacity is {1:{fm}} kN"
              .format(self.mode, self.F_t_Rd*1.0e-3, fm=fform))

        print("F_t_wc_Rd = {0:{fm}} kN (Column web in transverse tension)".format(self.F_t_wc_Rd*1.0e-3,
                                                                                  fm=fform))
        print("F_t_f_Rd = {0:{fm}} kN (Column flange in bending)".format(self.F_t_f_Rd*1.0e-3, fm=fform))
        print(
            "   F_T_1_Rd = {0:{fm}} kN\n   F_T_2_Rd = {1:{fm}} kN\n   F_T_3_Rd = {2:{fm}} kN"
            "\n   F_T_1_2_Rd = {3:{fm}} kN".format(
                self.F_T_Rdf[0]*1.0e-3, self.F_T_Rdf[1]*1.0e-3,
                self.F_T_Rdf[2]*1.0e-3, self.F_T_Rdf[3]*1.0e-3, fm=fform))
        print("   Breaking mode of flange in bending: {0}".format(self.mode_f))

        print("F_t_wb_Rd = {0:{fm}} kN (Beam web in tension)".format(self.F_t_wb_Rd*1.0e-3, fm=fform))
        print("F_t_p_Rd = {0:{fm}} kN (End-plate in bending)".format(self.F_t_p_Rd*1.0e-3, fm=fform))
        print(
            "   F_T_1_Rd = {0:{fm}} kN\n   F_T_2_Rd = {1:{fm}} kN\n   F_T_3_Rd = {2:{fm}} kN"
            "\n   F_T_1_2_Rd = {3:{fm}} kN".format(
                self.F_T_Rdp[0] * 1.0e-3, self.F_T_Rdp[1]*1.0e-3,
                self.F_T_Rdp[2] * 1.0e-3, self.F_T_Rdp[3]*1.0e-3, fm=fform))
        print("   Breaking mode of end-plate in bending: {0}\n".format(self.mode_p))


class ShearBlock:
    NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

    def __init__(self, h, b, t):
        print(self.NOT_IMPLEMENTED + 'Shear block!')
        self.h = h              # Height of block
        self.b = b              # Width of block
        self.t = t              # Thickness of block

        self.V_Rd = 0.0         # Shear capacity block


# Imported libraries
from classes import EndPlate
from timeit import default_timer as timer
import matplotlib.pyplot as plt
from write_excel import write_to_excel


plt.close("all")
plt.ion()

start = timer()

# Main program code
column_config = "continuous column"     # "continuous column" or "column end"
column_end_length = 0                   # Column length over highest point of end-plate [mm],
                                        # If column_config="continuous column" and given value column_end_length=0
                                        # value column_end_length=4500.0 is used in calculations
                                        # Not implemented yet: "none"(= beam to beam connection)
sway = "sway"
cost_function = 'Diaz'                  # Used cost function: 'Haapio' or 'Diaz'
test_method = "VT"                      # Test method for final inspection of connection
                                        # VT = Visual testing, UT = Ultrasonic testing, MT = Magnetic particle testing
paint = "alkyd"                         # Paint used in final painting of parts: alkyd, epoxy, polyurethane or acryl

Connection1 = EndPlate(column_config, column_end_length, sway, cost_function, test_method, paint)

# Inputs: L [mm]
# Inputs for standard size: size, material = "S355"
# Inputs for welded I-profile [mm]: h, b, t_f, t_w, r=a
# Additional inputs: N [kN], Vy [kN], Vz [kN], Mt [kNm], My [kNm], Mz [kNm]
# Strong direction in shear: Z, Strong direction under moment around Y
Connection1.add_column(size="HE 220 B", material="S275", L=4000)
# Inputs: N [kN], Vy [kN], Vz [kN], Mt [kNm], My [kNm], Mz [kNm]
# Strong direction in shear Z, and strong direction under moment around Y
Connection1.add_column_load(side='over')     # Loading over connection
Connection1.add_column_load(side='under')    # Loading under connection

# Inputs: beam_id = 0 or 1, L [mm]
# Inputs for standard size: beam_id, size, material = "S355"
# Inputs for welded I-profile [mm]: h, b, t_f, t_w, r=a
Connection1.add_beam(beam_id=0, size="IPE 500", material="S275", L=6000)
# Inputs: N [kN], Vy [kN], Vz [kN], Mt [kNm], My [kNm], Mz [kNm]
# Strong direction in shear: Z, Strong direction under moment around Y
Connection1.side[0].beam.loading()

b = Connection1.side[0].beam.b
# Inputs: beam_id, plate_width, plate_height, plate_thickness, upper_overhang, material = "S355"
# Additional inputs: lower_overhang
# Using parameter value equ=True forces beam to end-plate welds to be designed equivalent strength with beam
Connection1.add_end_plate(beam_id=0, b=b, t=12, upper_overhang=80, material="S275", equ=True)

Connection1.side[0].d = 20          # Nominal bolt size, overrides nominal bolt size given for individual row
Connection1.side[0].e = 40          # Edge distance, overrides edge distance given for individual row
# Inputs: bolt_z_coord, bolt_edge_dist, bolt_nominal_diameter_d, bolt_material = 8.8
# Additional inputs: d_washer=0.0, t_washer=0.0, thread_in_shear_plane=True, shear_row=False
Connection1.add_row(beam_id=0, z=140)
Connection1.add_row(beam_id=0, z=30)

try:
    Connection1.update()
except:
    Connection1.info(values=True, cost_info=True, figures=True)

end = timer()
print("Elapsed time: " + str(end - start) + " [s]\n")

#write_to_excel(cn=Connection1, file_name='testi', sheet_name='sheet1',
#               path='C:\\Users\LaHie1\Documents\Dippa\Excel taulukot\\')

Connection1.info(values=True, cost_info=False, figures=True)

print("END")
plt.show("front 1")

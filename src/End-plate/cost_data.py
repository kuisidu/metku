# Imported libraries

# File contains cost data of parts and materials used in manufacturing of end-plate connection

work_mins_in_a = 252.0*8.0*60.0                 # Work minutes in year [min]

c_land = 12.0                                   # Cost of the land for the building [e/m^2], Haapio, 3.3.3
c_construction = 900 + c_land                   # Real estate construction cost [e/m^2], Haapio, 4.3.3
RE_maintenance = 72.0                           # Real estate maintenance cost [e/(m^2*a)], Haapio, 4.3.6
c_En_unit = 0.1                                 # Unit cost of electricity [e/kWh], Haapio, 4.3.5

# Hot rolled plate

# Steel basic price [e/tonne]
c_SMBPl = 1169.0

# Steel grade add-on [e/tonne]
def c_SMG(grade):
    if grade == "S235":
        return 0.0
    elif grade == "S355":
        return 35.0
    else:
        return 35.0             # If plate material not S235 or S355 highest grade-add on used


# Thickness add-on [e/tonne]
# Thickness t in [mm]
def c_SMT(t):
    if 5.0 <= t < 6.0:
        return 140.0                    # [e/tonne]
    elif 6.0 <= t < 7.0:
        return 82.0
    elif 7.0 <= t < 8.0:
        return 35.0
    elif 8.0 <= t < 10.0:
        return 23.0
    elif 10.0 <= t < 12.0:
        return 12.0
    elif 12.0 <= t < 15.0:
        return 0.0
    elif 15.0 <= t < 20.0:
        return 0.0
    elif 20.0 <= t < 40.0:
        return 12.0
    elif 40.0 <= t < 50.0:
        return 12.0
    elif 50.0 <= t < 60.0:
        return 47.0
    elif 60.0 <= t < 80.0:
        return 58.0
    else:
        return 140.0                  # If thickness not in defined range


# Quantity add-on
# Mass in [kg]
def c_SMQ(mass):
    if mass < 5000:
        return 64.0                   # [-]
    elif mass < 10000:
        return 34.0
    elif mass < 15000:
        return 8.0
    else:
        return 0.0


# Quantity per thickness add-on
# Mass in [kg]
def c_SMTQ(mass):
    if mass < 2000:
        return 90.0
    elif mass < 3000:
        return 90.0
    elif mass < 5000:
        return 33.0
    else:
        return 0.0


# Ultrasonic inspection add-on (acc. EN 10160)
# Quality classes for flat product body (S0, S1, S2, S3) and and edges (E0, E1, E2, E3, E4)
def c_SMUT(S, E):
    if S == "S0" and E == "E0":
        return 44
    elif  S == "S1" and E == "E1":
        return 48
    elif  S == "S2" and E == "E2":
        return 69
    elif  S == "S3" and E == "E3":
        return 73
    elif  S == "S1" and E == "E2":
        return 76
    else:
        return 0


# Unit prices of profiles [e/kg]
# Hot rolled
def profile_price(profile, grade):
    if profile == "IPE" or profile == "HEA":
        if grade == "S235":
            return 1.5              # [e/kg]
        elif grade == "S355":
            return 1.6              # [e/kg]
        else:
            print("Unit price for grade " + grade + " not defined! Grade S355 unit price used!")
            return 1.6              # [e/kg]
    elif profile == "RHS":
        if grade == "S355":
            return 1.88
        else:
            print("Unit price for grade " + grade + " not defined! Grade S355 unit price used!")
            return 1.88  # [e/kg]
    else:
        print("Given profile " + profile + " not defined! Possible profiles IPE, HEA or RHS!")


# Unit price of bolts, DIN931, Partially threaded
# Hexacon 8.8 hot dip galvanised
# Table values {length [mm]: {nominal_bolt_size [mm]: unit_price [e/pcs]}}
bolt_unit_cost_part_thread = {20:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              25:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              30:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              35:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              40:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              45:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              50:  {12: 0.39,  16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              55:  {12: 0.50,  16: 0.78,  20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              60:  {12: 0.45,  16: 0.83,  20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              65:  {12: 0.56,  16: 0.89,  20: 1.57,  24: "N/A", 30: "N/A", 36: "N/A"},
                              70:  {12: 0.51,  16: 0.94,  20: 1.61,  24: "N/A", 30: "N/A", 36: "N/A"},
                              75:  {12: "N/A", 16: 1.00,  20: 1.70,  24: 2.96,  30: "N/A", 36: "N/A"},
                              80:  {12: 0.58,  16: 1.04,  20: 1.78,  24: 3.01,  30: "N/A", 36: "N/A"},
                              90:  {12: 0.64,  16: 1.16,  20: 1.96,  24: 3.32,  30: 6.99,  36: "N/A"},
                              100: {12: 0.71,  16: 1.26,  20: 2.11,  24: 3.61,  30: 7.40,  36: "N/A"},
                              110: {12: 0.76,  16: 1.41,  20: 2.32,  24: 3.89,  30: 7.77,  36: "N/A"},
                              120: {12: 0.83,  16: 1.53,  20: 2.50,  24: 4.16,  30: 8.14,  36: 17.46},
                              130: {12: 0.92,  16: 1.63,  20: 2.65,  24: 4.35,  30: 8.66,  36: 17.65},
                              140: {12: 0.99,  16: 1.73,  20: 2.83,  24: 4.63,  30: 9.38,  36: 18.24},
                              150: {12: 1.04,  16: 1.86,  20: 2.98,  24: 4.82,  30: 9.85,  36: "N/A"},
                              160: {12: 1.18,  16: 2.00,  20: 3.20,  24: 5.35,  30: 10.13, 36: 19.31},
                              180: {12: 1.25,  16: 2.22,  20: 3.61,  24: 5.91,  30: 11.02, 36: 20.67},
                              200: {12: 1.43,  16: 2.48,  20: 3.98,  24: 6.34,  30: 11.78, 36: 21.84},
                              220: {12: 2.06,  16: 3.15,  20: 6.94,  24: 7.45,  30: 12.47, 36: "N/A"},
                              240: {12: 3.28,  16: 3.38,  20: 7.35,  24: 10.11, 30: "N/A", 36: 24.57},
                              260: {12: 3.57,  16: 3.91,  20: 10.59, 24: 12.97, 30: 18.47, 36: "N/A"},
                              280: {12: "N/A", 16: 5.26,  20: "N/A", 24: 13.63, 30: "N/A", 36: 26.72},
                              300: {12: "N/A", 16: 5.69,  20: 11.43, 24: 14.07, 30: 20.07, 36: "N/A"},
                              'max': 26.72}

# Unit price of bolts, DIN933, Fully threaded
# Hexacon 8.8 hot dip galvanised
# Table values {length [mm]: {nominal_bolt_size [mm]: unit_price [e/pcs]}}
bolt_unit_cost_full_thread = {20:  {12: 0.26,  16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              25:  {12: 0.22,  16: 0.49,  20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              30:  {12: 0.24,  16: 0.52,  20: 0.93,  24: "N/A", 30: "N/A", 36: "N/A"},
                              35:  {12: 0.26,  16: 0.55,  20: 1.28,  24: "N/A", 30: "N/A", 36: "N/A"},
                              40:  {12: 0.30,  16: 0.58,  20: 1.01,  24: 2.58,  30: "N/A", 36: "N/A"},
                              45:  {12: 0.32,  16: 0.61,  20: 1.13,  24: "N/A", 30: "N/A", 36: "N/A"},
                              50:  {12: 0.36,  16: 0.65,  20: 1.14,  24: 2.78,  30: "N/A", 36: "N/A"},
                              55:  {12: "N/A", 16: 0.72,  20: 1.30,  24: 2.83,  30: "N/A", 36: "N/A"},
                              60:  {12: 0.39,  16: 0.73,  20: 1.36,  24: 2.91,  30: 5.73,  36: "N/A"},
                              65:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              70:  {12: 0.45,  16: 1.04,  20: 2.00,  24: 3.04,  30: 6.09,  36: 14.06},
                              75:  {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              80:  {12: 0.64,  16: 1.21,  20: 2.26,  24: 3.30,  30: 6.44,  36: "N/A"},
                              90:  {12: 0.84,  16: 1.51,  20: 2.32,  24: 3.73,  30: 6.70,  36: "N/A"},
                              100: {12: 0.92,  16: 2.40,  20: 2.39,  24: 3.85,  30: 7.88,  36: 16.08},
                              110: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: 8.76,  36: "N/A"},
                              120: {12: 1.27,  16: 2.54,  20: "N/A", 24: 4.80,  30: 12.60, 36: 19.58},
                              130: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              140: {12: 1.36,  16: 5.64,  20: "N/A", 24: "N/A", 30: 13.70, 36: "N/A"},
                              150: {12: "N/A", 16: "N/A", 20: 6.70,  24: 9.80,  30: "N/A", 36: "N/A"},
                              160: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              180: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: 27.56},
                              200: {12: "N/A", 16: "N/A", 20: 7.50,  24: "N/A", 30: 25.30, 36: "N/A"},
                              220: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              240: {12: "N/A", 16: "N/A", 20: "N/A", 24: 16.25, 30: "N/A", 36: "N/A"},
                              260: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              280: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              300: {12: "N/A", 16: "N/A", 20: "N/A", 24: "N/A", 30: "N/A", 36: "N/A"},
                              'max': 27.56}

# Unit price of nut, DIN 934-8/10
# Nut hot dip galvanised
# Table values {nominal_bolt_size [mm], unit_price [e/pcs]}
nut_unit_cost = {12: 0.10, 16: 0.27, 20: 0.54, 24: 0.98, 30: 2.13, 36: 4.45}

# Unit price of washer, DIN 7989
# Washer hot dip galvanised
# Table values {nominal_bolt_size [mm], unit_price [e/pcs]}
washer_unit_cost = {12: 0.30, 16: 0.46, 20: 0.65, 24: 0.82, 30: 1.41, 36: 1.76}

# Unit rate of drill bit
# A530 HSS TiN 4xd
# Haapio, table 8
drill_bit_unit_rate = {14: 0.09, 18: 0.14, 22: 0.16, 26: 0.26, 32: 0.37}

# Cost factors used in cost functions by Diaz et al.
Diaz = {'rho':      7820,        # Steel mass density [kg/m^3]
        'f_c':      1.03,        # Cutting factor that increases labor [-]
        'f_w':      1.4,         # Welding factor that increases labor [-]
        'k_c':      0.307,       # Unit cost of for cutting [e(min]
        'k_c_m_o':  0.0016,      # Unit cost of oxygen [e/l]
        'k_c_m_p':  0.0020,      # Unit cost of propane [e/l]
        'k_h':      0.323,       # Unit cost for hole forming [e/min]
        'k_p':      0.43,        # Unit cost for painting [e/min]
        'k_p_mt':   3.8,         # Unit cost of paint [e/l]
        'k_s':      0.4,         # Unit cost of steel [e/kg]
        'k_w':      0.123,       # Unit cost for welding [e/min]
        'k_w_mt':   1.4,         # Unit cost of weld [e/kg]
        'M_p':      0.15,        # Paint consumption [l/m^2]
        'T_c_ex':   2.0,         # Additional cutting time [min]
        't_ex':     1.4,         # Additional drilling path [cm]
        'T_h_ex':   11.9,        # Additional hole forming time [min]
        'T_p':      7.0,         # Time consumption for painting [min/m^2]
        'T_w_ex':   0.3}         # Additional welding time [min]

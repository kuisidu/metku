from anna_semi_rigid_frame import SemiRigidFrame

# Example 1: 96semi_rigid.pdf, page 7

# Frame properties

# All lengths in mm
storey_height = 3650
bay_length = 7300
storeys = 2
bays = 1
num_elements = 5

# Create semi-rigid frame
SRF = SemiRigidFrame(storeys=storeys,
                     bays=bays,
                     storey_height=storey_height,
                     bay_length=bay_length,
                     num_elements=num_elements)
# Add loads
load_id = 2
# Line loads N/mm
SRF.members["Beam0"].add_line_load(load_id, [-62.5, -62.5], "y")
SRF.members["Beam1"].add_line_load(load_id, [-38, -38], "y")
# Point loads
node1 = num_elements
node2 = num_elements*2
node3 = num_elements*3

#Loads were in KN they should be in N
load_vec1 = [36000,0,0]
load_vec2 = [18000,0,0]
factor = 1.0
SRF.add_point_load(load_id, node1, load_vec1, factor)
SRF.add_point_load(load_id, node2, load_vec2, factor)
# Draw frame

SRF.draw_frame()
SRF.optimize_frame(optimizer='slsqp',maxiter=50, swarmsize=150)

#l optimitzacio pso no funciona no troba starting point factible


"""
# Optimal solution from 96semi_rigid.pdf for example 1

# Column profiles
SRF.members["Column0"].profile = "IPE 400"
SRF.members["Column1"].profile = "IPE 400"
SRF.members["Column2"].profile = "IPE 400"
SRF.members["Column3"].profile = "IPE 400"

# Beam profiles and joint rigidities
SRF.members["Beam0"].profile = "IPE 500"
#SRF.members["Beam0"].alpha_to_Sj([0.59, 0.59])
SRF.members["Beam0"].Sj1 = 5.9 *1e4 # 10^4 kNm/rad
SRF.members["Beam0"].Sj2 = 5.9 *1e4 # 10^4 kNm/rad

SRF.members["Beam1"].profile = "IPE 400"
#SRF.members["Beam1"].alpha_to_Sj([0.66, 0.66])
SRF.members["Beam1"].Sj1 = 3.8 *1e4 # 10^4 kNm/rad
SRF.members["Beam1"].Sj2 = 3.8 *1e4 # 10^4 kNm/rad

SRF.generate_frame()
SRF.calculate()
SRF.draw_frame() 
"""
"""
# Optimized solution

# Column profiles
SRF.members["Column0"].profile = "IPE 360"
SRF.members["Column1"].profile = "IPE 330"
SRF.members["Column2"].profile = "IPE 360"
SRF.members["Column3"].profile = "IPE 330"

# Beam profiles and joint rigidities
SRF.members["Beam0"].profile = "IPE 360"
SRF.members["Beam0"].Sj1 = 8.9 *1e4 # 10^4 kNm/rad
SRF.members["Beam0"].Sj2 = 8.9 *1e4 # 10^4 kNm/rad

SRF.members["Beam1"].profile = "IPE 330"
SRF.members["Beam1"].Sj1 = 7.6 *1e4 # 10^4 kNm/rad
SRF.members["Beam1"].Sj2 = 7.6 *1e4 # 10^4 kNm/rad

SRF.generate_frame()
SRF.calculate()
SRF.draw_frame() 
"""

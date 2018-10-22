from semi_rigid_frame import SemiRigidFrame
# Example 2: 96semi_rigid.pdf, page 8

# Frame properties
storey_height = 3.65
bay_length = 6.1
storeys = 3
bays = 2
num_elements = 2

# Create semi-rigid frame
SRF = SemiRigidFrame(storeys=storeys,
                     bays=bays,
                     storey_height=storey_height,
                     bay_length=bay_length,
                     num_elements=num_elements)
# Add loads
load_id = 2
# Line loads

SRF.members["Beam0"].add_line_load(load_id, [-38, -38], "y")
SRF.members["Beam1"].add_line_load(load_id, [-38, -38], "y")
SRF.members["Beam2"].add_line_load(load_id, [-31, -31], "y")
SRF.members["Beam3"].add_line_load(load_id, [-38, -38], "y")
SRF.members["Beam4"].add_line_load(load_id, [-38, -38], "y")
SRF.members["Beam5"].add_line_load(load_id, [-31, -31], "y")

# Calculate node numbers and create load vectors
node1 =  num_elements
node2 = num_elements*2
node3 = num_elements*3
load_vec1 = [36.0,0,0]
load_vec2 = [36.0,0,0]
load_vec3 = [18.0,0,0]
factor = 1.0
# Point loads
SRF.add_point_load(load_id, node1, load_vec1, factor)
SRF.add_point_load(load_id, node2, load_vec2, factor)
SRF.add_point_load(load_id, node3, load_vec3, factor)

#Draw frame
SRF.draw_frame()
# Import dependencies



def create_structure(L, H0, H1, H2, dx, n):
    
    simple_truss = {
        "H0": H0,
        etc.

    
    truss = Truss2D(simple=simple_truss)
    
    simple_frame = [1, 1, H0, L]
    frame = Frame2D(simple=simple_frame, create_beams=False, supports='fixed') # Tuet, aina fixed?
    
    # Kuormat?
    pass
    
def create_continuous_variable_groups(structure, col_bounds, col_values):

    # HUOM! objectit oltava poikkileikkausolioita
    # mem.cross_section

    groups = []
    
    COL_group = {
        'name': 'Columns',
        'var_type': 'continuous',
        'values': col_values, # [300, 100, 10, 5],
        'bounds': col_bounds, # [[100, 800], [100, 300], [5, 50], [5, 50]],
        'properties': ['h', 'b', 'tf', 'tw'],
        'objects': structure.columns
    }
    
    groups.append(COL_group)
    
    # Ristikon osille samaan tapaan
    truss = structure.truss[0]
    top_chords = 
    
    return groups
    
def create_discrete_variable_groups(structure):
    # TODO
    pass


def create_constraint_groups():
    pass
    


if __name__ == '__main__':

    structure = create_structure()
    var_groups = create_continuous_variable_groups()
    
    problem = StructuralProblem(name="Example",
                                var_groups=var_groups,
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
                                }




import os
from frame2d.frame2d import SteelBeam

def to_robot_3D(filename,
                end_frame1,
                long_frame1,
                main_frame,
                num_frames,
                s,
                L,
                end_frame2=None,
                long_frame2=None,
                stiff_truss=None):
        """  Creates an Aurodesk Robot Structural Analysis .str -file
        
            Parameters
            ----------
            :param filename: Name of the created file
            :param num_frames: Number of frames 
            :param s: Spacing between created frames
            :param brace_profile: Profile for vertical braces
                
            :type filename: string
            :type num_frames: int
            :type s: float
            :type brace_profile: string
        """
        if not end_frame2:
            end_frame2 = end_frame1
        if not long_frame2:
            long_frame2 = long_frame1
        num_frames -= 1
        L1=L/2
        L2=L/2
        nodes = []
        elements = []
        profiles = []
        material = []
        releases = {}
        pointloads = []
        #vert_braces = []
        # 2*end_frames + 2*long_frames + 2 = 5
        for i in range(num_frames+7):
            if i == 0:
                frame = end_frame1                
            elif i == num_frames+2:
                frame = end_frame2
            elif i == num_frames+3:
                frame = long_frame1
            elif i == num_frames+4:
                frame = long_frame2
            elif i == num_frames+5:
                if stiff_truss:
                    frame = stiff_truss
                    # Omit eccentricity elements
                    frame.truss = False 
                else:
                    continue
            elif i == num_frames+6:
                if stiff_truss:
                    frame = stiff_truss
                    # Omit eccentricity elements
                    frame.truss = False 
                else:
                    continue
            else:
                frame = main_frame
            
            # Loop through members
            for member in frame.members.values():
                n1, n2 = member.coordinates
                if i == num_frames+3:
                    if isinstance(member, SteelBeam):
                        n1 = [0, n1[0], n1[1]]
                        n2 = [0, n2[0], n2[1]]
                    else:
                        continue
                elif i == num_frames+4:
                    if isinstance(member, SteelBeam):
                        n1 = [L, n1[0], n1[1]]
                        n2 = [L, n2[0], n2[1]]
                    else:
                        continue
               # Wind truss 
                elif i == num_frames+5:                   
                    if member.mtype=='web':
                        idx1 = int(abs(n1[0]-1e-3) // (main_frame.truss[0].L1))
                        idx2 = int(abs(n2[0]-1e-3)//  (main_frame.truss[0].L2))
                        idxT1 = int(idx1 // 2)
                        idxT2 = int(idx2 // 2)
                        idx1 -= 2*idxT1
                        idx2 -= 2*idxT2
                        if main_frame.truss:   
                            chord1 = main_frame.truss[idxT1].top_chords[idx1]
                            chord2 = main_frame.truss[idxT2].top_chords[idx2]
                        else:
                            beams = [b for b in main_frame.members.values() if b.mtype == "beam"]
                            chord1 = beams[idx1]
                            chord2 = beams[idx2]
                        n1 = [n1[0], 0, chord1.shape(n1[0])]
                        n2 = [n2[0], s, chord2.shape(n2[0])]                
                    else:
                        continue
                # Wind truss
                elif i == num_frames+6:
                    if member.mtype=='web':
                        idx1 = int(abs(n1[0]-1e-3) //  (main_frame.truss[0].L1))
                        idx2 = int(abs(n2[0]-1e-3)//  (main_frame.truss[0].L2))
                        idxT1 = int(idx1 // 2)
                        idxT2 = int(idx2 // 2)
                        idx1 -= 2*idxT1
                        idx2 -= 2*idxT2
                        if main_frame.truss:   
                            chord1 = main_frame.truss[idxT1].top_chords[idx1]
                            chord2 = main_frame.truss[idxT2].top_chords[idx2]
                        else:
                            beams = [b for b in main_frame.members.values() if b.mtype == "beam"]
                            chord1 = beams[idx1]
                            chord2 = beams[idx2]
                        n1 = [n1[0], s*(num_frames+2), chord1.shape(n1[0])]
                        n2 = [n2[0], s*(num_frames+1), chord2.shape(n2[0])]                                    
                    else:
                        continue
                else:
                    n1 = [n1[0], i*s, n1[1]]
                    n2 = [n2[0], i*s, n2[1]]
                    
                if n1 not in nodes:
                    nodes.append(n1)
                if n2 not in nodes:
                    nodes.append(n2)
                elem = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                elements.append(elem)
                splitted_val = member.profile.split(" ")
                profile_type = splitted_val[0]
                """
                # Add columns end node to vertical braces list
                if num_frames != 1 and member.mtype == 'top_chord':
                    if nodes.index(n1) +1 not in vert_braces:
                        vert_braces.append(nodes.index(n1) +1)
                    
                    if nodes.index(n2) +1 not in vert_braces:
                        vert_braces.append(nodes.index(n2) +1)
                """
                
                if profile_type == 'RHS':
                    profile = 'RRHS ' + splitted_val[1]
                elif profile_type == 'HE':
                    profile = splitted_val[0] + splitted_val[2] + splitted_val[1]
                elif profile_type == "SHS":
                    dims = splitted_val[1].split("X")
                    profile = 'RRHS ' + str(dims[0]) + 'X' + str(dims[0]) + 'X' + str(dims[1])
                else:
                    profile = member.profile
                # Webs
                if member.mtype == "web" or member.mtype == "beam" or\
                    member.mtype=="top_chord" or member.mtype == "bottom_chord":
                    Sj1 = min(1e10, member.Sj1)
                    Sj2 = min(1e10, member.Sj2)
                    releases[elements.index(elem) +1] = f'ORIgin RY Hy={Sj1}  END RY Hy={Sj2}'
    
                profiles.append(profile)
                material.append(member.material)
                
                # Eccentricity elements
                if len(member.ecc_coordinates):
                    for coords in member.ecc_coordinates:
                        n1, n2 = coords
                        if num_frames != 1:
                            n1 = [n1[0], i*s, n1[1]]
                            n2 = [n2[0], i*s, n2[1]]
                        if n1 not in nodes:
                            nodes.append(n1)
                        if n2 not in nodes:
                            nodes.append(n2)
                        elem = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                        elements.append(elem)
                        profiles.append("HEA 1000")
                        material.append("S355")
                        #if n1 < n2:
                        #    releases[elements.index(elem) +1] = 'END RY'
                        #else:
                        #    releases[elements.index(elem) +1] = 'ORIgin RY'
                        
                        
    
            
            if frame.truss:
                for truss in frame.truss:
                    for joint in truss.joints.values():
                        # Y joint
                        if len(joint.nodal_coordinates) == 2:
                            n1, n2 = sorted(joint.nodal_coordinates)
                            if num_frames != 1:
                                n1 = [n1[0], i*s, n1[1]]
                                n2 = [n2[0], i*s, n2[1]]
                            if n1 not in nodes:
                                nodes.append(n1)
                            if n2 not in nodes:
                                nodes.append(n2)
                            # Eccentricity element
                            elem = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                            elements.append(elem)
                            profiles.append("HEA 1000")
                            material.append("S355")
                        # K joint
                        elif len(joint.nodal_coordinates) == 5:
                            coords = sorted(joint.nodal_coordinates)
                            n1, n2, n3, n4, n5 = coords
                            if num_frames != 1:
                                n1 = [n1[0], i*s, n1[1]]
                                n2 = [n2[0], i*s, n2[1]]
                                n3 = [n3[0], i*s, n3[1]]
                                n4 = [n4[0], i*s, n4[1]]
                                n5 = [n5[0], i*s, n5[1]]
                            if n1 not in nodes:
                                nodes.append(n1)
                            if n2 not in nodes:
                                nodes.append(n2)
                            if n4 not in nodes:
                                nodes.append(n4)
                            if n5 not in nodes:
                                nodes.append(n5)
                            # Eccentricity elements
                            elem1 = [nodes.index(n1) + 1, nodes.index(n2) + 1]
                            elem2 = [nodes.index(n4) + 1, nodes.index(n5) + 1]
                            elements.append(elem1)
                            elements.append(elem2)
                            profiles.append("HEA 1000")
                            material.append("S355")
                            profiles.append("HEA 1000")
                            material.append("S355")
                            #releases[elements.index(elem1)+1] = "END RY"
                            #releases[elements.index(elem2)+1] = "END RY"
                    else:
                        pass
                    
                    
            for pointload in frame.point_loads.values():
                try:
                    n1 = pointload.coordinate
                    n1 = [n1[0], i*s, n1[1]]
                    idx = nodes.index(n1) + 1
                except ValueError:
                    n1 = pointload.coordinate
                    n1 = [n1[0], i*s, n1[1]]
                    nodes.append(n1)
                    idx = nodes.index(n1) +1
                FX, FZ, MY = pointload.v
                if FX:
                    pointloads.append(f'    {idx} FX={FX}')
                elif FX and FZ:
                    pointloads.append(f'    {idx} FX={FX}  FZ={FZ}')
                elif FX and MY:
                    pointloads.append(f'    {idx} FX={FX}  MY={MY}')
                elif FZ:
                    pointloads.append(f'    {idx} FZ={FZ}')
                elif FZ and MY:
                    pointloads.append(f'    {idx} FZ={FZ}  MY={MY}')
                elif MY:
                    pointloads.append(f'    {idx} MY={MY}')
                elif FX and FZ and MY:
                    pointloads.append(f'    {idx} FX={FX}  FZ={FZ}  MY={MY}')


        with  open(filename + '.str', 'w') as f:
            f.write("ROBOT97 \n")
            if num_frames != 1:
                f.write("FRAme SPAce \n")
            f.write("NUMbering DIScontinuous \n")
            f.write(f'NODes {len(nodes)}  ELEments {len(elements)} \n')
            f.write("UNIts \n")
            f.write("LENgth=mm	Force=N \n")
            f.write("NODes \n")
            for i, node in enumerate(nodes):
                if num_frames != 1:
                    f.write(f'{i+1}   {node[0]}    {node[1]}    {node[2]} \n')
                else:
                    f.write(f'{i+1}   {node[0]}    {node[1]} \n')

            f.write('ELEments \n')
            for i, element in enumerate(elements):
                f.write(f' {i+1}  {element[0]}    {element[1]} \n')

            f.write('PROperties \n')
            for i in range(len(elements)):
                f.write(f' "{material[i]}" \n')
                f.write(f' {i+1}  ')
                f.write(f' {profiles[i]}    \n')

            f.write("SUPports \n")
            # End frames
            if end_frame1.is_calculated:
                for sup in end_frame1.supports.values():
                    n1 = sup.coordinate
                    n1 = [n1[0], 0, n1[1]]
                    idx = nodes.index(n1) + 1
                    if sup.dofs[0] == [-1]:
                        f.write(f'{idx}  UX  \n')
                    elif sup.dofs[0:1] == [-1, -1]:
                        f.write(f'{idx}  UX UZ \n')
                    elif sup.dofs == [-1, -1, -1]:
                        f.write(f'{idx}  \n')
                    elif sup.dofs[1] == [-1]:
                        f.write(f'{idx}  UZ  \n')
                    elif sup.dofs[1:2] == [-1, -1]:
                        f.write(f'{idx}  UZ RY  \n')
                    else:
                        f.write(f'{idx}  RY \n')
            else:
                for sup in end_frame1.supports.values():
                        n1 = sup.coordinate
                        n1 = [n1[0], 0, n1[1]]
                        idx = nodes.index(n1) + 1
                        if sup.dofs == [1, 0, 0]:
                            f.write(f'{idx}  UX  \n')
                        elif sup.dofs == [1, 1, 0]:
                            f.write(f'{idx}  UX UZ \n')
                        elif sup.dofs == [1, 1, 1]:
                            f.write(f'{idx}  \n')
                        elif sup.dofs == [0, 1, 0]:
                            f.write(f'{idx}  UZ  \n')
                        elif sup.dofs == [0, 1, 1]:
                            f.write(f'{idx}  UZ RY  \n')
                        else:
                            f.write(f'{idx}  RY \n')
            if end_frame2.is_calculated:
                for sup in end_frame2.supports.values():
                    n1 = sup.coordinate
                    n1 = [n1[0], (num_frames+2)*s, n1[1]]
                    idx = nodes.index(n1) + 1
                    if sup.dofs[0] == [-1]:
                        f.write(f'{idx}  UX  \n')
                    elif sup.dofs[0:1] == [-1, -1]:
                        f.write(f'{idx}  UX UZ \n')
                    elif sup.dofs == [-1, -1, -1]:
                        f.write(f'{idx}  \n')
                    elif sup.dofs[1] == [-1]:
                        f.write(f'{idx}  UZ  \n')
                    elif sup.dofs[1:2] == [-1, -1]:
                        f.write(f'{idx}  UZ RY  \n')
                    else:
                        f.write(f'{idx}  RY \n')
            else:
                for sup in end_frame2.supports.values():
                        n1 = sup.coordinate
                        n1 = [n1[0], (num_frames+2)*s, n1[1]]
                        idx = nodes.index(n1) + 1
                        if sup.dofs == [1, 0, 0]:
                            f.write(f'{idx}  UX  \n')
                        elif sup.dofs == [1, 1, 0]:
                            f.write(f'{idx}  UX UZ \n')
                        elif sup.dofs == [1, 1, 1]:
                            f.write(f'{idx}  \n')
                        elif sup.dofs == [0, 1, 0]:
                            f.write(f'{idx}  UZ  \n')
                        elif sup.dofs == [0, 1, 1]:
                            f.write(f'{idx}  UZ RY  \n')
                        else:
                            f.write(f'{idx}  RY \n')
                            
            # Longitudal frames
            if long_frame1.is_calculated:
                for sup in long_frame1.supports.values():
                    n1 = sup.coordinate
                    n1 = [n1[1], n1[0], 0]
                    idx = nodes.index(n1) + 1
                    if sup.dofs[0] == [-1]:
                        f.write(f'{idx}  UX  \n')
                    elif sup.dofs[0:1] == [-1, -1]:
                        f.write(f'{idx}  UX UZ \n')
                    elif sup.dofs == [-1, -1, -1]:
                        f.write(f'{idx}  \n')
                    elif sup.dofs[1] == [-1]:
                        f.write(f'{idx}  UZ  \n')
                    elif sup.dofs[1:2] == [-1, -1]:
                        f.write(f'{idx}  UZ RY  \n')
                    else:
                        f.write(f'{idx}  RY \n')
            else:
                for sup in long_frame1.supports.values():
                        n1 = sup.coordinate
                        n1 = [n1[1], n1[0], 0]
                        idx = nodes.index(n1) + 1
                        if sup.dofs == [1, 0, 0]:
                            f.write(f'{idx}  UX  \n')
                        elif sup.dofs == [1, 1, 0]:
                            f.write(f'{idx}  UX UZ \n')
                        elif sup.dofs == [1, 1, 1]:
                            f.write(f'{idx}  \n')
                        elif sup.dofs == [0, 1, 0]:
                            f.write(f'{idx}  UZ  \n')
                        elif sup.dofs == [0, 1, 1]:
                            f.write(f'{idx}  UZ RY  \n')
                        else:
                            f.write(f'{idx}  RY \n')
                            
            if long_frame2.is_calculated:
                for sup in long_frame2.supports.values():
                    n1 = sup.coordinate
                    n1 = [L, n1[0], 0]
                    idx = nodes.index(n1) + 1
                    if sup.dofs[0] == [-1]:
                        f.write(f'{idx}  UX  \n')
                    elif sup.dofs[0:1] == [-1, -1]:
                        f.write(f'{idx}  UX UZ \n')
                    elif sup.dofs == [-1, -1, -1]:
                        f.write(f'{idx}  \n')
                    elif sup.dofs[1] == [-1]:
                        f.write(f'{idx}  UZ  \n')
                    elif sup.dofs[1:2] == [-1, -1]:
                        f.write(f'{idx}  UZ RY  \n')
                    else:
                        f.write(f'{idx}  RY \n')
            else:
                for sup in long_frame2.supports.values():
                        n1 = sup.coordinate
                        n1 = [L, n1[0], 0]
                        idx = nodes.index(n1) + 1
                        if sup.dofs == [1, 0, 0]:
                            f.write(f'{idx}  UX  \n')
                        elif sup.dofs == [1, 1, 0]:
                            f.write(f'{idx}  UX UZ \n')
                        elif sup.dofs == [1, 1, 1]:
                            f.write(f'{idx}  \n')
                        elif sup.dofs == [0, 1, 0]:
                            f.write(f'{idx}  UZ  \n')
                        elif sup.dofs == [0, 1, 1]:
                            f.write(f'{idx}  UZ RY  \n')
                        else:
                            f.write(f'{idx}  RY \n')
            
            f.write("RELeases \n")
            for elem in releases.keys():
                f.write(f'ELEments{elem} {releases[elem]} \n')
            
            f.write("LOAds \n")
            f.write("CASe # 1 LC1 \n")
            if len(main_frame.line_loads):
                f.write("ELEments \n")
                for i in range(num_frames):
                    for lineload in main_frame.line_loads.values():
                        n1, n2 = lineload.member.coordinates
                        if num_frames != 1:
                            n1 = [n1[0], i*s, n1[1]]
                            n2 = [n2[0], i*s, n2[1]]
                        n1_idx = nodes.index(n1) + 1
                        n2_idx = nodes.index(n2) + 1
                        idx = elements.index([n1_idx, n2_idx]) + 1
                        if lineload.direction == 'y':
                            dir = 'PZ'
                        else:
                            dir = 'PX'
                        q0, q1 = lineload.values
                        if q0 != q1:
                            f.write(f' {idx} X=0.0 {dir}={q0} TILl  ')
                            f.write(f'X=1.000  {dir}={q1}      RElative \n')
        
                        else:
                            f.write(f' {idx} {dir}={q0} \n')
            if len(pointloads):
                f.write("NODes \n")
                for val in pointloads:
                    f.write(val + "\n")
                    

            f.write("END")
        
        print(f'{filename}.str created to: \n{os.getcwd()}')
        
        
        
        
        
        

from frame2d.frame2d import Frame2D, FrameMember, FixedSupport, LineLoad

class Truss2D(Frame2D):
    
    def __init__(self, H, L, n, flip=False):
        super().__init__()
        self.H = H
        self.L = L
        self.n = n
        self.flip = flip
        self.joints = {}
        
    def add(self, this):
        super().add(this)
        
        if isinstance(this, TrussJoint):
            if this.jid not in self.joints:
                self.joints[this.jid] = this
        
    def generate_chords(self):
        self.add(FrameMember([[0,self.H],[self.L,self.H]], mtype="chord"))
        self.add(FrameMember([[0,0],[self.L,0]], mtype="chord"))
        
    def generate_webs(self):
        top_nodes = self.n // 2 + 1
        bottom_nodes = top_nodes - 1
        
        dL1 = self.L / bottom_nodes
        dL2 = self.L / bottom_nodes

        x1 = 0
        x2 = dL2/2
        for i in range(self.n):
            if self.flip:
                self.add(FrameMember([[x2,self.H],[x1,0]], mtype="web"))
            else:
                self.add(FrameMember([[x1,self.H],[x2,0]], mtype="web"))
            if i%2 == 0:            
                x1 += dL1
            else:
                x2 += dL2
                
        
        
class TrussJoint():
    def __init__(self, chord, loc, joint_type="N", g1=0.05, g2=0.05):

        self.chord = chord
        self.chord_elements = {}
        self.jid = None
        self.__loc = loc
        self.__coordinate = chord.local(loc)
        self.cnode = None # central node
        self.chord_nodes = []
        self.joint_type = joint_type
        self.__g1 = g1
        self.__g2 = g2
        self.nodal_coordinates = []
        self.nodes = {}
        self.elements = {}
        self.element_ids = []
        self.webs = {}
        self.is_generated = False
        # Index for eccentricity elements' material and section
        self.idx = None
        self.rhs_joint = None
                
    def add_node_coord(self, coord):
        coord = [round(c, PREC) for c in coord]
        if coord not in self.nodal_coordinates:
            self.nodal_coordinates.append(coord)
    
    def add_web(self, web, num_elements):
        self.webs[web.mem_id] = web
        #if self.is_generated:
        self.calc_nodal_coordinates()

    @property
    def loc(self):
        return round(self.__loc, PREC)
    
    @loc.setter
    def loc(self, val):
        self.__loc = round(val, PREC)

    @property
    def coordinate(self):     
        return [round(c, PREC) for c in self.chord.local(self.loc)]

    @coordinate.setter
    def coordinate(self, val):
        # (x, y) -coordinate
        if isinstance(val, list):            
            local = self.chord.global_coord(val)
            self.coordinate = local
        # local coordinate
        elif isinstance(val, float):
            self.loc = val
            self.calc_nodal_coordinates()
                        
    @property
    def g1(self):
        return self.__g1
    
    @g1.setter
    def g1(self, val):
        self.__g1 = val
        self.calc_nodal_coordinates()
        
    @property
    def g2(self):
        return self.__g2
    
    @g2.setter
    def g2(self, val):
        self.__g2 = val
        self.calc_nodal_coordinates()
               
    def calculate(self):
        """
        """
        if self.joint_type == "K":
            self.rhs_joint = RHSKGapJoint(self.chord.cross_section,
                                         [w.cross_section for w in self.webs.values()],
                                         [math.degrees(w.angle) for w in self.webs.values()],
                                         self.g1 * 1000) # m to mm
        elif self.joint_type == 'N':
            pass
        elif self.joint_type == 'Y':
            pass
        elif self.joint_type == 'KT':
            pass
        if self.rhs_joint:
            self.rhs_joint.V0 = max([max(V.shear_force) for V in self.chord_elements.values()])
            # Chord face failure
            r1 = self.rhs_joint.chord_face_failure() / 1000 # N to kN
            # Punching shear
            r2 = self.rhs_joint.punching_shear() / 1000 # N to kN
            # Brace failure
            r3 = self.rhs_joint.brace_failure() / 1000 # N to kN
            
            return [r1, r2, r3]
    
    def chord_face_failure(self):
        if self.rhs_joint:
            results = self.rhs_joint.chord_face_failure() 
            return results / 1000
    
    def punching_shear(self):
        if self.rhs_joint:
            results = self.rhs_joint.punching_shear() 
            return results / 1000
        
    def brace_failure(self):
        if self.rhs_joint:           
            results = self.rhs_joint.brace_failure()
            return results / 1000 # N to kN
    
    def chord_shear(self):
        if self.rhs_joint:
            self.rhs_joint.V0 = max([max(V.shear_force) for V in self.chord_elements.values()])
            results = self.rhs_joint.chord_shear()
            return results 
    
    def calc_nodal_coordinates(self):
        """
        Calculates eccentricity nodes coordinates and adds them to chord's
        node-dict, changes webs end coordinates and saves all calculated
        nodal coordinates to own nodal_coordinates -list
        
        coord1 : eccentricity coordinate on chord
        coord2 : web's new end coordinate
        v : chord's position vector
        u : vector perpendicular to v
        
        1: Calc ecc_y, distance from chord's center line to edge
        2: Calc ecc_x, distance from joint to web's edge
        3: Calc coord1, ecc_x * v
        4: Calc coord2, coord1 + ecc_y * u
        5: Move web to coord2
        6: Save coordinates 
        
        """
        # Chord's position vector
        v = np.array([math.cos(self.chord.angle), math.sin(self.chord.angle)])
        # Vector perpendicular to v
        u = np.array([-math.sin(self.chord.angle), math.cos(self.chord.angle)])
        # Remove previously calculated coordinates from chord
        for coord in self.nodal_coordinates:
            if coord in self.chord.added_coordinates:
                self.chord.added_coordinates.remove(coord)
            if coord in self.chord.nodal_coordinates and \
                coord not in self.chord.coordinates:
                self.chord.nodal_coordinates.remove(coord)
        # Initialize nodal coordinates as an empty array
        self.nodal_coordinates = []
        central_coord = np.asarray(self.coordinate)
        # Eccentricity along y-axis
        ecc_y = self.chord.h / 2000 # div by 2 and 1000 (mm to m)

        # Add joint's central coordinate to joint's and chord's list
        self.add_node_coord(list(central_coord))
        self.chord.add_node_coord(list([central_coord[0], self.chord.shape(central_coord[0])]))
        # If joint type is Y or T
        if len(self.webs) == 1:
            web = list(self.webs.values())[0]
            # angle between chord and web
            theta = abs(self.chord.angle - web.angle)
            # Eccentricity along x-axis
            ecc_x = abs(self.g1 * math.cos(self.chord.angle)
                    + web.h/2000 / math.sin(theta)) 

            ecc_x = round(ecc_x, PREC)
            # TOP CHORD
            if self.chord.mtype == 'top_chord':
                
                if self.loc == 0:
                    coord1 = list(central_coord + v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) - u *ecc_y)
                elif self.loc == 1:
                    coord1 = list(central_coord - v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) - u *ecc_y)
                else:
                    coord1 = list(central_coord)
                    coord2 = list(np.asarray(coord1) - u *ecc_y)

                self.add_node_coord(coord2)
                self.chord.add_node_coord(coord1)
                web.top_coord = coord2
                
            # BOTTOM CHORD
            else:
                if self.loc == 0:
                    coord1 = list(central_coord + v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) + u *ecc_y)
                elif self.loc == 1:
                    coord1 = list(central_coord - v*ecc_x)
                    self.coordinate = coord1
                    coord2 = list(np.asarray(coord1) + u *ecc_y)
                else:
                    coord1 = list(central_coord)
                    coord2 = list(np.asarray(coord1) + u *ecc_y)

                self.add_node_coord(coord2)
                self.chord.add_node_coord(coord1)
                web.bot_coord = coord2
            web.calc_nodal_coordinates()    
            
        # If joint type is K or KT
        elif len(self.webs) >= 2:
            # Iterate trhough every web member in the joint
            for web in self.webs.values():
                # angle between chord and web
                theta = abs(self.chord.angle - web.angle)
                # Eccentricity along x-axis
                ecc_x = abs(self.g1/2 * math.cos(self.chord.angle)
                        + web.h/2000 / math.sin(theta))      
                # If joint type is KT, there's two gap values g1 and g2
                if len(self.webs) == 3:
                    # Find the leftmost, rightmost and middle web with x-values
                    # leftmost is the one with smallest x (miv_val)
                    # rightmost is the one with highest x (max_val)
                    # middle is the one between those two
                    min_val = min([x.bot_coord[0] for x in self.webs.values()])
                    max_val = max ([x.bot_coord[0] for x in self.webs.values()])
                    middle_web = [x for x in self.webs.values() if x.bot_coord[0] != min_val and x.bot_coord[0] != max_val][0]
                    # eccentricity for the leftmost web
                    if web.bot_coord[0] == min_val:
                        ecc_x = self.g1 *math.cos(self.chord.angle) + web.h/2000 / math.sin(theta) + middle_web.h/2000
                    # eccentricity for the rightmost web
                    elif web.bot_coord[0] == max_val:
                        ecc_x = self.g2 *math.cos(self.chord.angle) + web.h/2000 / math.sin(theta) + middle_web.h/2000
                    # eccentricity for the middle web
                    else:
                        ecc_x = 0
                # TOP CHORD
                if self.chord.mtype == 'top_chord':
                    # if the web is the leftmost web
                    if web.bot_coord[0] == min([x.bot_coord[0] for x in self.webs.values()]):
                        ecc_x = -ecc_x                    
                    # Web's eccentricity node on chord's center line
                    coord1 = central_coord + ecc_x*v
                    # Check that coordinate is between chord's coordinates
                    if list(coord1)[0] < self.chord.coordinates[0][0] or\
                        list(coord1)[0] > self.chord.coordinates[1][0]:
                        coord1 = central_coord
                    # Web's eccentricity node on chord's edge
                    coord2 = coord1 - ecc_y*u
                    # Move chord's coordinate to newly calculated location
                    #coord2 = [round(c, 3) for c in coord2]
                    web.top_coord = list(coord2)
                # BOTTOM CHORD
                else:                    
                    # if the web is the leftmost web
                    if web.top_coord[0] == min([x.top_coord[0] for x in self.webs.values()]):
                        ecc_x = -ecc_x
                    # Web's eccentricity node on chord's center line
                    coord1 = central_coord + ecc_x*v
                    # Check that coordinate is between chord's coordinates
                    if list(coord1)[0] < self.chord.coordinates[0][0] or\
                        list(coord1)[0] > self.chord.coordinates[1][0]:
                        coord1 = central_coord
                    # Web's eccentricity node on chord's edge
                    coord2 = coord1 + ecc_y*u
                    # Move chord's coordinate to newly calculated location
                    #coord2 = [round(c, 3) for c in coord2]
                    web.bot_coord = list(coord2)                    
                # Calculate web's new nodal locations 
                web.calc_nodal_coordinates()
                # Add eccentricity nodes' coordinates to joint and chord
                # Node on the chord's center line needs to be exactly on the
                # center line, so we need to calculate the y-coord with
                # chord's shape function
                coord1 = list([coord1[0], self.chord.shape(coord1[0])])
                #coord1 = [round(c,3) for c in coord1]
                # Node coordinate on chord's center line
                self.add_node_coord(coord1)
                self.chord.add_node_coord(coord1)
                # Node coordinate on web's end
                self.add_node_coord(list(coord2))
                
        # Move existing nodes
        if len(self.nodes):
            for i, node in enumerate(self.nodes.values()):                
                node.x = self.nodal_coordinates[i]
          
        # Set joint type
        if len(self.nodal_coordinates) == 0:
            pass
        elif len(self.nodal_coordinates) == 2:
            self.joint_type = 'Y'
        elif len(self.nodal_coordinates) == 3:
            self.joint_type = 'Y'
        elif len(self.nodal_coordinates) == 4:
            self.joint_type = 'N'
        elif len(self.nodal_coordinates) == 5:
            self.joint_type = 'K'
            self.rhs_joint = RHSKGapJoint(self.chord.cross_section,
                                         [w.cross_section for w in self.webs.values()],
                                         [math.degrees(w.angle) for w in self.webs.values()],
                                         self.g1 * 1000) # m to mm
        elif len(self.nodal_coordinates) == 6:
            self.joint_type = 'KT'
        else:
            raise ValueError('Too many nodes for one joint')
            
            
    def round_coordinates(self, prec=PREC):       
        for i, coord in enumerate(self.nodal_coordinates):
            self.nodal_coordinates[i] = [round(c, prec) for c in coord]         
            
    def generate(self, fem_model):
        """
        """
        self.round_coordinates()
        self.generate_nodes(fem_model)
        self.add_section_and_material(fem_model)
        self.generate_eccentricity_elements(fem_model)
        self.get_chord_elements()
        
                
    def generate_nodes(self, fem_model):
        """
        Generates eccentricity nodes to chord
        :param fem_model: FrameFEM -instance, model where nodes are added
        """
        if not self.is_generated:
            self.is_generated = True
            for coord in self.nodal_coordinates:
                idx = fem_model.nodal_coords.index(coord)
                self.nodes[idx] = fem_model.nodes[idx]
                if coord == self.coordinate:
                    self.cnode = fem_model.nodes[idx]
                if coord in self.chord.nodal_coordinates:
                    self.chord_nodes.append(fem_model.nodes[idx])
        # sort chord nodes from leftmost to rightmost         
        self.chord_nodes = sorted(self.chord_nodes, key= lambda node: node.x[0])
                    
    def get_chord_elements(self):
        """
        Gets the joint's elements which are on chord
        """
        if len(self.chord_nodes) > 2:
            for i in range(2):
                nodes = self.chord_nodes[i:i+2]
                for key in self.chord.elements.keys():
                    element = self.chord.elements[key]
                    if nodes[0] in element.nodes and nodes[1] in element.nodes:
                        self.chord_elements[key] = element
        
        
    def add_section_and_material(self, fem_model):
        """
        Add's eccentricity elements section and material properties
        to the fem model.
        Eccentricity elements are modelled as infinitely rigid elements
        """
        self.idx = len(fem_model.sections)
        
        # Material(E, nu, rho)
        fem_model.add_material(210e3, 0.3, 7850e-9)
        # BeamSection(A, Iy)
        sect = BeamSection(38880*1e-6, 6299657109*1e-12)
        fem_model.add_section(sect)
        
        # USE CHORD'S PROPERTIES
        #self.idx = self.chord.mem_id
  
    
    def generate_eccentricity_elements(self, fem_model):
        """
        Generates eccentricity elements
        """
        # Create a list of node instances sorted by node's x-coordinate
        nodes = [x for _,x in sorted(zip(self.nodal_coordinates,
                                         list(self.nodes.values())))]
        chord_nodes = []
        web_nodes = []
        for node in nodes:
            if node in self.chord.nodes.values():
                chord_nodes.append(node)
            else:
                web_nodes.append(node)
                
        web_nodes = sorted(web_nodes, key=lambda node: node.x[0])
        chord_nodes = sorted(chord_nodes, key=lambda node: node.x[0])

        index = len(fem_model.elements)
        if self.joint_type == 'Y':
            n1, n2 = nodes
            self.elements[index] = EBBeam(n1,
                                        n2,
                                        fem_model.sections[self.idx], 
                                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            
        elif self.joint_type == 'N':
            wn1, wn2 = web_nodes
            cn1, cn2 = chord_nodes
            
            self.elements[index] = EBBeam(cn1,
                                         wn1,
                                        fem_model.sections[self.idx], 
                                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            # Chrod and web
            self.elements[index+1] = EBBeam(cn2,
                                            wn2,
                                            fem_model.sections[self.idx], 
                                            fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+1])
            self.element_ids.append(index+1)
            
            
        elif self.joint_type == 'K':
            cn1, cn2, cn3 = chord_nodes
            wn1, wn2 = web_nodes
            # Chrod and web
            self.elements[index] = EBBeam(cn1,
                                        wn1,
                                        fem_model.sections[self.idx], 
                                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            # Chrod and web
            self.elements[index+1] = EBBeam(cn3,
                                            wn2,
                                            fem_model.sections[self.idx], 
                                            fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+1])
            self.element_ids.append(index+1)
            
        elif self.joint_type == 'KT':
            cn1, cn2, cn3 = chord_nodes
            wn1, wn2, wn3 = web_nodes
            # Chrod and web
            self.elements[index] = EBBeam(cn1,
                                        wn1,
                                        fem_model.sections[self.idx], 
                                        fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index])
            self.element_ids.append(index)
            # Central and web
            self.elements[index+1] = EBBeam(cn2,
                                            wn2,
                                            fem_model.sections[self.idx], 
                                            fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+1])
            self.element_ids.append(index+1)
            # Chrod and web
            self.elements[index+2] = EBBeam(cn3,
                                            wn3,
                                            fem_model.sections[self.idx], 
                                            fem_model.materials[self.idx])
            # Add it to fem model
            fem_model.add_element(self.elements[index+2])
            self.element_ids.append(index+2)
            
        
    def plot(self, print_text=True, color=True):

        x, y = self.coordinate
        if color:
            if self.is_strong_enough:
                color = 'green'
            else:
                color = 'red'
        else:
            color = 'cyan'

        horzalign = 'center'
        vertalign = 'center'
        plt.scatter(x, y,s=100, c=color)
        plt.text(x, y, str(self.jid), horizontalalignment=horzalign,
                 verticalalignment=vertalign, fontsize=12)
        
    def plot_joint(self, length=0.3, show=True):
        """
        Plots joint
        """
        X, Y = self.coordinate  
        end_line = None
        col = None
        # CHORD
        # Chord's position vector
        v = np.array([math.cos(self.chord.angle), math.sin(self.chord.angle)])
        # Vector perpendicular to v
        u = np.array([-math.sin(self.chord.angle), math.cos(self.chord.angle)])
        if X-length*2 >= self.chord.local(0)[0]:
            start = X - length*2
        else:
            start = self.chord.coordinates[0][0]
            end_line = self.chord.coordinates[0][0]
        
        if X + length*2 <= self.chord.local(1)[0]:
            end = X + length*2
        else:
            end = self.chord.coordinates[1][0]  
            end_line = self.chord.coordinates[1][0]
            
        X0 = [start, X, end]
        Y0 = [Y-(X-start)*math.tan(self.chord.angle),
              Y,
              Y+(end-X)*math.tan(self.chord.angle)]
        X1 = [start, X, end]
        Y1 = [Y-(X-start)*math.tan(self.chord.angle) - self.chord.h/2000,
              Y-self.chord.h/2000,
              Y+(end-X)*math.tan(self.chord.angle)-self.chord.h/2000]
        X2 = [start, X,end]
        Y2 = [Y-(X-start)*math.tan(self.chord.angle)+self.chord.h/2000,
              Y+self.chord.h/2000,
              Y+(end-X)*math.tan(self.chord.angle)+self.chord.h/2000]
        
        plt.plot(X0, Y0, 'k--')
        plt.plot(X1, Y1, 'k')
        plt.plot(X2, Y2, 'k')
        if end_line:
            plt.plot([end_line, end_line],[Y-self.chord.h/2000, Y+self.chord.h/2000],'k')
            
        # COLUMN
        if start in self.chord.columns.keys():
            col = self.chord.columns[start]
            end = start - col.h/2000
            inner = start + col.h/2000
            Y1 = Y- length
            Y2 = Y + length
            if Y2 >= col.coordinates[1][1]:
                Y2 = col.coordinates[1][1]
            # outer
            plt.plot([end, end], [Y1, Y2], 'k')
            # mid
            plt.plot([start, start],[Y1, Y2], 'k--')
            # inner
            plt.plot([inner, inner], [Y1, Y2], 'k')
            
        elif end in self.chord.columns.keys():
            col = self.chord.columns[end]
            out = end + col.h/2000
            inner = end - col.h/2000
            Y1 = Y- length
            Y2 = Y + length
            if Y2 >= col.coordinates[1][1]:
                Y2 = col.coordinates[1][1]
            # outer
            plt.plot([out, out], [Y1, Y2], 'k')
            # mid
            plt.plot([end, end],[Y1, Y2], 'k--')
            # inner
            plt.plot([inner, inner], [Y1, Y2], 'k')
            
        
        # WEBS
        for web in self.webs.values():
            theta = abs(self.chord.angle - web.angle)
            if self.chord.mtype == "top_chord":
                coord = web.top_coord
                k = -1
            else:
                coord = web.bot_coord
                k = 1
            if web.angle < 0:
                X0 = [coord[0], coord[0] - k*length]
                Y0 = [coord[1], coord[1] - k*(math.tan(web.angle)*length)]
                X1 = [coord[0]+web.h/2000 / math.sin(theta) , coord[0]+web.h/2000 / math.sin(theta) - k*length]
                Y1 = [coord[1], coord[1] - k*(math.tan(web.angle)*length)]
                X2 = [coord[0]-web.h/2000 / math.sin(theta) , coord[0]-web.h/2000 / math.sin(theta) - k*length]
                Y2 = [coord[1], coord[1] - k*(math.tan(web.angle)*length)]
            
            elif math.degrees(web.angle) == 90:
                X0 = [coord[0], coord[0]]
                Y0 = [coord[1], coord[1] + k*2*length]
                X1 = [coord[0]+web.h/2000 , coord[0]+web.h/2000]
                Y1 = [coord[1], coord[1] + k*2*length]
                X2 = [coord[0]-web.h/2000 , coord[0]-web.h/2000]
                Y2 = [coord[1], coord[1] + k*2*length]
            else:
                X0 = [coord[0], coord[0] + k*length]
                Y0 = [coord[1], coord[1] + k*(math.tan(web.angle)*length)]
                X1 = [coord[0] + web.h/2000 / math.sin(theta) , coord[0]+web.h/2000 / math.sin(theta) + k*length]
                Y1 = [coord[1], coord[1] + k*(math.tan(web.angle)*length)]
                X2 = [coord[0] - web.h/2000 / math.sin(theta) , coord[0]-web.h/2000 / math.sin(theta) + k*length]
                Y2 = [coord[1], coord[1] + k*(math.tan(web.angle)*length)]

            plt.plot(X0, Y0, 'k--')
            plt.plot(X1, Y1, 'k')
            plt.plot(X2, Y2, 'k')
            
        # NODES
        for key in self.nodes.keys():
            node = self.nodes[key]
            x, y = node.x
            plt.scatter(x, y, s=100, c='pink')
            plt.text(x, y, str(key), horizontalalignment='center',
                     verticalalignment='center', fontsize=15)
            if x < X:
                x -= length
            elif x > X:
                x += length/1.5
            if y < Y:
                y -= length
            elif y > Y:
                y += length

            #plt.text(x, y+0.08, f'Fx {node.Fx:.0f} kN', fontsize=11)
            #plt.text(x, y+0.05, f'Fy {node.Fy:.0f} kN', fontsize=11)
            #plt.text(x, y+0.02, f'Mz {node.Mz:.0f} kNm', fontsize=11)
            
        # ECCENTRICITY ELEMENTS
        for elem in self.elements.values():
            n1, n2 = elem.nodes
            x1, y1 = n1.x
            x2, y2 = n2.x
            plt.plot([x1,x2], [y1, y2], c='b')
            
        plt.axis('equal')
        plt.title("Joint " + str(self.jid))
        if show:
            plt.show()

        
        
if __name__ == "__main__":
    
    H = 5
    L = 5
    n = 6
    
    truss = Truss2D(H, L, n, True)
    truss.generate_chords()
    truss.generate_webs()
    truss.add(FixedSupport([0,0]))
    truss.add(FixedSupport([L,0]))
    truss.generate()
    truss.plot(False)
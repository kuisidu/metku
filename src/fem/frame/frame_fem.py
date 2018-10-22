""" Class for frame analysis by finite element method 

    Written by Kristo Mela

    Attributes:
    
    nodes -- nodal coordinates
    elements -- finite element data
    properties -- element properties
    materials -- material data
    supports -- support conditions
    loads -- load data
    loadcases -- load case data
    dofs -- number of degrees of freedom
    
    Methods:
        
    nels -- number of elements
    nnodes -- number of nodes
    nload -- number of loads
    nloadcases -- number of load cases
    add_node -- insert new node
    add_material -- new material
    add_section -- new element profile
    add_element -- new element
    add_load -- new load
    add_loadcase -- new load case
    add_support -- new support
    nodal_dofs -- assign degrees of freedom to nodes
    element_node_indicies -- return array of node indices of an element
    global_stiffness_matrix -- assemble stiffness matrix
    global_load_vector -- assemble load vector
    linear_statics -- perform structural analysis for a single load case
    

    Uses classes:
    Material
    Section (BeamSection)
    Support
    LoadCase
    FEMNode
    Element 
        Rod
        EBBeam        
    Load
        PointLoad
        LineLoad

"""
    
import sys
import math
import numpy as np
#import scipy as sp
import scipy.sparse.linalg as sp
import matplotlib.pyplot as plt

class FrameFEM:

    def __init__(self):

        self.nodes = []
        self.nodal_coords = []
        self.elements = []
        self.sections = []
        self.materials = []
        self.supports = []
        self.loads = []
        self.loadcases = []
    
        # number of degrees of freedom
        self.dofs = 0

    def nels(self):
        """ Number of elements """
        return len(self.elements)
    
    def nnodes(self):
        """ Number of nodes """
        return len(self.nodes)

    def nsupp(self):
        """ Number of supported nodes """
        return len(self.supports)

    def nload(self):
        """ Number of loads """
        return len(self.loads)

    def nloadcases(self):
        """ Number of load cases """
        return len(self.loadcases)

    def add_node(self,x,y):
        newNode = FEMNode(x,y)
        self.nodal_coords.append([x, y])
        
        for i in range(self.nloadcases()):
            newNode.u = np.vstack((newNode.u,[0.0,0.0,0.0]))
            
        self.nodes.append(newNode)

    def add_material(self,E,nu,density):
        """ Adds new material """
        newMat = Material(E,nu,density)
        self.materials.append(newMat)

    def add_section(self,newSection):
        """ Add new cross section """
        self.sections.append(newSection)

    def add_element(self,new_elem):
        """ Add new element """
        self.elements.append(new_elem)

    def add_support(self,sid,nid,dof,val):
        """ Add support condition """
        new_support = Support(sid,self.nodes[nid],dof,val)
        self.supports.append(new_support)

    def add_load(self,load):
        self.loads.append(load)

    def add_loadcase(self,supp_id=1,load_id=2):
        """ Add load case 
            input:
                supp_id .. integer stating the support condition label
                load_id .. integer stating the load label
        """
        new_case = LoadCase(supp_id,load_id)
        self.loadcases.append(new_case)
        
        """ Add displacements to nodes for the new load case """
        for node in self.nodes:            
            if node.u.size == 0:
                node.u = [0.0,0.0,0.0]
            else:
                node.u = np.vstack((node.u,[0.0,0.0,0.0]))

    def nodal_dofs(self):
        """ Assign nodal degrees of freedom
            
            Numbering of dofs
        """

        """ HERE THE ELEMENTS CAN BE USED TO SET SOME OF THE 
            DOFS TO ZERO. FOR EXAMPLE, TRUSS ELEMENTS DO NOT HAVE
            ROTATIONAL DEGREES OF FREEDOM, AND IF SOME NODES
            ARE CONNECTED TO ONLY TRUSS ELEMENTS, THOSE NODES WILL NOT
            HAVE A ROTATION.
        """        
        for elem in self.elements:
            elem.init_dofs()
        
        """ First set the flag of supported dofs to 1 
            These will be zeroed later
        """
        for s in self.supports:
            for d in s.dof:
                s.node.dofs[d] = 1
                #self.nodes[s.nid].dofs[d] = 1
        
        """ Make running numbering of dofs.
            If the dof is marked as 1, then its dof is set to 0.
        """
        for n in self.nodes:
            for i in range(len(n.dofs)):                
                if n.dofs[i] > 0:
                    n.dofs[i] = -1
                else:
                    n.dofs[i] = self.dofs
                    self.dofs += 1
        
    def element_node_indices(self,e):
        """ Get node indices of an element 
            input: e .. element index
            output: ndx .. list of element node indices
        """
        ndx = []
        for n in self.elements[e].nodes:            
            ndx.append(self.nodes.index(n))
            
        return ndx
    
    def global_stiffness_matrix(self):
        """ Construct the global stiffness matrix """
        dofs = self.dofs
        
        """ initialize zero stiffness matrix 
            NOTE: for larger structures, K should probably be
            a sparse matrix
        """
        K = np.zeros((dofs,dofs))        
        
        for elem in self.elements:
            ke = elem.stiffness_matrix()
            # get element degrees of freedom
            # change the list to numpy array
            ve = np.array(elem.global_dofs())
            # find non-zero dofs
            nz = ve >= 0
            q = ve[nz]
            
            """ extract the submatrix of the element
                and add it to the global stiffness matrix
            """            
            K[np.ix_(q,q)] += ke[np.ix_(nz,nz)]
            
        return K
    
    def global_geometric_stiffness_matrix(self):
        """ Construct the global geometric stiffness matrix """
        dofs = self.dofs
        
        """ initialize zero geometric stiffness matrix 
            NOTE: for larger structures, KG should probably be
            a sparse matrix
        """
        KG = np.zeros((dofs,dofs))
        
        for elem in self.elements:
            kg = elem.geometric_stiffness_matrix()
            # get element degrees of freedom
            # change the list to numpy array
            ve = np.array(elem.global_dofs())
            # find non-zero dofs
            nz = ve >= 0
            q = ve[nz]
            
            """ extract the submatrix of the element
                and add it to the global stiffness matrix
            """            
            KG[np.ix_(q,q)] += kg[np.ix_(nz,nz)]
            
        return KG
            
    def global_load_vector(self,sid):
        """ Construct global load vector """
                        
        """ Initialize load vector """
        global_load = np.zeros(self.dofs)
        
        for load in self.loads:
            """ Get the loads to be added to the vector 
                and the global degrees of freedom of these loads
            """
            if load.sid == sid:
                v, vdofs = load.load_and_dofs()                
                global_load[vdofs] += v
                
        #print("FrameFEM  ", '\n', global_load)
        return global_load

    def linear_statics(self,lcase=0):
        """ Perform linear elastic static analysis for a given load case
        
        """
                
        K = self.global_stiffness_matrix()
        
        load_id = self.loadcases[lcase].load
        p = self.global_load_vector(load_id)

        """ Solve global displacement vector """
        u = np.linalg.solve(K,p)
        
                            
        """ Distribute displacements to nodes """
        for node in self.nodes:
            # Get nodal dofs
            dofs = np.array(node.dofs)
            
            # Find free dofs
            free_dofs = dofs>=0 
                        
            # Substitute free dofs from global displacement vector                        
            for i in range(len(free_dofs)):
                if free_dofs[i]:
                    node.u[i] = u[dofs[i]]
            
            """
            This substition did not work !!!
            node.u[free_dofs] = u[dofs[free_dofs]]
            """
        
        """ Calculate element internal forces """
        for el in self.elements:
            el.internal_forces()
        
        return u,K
    
    def linear_buckling(self,lcase=0, k=4):
        """ Perform linear buckling analysis for a given load case
        """
        
        """ Solve internal forces of members 
            and return stiffness matrix
        """
        u,K = self.linear_statics(lcase)
        KG = self.global_geometric_stiffness_matrix()        
        w,v = sp.eigsh(K,k=k,M=-KG,sigma=0.0,which='LA',mode='normal')
        #w,v = sp.eigvals(K,b=-KG)
        
        """ Distribute buckling displacements to nodes """
        for node in self.nodes:
            # Get nodal dofs
            dofs = np.array(node.dofs)
            
            # Find free dofs
            free_dofs = dofs>=0 
                        
            # Substitute free dofs from global displacement vector                        
            for i in range(len(free_dofs)):
                if free_dofs[i]:
                    node.v[i] = v[dofs[i]]
        

        return w,v
    
    
    def draw(self):
        """ draw nodes """
        for n in self.nodes:
            plt.plot(n.x[0],n.x[1],'ro')
        
        for i in range(self.nnodes()):
            plt.text(self.nodes[i].x[0],self.nodes[i].x[1],str(i))
        
        """ draw members """
        for i in range(self.nels()):
            #X = self.member_coord(i)
            X = self.elements[i].coord()
            plt.plot(X[:,0],X[:,1],'b')
            Xmid = X[0,:]+0.5*(X[1,:]-X[0,:])
            plt.text(Xmid[0],Xmid[1],str(i))

    def print_displacements(self):
        node_cnt = 0
        for node in self.nodes:
            print('*** Node {:d} ***'.format(node_cnt))
            node_cnt += 1
            print('u = {0:5.4f} mm'.format(node.u[0]))
            print('v = {0:5.4f} mm'.format(node.u[1]))
            print('phi = {0:5.4f} rad'.format(node.u[2]))

class Material:
    """ Class for material data
       
        For now, only linear elastic material is included.
        
        Attributes:
            young .. Young's modulus
            nu .. Poisson ratio
            density .. material density
    """

    def __init__(self,E,nu,rho):
        self.young = E
        self.nu = nu
        self.density = rho

class Section:
    """ Class for cross section properties
    
        Attributes:
            area .. cross-sectional area
            
        All sections must have 'area', but specific profiles
        can have other properties as well
    """

    def __init__(self,A):
        self.area = A

class BeamSection(Section):
    """ Class for beam element cross sections 
    
        Attributes:
            Iy .. second moment of area with respect to major axis
    """

    def __init__(self,A,Iy):
        Section.__init__(self,A)
        self.Iy = Iy


class Support:
    """ Supported nodes / degrees of freedom """

    def __init__(self,sid,node,dof,val=0.0):
        """ Input:
            sid -- label for supports
            node -- supported node
            dof -- supported degrees of freedom. Array with at most three
                   elements having values from 0 to 2 (see FEMNode for
                   ordering of the degrees of freedom)
            val -- value for supports (default = 0.0, but theoretically
                    the displacement could have a non-zero value
        """
        self.sid = sid
        self.node = node
        self.dof = dof
        self.val = val

class Load:
    """ General class for loads """
    
    def __init__(self,sid):
        """ Constructor """
        self.sid = sid

class PointLoad(Load):
    """ Class for point loads (forces and moments) """
    
    def __init__(self,sid,node,v,f):
        """ Input:
            sid -- load id
            nid -- node subjected to the load
            v -- load vector: v[0] = Fx, v[1] = Fy, v[2] = Mz
            f -- scaling factor
        """
        Load.__init__(self,sid)

        self.node = node
        self.v = v
        self.f = f

    def load_and_dofs(self):
        """ Returns the loads to be inserted to the global load vector
            and the corresponding global degrees of freedom
        """
        
        """ Load vector """
        F = self.f*np.array(self.v)
        
        """ Degrees of freedom """
        dofs = np.array(self.node.dofs)
        
        """ Find free degrees of freedom (those with numbering
                                          greater than -1) 
        """
        nzero_dofs = dofs >= 0
        
        """ Get free degrees of freedom and the corresponding part
            of the load vector
        """
        d = dofs[nzero_dofs]
        Fnz = F[nzero_dofs]
        
        return Fnz, d
                
class LineLoad(Load):
    """ Class for 2D line loads (forces and moments) """
    
    def __init__(self,sid,eid,xloc,qval,direction,coords=1):
        """ Input:
            sid -- load id
            eid -- element subjected to the load
            xloc -- starting and ending locations of the load (local coordinates)
                    xloc[0] .. staring coordinate 0 <= xloc[0] <= 1
                    xloc[1] .. ending coordinate 0 <= xloc[0] < xloc[1] <= 1
            qval -- value of the load at the coordinates xloc
            direction -- 0 .. x-axis
                         1 .. y-axis
            coords .. 0 .. local
                      1 .. global
        """
        Load.__init__(self,sid)

        self.elem = eid
        self.xloc = xloc
        self.qval = qval
        self.dir = direction
        self.coords = coords
        
    def load_and_dofs(self):
        """ Returns the loads to be inserted to the global load vector
            and the corresponding global degrees of freedom
        """ 
        
        """ Get load vector """        
        if self.dir == "x":
            q = np.array([self.qval[0],0,0])
        else:
            q = np.array([0,self.qval[0],0])
        
        """ Get equivalent nodal loads """        
        F = self.elem.equivalent_nodal_loads(q)
                
        
        """ Number of degrees of freedom """
        dofs = np.array(self.elem.global_dofs())
                        
        
        """ Find free degrees of freedom (those with numbering
                                          greater than -1) 
        """
        nzero_dofs = dofs >= 0
        
        """ Get free degrees of freedom and the corresponding part
            of the load vector
        """
        d = dofs[nzero_dofs]
        Fnz = F[nzero_dofs]
                                
        return Fnz, d
        
class LoadCase:
    """ Class for storing load case data """

    def __init__(self,support,load):
        self.support = support
        self.load = load

class FEMNode:
    """ Class for FEM nodes.
        
        Basically only for storing data
    """
    
    def __init__(self,x,y):
        """ x and y are the coordinates of the node 
            
            dofs -- numbers of nodal degrees of freedom
            u -- nodal displacements in the global coordinate system
            v -- nodal displacements for linear buckling 
        """
        self.x = np.array([x,y])
        """ Setting all dofs to 0 means that the dofs are free.
            The method nodal_dofs of class FrameFEM assigns the actual
            numbers to the dofs.
            
            dofs[0] -- ux
            dofs[1] -- uy
            dofs[2] -- rz
            
            For 3D frames, the order of dofs must be redefined
        """
        self.dofs = [0,0,0]
        # maybe u could also be a simple list?
        #self.u = np.array([0.0,0.0,0.0])
        """ Initialize displacements as an empty array
            NOTE: for multiple load cases, the dimension of u must be increased
            for each load case
        """
        self.u = np.array([])
        self.v = [0,0,0]

class Element:
    """ 1D finite elements: bars and beams
        Attributes:
            
        nodes -- array of FEMNode variables, member end nodes
        material -- material properties (Material)
        section -- cross-section properties (Section)
        axial_force -- axial force, list
        
    """

    def __init__(self,n1,n2,section,material):
        """ Input:
            n1 -- node 1 (reference to FEMNode variable)
            n2 -- node 2 (reference to FEMNode variable)
            section -- cross-section properties
            material -- Material type variable
                        
        """
        self.nodes = [n1,n2]
        self.material = material
        self.section = section
        self.axial_force = [0.0,0.0]
        self.floc = np.array([])

    def coord(self):
        """ Nodal coordinates of the element """
        X = np.array([self.nodes[0].x,self.nodes[1].x])
        return X

    def length(self):
        """ Member length """
        X = self.coord()
        L = np.linalg.norm(X[1,:]-X[0,:])
        return L
    
    def direction_cosines(self):
        """ Direction cosines """
        X = self.coord()
        dX = X[1,:]-X[0,:]
        L = self.length()
        c = dX/L
        return c

    def stiffness_matrix(self):
        return None
    
    def equivalent_nodal_loads(self,q):
        """ Compute equivalent nodal loads with line loads in
            vector q
        """
        return None
    
    def init_dofs(self):
        """ Initialize dofs 
            Does nothing for beam elements
        """  

    def nodal_displacements(self):
        """ Get nodal displacements of an element 
            Requires previously performed structural analysis such
            that nodal displacements are available.
        """
        
        return np.hstack((self.nodes[0].u,self.nodes[1].u))
    
    def local_displacements(self):
        """ Nodal displacements in local coordinates """
        T = self.transformation_matrix()
        q = self.nodal_displacements()
        return T.dot(q)



# For testing
if __name__ == "__main__":
    import sys
    
    from hollow_sections import SHS
    
    p = SHS(100,4)
    
    f = FrameFEM()
    
    # Create nodes   
    L = 5000.0
    h = 5000.0
    
    nbeam_els = 2
    dx = L/nbeam_els
    
    f.add_node(0.0,0.0)
    f.add_node(0.0,h)
    
    X = 0.0
    for i in range(0,nbeam_els):
        f.add_node(X+(i+1)*dx,h)
    
    #f.add_node(0.5*L,h)
    #f.add_node(L,h)
    f.add_node(L,0.0)
    
    #f.add_node(0.0,6000.0)
    #f.add_node(6000.0,6000.0)
    
    # Add material
    f.add_material(210e3,0.3,7850e-9)
    
    # Add section properties
    s = BeamSection(p.A,p.I[0])
    #s = BeamSection(100000,60000e4)
    f.add_section(s)
    
    # Add first element
    
    """
    for i in range(0,f.nnodes()-1):
        el = EBBeam(f.nodes[i],f.nodes[i+1],s,f.materials[0])  
        f.add_element(el)
    

    """    
    # column
    el = EBBeam(f.nodes[0],f.nodes[1],s,f.materials[0])
    f.add_element(el)
    # first beam element
    Sj1 = 0.0
    Sj2 = np.inf
    #rstiff = 30000
    el = EBSemiRigidBeam(f.nodes[1],f.nodes[2],s,f.materials[0],rot_stiff=[Sj1,Sj2])
    #el = EBBeam(f.nodes[1],f.nodes[2],s,f.materials[0])
    f.add_element(el)
    for i in range(0,nbeam_els-2):
        #el = EBSemiRigidBeam(f.nodes[i+2],f.nodes[i+3],s,f.materials[0],rot_stiff=[np.inf,np.inf])
        el = EBBeam(f.nodes[i+2],f.nodes[i+3],s,f.materials[0])
        f.add_element(el)
    
    
    el = EBSemiRigidBeam(f.nodes[-3],f.nodes[-2],s,f.materials[0],rot_stiff=[Sj2,Sj1])
    #el = EBBeam(f.nodes[-3],f.nodes[-2],s,f.materials[0])
    f.add_element(el)
    el = EBBeam(f.nodes[-2],f.nodes[-1],s,f.materials[0])
    f.add_element(el)
    
    #el = EBSemiRigidBeam(f.nodes[0],f.nodes[1],s,f.materials[0],[0.0e15,0.0e15])    
    #f.add_element(el)
    #el2 = EBBeam(f.nodes[1],f.nodes[2],s,f.materials[0])
    #f.add_element(el2)
    
    #el2 = EBSemiRigidBeam(f.nodes[1],f.nodes[2],s,f.materials[0],[1e4*1e3,kRigid])
    #f.add_element(el2)
    
    # Add supports
    f.add_support(1,0,[0,1,2],0.0)
    print(f.nnodes())
    f.add_support(1,f.nnodes()-1,[0,1,2],0.0)
    #f.add_support(1,0,[0,1,2],0.0)
    #f.add_support(1,2,[0,1,2],0.0)
    
    # Assign nodal dofs
    f.nodal_dofs()
    
    # Add load    
    """
    load1 = LineLoad(2,f.elements[1],[0,1],-5*np.array([1,1]),1)
    load2 = LineLoad(2,f.elements[2],[0,1],-5*np.array([1,1]),1)
    f.add_load(load1)
    f.add_load(load2)
    """
    nload = f.nodes[-(2+int(nbeam_els/2))]
    print(nload.x)
    load1 = PointLoad(2,f.nodes[-(2+int(nbeam_els/2))],[0.0,-25000.0,0.0],1.0)
    f.add_load(load1)
    #f.add_load(load2)
    
    f.draw()
    # Add point load
    """
    load1 = PointLoad(2,f.nodes[0],[0.0,0.0,100.0],1.0)
    f.add_load(load1)
    load2 = PointLoad(2,f.nodes[1],[0.0,0.0,-100.0],1.0)
    f.add_load(load2)
    """
    # Add load case
    f.add_loadcase(supp_id=1,load_id=2)
    
    p = f.global_load_vector(2)
    
    #print(p)
    """ Two-storey frame """
    
    
    
    """
    # Create nodes   
    f.add_node(0,0)
    f.add_node(1000.0,500.0)
    f.add_node(0,500.0)
    
    # Add material
    f.add_material(210e3,0.3,7850e-9)
    
    # Add section properties
    s = Section(p.A)
    f.add_section(s)
    
    # Add first element
    el = Rod(f.nodes[0],f.nodes[1],f.sections[0],f.materials[0])
    f.add_element(el)

    # Add second element
    el2 = Rod(f.nodes[2],f.nodes[1],f.sections[0],f.materials[0])
    f.add_element(el2)
    
    # Add supports
    f.add_support(1,0,[0,1],0)
    f.add_support(1,2,[0,1],0)
    
    # Assign nodal dofs
    f.nodal_dofs()
    
    # Add point load
    load1 = PointLoad(2,f.nodes[1],[-100.0,0,0],1.0)
    f.add_load(load1)
    
    # Add load case
    f.add_loadcase(supp_id=1,load_id=2)
    
    p = f.global_load_vector(2)
    """
    
    """
    k = f.elements[0].stiffness_matrix()
    print(k)

    sBeam = BeamSection(p.A,p.I[0])

    bel = EBBeam(f.nodes[0],f.nodes[1],sBeam,f.materials[0])

    f.add_element(bel)

    kbeam = bel.stiffness_matrix()
    print(kbeam)

    f.add_support(1,0,[0,1],0)

    f.nodal_dofs()

    print(f.element_node_indices(0))
    print(f.element_node_indices(1))
    
    print(f.nodes[0].dofs)
    print(f.nodes[1].dofs)
    
    for n in f.nodes:
        print(n.dofs)
    
    for e in f.elements:
        print(e.global_dofs())
    """
    
    K = f.global_stiffness_matrix()
    #print(K)
    #print(K.shape)
    
    u = f.linear_statics()
    
    
    
    K = f.global_stiffness_matrix()
    #print(K)
    #print(np.shape(K))
    
    #print(u)

    print("Results of analysis:")

    print("Nodal displacements:")
    node_cnt = 0
    for node in f.nodes:
        print('*** Node {:d} ***'.format(node_cnt))
        node_cnt += 1
        print('u = {0:5.4f} mm'.format(node.u[0]))
        print('v = {0:5.4f} mm'.format(node.u[1]))
        print('phi = {0:5.4f} rad'.format(node.u[2]))

    print("Member forces:")
    mem_cnt = 0
    for el in f.elements:
        print('*** Element {:d} ***'.format(mem_cnt))
        mem_cnt += 1
        print(el.floc)        
        print('Axial force (node 1) = {0:5.3f} kN'.format(el.axial_force[0]*1e-3))
        print('Axial force (node 2) = {0:5.3f} kN'.format(el.axial_force[1]*1e-3))
        print('Bending moment (node 1) = {0:5.3f} kNm'.format(el.bending_moment[0]*1e-6))
        print('Bending moment (node 2) = {0:5.3f} kNm'.format(el.bending_moment[1]*1e-6))
        print('Shear force (node 1) = {0:5.3f} kN'.format(el.shear_force[0]*1e-3))
        print('Shear force (node 2) = {0:5.3f} kN'.format(el.shear_force[1]*1e-3))
        #print(el.shear_force)
    
        
    #f.draw()
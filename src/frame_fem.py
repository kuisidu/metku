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
import matplotlib.pyplot as plt

# Constants
kRigid = 1e20 # articifial large stiffness for rigid joints
kHinge = 1e-3 # artificial small stiffness for hinges

class FrameFEM:

    def __init__(self):

        self.nodes = []
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
            # if node.u.size == 0:
            if len(node.u) == 0:
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
        return u
    
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
        #X = self.coord()
        L = math.sqrt(sum((self.nodes[0].x-self.nodes[1].x)**2))
        #L = np.linalg.norm(X[1,:]-X[0,:])
        return L
    
    def direction_cosines(self):
        """ Direction cosines """
        #X = self.coord()
        #dX = X[1,:]-X[0,:]
        dX = self.nodes[1].x-self.nodes[0].x
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

class Rod(Element):
    """ Rod element, carries only axial loads """

    def __init__(self,n1,n2,section,material):
        Element.__init__(self,n1,n2,section,material)


    def transformation_matrix(self):
        """ From local to global coordinates """
        c = self.direction_cosines()
        L1 = np.hstack( (c,np.zeros(2)) )
        L2 = np.hstack( (np.zeros(2),c) )
        L = np.vstack( (L1,L2))
        return L
    
    def global_dofs(self):
        """ Get numbering of element global degrees of freedom """        
        return self.nodes[0].dofs[:2] + self.nodes[1].dofs[:2]

    def init_dofs(self):
        """ For truss members there are no rotational degrees of freedom
            so these can be set to be neglected (=1)
        """
        for n in self.nodes:
            n.dofs[2] = 1
        
    def local_stiffness_matrix(self):
        """ Stiffness matrix in local coordinates """
        E = self.material.young
        A = self.section.area
        Le = self.length()
    
        """ Stiffness matrix in local coordinates """
        return E*A/Le*np.array([[1,-1],[-1,1]])        
        
    def stiffness_matrix(self):
        """ Compute the stiffness matrix """
        k0 = self.local_stiffness_matrix()
    
        """ Transformation matrix """
        L = self.transformation_matrix()
        
        """ Transform stiffness matrix to global coordinates """
        ke = L.transpose().dot(k0.dot(L))

        return ke

    def equivalent_nodal_loads(self,q):
        """ Equivalent nodal loads for load in vector q 
            Ignore bending moment
        """
        
        fglob = np.zeros(4)
        
        dc = self.direction_cosines()
        L = self.length()
        
        fglob[[0,2]] = 0.5*dc[1]*q[0]*L
        fglob[[1,3]] = 0.5*dc[0]*q[1]*L
             
        return fglob

    def internal_forces(self):
        """ Calculates internal forces (axial force) """
        q = self.local_displacements()
        E = self.material.young
        A = self.section.area
        L = self.length()
        self.axial_force[0] = E*A/L*(q[1]-q[0])

class EBBeam(Element):
    """ Euler-Bernoulli beam element 
        
        Attributes:
            bending_moment -- bending moment in the end nodes [list]
            shear_force -- shear force in the end nodes [list]
    """
    
    def __init__(self,n1,n2,section,material):
        Element.__init__(self,n1,n2,section,material)
        
        self.bending_moment = [0.0,0.0]
        self.shear_force = [0.0,0.0]        
    
    def transformation_matrix(self):
        """ From local to global coordinates """
        c = self.direction_cosines()
        """
        L = np.zeros((6,6))
        L[0,0:2] = c[0:2]
        L[1,0] = -L[0,1]
        L[1,1] = L[0,0]
        L[2,2] = 1.0
        L[3:6,3:6] = L[0:3,0:3]
        """
        L = np.array([[c[0],c[1],0,0,0,0],
                    [-c[1],c[0],0,0,0,0],
                    [0,0,1.0,0,0,0],
                    [0,0,0,c[0],c[1],0],
                    [0,0,0,-c[1],c[0],0],
                    [0,0,0,0,0,1.0]])

        return L
    
    def global_dofs(self):
        """ Get numbering of element globall degrees of freedom """        
        return  self.nodes[0].dofs + self.nodes[1].dofs
    
    def local_stiffness_matrix(self):
        """ Stiffness matrix in local coordinates """
        E = self.material.young
        A = self.section.area
        I1 = self.section.Iy
        Le = self.length()
        
        rodc = E*A/Le
        EI1 = E*I1
        """
        k0 = np.zeros((6,6))
        #
        k0[[0,3],[0,3]] = rodc
        k0[[0,3],[3,0]] = -rodc
        #
        bend1 = 12*EI1/Le**3
        k0[[1,4],[1,4]] = bend1
        k0[[1,4],[4,1]] = -k0[1,1]
        #
        bend2 = 6*EI1/Le**2
        k0[[1,2,1,5],[2,1,5,1]] = bend2
        #
        k0[[2,4,4,5],[4,2,5,4]] = -bend2
        #
        k0[[2,5],[5,2]] = 2*EI1/Le
        #
        k0[[2,5],[2,5]] = 4*EI1/Le
        
        """
        # This implementation is 20% faster to run
        k0 = np.array([[E*A/Le,0,0,-E*A/Le,0,0],
              [0,12*EI1/Le**3,6*EI1/Le**2 ,0,-12*EI1/Le**3,6*EI1/Le**2],             
              [0,6*EI1/Le**2, 4*EI1/Le,0,-6*EI1/Le**2,2*EI1/Le],             
              [-E*A/Le, 0,0,E*A/Le, 0,0],              
              [0,-12*EI1/Le**3,-6*EI1/Le**2,0,12*EI1/Le**3,-6*EI1/Le**2],             
              [0,6*EI1/Le**2,2*EI1/Le,0,-6*EI1/Le**2,4*EI1/Le]])
        
        return k0
    
    def stiffness_matrix(self):
        """ Compute the stiffness matrix """
        
        k0 = self.local_stiffness_matrix()                
        
        #k0 = CheckReleases[fem.elem[N],k0];
                        
        # local-global matrix
        L = self.transformation_matrix()
        
        # globaali elementin matriisi
        ke = L.transpose().dot(k0.dot(L))
        
        return ke

    def equivalent_nodal_loads(self,q):
        """ Equivalent nodal loads for load in vector q
        """
        
        floc = np.zeros(6)
        
        T = self.transformation_matrix()
        L = self.length()
        
        # Load vector transformed into element local coordinate system
        qloc = T[:3,:3].dot(q)
        
        
        # Construct nodal load in local coordinates
        # qloc[0,1,2] = axial, shear, moment of node 1
        # qloc[3,4,5] = axial, shear, moment ofnode 2
        # Axial force
        floc[[0,3]] = 0.5*qloc[0]*L
        # Shear force
        floc[[1,4]] = 0.5*qloc[1]*L
        # Moment
        floc[2] = qloc[1]*L**2/12.0
        floc[5] = -floc[2]
        
        # Store local loads
        self.floc = floc
        
        #print(self.floc)
        
        fglob = T.transpose().dot(floc)
        return fglob

    def shape_fun(self,x):
        """ Evaluate Hermitean shape functions at local coordinate x """
        L = self.length()
        N = np.zeros(4)
        N[0] = 1-3*x**2/L**2 + 2*x**3/L**3
        N[1] = x-2*x**2/L+x**3/L**2
        N[2] = 3*x**2/L**2-2*x**3/L**3
        N[3] = -x**2/L**2 + x**3/L**2
         
        return N

    def internal_forces(self):
        """ Calculate internal forces 
            NOTE: these internal forces do not take
            loads along the element into account!
            
            Works only for a single load case!
        """
        
        """ Get nodal displacements in local coordinates
            and multiply them with local element stiffness matrix
            to get internal forces in member's local coordinate system.
        """
        q = self.local_displacements()        
        ke = self.local_stiffness_matrix()
        
        #print(self.floc)
        #R = ke.dot(q) - self.floc
        
        if len(self.floc) > 0:
            R = ke.dot(q)-self.floc
        else:
            R = ke.dot(q)
                  
        #print(R)
        
        
        """ Any load on the element not acting on a node must be
            taken into account here. This requires the following steps:
            1. Identify loads acting on element
            2. Compute their equivalent nodal loads
            3. Modify internal forces
            
            Probably an easy solution is to include an
            Attribute floc for equivalent nodal loads. The nodal loads are
            saved there, when they are first calculated for analysis
        """        
        
        self.axial_force[:2] = [-R[0],R[3]]
        self.bending_moment[:2] = [-R[2],R[5]]
        self.shear_force[:2] = [R[1],-R[4]]
        
        

class EBSemiRigidBeam(Element):
    """ Euler-Bernoulli beam element with possibility to include
        semi-rigid joints. Based on Jalkanen (2004):
        JOUSTAVASTI TUETTU TASOPALKKIELEMENTTI, Journal of Structural Mechanics
        
        Attributes:
            bending_moment -- bending moment in the end nodes [list]
            shear_force -- shear force in the end nodes [list]
            rot_stiff -- rotational stiffness of nodes [list]
    """
    
    def __init__(self,n1,n2,section,material,rot_stiff=[np.inf,np.inf]):
        Element.__init__(self,n1,n2,section,material)
        self.bending_moment = [0.0,0.0]
        self.shear_force = [0.0,0.0]
        #self.floc = np.array([])
        for i in range(2):
            if rot_stiff[i] <= kHinge:
                rot_stiff[i] = kHinge
            elif rot_stiff[i] >= kRigid:
                rot_stiff[i] = kRigid
                        
        #print(rot_stiff)
        self.rot_stiff = rot_stiff
        #print(self.rot_stiff)
    
    def transformation_matrix(self):
        """ From local to global coordinates """
        c = self.direction_cosines()
        """
        L = np.zeros((6,6))
        L[0,0:2] = c[0:2]
        L[1,0] = -L[0,1]
        L[1,1] = L[0,0]
        L[2,2] = 1.0
        L[3:6,3:6] = L[0:3,0:3]
        """
        L = np.array([[c[0],c[1],0,0,0,0],
            [-c[1],c[0],0,0,0,0],
            [0,0,1.0,0,0,0],
            [0,0,0,c[0],c[1],0],
            [0,0,0,-c[1],c[0],0],
            [0,0,0,0,0,1.0]])
        return L
    
    def gamma(self):
        """ Eq. (12) of Jalkanen (2004) """
        
        gamma = []
        
        E = self.material.young
        I1 = self.section.Iy
        L = self.length()
        g0 = E*I1/L
        
        """
        for k in self.rot_stiff:
            # k = 0 means a hinged connection
            if k == 0.0:
                gamma.append(np.inf)
            else:                
                gamma.append(g0/k)
        """        
        for k in self.rot_stiff:
            """ k = 0 means a hinged connection """
            if k == 0.0:
                gamma.append(np.inf)
            elif k == kRigid:
                """ set gamma = 0.0 for rigid ends """
                gamma.append(0.0)
            else:
                gamma.append(g0/k)
        
        return gamma
    
    def delta(self):
        """ Eq. (15) of Jalkanen (2004) """
        gamma = self.gamma()
        delta = (1+4*gamma[0])*(1+4*gamma[1])-4*gamma[0]*gamma[1] 
        return delta
    
    def matrix_S(self):
        """ Eq. (16) of Jalkanen (2004) """
        #S = np.zeros((4,4))
        L = self.length()
        gamma = self.gamma()
            
        
        #S[1,0:4] = [-6*gamma[0]/L*(1+2*gamma[1]),-4*gamma[0]*(1+3*gamma[1]),6*gamma[0]/L*(1+2*gamma[1]),-2*gamma[0]]
        #S[3,0:4] = [-6*gamma[1]/L*(1+2*gamma[0]),-2*gamma[1],6*gamma[1]/L*(1+2*gamma[0]),-4*gamma[1]*(1+3*gamma[0])]
        
        
        S = np.array([[0,0,0,0],
                      [-6*gamma[0]/L*(1+2*gamma[1]),-4*gamma[0]*(1+3*gamma[1]),6*gamma[0]/L*(1+2*gamma[1]),-2*gamma[0]],
                      [0,0,0,0],
                      [-6*gamma[1]/L*(1+2*gamma[0]),-2*gamma[1],6*gamma[1]/L*(1+2*gamma[0]),-4*gamma[1]*(1+3*gamma[0])]])
        
        
        return S
    
    def global_dofs(self):
        """ Get numbering of element global degrees of freedom """        
        return  self.nodes[0].dofs + self.nodes[1].dofs
    
    def local_stiffness_matrix(self):        
        """ Stiffness matrix in local coordinates """        
        E = self.material.young
        A = self.section.area
        I1 = self.section.Iy
        Le = self.length()
              
        rodc = E*A/Le
        EI1 = E*I1
        bend1 = 12*EI1/Le**3
        bend2 = 6*EI1/Le**2
        
        #print(bend1,bend2)
        
        #k00 = np.zeros((4,4))
        k0 = np.zeros((6,6))
        """
        k00[[0,2],[0,2]] = bend1
        k00[[0,2],[2,0]] = -bend1
        
        
        k00[[0,1,0,3],[1,0,3,0]] = bend2          
        #
        k00[[1,2,2,3],[2,1,3,2]] = -bend2
        #
        k00[[1,3],[3,1]] = 2*EI1/Le
        #  
        k00[[1,3],[1,3]] = 4*EI1/Le
        """
        k00 = np.array([[12*EI1/Le**3,6*EI1/Le**2,-12*EI1/Le**3,6*EI1/Le**2],
                [6*EI1/Le**2,4*EI1/Le,-6*EI1/Le**2,2*EI1/Le],
                [-12*EI1/Le**3,-6*EI1/Le**2,12*EI1/Le**3,-6*EI1/Le**2],
                [6*EI1/Le**2,2*EI1/Le,-6*EI1/Le**2,4*EI1/Le]])
        
        
        delta = self.delta()
        S = self.matrix_S()
        
        #print(self.gamma())
        #print(delta)
        #print(k00)
        
        
        k1 = 1/delta*k00.dot(S)+1/delta*S.transpose().dot(k00)+1/delta**2*S.transpose().dot(k00).dot(S)
        
                        
        #s2 = S[1,:]        
        #s4 = S[3,:]        
        #k2 = self.rot_stiff[0]/delta**2*np.outer(s2,s2) + self.rot_stiff[1]/delta**2*np.outer(s4,s4)
        
        k2 = np.zeros((4,4))
        for i in range(0,2):
            """ if an end is rigid, then it does not contribute to k2,
                so only non-rigid ends count here
            """
            if self.rot_stiff[i] < kRigid:
                s = S[2*i+1,:]
                #print(s)
                k2 += self.rot_stiff[i]/delta**2*np.outer(s,s)
        
        
        #print("k00",k00[0][0])
        #print("k1", k1.shape)
        #print("k2", k2.shape)
        
        k = k00 + k1 + k2
        """
        k0[np.ix_([1,2,4,5],[1,2,4,5])] = k00+k1+k2        
        k0[[0,3],[0,3]] = rodc
        k0[[0,3],[3,0]] = -rodc
        
        """
        k0 = np.array([[E*A/Le,0,0,-E*A/Le,0,0],
              [0,k[0][0],k[0][1],0,k[0][2],k[0][3]],           
              [0,k[1][0],k[1][1],0,k[1][2],k[1][3]],          
              [-E*A/Le, 0,0,E*A/Le, 0,0],              
              [0,k[2][0],k[2][1],0,k[2][2],k[2][3]],      
              [0,k[3][0],k[3][1],0,k[3][2],k[3][3]]])
        
        #print(k0)
          
        return k0
    
    def stiffness_matrix(self):
        """ Compute the stiffness matrix """
        
        k0 = self.local_stiffness_matrix()                
        
        #k0 = CheckReleases[fem.elem[N],k0];
                        
        # local-global matrix
        L = self.transformation_matrix()
        
        # globaali elementin matriisi
        ke = L.transpose().dot(k0.dot(L))
        
        return ke

    def equivalent_nodal_loads(self,q):
        """ Equivalent nodal loads for load in vector q
        """
        
        #floc = np.zeros(6)
        #floc0 = np.zeros(6)
        
        T = self.transformation_matrix()
        L = self.length()
        
        # Load vector transformed into element local coordinate system
        qloc = T[:3,:3].dot(q)
        
        
        delta = self.delta()
        S = self.matrix_S()
        
        # Construct nodal load in local coordinates
        """
        floc[[0,3]] = 0.5*qloc[0]*L
        floc[[1,4]] = 0.5*qloc[1]*L
        floc[2] = qloc[1]*L**2/12.0
        floc[5] = -floc[2]
        """
        """
        # Construct nodal load in local coordinates
        floc0[[0,3]] = 0.5*qloc[0]*L
        floc0[[1,4]] = 0.5*qloc[1]*L
        floc0[2] = qloc[1]*L**2/12.0
        floc0[5] = -floc0[2]
        """
        
        floc0 = np.array([0.5*qloc[0]*L, 0.5*qloc[1]*L, qloc[1]*L**2/12.0,
                          0.5*qloc[0]*L, 0.5*qloc[1]*L, -qloc[1]*L**2/12.0])
        
        """
        floc[[0,3]] = 0.5*qloc[0]*L # axial nodal loads
        floc[1] = floc0[1] + S[1,0]/delta*floc0[2] + S[3,0]/delta*floc0[5]
        floc[2] = (1+S[1,1]/delta)*floc0[2] + S[3,1]/delta*floc0[5]
        floc[4] = S[1,2]/delta*floc0[2] + floc0[4] + S[3,2]/delta*floc0[5]
        floc[5] = S[1,3]/delta*floc0[2] + (1+S[3,3]/delta)*floc0[5]
        """
        
        floc = np.array([0.5*qloc[0]*L,
                    floc0[1] + S[1,0]/delta*floc0[2] + S[3,0]/delta*floc0[5],
                    (1+S[1,1]/delta)*floc0[2] + S[3,1]/delta*floc0[5],
                    0.5*qloc[0]*L,
                    S[1,2]/delta*floc0[2] + floc0[4] + S[3,2]/delta*floc0[5],
                    S[1,3]/delta*floc0[2] + (1+S[3,3]/delta)*floc0[5]])
        
        
        #print("floc", floc)
        self.floc = floc
        #print(q)
        #print("equivalent", self.floc)
        
        fglob = T.transpose().dot(floc)
        return fglob

    def internal_forces(self):
        """ Calculate internal forces 
            NOTE: these internal forces do not take
            loads along the element into account!
            
            Works only for a single load case!
        """
        q = self.local_displacements()        
        ke = self.local_stiffness_matrix()
        #print(ke)   
        
        #print("internal_forces", self.floc)
        
        if len(self.floc) > 0:
            R = ke.dot(q)-self.floc
        else:
            R = ke.dot(q)
        
        self.axial_force[:2] = [-R[0],R[3]]
        self.bending_moment[:2] = [-R[2],R[5]]
        self.shear_force[:2] = [R[1],-R[4]]

"""
# For testing
if __name__ == "__main__":
    import sys
    
    from hollow_sections import SHS
    
    p = SHS(100,4)
    
    f = FrameFEM()
    
    # Create nodes   
    f.add_node(0.0,0.0)
    f.add_node(3000.0,0.0)
    f.add_node(6000.0,0.0)
    
    #f.add_node(0.0,6000.0)
    #f.add_node(6000.0,6000.0)
    
    # Add material
    f.add_material(210e3,0.3,7850e-9)
    
    # Add section properties
    s = BeamSection(p.A,p.I[0])
    #s = BeamSection(100000,60000e4)
    f.add_section(s)
    
    # Add first element
    el = EBBeam(f.nodes[0],f.nodes[1],s,f.materials[0])  
    #el = EBSemiRigidBeam(f.nodes[0],f.nodes[1],s,f.materials[0],[0.0e15,0.0e15])    
    f.add_element(el)
    el2 = EBBeam(f.nodes[1],f.nodes[2],s,f.materials[0])  
    f.add_element(el2)
    
    #el2 = EBSemiRigidBeam(f.nodes[1],f.nodes[2],s,f.materials[0],[1e4*1e3,kRigid])
    #f.add_element(el2)
    
    # Add supports
    f.add_support(1,0,[0,1,2],0.0)
    #f.add_support(1,2,[1],0.0)
    #f.add_support(1,0,[0,1,2],0.0)
    #f.add_support(1,2,[0,1,2],0.0)
    
    # Assign nodal dofs
    f.nodal_dofs()
    
    # Add load    
    load1 = LineLoad(2,el,[0,1],-20*np.array([1,1]),1)
    load2 = LineLoad(2,el2,[0,1],-20*np.array([1,1]),1)
    f.add_load(load1)
    f.add_load(load2)
    
    f.draw()
    # Add point load

    load1 = PointLoad(2,f.nodes[0],[0.0,0.0,100.0],1.0)
    f.add_load(load1)
    load2 = PointLoad(2,f.nodes[1],[0.0,0.0,-100.0],1.0)
    f.add_load(load2)

    # Add load case
    f.add_loadcase(supp_id=1,load_id=2)
    
    p = f.global_load_vector(2)
    
    #print(p)
"""
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

    K = f.global_stiffness_matrix()
    #print(K)
    
    u = f.linear_statics()
    
    
    
    K = f.global_stiffness_matrix()
    #print(K)
    #print(np.shape(K))
    
    #print(u)

    print("Results of analysis:")

    for node in f.nodes:
        print(node.u)

    for el in f.elements:
        print(el.floc)
        print(el.axial_force)
        print(el.bending_moment)
        print(el.shear_force)
    
        
    #f.draw()

"""

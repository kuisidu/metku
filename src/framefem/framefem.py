
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


from abc import ABCMeta, abstractclassmethod
from scipy.sparse.linalg import eigsh


class FrameFEM:
    """ Class for frame analysis by finite element method

        Written by Kristo Mela

        Variables:
        -------------------
        :ivar nodes: FEMNodes
        :ivar nodal_coords: nodal coordinates [x, y]
        :ivar elements: Element -instances
        :ivar sections: Section -instances
        :ivar materials: Material -instances
        :ivar supports: Support -instances
        :ivar loads: Load -instances
        :ivar loadcases: LoadCase -instances
        :ivar dofs: number of degrees of freedom

        :vartype nodes: list
        :vartype nodal_coords: list
        :vartype elements: list
        :vartype sections: list
        :vartype materials: list
        :vartype supports: list
        :vartype loads: list
        :vartype loadcases: list
        :vartype dofs: int

    """

    def __init__(self):

        """ List containing all FEMNode -instances """
        self.nodes = []
        """ List of nodal coordinates """
        self.nodal_coords = []
        """ List of elements """
        self.elements = []
        """ List of sections """
        self.sections = []
        """ List of materials """
        self.materials = []
        """ List of supports """
        self.supports = []
        """ List of loads """
        self.loads = []
        """ List of load cases"""
        self.loadcases = []
        """ Number of degrees of freedom """
        self.dofs = 0
        """ Dimension of the problem (either 2 or 3) 
            Problem is 3D if any node as three coordinates
        """
        self.dim = 2
        
    def nels(self):
        """ Number of elements

        Returns:
        ---------
        :return: number of elements
        :rtype: int
        """
        return len(self.elements)

    def nnodes(self):
        """ Number of nodes

        Returns:
        --------
        :return: number of nodes
        :rtype: int

        """
        return len(self.nodes)

    def nsupp(self):
        """ Number of supported nodes

        Returns:
        --------
        :return: number of supported nodes
        :rtype: int

        """
        return len(self.supports)

    def nload(self):
        """ Number of loads

        Returns:
        --------
        :return: number of loads
        :rtype: int
        """
        return len(self.loads)

    def nloadcases(self):
        """ Number of load cases

        Returns:
        --------
        :return: number of load cases
        :rtype: int
        """
        return len(self.loadcases)

    def add_node(self, x, y, z=None):
        """ Adds node to fem model

        Parameters:
        -----------
        :param x: x-coordinate for node
        :param y: y-coordinate for node

        :type x: float
        :type y: float
        """

        newNode = FEMNode(len(self.nodes), x, y, z)
        if z is None:
            """ if z is 'None', there is no z-coordinate, and the problem
                is 2D. Otherwise, the value of z-coordinate is given, and the problem
                is 3D.
            """

            self.nodal_coords.append([x, y])            
            for i in range(self.nloadcases()):
                newNode.u = np.vstack((newNode.u, [0.0, 0.0, 0.0]))
        else:
            self.dim = 3
            self.nodal_coords.append([x, y, z])
            for i in range(self.nloadcases()):
                newNode.u = np.vstack((newNode.u, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

        self.nodes.append(newNode)
        return newNode

    def add_material(self, E, nu, density):
        """ Adds new material

        Parameters:
        -----------
        :param E: Young's modulus [MPa]
        :param nu: Poisson ratio
        :param density: density [kg / mm^3]

        :type E: float
        :type nu: float
        :type density: float
        """
        newMat = Material(E, nu, density)
        self.materials.append(newMat)
        return newMat

    def add_section(self, newSection):
        """ Add new cross-section

        Parameters:
        -----------
        :param newSection: section to be added to fem-model

        :type newSection: Section
        """
        self.sections.append(newSection)

    def add_element(self, new_elem):
        """ Add new element

        Parameters:
        -----------
        :param new_elem: element to be added

        :type new_elem: Element

        """
        self.elements.append(new_elem)

    def add_support(self, sid, nid, dof, val=0):
        """ Add support condition

        Parameters:
        -----------
        :param sid: support id
        :param nid: index of supported node
        :param dof: degrees of freedom
        :param val: value for support displacement

        :type sid: int
        :type nid: int
        :type dof: list
        :type val: float

        """
        new_support = Support(sid, self.nodes[nid], dof, val)
        self.supports.append(new_support)

    def add_load(self, load):
        """ Add load

        Parameters:
        ------------
        :param load: load to be added

        :type load: Load
        """
        self.loads.append(load)

    def add_loadcase(self, supp_id=1, load_id=2):
        """ Add load case

        Parameters:
        -----------
        :param supp_id: support id
        :param load_id: load id

        :type supp_id: int
        :type load_id: int
        """
        new_case = LoadCase(supp_id, load_id)
        self.loadcases.append(new_case)

        """ Add displacements to nodes for the new load case """
        for node in self.nodes:
            if len(node.u):
                # node.u = [0.0, 0.0, 0.0]
                # node.u = np.vstack((node.u, [0.0, 0.0, 0.0]))
                node.u = np.vstack(node.u,[0 for i in range(len(node.dofs))])
            else:
                node.u = [0 for i in range(len(node.dofs))]
                # node.u = np.vstack((node.u, [0.0, 0.0, 0.0]))

    def nodal_dofs(self):
        """ Assign nodal degrees of freedom
        """

        """ HERE THE ELEMENTS CAN BE USED TO SET SOME OF THE 
            DOFS TO ZERO. FOR EXAMPLE, TRUSS ELEMENTS DO NOT HAVE
            ROTATIONAL DEGREES OF FREEDOM, AND IF SOME NODES
            ARE CONNECTED TO ONLY TRUSS ELEMENTS, THOSE NODES WILL NOT
            HAVE A ROTATION.
        """
        for elem in self.elements:
            elem.init_dofs()

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



    def element_node_indices(self, e):
        """ Get node indices of an element

        Parameters:
        -----------
        :param e: element index

        :type e: int

        Returns:
        --------
        :return: list of element node indices
        :rtype: list
        """
        ndx = [self.nodes.index(n) for n in self.elements[e].nodes]
        # ndx = []
        # for n in self.elements[e].nodes:
        #    ndx.append(self.nodes.index(n))

        return ndx


    def global_stiffness_matrix(self):
        """ Constructs the global stiffness matrix

        Returns:
        --------
        :return: global stiffness matrix
        :rtype: np.array
        """

        """ initialize zero stiffness matrix 
            NOTE: for larger structures, K should probably be
            a sparse matrix
        """
        K = np.zeros((self.dofs, self.dofs))

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
            K[np.ix_(q, q)] += ke[np.ix_(nz, nz)]

        return K

    def global_geometric_stiffness_matrix(self):
        """ Constructs the global geometric stiffness matrix

        Returns:
        --------
        :return: global geometric stiffness matrix
        :rtype: np.array
        """
        dofs = self.dofs

        """ initialize zero geometric stiffness matrix 
            NOTE: for larger structures, KG should probably be
            a sparse matrix
        """
        KG = np.zeros((dofs, dofs))

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
            KG[np.ix_(q, q)] += kg[np.ix_(nz, nz)]

        return KG

    def global_load_vector(self, sid):
        """ Constructs global load vector

        Parameters:
        -----------
        :param sid: support id
        :type sid: int

        Returns:
        --------
        :return: global load vector
        :rtype: np.array
        """

        """ Initialize load vector """
        global_load = np.zeros(self.dofs)

        for load in self.loads:
            """ Get the loads to be added to the vector 
                and the global degrees of freedom of these loads
            """
            if load.sid == sid:
                v, vdofs = load.load_and_dofs()
                global_load[vdofs] += v

        # print("FrameFEM  ", '\n', global_load)
        return global_load

    def linear_statics(self, lcase=0):
        """ Perform linear elastic static analysis for a given load case

        Parameters:
        -----------
        :param lcase: index for load case

        :type lcase: int
        
        """
        if not self.dofs:
            self.nodal_dofs()

        K = self.global_stiffness_matrix()

        load_id = self.loadcases[lcase].load
        p = self.global_load_vector(load_id)
        
        
        """ Take supports into account
            The method is based on Filippa's Lecture Notes (3.5.2)
            
            The idea is that the rows and columns of the global stiffness
            matrix corresponding to the supported degrees of freedom are
            zeroed, and 1 is placed on the diagonal. Similarly, the
            elements of the load vector corresponding to the supported
            DOFs are zeroed.
        """
        supp_id = self.loadcases[lcase].support
        
        for supp in self.supports:
            if supp.sid == supp_id:                
                node_dofs = supp.node.dofs
                for dof in supp.dof:
                    # Get the global degree of freedom of the supported
                    # node's degree of freedom.
                    i = node_dofs[dof]
                    print(i)
                    if i > -1:
                        K[i,:] = 0
                        K[:,i] = 0
                        K[i,i] = 1
                        p[i] = 0
                        
        """ Solve global displacement vector """
        u = np.linalg.solve(K, p)
        
        
        """ Distribute displacements to nodes """
        for node in self.nodes:
            # Get nodal dofs
            dofs = np.array(node.dofs)

            # Find free dofs
            free_dofs = dofs >= 0

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


        return u, K

    def linear_buckling(self, lcase=0, k=4):
        """ Perform linear buckling analysis for a given load case

        Parameters:
        ------------
        :param lcase: load case id
        :param k: number of different buckling modes desired

        :type lcase: int
        :type k: int

        Returns:
        ---------
        :return: Array of k alpha_cr values and an vector representing buckling displacements

        :rtype: (np.array, np.array)
        """

        """ Solve internal forces of members 
            and return stiffness matrix
        """
        u, K = self.linear_statics(lcase)
        KG = self.global_geometric_stiffness_matrix()
        # v0 -- initial guess is set to np.ones -matrix,
        # this makes sure that w and v values doesn't change
        # after consecutive runs, default v0 is matrix with random values
        w, v = eigsh(K, k=k, M=-KG, sigma=0.0, which='LA', mode='normal', v0=np.ones(self.dofs))
        # w,v = sp.eigvals(K,b=-KG)

        """ Distribute buckling displacements to nodes """
        for node in self.nodes:
            # Get nodal dofs
            dofs = np.array(node.dofs)

            # Find free dofs
            free_dofs = dofs >= 0

            # Substitute free dofs from global displacement vector
                                    
            for i in range(len(free_dofs)):
                if free_dofs[i]:
                    node.v[i] = v[dofs[i]]
        
        return w, v

    def draw(self):
        """  Plots elements and nodes using matplotlib pyplot
        """

        fig = plt.figure()
        if self.dim == 3:
            ax = plt.axes(projection='3d')
        else:
            ax = plt.axes()

        """ draw nodes """
        for n in self.nodes:
            if self.dim == 2:
                ax.plot(n.coord[0], n.coord[1], 'ro')
            else:
                ax.scatter3D(n.coord[0], n.coord[1], n.coord[2],'ro')
                

        for i in range(self.nnodes()):
            if self.dim == 2:
                plt.text(self.nodes[i].coord[0], self.nodes[i].coord[1], str(i))
            #else:
                #plt.text(self.nodes[i].coord[0], self.nodes[i].coord[1], self.nodes[i].coord[2], str(i))

        """ draw members """
        for i in range(self.nels()):
            # X = self.member_coord(i)
            X = self.elements[i].coord()
            Xmid = X[0, :] + 0.5 * (X[1, :] - X[0, :])            
            
            if self.dim == 2:
                ax.plot(X[:, 0], X[:, 1], 'b')
                ax.text(Xmid[0], Xmid[1], str(i))
            else:
                ax.plot3D(X[:, 0], X[:, 1], X[:,2], 'b')
                ax.text(Xmid[0], Xmid[1], Xmid[2], str(i))

        plt.show()

    def print_displacements(self):
        """ Prints nodal displacements
        """
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


        Parameters:
        -----------
        :param E: Young's modulus [MPa]
        :param nu: Poisson's ratio
        :param rho: density [kg / mm^3]
        :param G: shear modulus. If it is not given, it is calculated from
                    E and nu

        :type E: float
        :type nu: float
        :type rho: float
        :type G: float

        Variables:
        ----------
        :ivar young: Young's modulus [MPa]
        :ivar nu: Poisson's ratio
        :ivar density: density [kg / mm^3]
        :ivar shear_modulus: Shear modulus [MPa]

        :vartype young: float
        :vartype nu: float
        :vartype density: float
        :vartype shear_modulus: float
    """

    def __init__(self, E, nu, rho, G=None):
        """ Young's modulus """
        self.young = E
        """ Poisson ratio """
        self.nu = nu
        """ Density """
        self.density = rho
        
        """ Shear modulus """
        if G == None:
            G = 0.5*E/(1+nu)
        
        self.shear_modulus = G
        


class Section:
    """ Class for cross section properties

        All sections must have 'area', but specific profiles
        can have other properties as well
    
        Parameters:
        -----------
        :param A: cross-sectional area [mm^2]

        :type A: float

        Variables:
        ----------
        :ivar area: cross-sectional area [mm^2]

        :vartype area: float

    """

    def __init__(self, A):
        self.area = A
        """ cross-sectional area [mm^2] """


class BeamSection(Section):
    """ Class for beam element cross sections 

        Parameters:
        -----------
        :param A: cross-sectional area [mm^2]
        :param Iy: second moment of area with respect to major axis [mm^4]
        :param Iz: second moment of area with respect to minor axis [mm^4]
        :param J: torsional constant [mm^4]

        :type A: float
        :type Iy: float
        :type Iz: float
        :type J: float

        Variables:
        ----------
        :ivar Iy: second moment of area [mm^4]

        :vartype Iy: float
    """

    def __init__(self, A, Iy, Iz=None, J=None):
        Section.__init__(self, A)
        """ Second moment of area [mm^4]"""
        self.Iy = Iy
        self.Iz = Iz
        """ Torsional constant """
        self.J = J
        


class Support:
    """ Supported nodes / degrees of freedom

        Parameters:
        -----------
        :param sid: support id
        :param node: supported node
        :param dof: degrees of freedom
        :param val: value for supports (default = 0.0, but theoretically
                    the displacement could have a non-zero value)

        :type sid: int
        :type node: FEMNode
        :type dof: list
        :type val: float

        Variables:
        -----------
        :ivar sid: support id
        :ivar node: supported node
        :ivar dof: degrees of freedom
        :ivar val: value for supports (default = 0.0, but theoretically
                    the displacement could have a non-zero value)

        :vqrtype sid: int
        :vqrtype node: FEMNode
        :vqrtype dof: list
        :vqrtype val: float

    """

    def __init__(self, sid, node, dof, val=0.0):
        self.sid = sid
        """ support id"""
        self.node = node
        """ supported node"""
        self.dof = dof
        """ degrees of freedom
            Alt. 1: dof = [1, 0,0 1, 0, 1] with ones for supported dofs
            Calling nodal_dofs() excludes the supported dofs
            
            Alt. 2: dof = [0, 3, 5] with indices of supported dofs
            Calling nodal_dofs() does not exclude the supported dofs but this
            is done in linear_statics()
        """
        
        #self.node.dofs = dof
        self.val = val
        """ value for supports"""


class Load:
    """ General class for loads

        Parameters:
        -----------
        :param sid: load id

        :type sid: int

        Variables:
        ----------
        :ivar sid: load id

        :vartype sid: int
    """

    def __init__(self, sid):
        self.sid = sid
        """ load id """

    @abstractclassmethod
    def load_and_dofs(self):
        """ Returns the loads to be inserted to the global load vector
            and the corresponding global degrees of freedom

            Returns:
            ---------
            :return: load to be inserted to the global load vector and
                    the corresponding global degrees of freedom
            :rtype: (int, int)
        """
        pass


class PointLoad(Load):
    """ Class for point loads (forces and moments)

        Parameters:
        ------------
        :param sid:  load id
        :param node: node subjected to the load
        :param v: load vector [Fx, Fy, Mz] [kN, kN, kNm]
        :param f: scaling factor

        :type sid: int
        :type node: FEMNode
        :type v: list
        :type f: float

        Variables:
        --------------
        :ivar sid:  load id
        :ivar node: node subjected to the load
        :ivar v: load vector [Fx, Fy, Mz] [kN, kN, kNm]
        :ivar f: scaling factor

        :vartype sid: int
        :vartype node: FEMNode
        :vartype v: list
        :vartype f: float

    """

    def __init__(self, sid, node, v, f):
        Load.__init__(self, sid)

        self.node = node
        """ node subjected to the load"""
        self.v = v
        """ load vector"""
        self.f = f
        """ scaling factor"""

    def load_and_dofs(self):
        """ Returns the loads to be inserted to the global load vector
            and the corresponding global degrees of freedom

            Returns:
            ---------
            :return: load to be inserted to the global load vector and
                    the corresponding global degrees of freedom
            :rtype: (int, int)

        """

        """ Load vector """
        F = self.f * np.array(self.v)

        
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
    """ Class for 2D line loads (forces and moments)

        Parameters:
        -----------
        :param sid: Load id
        :param eid: Element subjected to load
        :param xloc: Starting and ending locations of the load (local coordinates) [loc1, loc2]
        :param qval: Value of the load at the coordinates xloc [kN/m] [q1, q2]
        :param direction: Direction of load in global coordinates, 0: x, 1: y
        :param coords: (Optional, default=1) Coordinates used, 0: local, 1: global

        :type sid: int
        :type eid: Element
        :type xloc: list
        :type qval: list
        :type direction: int
        :type coords: int

    """

    def __init__(self, sid, eid, xloc, qval, direction, coords=1):

        Load.__init__(self, sid)

        self.elem = eid
        """ Elemnt subjected to load"""
        self.xloc = xloc
        """ Starting and ending locations of the load in local coordinates"""
        self.qval = qval
        """ Values of the load at the coordinates xloc"""
        self.dir = direction
        """ Direction of the load, 0: x, 1: y"""
        self.coords = coords
        """ Coordinates used, 0: local, 1: global (Default=1)"""

    def load_and_dofs(self):
        """ Returns the loads to be inserted to the global load vector
            and the corresponding global degrees of freedom

            Returns:
            ---------

            :return: load to be inserted to the global load vector and
                    the corresponding global degrees of freedom
            :rtype: (int, int)

        """

        # Get load vector
        if self.dir == "x":
            q = np.array([self.qval[0], 0, 0])
        else:
            q = np.array([0, self.qval[0], 0])

        # Get equivalent nodal loads
        F = self.elem.equivalent_nodal_loads(q)

        # Number of degrees of freedom
        dofs = np.array(self.elem.global_dofs())

        # Find free degrees of freedom (those with numbering greater than -1)
        nzero_dofs = dofs >= 0

        # Get free degrees of freedom and the corresponding part
        #    of the load vector

        d = dofs[nzero_dofs]
        Fnz = F[nzero_dofs]

        return Fnz, d


class LoadCase:
    """ Class for storing load case data

        Parameters:
        -----------
        :param support: Support id, supports having this id are added to the calculation model
        :param load: Load id, loads having this id are added to the calculation model

        :type support: int
        :type load: int

        Variables:
        -----------
        :ivar support: Supports id, supports having this id are added to the calculation model
        :ivar load: Loads id, loads having this id are added to the calculation model

        :vartype support: int
        :vartpe load: int

    """

    def __init__(self, support, load):
        self.support = support
        self.load = load


class FEMNode:
    """ Class for FEM nodes.

        Parameters:
        -----------
        :param nid: node id (integer)
        :param x: x-coordinate of the node
        :param y: y-coordinate of the node
        :param z: z-coordinate of the node (optional)

        :type x: float
        :type y: float
        :type z: float

        Variables:
        ----------
        :ivar x: Coordinates of the node
        :ivar dofs: Degrees of freedom for node [0, 1, 2] [Ux, Uy, Rz]
                In 3D: [Ux, Uy, Uz, Rx, Ry, Rz]
        :ivar u: Nodal displacements due to linear static analysis
        :ivar v: Nodal displacements due to linear buckling analysis

        :vartype x: np.array
        :vartype dofs: list
        :vartype u: np.array
        :vartype v: np.array


    """

    def __init__(self, nid, x, y, z=None):
        """ Constructor
            If z-coordinate is given, then the problem becomes immediately 3D
        """

        """ Node id """
        self.nid = nid
        

        """ Node's degrees of freedom """
        # NOTE: for multiple load cases, the dimension of u must be increased for each load case
        self.u = np.array([])

        """ Nodal coordinates"""
        if z is not None:
            self.coord = np.array([x, y, z])
            # Degrees of freedom (integer values)
            self.dofs = [0, 0, 0, 0, 0, 0]
            
            # Nodal displacements
            # [Ux, Uy, Uz, Rx, Ry, Rz]
            self.v = [0, 0, 0, 0, 0, 0]
        else:
            self.coord = np.array([x, y])
            self.dofs = [0, 0, 0]
            # [Ux, Uy, Uz]
            self.v = [0,0,0]
        
            
        """ Nodal displacement vector (linear buckling)"""
        # List of FrameMember-type objects that are connected to this node
        self.parents = []

    @property
    def x(self):
        return self.coord[0]

    @x.setter
    def x(self, val):
        self.coord[0] = val

    @property
    def y(self):
        return self.coord[1]

    @y.setter
    def y(self, val):
        self.coord[1] = val

    @property
    def z(self):
        return self.coord[2]

    @z.setter
    def z(self, val):
        self.coord[2] = val

class Element(metaclass=ABCMeta):
    """ Class for 1D finite elements: bars and beams

        Parameters:
        -----------
        :param n1: Start node
        :param n2: End node
        :param section: Element's section
        :param material: Element's material

        :type n1: FEMNode
        :type n2: FEMNode
        :type section: Section
        :type material: Material

        Variables:
        -----------
        :ivar nodes: Element's nodes
        :ivar material: Element's material
        :ivar section: Element's section
        :ivar axial_force: Axial force on element's nodes
        :ivar shear_force: Shear force on element's nodes
        :ivar bending_moment: Bending moment on element's nodes
        :ivar floc: Local loads

        :vartype nodes: list
        :vartype material: Material
        :vartype section: Section
        :vartype axial_force: list
        :vartype shear_force: list
        :vartype bending_moment: list
        :vartype floc: np.array

    """

    def __init__(self, n1, n2, section, material):
        self.nodes = [n1, n2]
        self.material = material
        self.section = section
        self.axial_force = [0.0, 0.0]
        self.shear_force = [0.0, 0.0]
        self.bending_moment = [0.0, 0.0]
        self.floc = np.array([])

    def coord(self):
        """ Nodal coordinates of the element

            Returns:
            --------
            :return: Nodal coordinates of the element
            :rtype: np.array

        """
        X = np.array([self.nodes[0].coord, self.nodes[1].coord])
        return X

    def direction_cosines(self):
        """ Calculates element's direction cosines

            Returns:
            --------
            :return: Element's direction cosines
            :rtype: np.array

        """
        X = self.coord()
        dX = X[1, :] - X[0, :]
        L = self.length()
        c = dX / L
        return c

    @abstractclassmethod
    def equivalent_nodal_loads(self, q):
        """ Compute equivalent nodal loads with line loads in
            vector q
        """
        pass

    def length(self):
        """ Member length

            Returns:
            --------
            :return: length of the element
            :rtype: float

        """
        X = self.coord()
        L = np.linalg.norm(X[1] - X[0])
        
        
        return round(L, 3)

    @abstractclassmethod
    def local_stiffness_matrix(self):
        """ Generates stiffness matrix in local coordinates

            Returns:
            --------
            :return: element's stiffness matrix in local coordinates
            :rtype: np.array

        """
        pass

    @abstractclassmethod
    def local_geometric_stiffness_matrix(self):
        """ Geometric stiffness matrix in local coordinates

            Returns:
            ---------
            :return: element's geometric stiffness matrix in local coordinates
            :rtype: np.array
        """

    @abstractclassmethod
    def stiffness_matrix(self):
        """ Computes element's global stiffness matrix

            Returns:
            ---------
            :return: element's stiffness matrix in global coordinates
            :rtype: np.array

        """
        pass

    def init_dofs(self):
        """ Initialize dofs 
            Does nothing for beam elements
        """
        pass

    @abstractclassmethod
    def global_dofs(self):
        """ Get numbering of element's global degrees of freedom

            Returns:
            --------
            :return: list of end nodes dofs added together
            :rtype: list
        """
        pass

    def nodal_displacements(self):
        """ Get nodal displacements of an element 
            Requires previously performed structural analysis such
            that nodal displacements are available.

            Returns:
            --------
            :return: Element's nodal displacements
            :rtype: np.array
        """
        return np.concatenate((self.nodes[0].u, self.nodes[1].u))

    def local_displacements(self):
        """ Nodal displacements in local coordinates

            Returns:
            ---------
            :return: Nodal displacements in local coordinates
            :rtype: np.array

        """
        T = self.transformation_matrix()
        q = self.nodal_displacements()
        
        #print(T)
        #print(q)
        return T.dot(q)

    @abstractclassmethod
    def transformation_matrix(self):
        """ Calculates transformation matrix from local to global coordinates

            Returns:
            --------
            :return: element's transformation matrix
            :rtype: np.array

         """
        pass

    @abstractclassmethod
    def internal_forces(self):
        """ Calculates internal forces and saves them to element's attributes

            axial_force = [n1, n2]-- axial force on node1 and node2
            shear_force = [n1, n2] -- shear force on node1 and node2
            bending_moment = [n1, n2] -- bending moment on node1 and node2
        """
        pass

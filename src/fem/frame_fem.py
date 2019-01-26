import sys
import math
import numpy as np
# import scipy as sp
import scipy.sparse.linalg as sp
import matplotlib.pyplot as plt


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

        Uses classes:
        -------------
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
    def __init__(self):

        self.nodes = []
        """ List containing all FEMNode -instances """
        self.nodal_coords = []
        """ List of nodal coordinates """
        self.elements = []
        """ List of elements """
        self.sections = []
        """ List of sections """
        self.materials = []
        """ List of materials """
        self.supports = []
        """ List of supports """
        self.loads = []
        """ List of loads """
        self.loadcases = []
        """ List of load cases"""
        self.dofs = 0
        """ Number of degrees of freedom """

    def nels(self):
        """
        Returns number of elements

        Retrurns:
        ---------
        :return: number of elements
        :rtype: int
        """
        return len(self.elements)

    def nnodes(self):
        """
        Number of nodes

        Returns:
        --------
        :return: number of nodes
        :rtype: int

        """
        return len(self.nodes)

    def nsupp(self):
        """
        Number of supported nodes

        Returns:
        --------
        :return: number of supported nodes
        :rtype: int

        """
        return len(self.supports)

    def nload(self):
        """
        Number of loads

        Returns:
        --------
        :return: number of loads
        :rtype: int
        """
        return len(self.loads)

    def nloadcases(self):
        """
        Number of load cases

        Returns:
        --------
        :return: number of load cases
        :rtype: int
        """
        return len(self.loadcases)

    def add_node(self, x, y):
        """
        Adds node to fem model

        Parameters:
        -----------
        :param x: x-coordinate for node
        :param y: y-coordinate for node

        :type x: float
        :type y: float
        """
        newNode = FEMNode(x, y)
        self.nodal_coords.append([x, y])

        for i in range(self.nloadcases()):
            newNode.u = np.vstack((newNode.u, [0.0, 0.0, 0.0]))

        self.nodes.append(newNode)

    def add_material(self, E, nu, density):
        """
        Adds new material

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

    def add_section(self, newSection):
        """
        Add new cross-section

        Parameters:
        -----------
        :param newSection: section to be added to fem-model

        :type newSection: Section
        """
        self.sections.append(newSection)

    def add_element(self, new_elem):
        """
        Add new element

        Parameters:
        -----------
        :param new_elem: element to be added

        :type new_elem: Element

        """
        self.elements.append(new_elem)

    def add_support(self, sid, nid, dof, val):
        """
        Add support condition

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
        """
        Add load

        Parameters:
        ------------
        :param load: load to be added

        :type load: Load
        """
        self.loads.append(load)

    def add_loadcase(self, supp_id=1, load_id=2):
        """
        Add load case

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
            if node.u.size == 0:
                node.u = [0.0, 0.0, 0.0]
            else:
                node.u = np.vstack((node.u, [0.0, 0.0, 0.0]))

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
                # self.nodes[s.nid].dofs[d] = 1

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
        """
        Get node indices of an element

        Parameters:
        -----------
        :param e: element index

        :type e: int

        Returns:
        --------
        :return: list of element node indices
        :rtype: list
        """
        ndx = []
        for n in self.elements[e].nodes:
            ndx.append(self.nodes.index(n))

        return ndx

    def global_stiffness_matrix(self):
        """
        Constructs the global stiffness matrix

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
        """
        Constructs the global geometric stiffness matrix

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
        """
        Constructs global load vector

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
        """
        Perform linear elastic static analysis for a given load case

        Parameters:
        -----------
        :param lcase: index for load case

        :type lcase: int
        
        """

        K = self.global_stiffness_matrix()

        load_id = self.loadcases[lcase].load
        p = self.global_load_vector(load_id)

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
        """
        Perform linear buckling analysis for a given load case

        Parameters:
        ------------
        :param lcase: load case id
        :param k: number of different buckling modes desired

        :type lcase: int
        :type k: int

        Returns:
        ---------
        :return: Array of k alpha_cr values and  an array representing buckling displacements

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
        w, v = sp.eigsh(K, k=k, M=-KG, sigma=0.0, which='LA', mode='normal', v0=np.ones(self.dofs))
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
        """
        Plots elements and nodes using matplotlib pyplot
        """

        """ draw nodes """
        for n in self.nodes:
            plt.plot(n.x[0], n.x[1], 'ro')

        for i in range(self.nnodes()):
            plt.text(self.nodes[i].x[0], self.nodes[i].x[1], str(i))

        """ draw members """
        for i in range(self.nels()):
            # X = self.member_coord(i)
            X = self.elements[i].coord()
            plt.plot(X[:, 0], X[:, 1], 'b')
            Xmid = X[0, :] + 0.5 * (X[1, :] - X[0, :])
            plt.text(Xmid[0], Xmid[1], str(i))

        plt.show()

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


        Parameters:
        -----------
        :param E: Young's modulus [Mpa]
        :param nu: Poisson's ratio
        :param rho: density [kg / mm^3]

        :type E: float
        :type nu: float
        :type rho: float

        Variables:
        ----------
        :ivar young: Young's modulus [Mpa]
        :ivar nu: Poisson's ratio
        :ivar density: density [kg / mm^3]

        :vartype young: float
        :vartype nu: float
        :vartype density: float
    """

    def __init__(self, E, nu, rho):


        self.young = E
        """ Young's modulus """
        self.nu = nu
        """ Poisson ratio """
        self.density = rho
        """ Density """


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
        :param Iy: second moment of area [mm^4]

        :type A: float
        :type Iy: float

        Variables:
        ----------
        :ivar Iy: second moment of area [mm^4]

        :vartype Iy: float
    """

    def __init__(self, A, Iy):
        Section.__init__(self, A)
        self.Iy = Iy
        """ Second moment of area [mm^4]"""


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
        """ degrees of freedom"""
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







    """

    def __init__(self, sid, eid, xloc, qval, direction, coords=1):
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
        Load.__init__(self, sid)

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
            q = np.array([self.qval[0], 0, 0])
        else:
            q = np.array([0, self.qval[0], 0])

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

    def __init__(self, support, load):
        self.support = support
        self.load = load


class FEMNode:
    """ Class for FEM nodes.
        
        Basically only for storing data
    """

    def __init__(self, x, y):
        """ x and y are the coordinates of the node 
            
            dofs -- numbers of nodal degrees of freedom
            u -- nodal displacements in the global coordinate system
            v -- nodal displacements for linear buckling 
        """
        self.x = np.array([x, y])
        """ Setting all dofs to 0 means that the dofs are free.
            The method nodal_dofs of class FrameFEM assigns the actual
            numbers to the dofs.
            
            dofs[0] -- ux
            dofs[1] -- uy
            dofs[2] -- rz
            
            For 3D frames, the order of dofs must be redefined
        """
        self.dofs = [0, 0, 0]
        # maybe u could also be a simple list?
        # self.u = np.array([0.0,0.0,0.0])
        """ Initialize displacements as an empty array
            NOTE: for multiple load cases, the dimension of u must be increased
            for each load case
        """
        self.u = np.array([])
        self.v = [0, 0, 0]
        # List of FrameMember-type objects that are connected to this node
        self.parents = []

        # Forces
        self.Fx = 0
        self.Fy = 0
        self.Mz = 0


class Element:
    """ 1D finite elements: bars and beams
        Attributes:
            
        nodes -- array of FEMNode variables, member end nodes
        material -- material properties (Material)
        section -- cross-section properties (Section)
        axial_force -- axial force, list
        
    """

    def __init__(self, n1, n2, section, material):
        """ Input:
            n1 -- node 1 (reference to FEMNode variable)
            n2 -- node 2 (reference to FEMNode variable)
            section -- cross-section properties
            material -- Material type variable
                        
        """
        self.nodes = [n1, n2]
        self.material = material
        self.section = section
        self.axial_force = [0.0, 0.0]
        self.floc = np.array([])

    def coord(self):
        """ Nodal coordinates of the element """
        X = np.array([self.nodes[0].x, self.nodes[1].x])
        return X

    def length(self):
        """ Member length """
        X = self.coord()
        L = np.linalg.norm(X[1, :] - X[0, :])
        return L

    def direction_cosines(self):
        """ Direction cosines """
        X = self.coord()
        dX = X[1, :] - X[0, :]
        L = self.length()
        c = dX / L
        return c

    def stiffness_matrix(self):
        return None

    def equivalent_nodal_loads(self, q):
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

        return np.hstack((self.nodes[0].u, self.nodes[1].u))

    def local_displacements(self):
        """ Nodal displacements in local coordinates """
        T = self.transformation_matrix()
        q = self.nodal_displacements()
        return T.dot(q)

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

        if self.floc.size > 0:
            R = ke.dot(q) - self.floc
        else:
            R = ke.dot(q)

        """ Any load on the element not acting on a node must be
            taken into account here. This requires the following steps:
            1. Identify loads acting on element
            2. Compute their equivalent nodal loads
            3. Modify internal forces
            
            Probably an easy solution is to include an
            Attribute floc for equivalent nodal loads. The nodal loads are
            saved there, when they are first calculated for analysis
        """

        self.axial_force[:2] = [-R[0], R[3]]
        self.bending_moment[:2] = [-R[2], R[5]]
        self.shear_force[:2] = [R[1], -R[4]]
        # Save force and moment values to nodes
        # This may cause issues when the same node is connected to elements
        # that are in different member
        self.nodes[0].Fx, self.nodes[1].Fx = self.axial_force
        self.nodes[0].Fy, self.nodes[1].Fy = self.shear_force
        self.nodes[0].Mz, self.nodes[1].Mz = self.bending_moment

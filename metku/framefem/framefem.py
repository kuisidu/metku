# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy as sp
import time
from copy import deepcopy


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
        :vartype loadcases: dict
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
        """ Dict of load cases
            The key is the load id of the corresponding load case
            Value is a LoadCase object
        """
        self.loadcases = {}
        """ Number of degrees of freedom """
        self.dofs = 0
        """ Dimension of the problem (either 2 or 3) 
            Problem is 3D if any node as three coordinates
        """
        self.dim = 2
        """ Supported degrees of freedom """
        self.__supp_dofs = []
        self.supp_dofs = []
        
        self.load_factors = []
        self.buckling_modes = []
        
    def __repr__(self):
        
        s = f"Frame {self.dim}D model: {self.nels()} elements, {self.nnodes()} nodes"
        
        return s
        
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
    
    @property
    def supp_dofs(self):
        """ Get supported degrees of freedom """
        
        if len(self.__supp_dofs) == 0:
            self.supported_dofs()
        
        return self.__supp_dofs
    
    @supp_dofs.setter
    def supp_dofs(self,val):
        """ Set supported degrees of freedom """        
        
        self.__supp_dofs = val
        
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
            for load_id in self.loadcases.keys():   
                newNode.u[load_id] = np.zeros(len(newNode.dofs))
                #newNode.u[load_id] = [0.0, 0.0, 0.0]
        else:
            self.dim = 3
            self.nodal_coords.append([x, y, z])
            #for i in range(self.nloadcases()):
            #    newNode.u = np.vstack((newNode.u, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

        self.nodes.append(newNode)
        return newNode

    @property
    def nfree_dofs(self):
        """ number of free degrees of freedom """
        if len(self.supp_dofs) == 0:
            self.supported_dofs()
        
        nsupp_dofs = len(self.supp_dofs)
        
        return self.dofs-nsupp_dofs
    
    def element_free_dofs(self,el):
        """ Returns free dofs of an element """

        if isinstance(el,int):
            el = self.elements[el]
        
        el_dofs = el.global_dofs()
        supp_dofs = self.supp_dofs

        return [a for a in el_dofs if not a in supp_dofs]
        
        
    def add_material(self, E, nu, density, G=None):
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
        newMat = Material(E, density, nu=nu, G=G)
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
        
        for lcase in self.loadcases.keys():
            new_elem.add_loadcase(lcase)
        
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
        """ Adds new load case, if load id is same replaces old one

        Parameters:
        -----------
        :param supp_id: support id
        :param load_id: load id

        :type supp_id: int
        :type load_id: int
        """
        #print('Create Load case.')
        new_case = LoadCase(supp_id, load_id)
        self.loadcases[load_id] = new_case

        """ Add displacements to nodes for the new load case """
        ncases = self.nloadcases()
        for node in self.nodes:
            node.u[load_id] = np.zeros(len(node.dofs))
            """
            if ncases == 1:
                print(node.u[0])
                node.u[0] = np.zeros(len(node.dofs))
            else:
                node.u.append(np.zeros(len(node.dofs)))
            """
            """    
            if len(node.u) > 0:
                # node.u = [0.0, 0.0, 0.0]
                # node.u = np.vstack((node.u, [0.0, 0.0, 0.0]))
                node.u = np.vstack(node.u,[0 for i in range(len(node.dofs))])
            else:
                node.u = [0 for i in range(len(node.dofs))]
                # node.u = np.vstack((node.u, [0.0, 0.0, 0.0]))
            """        
        for el in self.elements:
            el.add_loadcase(load_id)
            
    def add_release(self,eid,dofs):
        """ Adds releases to an element
            input: eid .. element number
                   dofs .. list of degrees of freedom for the given element
                           for which the corresponding internal forces
                           are released (set to zero)
        """
        print("vapautus")
        self.elements[eid].releases = dofs
    
    def nodal_dofs(self):
        """ Assign nodal degrees of freedom
        """

        """ HERE THE ELEMENTS CAN BE USED TO SET SOME OF THE 
            DOFS TO ZERO. FOR EXAMPLE, TRUSS ELEMENTS DO NOT HAVE
            ROTATIONAL DEGREES OF FREEDOM, AND IF SOME NODES
            ARE CONNECTED TO ONLY TRUSS ELEMENTS, THOSE NODES WILL NOT
            HAVE A ROTATION.
        """
        self.dofs = 0
        
        for elem in self.elements:
            elem.init_dofs()

        """ Make running numbering of dofs.
            If the dof is marked as 1, then its dof is set to 0.
        """
        for n in self.nodes:
            for i in range(len(n.dofs)):
                #if n.dofs[i] > 0:
                """if not n.dofs[i]:
                    print('Not dof')
                    n.dofs[i] = -1
                else:
                """
                if n.dofs[i] >= 0:
                    #print('Dof')
                    n.dofs[i] = self.dofs
                    self.dofs += 1
            #print(n.dofs)



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

    def global_free_dofs(self):
        """ List of free degrees of freedom """
        
        return list(set(np.arange(self.dofs))-set(self.supp_dofs))

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
            q = ve[nz].astype(int)
            """ extract the submatrix of the element
                and add it to the global stiffness matrix
            """
            # a = np.arange(len(q))
            # row = np.vstack(a[nz])
            # col = a[nz]
            # K[np.vstack(q), q] += ke[row, col]
            #print(ke[np.ix_(nz, nz)])
            #print(ve,q)
            K[np.ix_(q, q)] += ke[np.ix_(nz, nz)]

        return K

    def global_geometric_stiffness_matrix(self,lcase=0):
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
            kg = elem.geometric_stiffness_matrix(lcase)
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
                #print('Calculate nodal forces')
                v, vdofs = load.load_and_dofs()
                #print(v,vdofs)
                global_load[vdofs.astype(int)] += v

        # print("FrameFEM  ", '\n', global_load)
        return global_load
    
    def supported_dofs(self):
        """ Determine supported dofs. 
        
            Assumption: there is only one support case, so no need to
            distinguish between different support IDs.
        """
        supp_dofs = []
        for supp in self.supports:
            supp_dofs += list(supp.node.dofs[supp.dof])      
        
        self.supp_dofs = supp_dofs

    def linear_statics(self, lcase=0, support_method='ZERO'):
        """ Perform linear elastic static analysis for a given load case

        Parameters:
        -----------
        :param lcase: index for load case
        :param support_method: how supported DOFS are treated
                              'ZERO' .. columns and rows corresponding to
                              supported DOFS are zeroed and a 1 is put
                              to the diagonal
                              
                              'REM' .. columns and rows corresponding to the
                              supported DOFS are removed.

        :type lcase: int
        
        """
        if not self.dofs:
            self.nodal_dofs()

        #start = time.time()
        K = self.global_stiffness_matrix()
        
        #print(lcase)
        #print(self.loadcases)
        load_id = self.loadcases[lcase].load        
        p = self.global_load_vector(load_id)
        
       
        #print(p)
        """ Take supports into account
            The method is based on Filippa's Lecture Notes (3.5.2)
            
            The idea is that the rows and columns of the global stiffness
            matrix corresponding to the supported degrees of freedom are
            zeroed, and 1 is placed on the diagonal. Similarly, the
            elements of the load vector corresponding to the supported
            DOFs are zeroed.
        """
        supp_id = self.loadcases[lcase].support
        
        rem_dofs = []
        """ Vector of global unsupported dofs """
        glob_dofs = np.arange(self.dofs)
        
        
        for supp in self.supports:
            if supp.sid == supp_id:                
                node_dofs = supp.node.dofs
                for dof in supp.dof:
                    # Get the global degree of freedom of the supported
                    # node's degree of freedom.
                    i = node_dofs[dof]
                    #  print(i)
                    if i > -1:
                        if support_method == 'ZERO':
                            K[i,:] = 0
                            K[:,i] = 0
                            K[i,i] = 1
                            p[i] = 0
                        elif support_method == 'REM':
                            rem_dofs.append(i)
        
        self.supp_dofs = rem_dofs
        
        
        if support_method == 'REM':
            rem_dofs = np.array(rem_dofs,dtype=int)
            
            K = np.delete(K,rem_dofs,0)
            K = np.delete(K,rem_dofs,1)
            p = np.delete(p,rem_dofs)
            glob_dofs = np.delete(glob_dofs,rem_dofs)
                        
            #print(rem_dofs)
            #print(K)
            #print(p)
            self.K = K
            self.p = p
            #print(glob_dofs)
            
                        
        """ Solve global displacement vector """

        #Kc, low = sp.linalg.cho_factor(K)
        #uc = sp.linalg.cho_solve((Kc, low), p)

        # print(f'K shape {K.shape}')
        # print(f'K {K}')
        # print(f'p {p}')
        u = np.linalg.solve(K, p)
        self.u = u
        #print(u)

        """ Substitute obtained displacements to nodes """
        if support_method == 'ZERO':
            """ Distribute displacements to nodes """
            for node in self.nodes:
                # Get nodal dofs
                dofs = np.array(node.dofs)
    
                # Find free dofs
                free_dofs = dofs >= 0
                # Substitute free dofs from global displacement vector                        
                for i, free_dof in enumerate(free_dofs):
                    if free_dof:
                        node.u[lcase][i] = u[dofs[i]]

            """
            This substition did not work !!!
            node.u[free_dofs] = u[dofs[free_dofs]]
            """
        elif support_method == 'REM':
            for ui, d in zip(u,glob_dofs):                
                for node in self.nodes:
                    try:
                        #node.u[lcase][node.dofs.index(d)] = np.array(ui)
                        node.u[lcase][np.where(node.dofs==d)] = np.array(ui)
                    except ValueError:
                        pass
                    
                    


        """ Calculate element internal forces """
        for i, el in enumerate(self.elements):            
            el.internal_forces(lcase=lcase)
        #end = time.time()
        # print("FRAMEFEM TIME: ", end - start)
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
        
        TODO!!!
            the global geometric stiffness matrix does not include any fixing of
            supported dofs now.
            
            This is because the logic for fixing supported dofs in linear statics was
            changed from row and column elimination to setting elements to zero or one.
            
            Solution:
                when using linear buckling analysis, always use "REM" as method of
                support handling. Then, remove the rows and columns of the global
                geometric stiffness matrix as well.
        """

        """ Solve internal forces of members 
            and return stiffness matrix
        """
        u, K = self.linear_statics(lcase,support_method="REM")
        KG = self.global_geometric_stiffness_matrix(lcase)
        
        """
        for supp in self.supports:
            if supp.sid == supp_id:                
                node_dofs = supp.node.dofs
                for dof in supp.dof:
                    # Get the global degree of freedom of the supported
                    # node's degree of freedom.
                    i = node_dofs[dof]
                    #  print(i)
                    if i > -1:
                        if support_method == 'ZERO':
                            KG[i,:] = 0
                            KG[:,i] = 0
                            KG[i,i] = 1                            
                        elif support_method == 'REM':
                            rem_dofs.append(i)
        
        if support_method == 'REM':
            KG = np.delete(K,rem_dofs,0)
            KG = np.delete(K,rem_dofs,1)            
            glob_dofs = np.delete(glob_dofs,rem_dofs)
        """
        KG = np.delete(KG,self.supp_dofs,0)
        KG = np.delete(KG,self.supp_dofs,1) 
        
        
        # v0 -- initial guess is set to np.ones -matrix,
        # this makes sure that w and v values doesn't change
        # after consecutive runs, default v0 is matrix with random values
        (w, buckling_modes) = eigsh(K, k=k, M=-KG, sigma=1.0, which='LA', mode='buckling')
        #w = sp.linalg.eigvals(K,b=-KG)
        #(w,buckling_modes) = sp.linalg.eig(K,b=-KG)
        #(w,buckling_modes) = sp.sparse.linalg.eigs(K,k=6,M=-KG,which='SM',sigma=0)

        self.load_factors = w
        self.buckling_modes = buckling_modes

        """ Distribute buckling displacements to nodes """
        
        """
        for node in self.nodes:
               
            # Get nodal dofs
            
            dofs = np.array(node.dofs)
            
            # Find free dofs
            free_dofs = dofs >= 0
    
            # Substitute free dofs from global displacement vector
                                        
            for i in range(len(free_dofs)):
                if free_dofs[i]:
                    node.v[i] = v[dofs[i]]
        """
        
        for node in self.nodes:
            node.v = np.zeros((buckling_modes.shape[1],len(node.v)))
        
        
        glob_dofs = self.global_free_dofs()
        
        for i in range(buckling_modes.shape[1]):
            v = buckling_modes[:,i]
            for vi, d in zip(v,glob_dofs):                
                for node in self.nodes:
                    try:
                        node.v[i][list(node.dofs).index(d)] = vi
                    except ValueError:
                        pass
        
        return w, buckling_modes, KG

    def statics_matrix(self):
        """ Creates the statix matrix. This applicable for instances that
            only have Rod elements (extensible to EBBeam elements as well).
            
            Primary use of this method is for optimization of trusses (and frames)
            using the SAND approach
        """
        
        B = np.zeros([self.dofs,self.nels()])
        
        e = np.array([-1,1])
        for i, el in enumerate(self.elements):
            b = el.transformation_matrix().transpose().dot(e)
            d = el.global_dofs()
            B[np.ix_(d),i] = b
                        
        
        return np.delete(B,self.supp_dofs,0)
        
        #return B
        

    #def draw(self, deformed=False, buckling_mode=None, scale=1.0, axes=None):
    def draw(self, axes=None,**kwargs):
        #deformed=False, buckling_mode=None, scale=1.0, axes=None):
        """  Plots elements and nodes using matplotlib pyplot
        """

        #fig = plt.figure()
        if self.dim == 3:
            #ax = plt.axes(projection='3d')
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            if axes is None:
                fig, ax = plt.subplots(1)
            else:
                ax = axes

        """ draw nodes """
        
        for i, node in enumerate(self.nodes):
            node.plot(print_text=str(i),axes=ax)
        
        """ draw elements """
        if 'deformed' in kwargs:
            deformed = True
        else:
            deformed = False
        
        if 'buckling_mode' in kwargs:
            buckling_mode = True
            buck_mode = kwargs['buckling_mode']
        else:
            buckling_mode = False
        
        if 'scale' in kwargs:
            scale = kwargs['scale']
        else:
            scale = 1.0
        
        """ Draw elements """
        for i, el in enumerate(self.elements):
            """ Draw deformed shape is requested.
                the value of 'deformed' is the load case to be plotted
            """
            if 'deformed' in kwargs:
                lcase = kwargs['deformed']
                
                el.plot(print_text=str(i),axes=ax,deformed=lcase,scale=scale)
            else:
                el.plot(print_text=str(i),axes=ax)
                
        el_col = 'k'
        if deformed or buckling_mode is not None:
            lstyle = '--' #+ el_col
        else:
            lstyle = '-' + el_col
        
        # Harmaa: color = (0.6,0.6,0.6)
        if deformed:
            """
            for n in self.nodes:
                if self.dim == 2:
                    ax.plot(n.coord[0]+scale*n.u[0], n.coord[1]+scale*n.u[1], 'bo')
                else:
                    ax.scatter3D(n.coord[0], n.coord[1], n.coord[2],'ro')
            """
            lstyle = '-k'
            for el in self.elements:
                X = np.zeros((2,2))
                for i, n in enumerate(el.nodes):
                    X[i,:] = (n.coord + scale * np.array(n.u[:2]))
                    
                if self.dim == 2:
                    ax.plot(X[:, 0], X[:, 1], lstyle)                  
                else:
                    ax.plot3D(X[:, 0], X[:, 1], X[:,2], lstyle)
                
        
        if buckling_mode:
            lstyle = '-k'
            for el in self.elements:
                #print(el.nodes)
                X = np.zeros((2,2))
                for i, n in enumerate(el.nodes):
                    #print(i, buckling_mode, n.v[buck_mode][:2])
                    X[i,:] = (n.coord + scale * n.v[buck_mode][:2])
                    
                if self.dim == 2:
                    ax.plot(X[:, 0], X[:, 1], lstyle)                  
                else:
                    ax.plot3D(X[:, 0], X[:, 1], X[:,2], lstyle)
        
        if self.dim == 2:
            ax.set_aspect('equal')
            
        if self.dim == 3:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
        
        """
        if show:
            plt.show()
        """
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

    def __init__(self, E, rho, nu=None, G=None):
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
        self.A = A
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
        
        self.I = [0,0]
        self.I[0] = Iy
        
        if Iz is not None:
            self.I[1] = Iz
        
        #self.Iy = Iy
        #self.Iz = Iz
        """ Torsional constant """
        self.J = J
    
    @property
    def Iy(self):
        
        return self.I[0]

    @property
    def Iz(self):
        
        return self.I[1]

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

    def __init__(self, sid, eid, xloc, qval, direction, coords='global'):

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
        if isinstance(self.elem,list):
            print(self.elem)
        F = self.elem.equivalent_nodal_loads(self)
        
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
        #self.u = np.array([])
        self.u = {}
        
        """ Nodal coordinates"""
        if z is not None:
            self.coord = np.array([x, y, z])
            # Degrees of freedom (integer values)
            #self.dofs = [0, 0, 0, 0, 0, 0]
            self.dofs = -np.ones(6)
            # Nodal displacements
            # [Ux, Uy, Uz, Rx, Ry, Rz]
            self.v = np.array([0, 0, 0, 0, 0, 0])
        else:
            self.coord = np.array([x, y])
            self.dofs = np.array([-1, -1, -1])
            # [Ux, Uy, Uz]
            self.v = np.array([0,0,0])
        
            
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
        if len(self.coord) == 2:
            return None
        return self.coord[2]

    @z.setter
    def z(self, val):
        self.coord[2] = val

    def __repr__(self):
        return f'FEMNode(nid={self.nid}, x={self.x}, y={self.y}, z={self.z})'

    def plot(self, print_text='', c='k', axes=None, style='.r'):

        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes

        if self.z is None:
            #ax.plot(X[:, 0], X[:, 1], lstyle, color=(0.6,0.6,0.6))
            ax.plot(self.x, self.y, style)
            ax.text(self.x, self.y, print_text)
        else:
            ax.scatter3D(self.x, self.y, self.z, style)
            ax.text(self.x, self.y, self.z, print_text)

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
        :vartype fint: dict .. internal forces
        :vartype axial_force: list
        :vartype shear_force: list
        :vartype bending_moment: list
        :vartype floc: np.array

    """

    def __init__(self, n1, n2, section, material):
        self.nodes = [n1, n2]
        self.material = material        
        self.section = section
        self.fint = {}
        #self.axial_force = [0.0, 0.0]
        #self.shear_force = [0.0, 0.0]
        #self.bending_moment = [0.0, 0.0]
        self.floc = {} #np.array([])
        
        """ Releases is a list of degrees of freedom for those
            internal forces that should be set to zero
        """
        self.releases = []
        
        """ Dimension of the element """
        self.dim = 2
        
        #self.x = self.coord()
        #self.len = self.length()
        #self.dcos = self.direction_cosines()

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
        
        #return (self.x[1]-self.x[0])/self.len

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
        
        #X = self.__x
        X = self.coord()
        L = np.linalg.norm(X[1] - X[0])        
        return round(L, 3)
        #return round(np.linalg.norm(self.x[1]-self.x[0]),3)

    @abstractclassmethod
    def local_stiffness_matrix(self, E, A, I, L):
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

    def nodal_displacements(self,lcase=0):
        """ Get nodal displacements of an element in global coordinates
            Requires previously performed structural analysis such
            that nodal displacements are available.

            Returns:
            --------
            :return: Element's nodal displacements
            :rtype: np.array
        """
        return np.concatenate((self.nodes[0].u[lcase],
                               self.nodes[1].u[lcase]))

    def local_displacements(self, lcase=0):
        """ Nodal displacements in local coordinates

            Returns:
            ---------
            :return: Nodal displacements in local coordinates
            :rtype: np.array

        """
        T = self.transformation_matrix()
        q = self.nodal_displacements(lcase)

        return T.dot(q)

    def add_loadcase(self,lcase=0):
        """ Creates a place holder for internal forces 
            fx .. axial force
            fy .. shear force with respect to horizontal axis
            fz .. shear force with respect to vertical axis
            mx .. torsion
            my .. bending moment with respect to horizontal axis
            mz .. bending moment with respect to vertical axis
            
            In 2d, the following internal forces are present
            fx, fz, my
        """
        
        self.fint[lcase] = {'fx': [0,0], 'fy': [0,0], 'fz': [0,0],\
                            'mx': [0,0], 'my': [0,0], 'mz': [0,0]}
        
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
    
    def plot(self, print_text='', c='k', axes=None, lstyle='-',
             deformed=None,buckling_mode=None,scale=1.0):

        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes

        X = self.coord()
        
        # Plot members
        #ax.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], c)
        # Plot text

        Xmid = X[0, :] + 0.5 * (X[1, :] - X[0, :])            
            
        if self.dim == 2:
            #ax.plot(X[:, 0], X[:, 1], lstyle, color=(0.6,0.6,0.6))
            if deformed is not None:                
                lstyle = '--'
                dstyle = '-'
                c = (0.6,0.6,0.6)
                Xd = deepcopy(X)
                u = np.array([n.u[deformed][:2] for n in self.nodes])
                
                Xd += scale*u
                
                ax.plot(Xd[:, 0], Xd[:, 1], dstyle, color='k')
                
            # Plot initial element
            ax.plot(X[:, 0], X[:, 1], lstyle, color=c)
                
            ax.text(Xmid[0], Xmid[1], print_text)
        else:
            ax.plot3D(X[:, 0], X[:, 1], X[:,2], lstyle)
            ax.text(Xmid[0], Xmid[1], Xmid[2], print_text)

        """
        if self.mtype == 'beam':
            horzalign = 'center'
            vertalign = 'bottom'

        elif self.mtype == 'column':
            horzalign = 'right'
            vertalign = 'center'

        else:
            horzalign = 'center'
            vertalign = 'center'

        x, y = self.to_global(0.3) - self.perpendicular * 50
        rot = np.degrees(self.angle)

        if print_text:
            ax.text(x, y, str(self.mem_id) + ": " + str(self.cross_section),
                     rotation=rot, horizontalalignment=horzalign,
                     verticalalignment=vertalign)
        else:
            ax.text(x, y, str(self.mem_id),
                     rotation=rot, horizontalalignment=horzalign,
                     verticalalignment=vertalign)
        """
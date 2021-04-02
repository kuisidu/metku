# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:14:54 2018

Euler-Bernoulli beam element in 2D

@author: kmela
"""

import numpy as np
from functools import lru_cache

CACHE_BOUND = 2**10

try:
    from metku.framefem import Element, LineLoad, FEMNode
except:
    from framefem import Element, LineLoad


class EBBeam(Element):
    """ Euler-Bernoulli beam element

        Parameters:
        -----------
        :param n1: start node
        :param n2: end node
        :param section: element's section properties
        :param material: element's material properties

        :type n1: FEMNode
        :type n2: FEMNode
        :type section: Section
        :type material: Material

        Variables:
        ----------
        :ivar bending_moment: bending moment values at element's ends
        :ivar shear_force: shear force values at element's ends

        :vartype bending_moment: list
        :vartype shear_force: list
    """

    def __init__(self, n1, n2, section, material):
        Element.__init__(self, n1, n2, section, material)

        #self.bending_moment = [0.0, 0.0]
        #self.shear_force = [0.0, 0.0]
        
        self.Krel = {'K11':None, 'K12':None, 'K21':None, 'K22':None}


    def transformation_matrix(self):
        """ Calculates transformation matrix from local to global coordinates
            Chandrupatla 3rd. Edition, p. 249
            Returns:
            --------
            :return: element's transformation matrix
            :rtype: np.array

         """
        c = self.direction_cosines()
        L = np.zeros((6, 6))
        L[0, 0:2] = c[0:2]
        L[1, 0] = -L[0, 1]
        L[1, 1] = L[0, 0]
        L[2, 2] = 1.0
        L[3:6, 3:6] = L[0:3, 0:3]
        return L


    def init_dofs(self):
        """ For beam members all degrees of freedom
            are in play.
        """
        for n in self.nodes:
            n.dofs[:3] = 0
                        

    def global_dofs(self):
        """ Get numbering of element's global degrees of freedom

            Returns:
            --------
            :return: list of end nodes dofs added together
            :rtype: list
        """
        return np.append(self.nodes[0].dofs,self.nodes[1].dofs)
        #return self.nodes[0].dofs + self.nodes[1].dofs

    @lru_cache(CACHE_BOUND)
    def local_stiffness_matrix(self, E, A, I1, Le):
        """ Generates stiffness matrix in local coordinates


            Returns:
            --------
            :return: element's stiffness matrix in local coordinates
            :rtype: np.array

        """
        rodc = E * A / Le
        EI1 = E * I1
        k0 = np.zeros((6, 6))
        #
        k0[[0, 3], [0, 3]] = rodc
        k0[[0, 3], [3, 0]] = -rodc
        #
        bend1 = 12 * EI1 / Le ** 3
        k0[[1, 4], [1, 4]] = bend1
        k0[[1, 4], [4, 1]] = -k0[1, 1]
        #
        bend2 = 6 * EI1 / Le ** 2
        k0[[1, 2, 1, 5], [2, 1, 5, 1]] = bend2
        #
        k0[[2, 4, 4, 5], [4, 2, 5, 4]] = -bend2
        #
        k0[[2, 5], [5, 2]] = 2 * EI1 / Le
        #
        k0[[2, 5], [2, 5]] = 4 * EI1 / Le
        
        if len(self.releases) > 0:
            """ There are releases in the element, so the stiffness matrix
                is modified.
            """
            rel = self.releases
            # nrel is the list of non-released forces
            nrel = np.setdiff1d(np.arange(6),rel)
            K22 = k0[np.ix_(rel,rel)]
            K11 = k0[np.ix_(nrel,nrel)]
            K12 = k0[np.ix_(nrel,rel)]
            K21 = k0[np.ix_(rel,nrel)]
            
            self.Krel['K11'] = K11
            self.Krel['K12'] = K12
            self.Krel['K21'] = K21
            self.Krel['K22'] = K22
            
            """
            print(rel,nrel)
            print('k0')
            print(k0)
            print('K22')
            print(K22)
            print('K11')
            print(K11)
            print('K12')
            print(K12)
            print('K21')
            print(K21)
            """
            #K221K21inv = np.linalg.inv(K22).dot(K21)
            kc = K11 - K12.dot(np.linalg.inv(K22).dot(K21))
            
            """
            print('kc = ')
            print(kc)
            print(kc.shape)
            """
            k0[np.ix_(nrel,nrel)] = kc
            k0[rel,:] = 0
            k0[:,rel] = 0
            #print(k0)

        return k0

    def stiffness_matrix(self):
        """ Computes element's global stiffness matrix

            Returns:
            ---------
            :return: element's stiffness matrix in global coordinates
            :rtype: np.array

        """
        E = self.material.young
        A = self.section.A
        I1 = self.section.I[0]
        Le = self.length()
        k0 = self.local_stiffness_matrix(E, A, I1, Le)

        # k0 = CheckReleases[fem.elem[N],k0];

        # local-global matrix
        L = self.transformation_matrix()

        # globaali elementin matriisi
        ke = L.transpose().dot(k0.dot(L))

        return ke

    def local_geometric_stiffness_matrix(self,lcase=0):
        """ Geometric stiffness matrix in local coordinates
            source: Cook et. al 1989, Section 14.2

            Returns:
            ---------
            :return: element's geometric stiffness matrix in local coordinates
            :rtype: np.array

        """
        #P = self.axial_force[1]
        P = self.fint[lcase]['fx'][1]
        Le = self.length()

        g0 = P / 30 / Le
        keg0 = g0 * np.array([[36, 3 * Le, -36, 3 * Le],
                              [3 * Le, 4 * Le ** 2, -3 * Le, -Le ** 2],
                              [-36, -3 * Le, 36, -3 * Le],
                              [3 * Le, -Le ** 2, -3 * Le, 4 * Le ** 2]])
        kg0 = np.zeros([6, 6])
        q = [1, 2, 4, 5]
        kg0[np.ix_(q, q)] = keg0

        return kg0

    def geometric_stiffness_matrix(self,lcase=0):
        """ Geometric stiffness matrix in global coordinates \n
            From: Cook et. al (1989), Section 14.2

            Returns:
            ---------
            :return: element's geometric stiffness matrix in global coordinates
            :rtype: np.array

        """

        # local-global transformation matrix
        L = self.transformation_matrix()

        kG0 = self.local_geometric_stiffness_matrix(lcase)

        kG = L.transpose().dot(kG0.dot(L))

        return kG

    def equivalent_nodal_loads(self, load):
        """ Equivalent nodal loads for a load 
        
            'load' is a Load type object. By default, it is assumed that
            'load' is LineLoad

            Returns:
            ---------
            :return: Equivalent nodal loads for load in vector q
            :rtype: np.array
        """

        floc = np.zeros(6)

        T = self.transformation_matrix()
        L = self.length()

        if isinstance(load,LineLoad):
            q1 = load.qval[0]
            q2 = load.qval[1]
            #print(q1,q2)
            if load.coords == 'local':
                #print('local coordinates')
                """ Load is given in local coordinates """
                if load.dir == "x":
                    """ Load is in the axial direction """
                elif load.dir == 'y':
                    """ Load is perpendicular to the member """
                    floc[1] = 7/20*q1*L + 3/20*q2*L
                    floc[4] = 3/20*q1*L + 7/20*q2*L                    
                    floc[2] = L**2*(1/20*q1+1/30*q2)
                    floc[5] = -L**2*(1/30*q1+1/20*q2)
                #print(floc)
            else:
                """ Load is given in global coordinates: it has to be transformed
                    first into local coordinates.
                """
                if load.dir == "x":
                    #q = np.array([q1, 0, 0])
                    q = np.array([1, 0, 0])
                else:
                    #q = np.array([0, q1, 0])
                    q = np.array([0, 1, 0])
                # Load vector transformed into element local coordinate system
                qloc = T[:3, :3].dot(q)
                
                """ Loads perpendicular to the axis of the element """
                qy1 = q1*qloc[1]
                qy2 = q2*qloc[1]
                floc[1] = 7/20*qy1*L + 3/20*qy2*L
                floc[4] = 3/20*qy1*L + 7/20*qy2*L
                floc[2] = L**2*(1/20*qy1+1/30*qy2)
                floc[5] = -L**2*(1/30*qy1+1/20*qy2)
                
                """ Loads parallel to the axis of the element """
                qx1 = q1*qloc[0]
                qx2 = q2*qloc[0]
                floc[0] = 0.5*L*(0.5*(qx1+qx2)-1/6*(qx2-qx1))
                floc[3] = 0.5*L*(0.5*(qx1+qx2)+1/6*(qx2-qx1))
                
                # Construct nodal load in local coordinates
                # qloc[0,1,2] = axial, shear, moment of node 1
                # qloc[3,4,5] = axial, shear, moment ofnode 2
                # Axial force
                #floc[[0, 3]] = 0.5 * qloc[0] * L
                # Shear force
                #floc[[1, 4]] = 0.5 * qloc[1] * L
                # Moment
                #floc[2] = qloc[1] * L ** 2 / 12.0
                #floc[5] = -floc[2]

                #print(floc)

        if len(self.releases) > 0:
            """ There are releases in the element, so the local load vector
                is modified.
            """
            rel = self.releases
            # nrel is the list of non-released forces
            nrel = np.setdiff1d(np.arange(6),rel)            
            
            #print(np.linalg.inv(self.Krel['K22']).dot(floc[rel]))
            floc[nrel] -= self.Krel['K12'].dot(np.linalg.inv(self.Krel['K22']).dot(floc[rel]))
            floc[rel] = 0            

        # Store local loads
        self.floc[load.sid] = floc

        # print(self.floc)

        fglob = T.transpose().dot(floc)
        
        
        return fglob

    def shape_fun(self, x):
        """ Evaluate Hermitean shape functions at local coordinate x

            Returns:
            ---------
            :return: Hermitean shape functions at local coordinate x
            :rtype: np.array


        """
        L = self.length()
        N = np.zeros(4)
        N[0] = 1 - 3 * x ** 2 / L ** 2 + 2 * x ** 3 / L ** 3
        N[1] = x - 2 * x ** 2 / L + x ** 3 / L ** 2
        N[2] = 3 * x ** 2 / L ** 2 - 2 * x ** 3 / L ** 3
        N[3] = -x ** 2 / L ** 2 + x ** 3 / L ** 2

        return N

    def internal_forces(self,lcase=0):
        """ Calculate internal forces
            NOTE: these internal forces do not take
            loads along the element into account!

            Works only for a single load case!
        """

        """ Get nodal displacements in local coordinates
            and multiply them with local element stiffness matrix
            to get internal forces in member's local coordinate system.
        """
        q = self.local_displacements(lcase)
        # ke = self.k0
        E = self.material.young
        A = self.section.A
        I1 = self.section.I[0]
        Le = self.length()
        ke = self.local_stiffness_matrix(E, A, I1, Le)

        try:            
            R = ke.dot(q) - self.floc[lcase]
        except:
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

        """
        self.axial_force[:2] = [-R[0], R[3]]
        self.bending_moment[:2] = [-R[2], R[5]]
        self.shear_force[:2] = [R[1], -R[4]]
        """
        self.fint[lcase]['fx'] = [-R[0], R[3]]
        self.fint[lcase]['my'] = [-R[2], R[5]]
        self.fint[lcase]['fz'] = [R[1], -R[4]]
        
class EBBeam3D(Element):
    """ Euler-Bernoulli beam element in 3D

        Parameters:
        -----------
        :param n1: start node
        :param n2: end node
        :param nref: reference node for orienting the element
        :param section: element's section properties
        :param material: element's material properties

        :type n1: FEMNode
        :type n2: FEMNode
        :type section: Section
        :type material: Material

        Variables:
        ----------
        :ivar bending_moment: bending moment values at element's ends
        :ivar shear_force: shear force values at element's ends

        :vartype bending_moment: list
        :vartype shear_force: list
        
        Degrees of freedom
        [0] .. ux1
        [1] .. uy1
        [2] .. uz1
        [3] .. rx1
        [4] .. ry1
        [5] .. rz1
        
        [6] .. ux2
        [7] .. uy2
        [8] .. uz2
        [9] .. rx2
        [10] .. ry2
        [11] .. rz2
        
        TODO:
            equivalent_nodal_loads
            internal_forces
            test all functions
    """
    
    def __init__(self, n1, n2, section, material, nref = None):
        Element.__init__(self, n1, n2, section, material)

        self.dim = 3
        
        """ Reference node for orienting the member
            This node must not be on the line passing through nodes n1 and n2.
        """
        if nref is None:
            """ Calculate third point separately 
                The 2D case is in xz-plane, so by default, the y
                axis is the out-of-plane axis. The member is oriented
                such that in the plane case, the local and global y axes
                coincide.
            """
            if np.linalg.norm(np.cross(self.direction_vector(),[0,1,0])) < 1e-6:
            #if abs(n2.coord[1]-n1.coord[1]) < 1e-6:
                """ local x axis coincides with the global
                    y axis, so the reference point is obtained by stepping
                    in the z direction
                """
                nref = FEMNode(-1,n2.coord[0],n2.coord[1],n2.coord[2]+1)
            elif np.linalg.norm(np.cross(self.direction_vector(),[0,0,1])) < 1e-6:
                """ member is vertial, so local x axis coincides with
                    z-axis                
                """
                nref = FEMNode(-1,n2.coord[0]+1,n2.coord[1],n2.coord[2])
            else:
                nref = FEMNode(-1,n2.coord[0],n2.coord[1],n2.coord[2]+1)
        
        self.nref = nref
        
    def init_dofs(self):
        """ For beam members all degrees of freedom
            are in play.
            
        """
        for n in self.nodes:
            n.dofs[:6] = 0
    
    def direction_vector(self):
        
        return super().direction_cosines()
            
    def direction_cosines(self):
        """ Calculates direction cosines for a 3D-beam element """
        
        l0 = np.zeros((3,3))
        
        X = self.coord()
        L = self.length()
        dX = X[1, :] - X[0, :]        
        # direction cosines along the local x axis
        l0[0,:] = dX / L
        """
        l = l0[0,0]
        m = l0[0,1]
        n = l0[0,2]
        D = np.sqrt(l**2+n**2)
        l0[2,:] = 1/D*np.array([n,0,-l])
        l0[1,:] = 1/D*np.array([-m*n,n*l+n**2,-m*n])
        """
        # calculate direction cosines along the local y axis
        # coordinates of the reference node
        X3 = self.nref.coord
        dX13 = X3-X[0,:]
        l13 = np.linalg.norm(dX13)
        
        # unit vector between point 1 and reference point
        V13 = dX13/l13
        VZ0 = np.cross(V13,l0[0,:])        
        l0[1,:] = VZ0/np.linalg.norm(VZ0)
        
        l0[2,:] = np.cross(l0[0,:],l0[1,:])
        
        return l0
        
    def transformation_matrix(self):
        """ Calculates transformation matrix from local to global coordinates

            Returns:
            --------
            :return: element's transformation matrix
            :rtype: np.array

         """
        l0 = self.direction_cosines()
        L = np.zeros((12, 12))
        L[0:3, 0:3] = l0
        L[3:6, 3:6] = l0
        L[6:9, 6:9] = l0
        L[9:12, 9:12] = l0
        return L

    def global_dofs(self):
        """ Get numbering of element's global degrees of freedom

            Returns:
            --------
            :return: list of end nodes dofs added together
            :rtype: list
        """
        return np.append(self.nodes[0].dofs,self.nodes[1].dofs)
        #return self.nodes[0].dofs + self.nodes[1].dofs
    
    
    def local_stiffness_matrix(self):
        """ Generates stiffness matrix in local coordinates


            Returns:
            --------
            :return: element's stiffness matrix in local coordinates
            :rtype: np.array

        """
        E = self.material.young
        G = self.material.shear_modulus
        A = self.section.A
        I1 = self.section.Iy
        I2 = self.section.Iz
        J = self.section.It
        
        Le = self.length()

        rodc = E * A / Le
        EI1 = E * I1
        EI2 = E * I2
        torc = G*J/Le
        k0 = np.zeros((12, 12))
        # Axial stiffness
        k0[[0, 6], [0, 6]] = rodc
        k0[[0, 6], [6, 0]] = -rodc
        # Torsional stiffness
        k0[[3, 9], [3, 9]] = torc
        k0[[3, 9], [9, 3]] = -torc
        # Bending with respect to y axis
        bend1y = 12 * EI2 / Le ** 3
        k0[[1, 7], [1, 7]] = bend1y
        k0[[1, 7], [7, 1]] = -bend1y
        #
        bend2y = 6 * EI2 / Le ** 2
        k0[[1, 1, 5, 11], [5, 11, 1, 1]] = bend2y
        #
        k0[[5, 7, 7, 11], [7, 5, 11, 7]] = -bend2y
        #
        k0[[5, 11], [11, 5]] = 2 * EI2 / Le
        #
        k0[[5, 11], [5, 11]] = 4 * EI2 / Le
        
        bend1z = 12 * EI1 / Le**3
        k0[[2, 8], [2, 8]] = bend1z
        k0[[8, 2], [2, 8]] = -bend1z

        bend2z = 6 * EI1 / Le ** 2
        k0[[8, 4, 8, 10], [4, 8, 10, 8]] = bend2z
        #
        k0[[2, 10, 2, 4], [10, 2, 4, 2]] = -bend2z

        k0[[4, 10], [10, 4]] = 2 * EI1 / Le
        #
        k0[[4, 10], [4, 10]] = 4 * EI1 / Le

        return k0
    
    def stiffness_matrix(self):
        """ Computes element's global stiffness matrix

            Returns:
            ---------
            :return: element's stiffness matrix in global coordinates
            :rtype: np.array

        """
        k0 = self.local_stiffness_matrix()

        # k0 = CheckReleases[fem.elem[N],k0];

        # local-global matrix
        L = self.transformation_matrix()

        # globaali elementin matriisi
        ke = L.transpose().dot(k0.dot(L))

        return ke
    
    def equivalent_nodal_loads(self, load):
        """ Equivalent nodal loads for a load 
        
            'load' is a Load type object. By default, it is assumed that
            'load' is LineLoad

            Returns:
            ---------
            :return: Equivalent nodal loads for load in vector q
            :rtype: np.array
        """

        """ Equivalent nodal loads in local coordinates.
            Six degrees of freedom per node.
            
            Degrees of freedom
            [0] .. ux1
            [1] .. uy1
            [2] .. uz1
            [3] .. rx1
            [4] .. ry1
            [5] .. rz1
            
            [6] .. ux2
            [7] .. uy2
            [8] .. uz2
            [9] .. rx2
            [10] .. ry2
            [11] .. rz2
        """
        floc = np.zeros(12)

        T = self.transformation_matrix()
        L = self.length()

        if isinstance(load,LineLoad):
            q1 = load.qval[0]
            q2 = load.qval[1]
            #print(q1,q2)
            if load.coords == 'local':
                #print('local coordinates')
                """ Load is given in local coordinates """
                if load.dir == "x":
                    """ Load is in the axial direction """
                elif load.dir == 'y':
                    """ Load is perpendicular to the member """
                    floc[1] = 7/20*q1*L + 3/20*q2*L
                    floc[7] = 3/20*q1*L + 7/20*q2*L                    
                    floc[5] = L**2*(1/20*q1+1/30*q2)
                    floc[11] = -L**2*(1/30*q1+1/20*q2)
                elif load.dir == 'z':
                    floc[2] = 7/20*q1*L + 3/20*q2*L
                    floc[8] = 3/20*q1*L + 7/20*q2*L                    
                    floc[4] = L**2*(1/20*q1+1/30*q2)
                    floc[10] = -L**2*(1/30*q1+1/20*q2)
                #print(floc)
            else:
                """ Load is given in global coordinates: it has to be transformed
                    first into local coordinates.
                """
                if load.dir == "x":                    
                    q = np.array([1, 0, 0,0,0,0])
                elif load.dir == 'y':                    
                    q = np.array([0, 1, 0,0,0,0])
                elif load.dir == 'z':                    
                    q = np.array([0, 0, 1,0,0,0])
                # Load vector transformed into element local coordinate system
                qloc = T[:3, :3].dot(q)
                
                """ Loads perpendicular to the axis of the element """
                qy1 = q1*qloc[1]
                qy2 = q2*qloc[1]
                floc[1] = 7/20*qy1*L + 3/20*qy2*L
                floc[4] = 3/20*qy1*L + 7/20*qy2*L
                floc[2] = L**2*(1/20*qy1+1/30*qy2)
                floc[5] = -L**2*(1/30*qy1+1/20*qy2)
                
                """ Loads parallel to the axis of the element """
                qx1 = q1*qloc[0]
                qx2 = q2*qloc[0]
                floc[0] = 0.5*L*(0.5*(qx1+qx2)-1/6*(qx2-qx1))
                floc[3] = 0.5*L*(0.5*(qx1+qx2)+1/6*(qx2-qx1))
                
                # Construct nodal load in local coordinates
                # qloc[0,1,2] = axial, shear, moment of node 1
                # qloc[3,4,5] = axial, shear, moment ofnode 2
                # Axial force
                #floc[[0, 3]] = 0.5 * qloc[0] * L
                # Shear force
                #floc[[1, 4]] = 0.5 * qloc[1] * L
                # Moment
                #floc[2] = qloc[1] * L ** 2 / 12.0
                #floc[5] = -floc[2]

                #print(floc)

        if len(self.releases) > 0:
            """ There are releases in the element, so the local load vector
                is modified.
            """
            rel = self.releases
            # nrel is the list of non-released forces
            nrel = np.setdiff1d(np.arange(12,rel))          
            
            #print(np.linalg.inv(self.Krel['K22']).dot(floc[rel]))
            floc[nrel] -= self.Krel['K12'].dot(np.linalg.inv(self.Krel['K22']).dot(floc[rel]))
            floc[rel] = 0            

        # Store local loads
        self.floc[load.sid] = floc

        # print(self.floc)

        fglob = T.transpose().dot(floc)
        
        
        return fglob

    def local_geometric_stiffness_matrix(self,lcase=0):
        """ Geometric stiffness matrix in local coordinates
            source: Cook et. al 1989, Section 14.2

            Returns:
            ---------
            :return: element's geometric stiffness matrix in local coordinates
            :rtype: np.array

        """
        #P = self.axial_force[1]
        P = self.fint[lcase]['fx'][1]
        Le = self.length()

        g0 = P / 30 / Le
        keg0 = g0 * np.array([[36, 3 * Le, -36, 3 * Le],
                              [3 * Le, 4 * Le ** 2, -3 * Le, -Le ** 2],
                              [-36, -3 * Le, 36, -3 * Le],
                              [3 * Le, -Le ** 2, -3 * Le, 4 * Le ** 2]])
        kg0 = np.zeros([6, 6])
        q = [1, 2, 4, 5]
        kg0[np.ix_(q, q)] = keg0

        return kg0

    def internal_forces(self,lcase=0):
            """ Calculate internal forces
    
            """
    
            """ Get nodal displacements in local coordinates
                and multiply them with local element stiffness matrix
                to get internal forces in member's local coordinate system.
            """
            q = self.local_displacements(lcase)
            ke = self.local_stiffness_matrix()
    
            print('q=')
            print(q)
    
            try:            
                R = ke.dot(q) - self.floc[lcase]
            except:
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
    
            """
            self.axial_force[:2] = [-R[0], R[3]]
            self.bending_moment[:2] = [-R[2], R[5]]
            self.shear_force[:2] = [R[1], -R[4]]
            Degrees of freedom
                [0] .. ux1
                [1] .. uy1
                [2] .. uz1
                [3] .. rx1
                [4] .. ry1
                [5] .. rz1
                
                [6] .. ux2
                [7] .. uy2
                [8] .. uz2
                [9] .. rx2
                [10] .. ry2
                [11] .. rz2
            """
            self.fint[lcase]['fx'] = [-R[0], R[6]]
            self.fint[lcase]['fy'] = [-R[1], R[7]]
            self.fint[lcase]['fz'] = [R[2], -R[8]]
            self.fint[lcase]['mx'] = [-R[3], R[9]]
            self.fint[lcase]['my'] = [-R[4], R[10]]
            self.fint[lcase]['mz'] = [R[5], -R[11]]
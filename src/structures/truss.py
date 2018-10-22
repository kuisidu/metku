""" General class for truss structures

    Classes:
        Truss -- general truss structure
        TrussMember -- general truss member
        TrussNode -- general truss node
   
    Written by Kristo Mela 18.4.2017
"""
import math
import numpy as np
import matplotlib.pyplot as plt

class Truss:
    """ General truss structure
        Attributes:
            members -- list of truss members
            nodes -- list of truss nodes
            ID -- table for numbering degrees of freedom
            loads -- array of nodal loads
            loadcases -- array of load case data
            
        Methods:
            add_member -- adds a new member
            add_node -- adds a new node
            nnodes -- get number of nodes
            nmembers -- get number of members
            nlc -- get number of load cases
            member_coord -- get nodal coordinates of a member
            member_length -- compute the length of a member
            member_direction_vector -- determine the direction vector of
                                       a member
            common_node -- determine the common node shared by two members
            angle_between_members -- determine the angle between two members
            members_at_node -- find members connected to a node
        """
    
    
    def __init__(self):
        """ 
        
        x -- 3 x nN matrix of nodal coordinates
        members   % array of truss members
                  % each element is a struct with the following fields:
                  % .n -- 2 x 1 vector with node indices of the members
                  % .profile -- member profile
        ID        % ID table of degrees of freedom
        loads     % array of nodal loads
        loadcases % array of load case data
        """
                
        self.members = {}
        self.nodes = {}
        self.ID = []
        self.load = []
        self.loadcases = []
        self.dim = 2
    
    def add_member(self,n1,n2,profile,label="Member"):
        """ Adds a new member to the truss
            input: 
                n1 .. end node 1
                n2 .. end node 2
                the nodes are "TrussNode" objects
                   
        """
        #new_member = TrussMember([n1,n2],p,type)
        new_member = TrussMember([self.nodes[n1],self.nodes[n2]],profile,label)
        self.members.append(new_member)
    
    def add_node(self,x,y,z=0,label="Node"):
        """ Adds a node to the truss """
        
        """ 
            The dimension of the truss (2D or 3D)
            is determined by the number of coordinates
        """
        if z != 0:
            self.dim = 3
            
        new_node = TrussNode(x,y,z,label)
        self.nodes.append(new_node)
    
    def nnodes(self):
        """ Number of nodes """
        return len(self.nodes)
    
    def nmembers(self):
        """ Number of members """
        return len(self.members)
    
    def nlc(self):
        """ Number of load cases """
        return len(self.loadcases)

    def common_node(self,m1,m2):
        """ Find the common node of members m1 and m2 
            Input:
                m1 .. index of member 1
                m2 .. index of member 2        
        """
        n1 = self.members(m1).nodes
        n2 = self.members(m2).nodes
        if n1[0] in n2:
            ncommon = n1[0]
        else:
            ncommon = n1[1]
            
        return ncommon

    def angle_between_members(self,m1,m2):
        """ Angle (in degrees) between two members 
            Input:
                m1 .. index of member 1
                m2 .. index of member 2        
        """
        
        # it is assumed that the members share a node
        N1 = self.members[m1].nodes
        N2 = self.members[m2].nodes
        
        # find common node
        if N1[0] in N2:
            Ncommon = N1[0]
            # get the nodes of each member that are not shared
            N1unique = N1[1]
            if N2[0] == Ncommon:
                N2unique = N2[1]
            else:
                N2unique = N2[0]
                            
        elif N1[1] in N2:
            Ncommon = N1[1]
            N1unique = N1[0]
            if N2[0] == Ncommon:
                N2unique = N2[1]
            else:
                N2unique = N2[0]
        else:
            Ncommon = None
                                            
        # coordinates of the common node
        Xcommon = self.nodes[Ncommon].x
        # Coordinates of the unique nodes
        X1 = self.nodes[N1unique].x
        X2 = self.nodes[N2unique].x
                    
        v1 = X1-Xcommon
        r1 = v1/np.linalg.norm(v1)
        v2 = X2-Xcommon
        r2 = v2/np.linalg.norm(v2)

        cd = np.dot(r1,r2)
        a = math.degrees(math.acos(cd))
        
        return a
    
    def members_at_node(self,n):
        """ Find members connected to node 'n' """
        nid = id(self.nodes[n])
    
        mems = []
    
        for i in range(self.nmembers()):
            if id(self.members[i].nodes[0]) == nid or \
                id(self.members[i].nodes[1]) == nid:
                mems.append(i)

        return mems
    
    def draw(self):
        """ draw nodes """
        for n in self.nodes:
            plt.plot(n.x[0],n.x[1],'ro')
        
        for i in range(self.nnodes()):
            plt.text(self.nodes[i].x[0],self.nodes[i].x[1],str(i))
        
        """ draw members """
        for i in range(self.nmembers()):
            #X = self.member_coord(i)
            X = self.members[i].coord()
            plt.plot(X[:,0],X[:,1],'b')
            Xmid = X[0,:]+0.5*(X[1,:]-X[0,:])
            plt.text(Xmid[0],Xmid[1],str(i))
        #plt.axis([0, 6, 0, 20])
        
        plt.show()
    
class TrussMember:
    """ Class for truss members.
    
        Basically only for storing data
    """
    
    def __init__(self,n,profile,label="Member"):
        """ Truss member
        
            n -- array of TrussNode objects
            profile -- cross-section of the member
            label -- member label
        """
        self.nodes = n
        self.profile = profile
        self.label = label

    def coord(self):
        """ Nodal coordinates of member m
        
            m -- member index
            xorig -- True: use original nodal locations
                     False: use modified nodal locastions
        """
        X = np.array([self.nodes[0].x,self.nodes[1].x])
        return X

    def length(self):
        X = self.coord()
        L = np.linalg.norm(X[1,:]-X[0,:])
        return L

    def direction_vector(self):
        """ Unit direction vector of member
        
        xorig -- True: use original nodal locations
        False: use modified nodal locastions
        """
        X = self.coord()
        v0 = X[1,:]-X[0,:]
        v = v0/np.linalg.norm(v0)
        return v

class TrussNode:
    """ Class for truss nodes.
    
        Basically only for storing data
    """
    
    def __init__(self,x,y,z=0,label="Node"):
        """ input: x, y and z are the coordinates of the node 
            label .. unique nodel label
        
        """
        self.x = np.array([x,y,z])
        self.label = label
        #self.xmod = np.array([x,y,z])

""" Testing function """
if __name__ == "__main__":
    import sys
    
    from hollow_sections import SHS
    
    t = Truss()
    t.add_node(0,0)
    t.add_node(1,1)
    t.add_node(2,0)
    t.add_node(0,1)
    t.add_node(3,1)
    
    p = SHS(100,5)
    
    t.add_member(0,1,p)
    t.add_member(0,2,p)
    t.add_member(3,1,p)
    t.add_member(2,1,p)
    t.add_member(2,4,p)
    t.add_member(1,4,p)
    
    for m in t.members:
        print(m.coord())
        print(m.length())
    
    print(t.members_at_node(2))
    
    t.draw()

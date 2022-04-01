# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 20:25:02 2021

Class for simple (pin-jointed) trusses. Primarily intended for optimization

@author: kmela
"""

import numpy as np
from timeit import default_timer as timer

from raami import Raami
from frame_node import FrameNode
from frame_member import FrameMember, SteelFrameMember
from frame_loads import PointLoad
from frame_supports import XYHingedSupport

class SimpleTruss(Raami):
    """ Class for simple trusses. Essentially, the structure is a collection of
        pin-jointed members, loaded at the nodes. Consequently, the members
        carry only axial forces.
        
        This class is mainly used for truss optimization problems.
        
        Methods:
        ------------
        add
        generate
        calculate
        plot
    """
    
    def __init__(self,nodes=None,members=None,profiles=None,name="Raami Truss"):
        """ Constructor
        

        Parameters:
        ------------
        :param nodes: list of nodal coordinates [[x1,y1],[x2,y2],...,[xn,yn]]
        :param members: list of member nodes [[e11,e12],[e21,e22],...,[em1,em2]]
            where eij is a node index in the list 'nodes'
        :param profiles: list of member profiles, corresponding to the members

        :type nodes: list of lists
        :type members: list of lists
        :type profile: list


        Variables:
        ----------
        :ivar nodes:
        
        """
        super().__init__(name)
    
        self._nodes = nodes
        self._members = members
        self._profiles = profiles
        
        """ Create truss nodes """
        for node in nodes:
            self.add(FrameNode(node))
        
        """ Create truss members """
        for mem, prof in zip(members,profiles):
            mem_nodes = [self.nodes[mem[0]],self.nodes[mem[1]]]
            self.add(SimpleTrussMember(mem_nodes,prof))
    
    def add(self, this):
        """
        Adds object to truss
        :param this: object to be added
        :return:
        """
        # Simple truss member
        if isinstance(this, SimpleTrussMember):
            # Give member a unique id
            # id is used in creating elements with wanted cross-sectional properties
            this.mem_id = int(len(self.members))

            # Create a link between the FrameMember object and the frame
            this.frame = self

            # Generate member coordinates
            #this.calc_nodal_coordinates()

            self.members[this.mem_id] = this
        else:
            super().add(this)
    #def generate(self):
    #    """ Generates the finite element model """ 
        
    
    
class SimpleTrussMember(SteelFrameMember):
    """ Class for simple truss members. These are a subclass of FrameMember,
        but with less functionalities, because only axial forces are carried.
        
        On the other hand, the profile can be generic, i.e. only 'A' and optionally 'I'
        can be given.
    """
    
    def __init__(self, nodes, section, mtype='bar', mem_id="", nel=1):

        super().__init__(nodes, section, mtype, mem_id, nel)

def three_bar_truss(L=1000,F1=100,F2=100,nlc=1):
    """ three bar truss example """
    from sections.steel.RHS import RHS
    #from metku.sections.cross_sections import CrossSection

    # Nodal coordinates    
    X = [[0,0],[-L,L],[0,L],[np.sqrt(3)*L,L]]
    
    """ Member connectivities. For example, [0,1] means that the member
        end nodes are the 0th and 1st nodes in the list.
    """
    mems = [[0,1],[0,2],[0,3]]
    
    #profs = [SHS(50,5),SHS(50,5),SHS(50,5)]
    """ Cross-sections
        Here, a name for the cross-section has to be given. This has to do
        with the way FrameMember has been implemented.
        
        Later, the profile data can be changed manually, e.g. by modifying
        the cross-sectional area 'A' directly.
    """
    profs = [RHS(50,50,3),RHS(50,50,3),RHS(50,50,3)]
    #profs = [CrossSection(A=100),CrossSection(A=100),CrossSection(A=100)]

    # Generate truss
    t = SimpleTruss(X,mems,profs)
    
    # Add hinged supports
    t.add(XYHingedSupport(X[1]))
    t.add(XYHingedSupport(X[2]))
    t.add(XYHingedSupport(X[3]))
    # Add point load
    # Voima F1 on pystysuuntainen ja F2 vaakasuuntainen
    if nlc == 1:
        p = PointLoad(X[0],[F2,F1,0])
        t.add(p)
    else:
        p1 = PointLoad(X[0],[F2,0,0],load_id=2)
        t.add(p1)
        p2 = PointLoad(X[0],[0,F1,0],load_id=3)
        t.add(p2)
    
    #t.plot()
    
    # Generate calculation model
    t.generate()
    # Calculate the responses
    t.calculate(support_method='REM')
    
    return t
        
def ten_bar_truss(L,F):
    """ Classic 10-bar truss """
    from sections.steel.RHS import RHS
    
    nodes = []
    nodes.append([2 * L, L])
    nodes.append([2 * L, 0])
    nodes.append([L, L])
    nodes.append([L, 0])
    nodes.append([0, 0])
    nodes.append([0, L])
    
    mems = [[4,2],[2,0],[5,3],[3,1],[2,3],[0,1],[4,3],[5,2],[2,1],[3,0]]
    
    profs = []
    for i in range(10):
        profs.append(RHS(50,50,3))
    
    t = SimpleTruss(nodes,mems,profs,"Ten Bar Truss")
    
    t.add(PointLoad(t.nodes[1],[0,-F,0]))
    t.add(PointLoad(t.nodes[3],[0,-F,0]))
    
    t.add(XYHingedSupport(t.nodes[4]))
    t.add(XYHingedSupport(t.nodes[5]))
    
    #t.generate()
    #t.calculate()
    
    #t.plot()
    return t

def fifteen_bar_truss(L=360,h=120,P=10000):
    """ by default, the units are imperial """
    
    from RHS import RHS
    
    nodes = []
    nodes.append([0, h])
    nodes.append([L/3, h])
    nodes.append([2*L/3, h])
    nodes.append([L, h])
    nodes.append([0, 0])
    nodes.append([L/3, 0])
    nodes.append([2*L/3, 0])
    nodes.append([L, 0])
    
    mems = [[0,1],[1,2],[2,3],[4,5],[5,6],[6,7],
            [5,1],[6,2],[7,3],
            [0,5],[4,1],
            [1,6],[5,2],
            [2,7],[6,3]]
    
    profs = []
    for i in range(len(mems)):
        profs.append(RHS(50,50,3))
    
    t = SimpleTruss(nodes,mems,profs)
        
    t.add(PointLoad(nodes[7],[0,-P,0]))
            
    t.add(XYHingedSupport(nodes[0]))
    t.add(XYHingedSupport(nodes[4]))
    
    t.generate()
    
    #print(type(t.point_loads['PointLoad'].coordinate))
    
    t.calculate(load_id='all',support_method="REM")
    
    t.plot()
    
    return t
    
if __name__ == '__main__':
    
    import cProfile    
    #t = three_bar_truss(L=3000,F1=-200e3,F2=-250e3,nlc=1)
    t = ten_bar_truss(L=3000,F=200e3)
    t.plot()
    #t = fifteen_bar_truss()
    
    
    t.generate_fem()
    t.structural_analysis(load_id=t.load_ids[0],support_method="REM")
    #t.generate()
    #load_id = t.load_ids[0]
    #start = timer()    
    #t.f.linear_statics(lcase=load_id,support_method="REM")
    #t.calculate(load_id="all",support_method="REM")
    #end = timer()
    #print(end-start)
    #cProfile.run('t.f.linear_statics(lcase=load_id,support_method="REM")')
    #cProfile.run('t.calculate(load_id="all",support_method="REM")')
        
       
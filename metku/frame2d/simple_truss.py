# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 15:02:50 2020

@author: kmela
"""

try:
    from metku.frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, \
        Support, Hinge, XYHingedSupport, PREC
    from metku.framefem.elements import EBBeam, EBSemiRigidBeam
    # from fem.elements.eb_semi_rigid_beam import EBSemiRigidBeam
    from metku.sections.steel import HEA

    from metku.framefem import FrameFEM, BeamSection
    from metku.structures.steel.steel_member import SteelMember
    from metku.eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint, \
        RHSYJoint
except:
    from frame2d.frame2d import Frame2D, FrameMember, PointLoad, LineLoad, \
        Support, Hinge, XYHingedSupport, PREC
    from framefem.elements import EBBeam, EBSemiRigidBeam
    # from fem.elements.eb_semi_rigid_beam import EBSemiRigidBeam
    from sections.steel import HEA
    from sections.steel.RHS import SHS

    from framefem import FrameFEM, BeamSection
    from structures.steel.steel_member import SteelMember
    from eurocodes.en1993.en1993_1_8.rhs_joints import RHSKGapJoint, RHSYJoint
    
import numpy as np

class SimpleTruss(Frame2D):
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
    
    def __init__(self,nodes=None,members=None,profiles=None):
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
        super().__init__(num_elements=1)
    
        self._nodes = nodes
        self._members = members
        self._profiles = profiles
        
        """ Create truss members """
        for mem, prof in zip(members,profiles):
            mem_nodes = [nodes[mem[0]],nodes[mem[1]]]
            self.add(SimpleTrussMember(mem_nodes,profile=prof))
    
    def add(self, this):
        """
        Adds object to truss
        :param this: object to be added
        :return:
        """
        # Simple truss member
        if isinstance(this, SimpleTrussMember):
            this.mem_id = int(len(self.members))
            self.members[this.mem_id] = this
            this.calc_nodal_coordinates(self.num_elements)            
        else:
            super().add(this)
    #def generate(self):
    #    """ Generates the finite element model """ 
        
    
    
class SimpleTrussMember(FrameMember):
    """ Class for simple truss members. These are a subclass of FrameMember,
        but with less functionalities, because only axial forces are carried.
        
        On the other hand, the profile can be generic, i.e. only 'A' and optionally 'I'
        can be given.
    """
    
    def __init__(self,
                 coordinates,
                 mem_id="",
                 profile="SHS 100x100x5.0",
                 material="S355",
                 mtype='bar',
                 num_elements=1):

        super().__init__(coordinates, mem_id, profile, material, num_elements,mtype=mtype)

def three_bar_truss(L=1000):
    """ three bar truss example """
    
    #from metku.sections.cross_sections import CrossSection

    # Nodal coordinates    
    X = [[0,0],[-L,L],[0,L],[round(np.sqrt(3)*L,PREC),L]]
    
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
    profs = ['SHS 50x50x5.0','SHS 50x50x5.0','SHS 50x50x5.0']
    #profs = [CrossSection(A=100),CrossSection(A=100),CrossSection(A=100)]

    # Generate truss
    t = SimpleTruss(X,mems,profs)
    
    # Add hinged supports
    t.add(XYHingedSupport(X[1]))
    t.add(XYHingedSupport(X[2]))
    t.add(XYHingedSupport(X[3]))
    # Add point load
    p = PointLoad(X[0],[-100,-100])
    t.add(PointLoad(X[0],[-100,-100,0]))
    
    # Generate calculation model
    t.generate()
    # Calculate the responses
    t.calculate()
    
    t.plot()
    
    return t
        
if __name__ == '__main__':
    t = three_bar_truss()
        
        
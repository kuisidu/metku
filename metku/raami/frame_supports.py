# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 15:03:22 2021

Frame Supports

@author: kmela
"""

# --------------------------- SUPPORT CLASSES ----------------------------
class Support:
    def __init__(self, node, dofs, supp_id=1):
        """
        

        Parameters
        ----------
        node : FrameNode
            Node to be supported.
        dofs : List
            List of supported degrees of freedom.
        supp_id : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        None.

        """
        """ Constructor
            coordinate .. nodal coordinates [list]
            dofs .. list of supported degrees of freedom
            supp_id .. each support has an ID that relates it to a given
                    support system. NOTE! In 'calculate' method, supp_id = 1 is always assumed!
        """
        self.node = node
        self.supp_id = supp_id
        self.dofs = dofs

    def add_support(self, fem_model):
        """ Adds support to the given FEM model """        
        try: 
            idx = fem_model.nodes.index(self.node.fem_node)    
            fem_model.add_support(self.supp_id, idx, self.dofs, val=0.0)
        except ValueError:
            print("Setting of support not successful.")


class FixedSupport(Support):
    '''
    In ROBOT this setting means that UX, UZ and RY check marks are present in the support dialog
    '''
    def __init__(self, node, supp_id=1):
        # super().__init__(coordinate, [1, 1, 1], supp_id)
        super().__init__(node, [0, 1, 2], supp_id)


class XHingedSupport(Support):
    '''
    In ROBOT this setting means that UX direction check marks is present in the support dialog, but UZ and RY are empty
    '''
    def __init__(self, node, supp_id=1):
        super().__init__(node, [0], supp_id)


class YHingedSupport(Support):
    '''
    In ROBOT this setting means that UZ direction check mark is present in the support dialog, but UX and RY are empty
    '''
    def __init__(self, node, supp_id=1):
        super().__init__(node, [1], supp_id)


class XYHingedSupport(Support):
    '''
    In ROBOT this setting means that UX and UZ direction check marks are present in the support dialog, but RY is empty
    '''
    def __init__(self, node, supp_id=1):
        super().__init__(node, [0, 1], supp_id)

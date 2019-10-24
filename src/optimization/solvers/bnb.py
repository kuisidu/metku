# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 19:29:05 2019

Branch and bound algorithm for optimization

Data structures:
    
    Class BnB .. for running the algorithm
    Class BnBNode .. contains information about branching at given node


@author: kmela
"""

from copy import deepcopy
import numpy as np
from treelib import Tree, Node

try:
    from src.optimization.solvers.optsolver import OptSolver
except:    
    from optimization.solvers.optsolver import OptSolver


CON_TOL = 1e-4

# A class that represents an individual node in a 
# Binary Tree 
class BnBNode: 
    def __init__(self,flb=-np.inf,branched_var=None,lim=None,sense=None): 
        """ Constructor 
            input:
                flb .. lower bound for the objective function
                branched_vars .. list of variables that are branched in this node
                lims .. branching limits for the variables
                sense .. branching direction ('<' or '>')
                
            parameters
                xR .. solution of the relaxation
                lb .. f(xR)
        """
        self.xR = []
        self.lb = flb
        self.var = branched_var
        self.lim = lim
        self.sense = sense
        
    def level(self):
        """ Level of the node: number of branchings """
        return len(self.var)
        
    def add_var_branch(self,var,lim,sense):
        """ Adds variable branching """
        
        self.var.append(var)
        self.lim.append(lim)
        self.sense.append(sense)
        
    def is_feasible(self):
        """ Checks if the lower bounding problem of 
            the node is feasible
        """
        if self.lb > -np.inf:
            return True
        else:
            return False
    
    def info(self):
        """ Prints relevant node data """        
        
        #print(self.sense,self.lims,self.vars)
        print('** BnB Node **')
        print('Branched variable:')
        if self.var is not None:
            print('{0:s} {1:s} {2:4.2f}'.format(self.var.name,self.sense,self.lim))
        #for i, var in enumerate(self.vars):
        #    print('{0:s} {1:s} {2:4.2f}'.format(var.name,self.sense[i],self.lims[i]))
        
    

class BnB(OptSolver):
    """ Class for Branch and Bound tree
        Includes functionalities for running branch and bound algorithm
        
        Data structures:
            Tree (from treelib): branching tree. Each node contains information
            of one branching. When a lower bounding problem is solved at a node,
            the tree is traversed until root to find the branching variables,
            limits and branching directions.
            
            
        
    """
    
    def __init__(self,problem,lb_solver=None):
        """ Constructor
            input:
                problem .. Optimization problem class object that
                            contains the problem data.
                            
                lb_solver .. OptSolver class object that implements the
                             solution method for the continuous relaxation
                            
            attributes:
                nodes .. list of nodes to be explored. Will be dynamically
                         updated during the iterations.
                max_iters .. maximum number of iterations allowed
                xbest .. variable values of the incumbent
                fbest .. objective function value of the incumbent
                
                
            self.constr_vals = np.array([])
            self.X = np.array([])
            self.problem = None
            self.fvals = []
            self.xvals = []
            self.best_f = np.inf
            self.best_x = None
        """
        super().__init__()
        
        self.created_nodes = 0
        self.problem = problem
        self.lb_solver = lb_solver
        # Initialize empty list of nodes
        # These are node identifier values correspondng to tree nodes
        self.nodes = []
        
        self.max_iters = 100
        
        self.orig_lb = []
        self.orig_ub = []
        
        self.set_original_bounds()
        
        # treelib Tree object that contains information about the 
        # branchings
        self.tree = Tree()
        #self.best_x = []
        #self.fbest = np.inf
        
    def add_node(self,node,tag=None,identifier=None,parent=None):
        
        #self.nodes.append(node)
        self.tree.create_node(tag,identifier,parent,data=node)
        
        self.nodes.append(identifier)
        
    def set_original_bounds(self):
        """ Stores the original variable bounds """
        for var in self.problem.vars:
            self.orig_lb.append(deepcopy(var.lb))
            self.orig_ub.append(deepcopy(var.ub))
            
    def set_variable_bounds(self):
        """ Stores the original variable bounds """
        for i, var in enumerate(self.problem.vars):            
            var.lb = deepcopy(self.orig_lb[i])
            var.ub = deepcopy(self.orig_ub[i])            
    
    def max_level_nodes(self):
        """ Finds the maximum level in the list of nodes, and returns
            the corresponding nodes
        """
        level_max = -1
        max_nodes = []    
        
        for node in self.nodes:            
            if node.level() > level_max:
                """ If the level of current node is larger than the current
                    maximum, set the level of the node as the maximum and
                    set the list of nodes to the current node
                """
                level_max = node.level()
                max_nodes = [node]
            elif node.level() == level_max:
                """ If the level of the current node equals the current maximum,
                    add the node to the list of nodes
                """
                max_nodes.append(node)
            
        return level_max, max_nodes
    
    def min_level_nodes(self):
        """ Finds the minimum level in the list of nodes, and returns
            the corresponding nodes
        """
        level_min = np.inf
        min_nodes = []
        
        for node in self.nodes:          
            if node.level() < level_min:
                """ If the level of current node is smaller than the current
                    minimum, set the level of the node as the minimum and
                    set the list of nodes to the current node
                """
                level_min = node.level()
                min_nodes = [node]
            elif node.level() == level_min:
                """ If the level of the current node equals the current minimum,
                    add the node to the list of nodes
                """
                min_nodes.append(node)
            
        return level_min, min_nodes
    
    def select_node(self,strategy="depth_first"):
        """ Selects node 
            input:
                strategy .. 
                'depth_first': choose the node from among those with lowest level
                'breath_first: chooses the node from among those with highest level
            
            output:
                treelib.node object
        """
                 
    
        if strategy == 'depth_first':
            """ Find the deepest nodes """
            depth = self.tree.depth()        
            nodes = list(self.tree.filter_nodes(lambda x: self.tree.depth(x)==depth))
            #print(nodes)
            #level, nodes = self.max_level_nodes()
            
        elif strategy == 'breath_first':
            """ Find the nodes with smallest level """
            level, nodes = self.min_level_nodes()
        
        
        print('Node selection:')
        print('Strategy: {0}'.format(strategy))        
        
        """ From among the set of nodes chosen according to the strategy,
            choose the one with the smallest lower bound of the objective
            function
        """
        fmin = np.inf
        for node in nodes:
            if node.data.lb < fmin:
                fmin = node.data.lb
                selected_node = node
        
        
        self.nodes.remove(selected_node.identifier)
        
        selected_node.data.info()
        
        return selected_node
    
    def lower_bound(self,node):
        """ Finds a lower bound for the problem at the 'node' by
            solving a relaxation of the problem and imposing the
            branching limits given by the node.
        """
        
        """ Set the variable bounds to the problem corresponding
            to the branches in the node
        """
        
        """ Set variable bounds to their original values """
        self.set_variable_bounds()
        
        for var in self.problem.vars:
            print("{0} <= {1:s} <= {2}".format(var.lb,var.name,var.ub))
        
        #subtree = self.tree.expand_tree(node.identifier,reverse=True)
        
        
        
        #print(subtree)
        
        #for node in subtree:
        for sub_node in list(self.tree.rsearch(node.identifier)):
            if not self.tree[sub_node].is_root():
                bnode = self.tree[sub_node]            
                if bnode.data.sense == '<':
                    bnode.data.var.ub = min(bnode.data.var.ub,bnode.data.lim)
                else:
                    bnode.data.var.lb = max(bnode.data.var.lb,bnode.data.lim)                
        
        for var in self.problem.vars:
            print("{0} <= {1:s} <= {2}".format(var.lb,var.name,var.ub))
        
        """
        for i, var in enumerate(node.vars):
            if node.sense[i] == '<':
                var.ub = min(var.ub,node.lims[i])
            else:
                var.lb = max(var.lb,node.lims[i])
        """
                
        """ Run 'solver'. 
            Solver should be LP or NLP solver
            Examples:
                a) few iterations of SLP
                b) NLP method to optimality -> heuristic, if no convexification is done
    
            If problem is infeasible, the node can be fathomed.
            If the problem is feasible, update the lower bound and set xR
        """
        flb, xlb, nit = self.lb_solver.solve(self.problem)
        if self.lb_solver.result.success == True:
            self.tree[node.identifier].data.lb = flb
            self.tree[node.identifier].data.xR = xlb   
        elif self.lb_solver.feasible:
            self.tree[node.identifier].data.lb = flb
            self.tree[node.identifier].data.xR = xlb 
            
        
        return self.lb_solver.feasible
        
    
    def upper_bound(self,node):
        """ Try to find feasible solution using a heuristic etc. """
        pass
    
    
    def variables_at_discrete_values(self,node):
        """ Checks, whether all discrete variables have discrete
            values at the solution of the relaxation at given node.
        """ 
        xR = node.xR
            
        print(xR)
        res = []
        
        for i, var in enumerate(self.problem.vars):
            res.append(var.is_allowed_value(xR[i]))
            
        
        return res
        
            
    def branch(self,node,strategy='max_int_violation'):
        """ Branching of a node """
                
        xR = node.data.xR
        #print(xR)
        
        if strategy == 'max_int_violation':
            max_violation = 0
            branching_var = None
            for i, var in enumerate(self.problem.vars):
                violation_i = var.discrete_violation(xR[i])
                if violation_i > max_violation:
                    max_violation = violation_i
                    branching_var = var                    
                    xr = xR[i] 
        
        """ Create branching values """
        low_branch_value = branching_var.smaller_discrete_value(xr)
        high_branch_value = branching_var.larger_discrete_value(xr)
        
        """ Create new nodes """
        new_node1 = BnBNode(flb=deepcopy(node.data.lb),branched_var=branching_var,lim=low_branch_value,sense='<')
        tag1 = branching_var.name + ' <= ' + str(low_branch_value)
        id1 = self.created_nodes + 1
        
        new_node2 = BnBNode(flb=deepcopy(node.data.lb),branched_var=branching_var,lim=high_branch_value,sense='>')
        tag2 = branching_var.name + ' >= ' + str(high_branch_value)
        #id2 = str(self.tree.size()+1)
        id2 = self.created_nodes + 2
        self.created_nodes += 2

        self.add_node(new_node1,tag1,id1,parent=node)
        self.add_node(new_node2,tag2,id2,parent=node)
        
        self.tree.show()
        
    
    def pre_process(self,node):
        """ Pre-processing of a node """
        pass
    
    def post_process(self,node):
        """ Post-processing of a node """
        pass    
    
    def solve(self):
        """ Runs the Branch and Bound algorithm """
        
        iteration = 0
        
        # Initialize tree
        root = BnBNode()
        self.add_node(root,"Root","root")
        
        self.tree.show()        
        
        while iteration <= self.max_iters or not self.nodes:
            
            print("Begin iteration {0}:".format(iteration))
            print("Unexplored nodes: {0}".format(len(self.nodes)))
            """ Choose node to be solved
                This is treelib.Node object
            """
            node = self.select_node(strategy="depth_first")            
            
            # Preprocessing:
            # Tighten bounds of the problem
            self.pre_process(node)
            
            # Find lower bound for the node
            feasible_lower_bound = self.lower_bound(node)
            
            if feasible_lower_bound:
                print("Lower bounding problem is feasible:")
                print("Lower bounding obtained: {0:4.2f} (f_incumbent = {1:4.2f})".format(node.data.lb,self.best_f))                
            else:
                print("Lower bounding problem is infeasible.")
            # Try to find a feasible solution in the node, to be
            # used as an upper bound for the problem.
            self.upper_bound(node)
            
            # Tighten the bounds based on the solution found
            self.post_process(node)
            
            """ 3. The node can be fathomed, if 
                a) lower bound is greater than best known feasible solution
                b) relaxed problem is infeasible
            """             
            #if node.data.is_feasible():
            if feasible_lower_bound:
                disc = self.variables_at_discrete_values(node.data)
                # Here, the relaxation is feasible
                if node.data.lb > self.best_f:
                    # Node can be fathomed
                    pass
                elif all(disc):
                    """ If all discrete variables have discrete values                    
                        in the solution xR, the node can be fathomed
                        If the node provides an improvement to the
                        best known solution, make it as the incumbent.
                    """
                    if node.data.lb < self.best_f:
                        self.best_f = node.lb
                        self.best_x = node.xR
                else:
                    """ In this case, some of the discrete variables
                        have non-integer/non-discrete values at xR.
                        Then, the node is branched.
                    """
                    print("Branching")
                    self.branch(node)
            else:
                """ Remove node from branching tree """                
                self.tree.remove_node(node.identifier)
                    
            
            iteration += 1
        
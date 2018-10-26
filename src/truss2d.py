# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:11:35 2018

@author: huuskoj
"""
from frame2d.frame2d import Frame2D, FrameMember
from fem.frame.frame_fem import FrameFEM

import math
import numpy as np
import matplotlib.pyplot as plt

class Truss2D(Frame2D):
    
    def __init__(self):
        super().__init__()
        self.top_chord = None
        self.bottom_chord = None
        self.webs = {}
        self.f = FrameFEM()
   
    def add(self, this):
        
        if isinstance(this, TopChord):
            self.top_chord = this
            self.members[len(self.members)] = this
        elif isinstance(this, BottomChord):
            self.bottom_chord = this   
            self.members[len(self.members)] = this
        elif isinstance(this, TrussWeb):
            this.mem_id = len(self.webs)
            self.webs[this.mem_id] = this
            self.members[len(self.members)] = this
            
            
            
class TrussMember(FrameMember):
    def __init__(self, coordinates, alpha1=25, alpha2=25, mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)        
    
    @property
    def angle(self):
        """
        Returns member's angle in radians
        """
        x0, y0 = self.coordinates[0]
        x1, y1 = self.coordinates[1]
        angle = math.atan((y1-y0)/(x1-x0))
        return angle
    
    def calc_nodal_coordinates(self, num_nodes):
        """
            Calculates node locations along member
            Coordinates used are global coordinates
            Adds locations to a list where they can be accessed later
            :param num_elements: int, number of elements to be created
        """
        self.nodal_coordinates = []
        start_node = self.coordinates[0]
        end_node = self.coordinates[1]

        x0 = start_node[0]
        y0 = start_node[1]
        Lx = end_node[0] - start_node[0]
        Ly = end_node[1] - start_node[1]
        dx = Lx / num_nodes
        dy = Ly / num_nodes
        for i in range(num_nodes):
            x = (i) * dx + x0
            y = (i) * dy + y0
            loc = i / num_nodes
            node = [x, y]
            if node not in self.nodal_coordinates:
                self.nodal_coordinates.append(node)
            if loc not in self.loc:
                self.loc.append(loc)
        if [end_node[0], end_node[1]] not in self.nodal_coordinates:
            self.nodal_coordinates.append([end_node[0], end_node[1]])
        if 1 not in self.loc:
            self.loc.append(1)
            
        if self.is_generated:
            for i, node in enumerate(self.nodes.values()):
                node.x = np.array(self.nodal_coordinates[i])
    
    def plot(self, print_text, c):

        X = self.coordinates
        if c:
            if self.is_strong_enough:
                color = 'green'
            else:
                color = 'red'
        color = 'k'
        # Plot members
        plt.plot([X[0][0], X[1][0]], [X[0][1], X[1][1]], color)
        # Plot text
        if print_text:
            # Calculate text location
            delta_y = X[1][1] - X[0][1]
            delta_x = X[1][0] - X[0][0]
            if delta_x == 0:
                rot = 90
            elif delta_y == 0:
                rot = 0
            else:
                rot = math.degrees(math.atan(delta_y / delta_x))
            
            x = (X[0][0] + X[1][0]) / 2
            y = (X[1][1] + X[0][1]) / 2
            horzalign = 'center'
            vertalign = 'bottom'
         
            plt.text(x, y, str(self.mem_id) + ": " + self.profile,
                     rotation=rot, horizontalalignment=horzalign,
                     verticalalignment=vertalign)
        
            
class TopChord(TrussMember):
    def __init__(self, coordinates, mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'top_chord'


class BottomChord(TrussMember):
    def __init__(self, coordinates,  mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'bottom_chord'

class TrussWeb(TrussMember):
    def __init__(self, coordinates, alpha1=0, alpha2=0, mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)
        self.mtype = 'web'
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        
        
class TrussJoint():
    def __init__(self, chord, coordinate, joint_type="N", g1=0.1, g2=0.1):
        
        self.chord = chord
        self.coordinate = coordinate
        self.joint_type = joint_type
        self.g1 = g1
        self.g2 = g2
        self.chord = None
        self.webs = {}
    

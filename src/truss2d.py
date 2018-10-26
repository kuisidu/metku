# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 10:11:35 2018

@author: huuskoj
"""
from frame2d.frame2d import Frame2D, FrameMember

class Truss2D(Frame2D):
    
    def __init__(self):
        self.top_chord = None
        self.bottom_chord = None
        self.webs = {}
        self.members = {}
    
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
    pass
            
class TopChord(TrussMember):
    def __init__(self, coordinates, alpha1=25, alpha2=25, mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'top_chord'


class BottomChord(TrussMember):
    def __init__(self, coordinates, alpha1=25, alpha2=25, mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'bottom_chord'

class TrussWeb(TrussMember):
    def __init__(self, coordinates, alpha1=25, alpha2=25, mem_id="",
                 profile="SHS 50x5", material="S355"):
        super().__init__(coordinates, mem_id, profile, material)

        self.mtype = 'web'
        self.alpha1 = alpha1
        self.alpha2 = alpha2

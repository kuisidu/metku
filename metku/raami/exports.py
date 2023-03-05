# -*- coding: utf-8 -*-
# Copyright 2022 Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Author(s): Kristo Mela
"""
Created on Tue Apr 19 10:14:37 2022

Raami functions for exporting to various programs

@author: kmela
"""

def write_elset(file,elset_name,els,instance=None):
    
    header = '*Elset, elset=' + elset_name
    
    if instance is None:
        header += '\n'
    else:
        header += ', instance=' + instance + '\n'
    
    file.write(header)
    i = 0
    while i < len(els):
        file.write(', '.join(str(r) for r in els[i:i+16]))
        file.write('\n')
        i += 16

class AbaqusOptions:
    """ Storage for ABAQUS export options """
    
    def __init__(self,**kwargs):
        """
        

        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        self.x_monitor = 0
        self.n_monitored = 1
        self.mpc = []   # multi-point constraints 
        self.ecc_elements = [] # list of eccentricity elements
        self.top_gap_elements = [] # gap elements of the top chord
        self.bottom_gap_elements = [] # gap elements of bottom chord
        
        self.elsets = {} # *Elset groups
        
        self.load_elsets = {}
        
        # Individual members to be included as element sets
        self.included_members = []
        
        attr = self.__dict__.keys()
        
        for key, val in kwargs.items():
            if key in attr:
                setattr(self,key,val)
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 10:14:37 2022

Raami functions for exporting to various programs

@author: kmela
"""


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
        self.mpc = []
        self.ecc_elements = []
        self.top_gap_elements = []
        self.bottom_gap_elements = []
        
        self.elsets = {}
        
        attr = self.__dict__.keys()
        
        for key, val in kwargs.items():
            if key in attr:
                setattr(self,key,val)
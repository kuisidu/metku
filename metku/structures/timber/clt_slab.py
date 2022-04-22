# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 21:54:04 2021

@author: kmela
"""

from metku.materials.timber_data import Timber, T

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.path as mpath

class CLTLayer:
    """ Class for CLT slab layers """
    
    def __init__(self,thickness=20,material=Timber(T.C24),orientation=0):
        """
        

        Parameters
        ----------
        thickness : float, optional
            Thickness of the layer. The default is 20.
        material : Timber, optional
            Material propertyies (see timber_data.py). The default is C24.
        orientation : float, optional
            Orientation of the grains. 0 means grains are long the main axis
            and 90 that they are perpendicular to the main axis. The default is 0.

        Returns
        -------
        None.

        """
        
        self.t = thickness
        self.mat = material
        self.orientation = orientation

        self.bx = 1000
        self.by = 1000
    
    def __repr__(self):
        
        return f"({self.t}/{self.mat}/{self.orientation})"
    
    def Ay(self):
        """ Area of the layer with the face in the y-z coordinates """
        
        return self.bx*self.t
    
    def Az(self):
        """ Area of the layer with the face in the x-z coordinates """
        
        return self.by*self.t
    
    def Iy(self):
        """ Second moment of area with respect of the centroid of
            the layer. x-axis is along the span of the slab, and y-axis
            is the horizontal axis
        """
        
        return self.bx*self.t**3/12

    def Iz(self):
        """ Second moment of area with respect of the centroid of
            the layer. x-axis is along the span of the slab, and y-axis
            is the horizontal axis
        """
        
        return self.by*self.t**3/12
    
    def EIy(self):
        """ Bending stiffness of the layer with respect to its
            neutral axis.
        """
        
        if self.orientation == 0:
            E = self.mat.E0mean
        else:
            E = self.mat.E90mean
        
        return E*self.Iy()

    def EIz(self):
        """ Bending stiffness of the layer with respect to its
            neutral axis.
        """
        
        if self.orientation == 0:
            E = self.mat.E90mean
        else:
            E = self.mat.E0mean
        
        return E*self.Iz()


class CLTSlab:
    """ Class for CLT slabs """
    
    def __init__(self,layers,width=1000,length=5000):
        """
        

        Parameters
        ----------
        layers : list of CLTLayer type objects or dict with lists for
                 the following keys: "thickness", "material", "orienatation"
            Layers of the slab.
        width : float
            width of the slab
        length : float
            length of the slab

        Returns
        -------
        None.

        """
        
        self.b = width
        self.L = length
        
        if isinstance(layers,list):
            self.layers = layers
        elif isinstance(layers,dict):
            thickness = layers['thickness']
            materials = layers['material']
            orientations = layers['orientation']
            
            self.layers = []
            
            for t, mat, orient in zip(thickness,materials,orientations):
                self.layers.append(CLTLayer(t,mat,orient))
        
        # Set layer width and length
        for layer in self.layers:
            layer.bx = self.b
            layer.by = self.L
    
    def __repr__(self):
        
        s = 'CLT Slab: '
        
        for i, layer in enumerate(self.layers):
            s += f'{i}: ' + layer.__repr__()
            
            if i < len(self.layers)-1:
                s += ' / '
        
        return s
    
    @property
    def hclt(self):
        """ Total height of the slab """
        
        return sum([layer.t for layer in self.layers])
    
    @property
    def zgc(self):
        """ Location of the centroid in thickness direction """
        return 0.5*self.hclt
    
    @property
    def Ay(self):
        """ Cross-sectional area in y-z-coordinates """
        
        return self.width*self.hclt
    
    @property
    def Az(self):
        """ Cross-sectional area in x-z-coordinates """
        
        return self.length*self.hclt
    
    def EIy(self):
        """ Bending stiffness of the CLT slab with respect to y-axis """
        
        EIy0 = sum([layer.EIy() for layer in self.layers])
        

    def draw(self,axes=None,origin=[0,0],clt_width=0):
        
        if axes is None:
            fig, ax = plt.subplots(1)
        else:
            ax = axes
        
        if clt_width == 0:
            clt_width = self.by
        
        h = self.hclt
        for layer in self.layers:
            if layer.orientation == 0:
                layer_hatch = ''
            else:
                layer_hatch = 'x'
            new_patch = patches.Rectangle((origin[0], origin[1]+h-layer.t),
                                       width=clt_width, height=layer.t,
                                       fill=False, hatch=layer_hatch)            
            ax.add_patch(new_patch)
            h -= layer.t
        
        
        ax.set_xlim(origin[0], 100)
        ax.set_ylim(0, self.hclt)
        
        
        ax.set_aspect('equal')
        
        plt.show()
        return ax

if __name__ == '__main__':
    
    layers = {'thickness':[20,20,20],
              'material':[Timber(T.C24),Timber(T.C24),Timber(T.C24)],
              'orientation':[0,90,0]}
    
    s = CLTSlab(layers)

    s.draw()
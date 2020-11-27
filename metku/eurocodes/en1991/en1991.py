# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 08:20:46 2020

Data for loads and load combinations

EN 1990 and EN 1991

@author: kmela
"""

class Load:
    """ General class for loads """
    
    def __init__(self,value=1.0,comb_factors=[1.0,1.0,1.0]):
        """ Constructor """
        
        self.value = value
        self.comb_factors = comb_factors
        
    def combined(self,combi=0):
        """ Returns the combined value
            :param combi: integer (0, 1, 2) for choosing the combination factor
        """
        return self.value*self.comb_factors[combi]

class ImposedLoad(Load):
    """ Class for imposed (live) loads """
    
    def __init__(self,value_ground=2.0,comb_factors=[0.7,0.5,0.2]):
        """ Constructor
        
            :param value_ground: value of snow load on ground [kN/m2]
        """
        
        super().__init__(value_ground,comb_factors)

class SnowLoad(Load):
    """ Class for snow loads """
    
    def __init__(self,value_ground=2.0,comb_factors=[0.7,0.5,0.2]):
        """ Constructor
        
            :param value_ground: value of snow load on ground [kN/m2]
        """
        
        if value_ground < 2.75:
            comb_factors = [0.7,0.4,0.2]
        else:
            comb_factors = [0.7,0.5,0.2]
            
        super().__init__(value_ground,comb_factors)
        
    def __repr__(self):
        
        return "Snow load: basic value {0:3.1f} kN/m2".format(self.value)
    
class WindLoad(Load):
    """ Class for wind loads """
    
    def __init__(self,value=0.6,comb_factors=[0.6,0.2,0.0]):
        """ Constructor
        
            :param value: basic value of wind load [kN/m2]
        """
            
        super().__init__(value,comb_factors)
        
    def __repr__(self):
        
        return "Wind load: basic value {0:3.1f} kN/m2".format(self.value)


if __name__ == "__main__":
    
    s = SnowLoad(2.5)
    w = WindLoad(0.6)
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-8 Rules for connections

@author: kmela
"""

import math
import numpy as np

from eurocodes.en1993.constants import gammaM5

class RHSJoint:
    """ Class for joints between RHS members """
    
    def __init__(self,chord_profile):
        """ Constructor
        """
        self.chord = chord_profile        
        
    def h0(self):
        return self.chord.H
        
    def b0(self):
        return self.chord.B
        
    def t0(self):
        return self.chord.T
    
    def r0(self):
        return self.chord.R
    
    def fy0(self):
        return self.chord.fy
    
    def N0(self):
        return self.chord.NEd
    
    def M0(self):
        return self.chord.MEd
        
    def gamma(self):
        return 0.5*self.h0()/self.t0()
        
    def b_straigth(self):
        """ Length of the straight part of the chord profile """
        return self.b0()-2*self.r0();
        
    def beta(self):
        """ Abstract method for beta that depends on the joint type """
        return None
    
    def eval_KN(self):
        """ Evaluate stress function """
        # by default, assume the chord face is in tension: kn = 1.0
        kn = 1.0
        n = self.eval_N()
        b = self.beta()
        if n > 0.0: # compression                
            kn = min(1.3-0.4*n/b,1.0)
        
        return kn
        
    def strength_reduction(self):
        """ Reduce the strength of joint for steel greater than S355
        """
        r = 1.0
        if self.fy0() > 355:
            r = 0.9
        elif self.fy0() > 460:
            r = 0.8
        return r
        
    def eval_N(self):
        """ Evaluate chord stress function """
        # s0Ed > 0 for compression
        s0Ed = self.chord_stress()
        n = s0Ed/self.fy0()/gammaM5
        return n
    
    
    def chord_stress(self):
        """ Chord stress ratio """
        A0 = self.chord.A
        Wel0 = self.chord.Wel[0];                 
        # s0Ed is positive, if the stress in the chord
        # is compressive.
        # CAREFUL WITH UNITS HERE!
        s0Ed = -self.N0()*1e3/A0+self.M0()*1e6/Wel0
        
        return s0Ed

        
    def new_evalN(self):
        """ Evaluate stress ratio as proposed """
        # s0Ed > 0 for compression
        s0Ed = self.new_chord_axial_stress()
        # n < 0 for compression
        n = -s0Ed/self.fy0()/gammaM5
        
        return n
        
    def new_chord_axial_stress(self):
        """ Axial stress in the chord face at the joint
            Based on the new proposal, where plastic bending resistance is 
            used instead of elastic
        """
        A0 = self.chord.A
        Wpl0 = self.chord.Wpl[1]
        # s0Ed < 0 for compression
        s0Ed = -self.N0()*1e3/A0+self.M0()*1e6/Wpl0
        return s0Ed



    
class RHSKGapJoint(RHSJoint):
    """ For K and N Gap joints between RHS members """
    
    def __init__(self,chord,braces,angles,gap=20):
        """ Constructor
        """
        RHSJoint.__init__(self,chord)
        self.braces = braces
        self.gap = gap
        self.angles = angles
        
    def h(self):
        """ Brace heights """
        return [self.braces[0].H,self.braces[1].H]
    
    def b(self):
        """ Brace widths """
        return [self.braces[0].B,self.braces[1].B]
    
    def t(self):
        """ Brace thickness """
        return [self.braces[0].T,self.braces[1].T]
    
    def fy(self):
        """ Brace yield strength """
        return [self.braces[0].fy,self.braces[1].fy]
    
    def beta(self):
        return sum([self.b(),self.h()])/4/self.b0()
    
    
    def eccentricity(self):
        t = self.angles
        h = self.h()
        e0 = 0.5*h[0]/math.sin(math.radians(t[0])) + 0.5*h[1]/math.sin(math.radians(t[1])) + self.gap
        e = e0*math.sin(math.radians(t[0]))*math.sin(math.radians(t[1]))/math.sin(math.radians(sum(t)))-0.5*self.h0()
        return e
    
    # Resistance formulas
    def chord_face_failure(self):
                        
        b = self.Beta()
        g = self.gamma()                           
        kn = self.eval_KN()
        r = self.strength_reduction()
        fy0 =  self.fy0()
        t0 = self.t0()
        NRd = []
        NRd0 = r*8.9*kn*fy0*t0**2*math.sqrt(g)*b/gammaM5
        for angle in self.angles:
            s = math.sin(math.radians(angle))            
            NRd.append(NRd0/s)
            
        return NRd
        
    def chord_shear(self):
        Aeta = self.a_eta()
        # sin of brace angles
        s = np.sin(np.radians(np.array(self.angles)))
        fy0 = self.fy0()
        # braces
        NiRd = fy0*Aeta/math.sqrt(3)/s/gammaM5
        # chord                                         
        N0Rd = fy0*((self.chord.A-Aeta)+Aeta*math.sqrt(1-(self.V0()*1e3/self.chord.VplRd)**2))/gammaM5
            
        r = self.strength_reduction()
        NiRd = r*NiRd
        N0Rd = r*N0Rd
        return NiRd, N0Rd
        
    def a_eta(self):
        a = math.sqrt(1/(1+4*self.gap**2/3/self.t0()**2))
        return (2*self.h0()+a*self.b0())*self.t0()
        
    def brace_failure(self):
        r = self.strength_reduction()
        b = np.array(self.b())
        h = np.array(self.h())
        t = np.array(self.t())
        fy = np.array(self.fy())
        fy0 = self.fy0()
        t0 = self.t0() 
        beff = min(b,10/(self.b0()/t0)*fy0*t0/(fy*t)*b)
            
        NiRd = r*self.fy*t*(2*h-4*t+b+beff)/gammaM5
        return NiRd
        
    def punching_shear(self):
        b = np.array(self.b())
        fy0 = self.fy0()
        t0 = self.t0()
        bep = min(b,10/(self.b0()/t0)*b)
        r = self.strength_reduction()
        
        s = np.sin(np.radians(np.array(self.angles)))
            
        NRd = r*fy0*t0/math.sqrt(3)/s*(2*self.h()/s+b+bep)/gammaM5
        return NRd

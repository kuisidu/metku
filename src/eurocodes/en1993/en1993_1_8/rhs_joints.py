# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-8 Rules for connections

@author: kmela
"""

import math
import numpy as np

try:
    from src.eurocodes.en1993.constants import gammaM5, E
except:
    from eurocodes.en1993.constants import gammaM5, E

class RHSJoint:
    """ Class for joints between RHS members """
    
    def __init__(self,chord_profile, N0=None, V0=None, M0=None):
        """ Constructor
        """
        self.chord = chord_profile
        if N0 is None:
            N0 = chord_profile.Ned
        self.N0 = N0

        if V0 is None:
            V0 = chord_profile.Ved
        self.V0 = V0

        if M0 is None:
            M0 = chord_profile.Med
        self.M0 = M0
        
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

    def gamma(self):
        return 0.5*self.h0()/self.t0()
        
    def b_straigth(self):
        """ Length of the straight part of the chord profile """
        return self.b0()-2*self.r0()
        
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
            kn = min(1.3-0.4*n/self.beta(), 1.0)
        
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
        s0Ed = -self.N0/A0+self.M0/Wel0
        
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
        s0Ed = -self.N0/A0+self.M0/Wpl0
        return s0Ed



    
class RHSKGapJoint(RHSJoint):
    """ For K and N Gap joints between RHS members """
    
    def __init__(self,chord,braces,angles,gap=20, **kwargs):
        """ Constructor
        """
        RHSJoint.__init__(self,chord, **kwargs)
        self.braces = braces
        self.gap = gap
        self.angles = angles
        
    def h(self):
        """ Brace heights """
        return np.asarray([self.braces[0].H,self.braces[1].H])
    
    def b(self):
        """ Brace widths """
        return np.asarray([self.braces[0].B,self.braces[1].B])
    
    def t(self):
        """ Brace thickness """
        return np.asarray([self.braces[0].T,self.braces[1].T])
    
    def fy(self):
        """ Brace yield strength """
        return np.asarray([self.braces[0].fy,self.braces[1].fy])
    
    def beta(self):
        return (sum(self.b())+ sum(self.h())) / (4*self.b0())
    
    
    def eccentricity(self):
        t = self.angles
        h = self.h()
        e0 = 0.5*h[0]/math.sin(math.radians(t[0])) + 0.5*h[1]/math.sin(math.radians(t[1])) + self.gap
        e = e0*math.sin(math.radians(t[0]))*math.sin(math.radians(t[1]))/math.sin(math.radians(sum(t)))-0.5*self.h0()
        return e
    
    # Resistance formulas
    def chord_face_failure(self):
                        
        b = self.beta()
        g = self.gamma()                           
        kn = self.eval_KN()
        r = self.strength_reduction()
        fy0 = self.fy0()
        t0 = self.t0()
        NRd = []
        NRd0 = r*8.9*kn*fy0*t0**2*math.sqrt(g)*b/gammaM5
        for angle in self.angles:
            s = math.sin(math.radians(angle))            
            NRd.append(NRd0/s)
            
        return np.asarray(NRd)
        
    def chord_shear(self):
        Aeta = self.a_eta()
        # sin of brace angles
        s = np.sin(np.radians(np.array(self.angles)))
        fy0 = self.fy0()
        # braces
        NiRd = fy0*Aeta/math.sqrt(3)/s/gammaM5
        # chord
        if self.V0/self.chord.shear_force_resistance() > 1:
            # print("Chord is not strong enough for shear forces")
            N0Rd = fy0*((self.chord.A-Aeta))/gammaM5

        else:
            N0Rd = fy0*((self.chord.A-Aeta)+Aeta*math.sqrt(1-(self.V0/self.chord.shear_force_resistance())**2))/gammaM5
            
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
        beff = 10/(self.b0()/t0)*fy0*t0/(fy*t)*b
        beff_arr = np.array([min(x,beff[i]) for i, x in enumerate(b)])
            
        NiRd = r*self.fy()*t*(2*h-4*t+b+beff_arr)/gammaM5
        return NiRd
        
    def punching_shear(self):
        b = self.b()
        fy0 = self.fy0()
        t0 = self.t0()
        bep = 10/(self.b0()/t0)*b 
        r = self.strength_reduction()
        
        s = np.sin(np.radians(np.array(self.angles)))
            
        NRd = r*fy0*t0/math.sqrt(3)/s*(2*self.h()/s+b+bep)/gammaM5
        return NRd


class RHSYJoint(RHSJoint):
    
    def __init__(self, chord_profile, brace, angle, **kwargs):
        super().__init__(chord_profile, **kwargs)
        self.brace = brace
        self.angle = angle
        
        
    def h(self):
        """ Brace heights """
        return self.brace.H
    
    def b(self):
        """ Brace widths """
        return self.brace.B
    
    def t(self):
        """ Brace thickness """
        return self.brace.T
    
    def fy(self):
        """ Brace yield strength """
        return self.brace.fy
    
    def beta(self):
        return self.b() / self.b0()
        
        
    def chord_face_failure(self):
                        
        b = self.beta()
        if b >= 1:
            b = 0.99
        kn = self.eval_KN()
        n = self.eval_N()
        r = self.strength_reduction()
        fy0 = self.fy0()
        t0 = self.t0()
        NRd = r**kn*fy0*t0**2 / ((1-b)*np.sin(np.radians(self.angle)))
        NRd = NRd * (2*n / np.sin(np.radians(self.angle)) + 4*np.sqrt(1-b))
        NRd = NRd/ gammaM5

        return NRd
    
    
    def slend(self):

        slend = 3.46*(self.h0()/(self.t0() -2) *np.sqrt(1/np.sin(np.radians(self.angle))))
        slend = slend / (np.pi * np.sqrt(E / self.fy0()))
        return slend
   
    def chord_web_buckling(self):
        alpha = self.chord.imp_factor[0]
        slend = self.slend()
        phi = 0.5*(1 + alpha*(slend - 0.2) + slend **2)
        chi = 1 / (phi + np.sqrt(phi**2 - slend**2))
        fb = chi * self.fy0()
        kn = self.eval_KN()
        r = self.strength_reduction()
        t0 = self.t0()
        h = self.h()
        NRd = r**kn*fb*t0 / (np.sin(np.radians(self.angle)))
        NRd = NRd * (2*h / np.sin(np.radians(self.angle)) + 10*t0)
        NRd = NRd / gammaM5
        return NRd
    
    def brace_failure(self):
        r = self.strength_reduction()
        b = self.b()
        h = self.h()
        t = self.t()
        fy = self.fy()
        fy0 = self.fy0()
        t0 = self.t0() 
        beff = 10/(self.b0()/t0)*fy0*t0/(fy*t)*b
            
        NiRd = r*self.fy()*t*(2*h-4*t+2*beff)/gammaM5
        return NiRd
    
    
        
        
    
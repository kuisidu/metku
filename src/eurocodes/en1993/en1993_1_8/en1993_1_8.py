# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-8 Rules for connections

@author: kmela
"""

import math

from eurocodes.en1993.constants import gammaM5

class RHSJoint:
    """ Class for joints between RHS members """
    properties
        chord
        N0
        Np
        M0
        V0
        fy0
    end
    
    properties (Constant = true)
        E = 210 % Young's Modulus [GPa]
        nu = 0.30 % Poisson's ratio
        G = 81 % Shear modulus [GPa]
        gM5 = 1.0;
    end
    
    properties (Dependent = true)
        Gamma
        H0
        B0
        T0
        %fy0
    end
    
    methods (Abstract)
        b = Beta(obj);
        eta = Eta(obj);        
    end
    
    methods
        function obj = RHSJoint(chordProfile,fy0,varargin)
            nargs = length(varargin);
            N0 = 0;
            M0 = 0;
            if nargs > 0
                N0 = varargin{1};
            end
            
            if nargs > 1
                M0 = varargin{2};
            end            
            obj.chord = chordProfile;  
            obj.fy0 = fy0;
            obj.N0 = N0;
            obj.M0 = M0;
        end                              
        
        function h0 = get.H0(obj)
            h0 = obj.chord.H;
        end
        
        function b0 = get.B0(obj)
            b0 = obj.chord.B;
        end
        
        function t0 = get.T0(obj)
            t0 = obj.chord.T;
        end
        
%         function fy0 = get.fy0(obj)
%             fy0 = obj.chord.fy;
%         end
        
        function g = get.Gamma(obj)
            g = obj.H0/2/obj.T0;
        end
        
        function bbar = BStraigth(obj)            
            bbar = obj.chord.B-2*obj.chord.R;
        end
        
        function kn = EvalKN(obj)
            n = obj.EvalN;
            b = obj.Beta;
            if n > 0 % compression                
                kn = min(1.3-0.4*n/b,1.0);
            else % tension
                kn = 1.0;
            end
        end
        
        function r = StrengthReduction(obj)
            r = 1.0;
            if obj.fy0 > 355
                r = 0.9;
            elseif obj.fy0 > 460
                r = 0.8;            
            end
        end
        
        function n = EvalN(obj)
            % s0Ed > 0 for compression
            s0Ed = obj.ChordStress;
            n = s0Ed/obj.fy0/obj.gM5;
        end
        
        function s0Ed = ChordStress(obj)
            A0 = obj.chord.A;
            Wel0 = obj.chord.Wel(1);                 
            % s0Ed is positive, if the stress in the chord
            % is compressive.
            s0Ed = -obj.N0*1e3/A0+obj.M0*1e6./Wel0;
            
            %Wpl0 = obj.chord.Wpl(1)           
            %s0Ed = -obj.N0*1e3/A0+obj.M0*1e6./Wpl0;
        end
        
        function n = NewEvalN(obj)
            % s0Ed > 0 for compression
            s0Ed = obj.NewChordAxialStress;
            % n < 0 for compression
            n = -s0Ed/obj.fy0/obj.gM5;
        end
        
        % Axial stress in the chord face at the joint
        % Based on the new proposal, where plastic bending resistance is 
        % used instead of elastic
        function s0Ed = NewChordAxialStress(obj)
            A0 = obj.chord.A;                        
            Wpl0 = obj.chord.Wpl(1);
            % s0Ed < 0 for compression
            s0Ed = -obj.N0*1e3/A0+obj.M0*1e6./Wpl0;
        end
    end
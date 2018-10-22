# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 16:06:13 2018

Cold-formed profiles

@author: kmela
"""

import numpy as np
import matplotlib.pyplot as plt
from open_prof import OpenProf

class ZProf:
    """ ZProf Class for open Z-sections
        Main purpose is to calculate cross-sectional properties as
        described in EN 1993-1-3, Annex C.    
    """
    """    
        properties (Depent = true)
            hWeb % Web height
            bTop % Width of Top flange
            bBottom % Width of Bottom flange
            gr % length needed for notional widths
            rm % corner rounding to the midline
        
    """
    def ZProf(self,Tnom,H,A,B,CA,CB,R=0,fyb=350,fub=420,Tcoat=0.04):
        """ Constructor 
            input:
                Tnom .. nominal wall thickness
                H .. total height
                A .. width of top flange
                B .. width of bottom flange
                CA .. length of top fold
                CB .. length of bottom fold
                R .. rounding
                fyb = 350 .. basic yield strength [MPa]
                fub = 420 .. basic ultimate strength [MPa]
                Tcoat = 0.04 .. coating [mm]        
        """
        
        """
            %rounding_method = 1; % corner node at intersection of midlines
            %rounding_method = 2; % corner by two elements
            rounding_method = 3; % corner by three elements
        """
        self.tnom = Tnom
        self.h = H
        self.a = A
        self.b = B
        self.ca = CA
        self.cb = CB            
        self.fyb = fy
        self.fub = fub
        self.tcoat = Tcoat
        self.tcore = Tnom-Tcoat
            
        if R == 0:
            """ According to EN 10162, Table A.2;
                S320GD+Z and S350GD+Z
            """
            self.r = 1.5*T
        else:
            self.r = R
            
        # store profile data here (instance of OpenProf class)
        self.prof = None
            
        """ Create nodes.
            a) The origin is in the intersection of
               midlines of bottom flange and web.
            b) corner roundings are taken into account by the "German
               approach"
               The nodes are numbered starting from the tip of
               the fold of the bottom flange. If CB = 0, there is
               no fold
            switch rounding_method
                case 1
                    y = [B-Tnom; 0; 0; -(A-Tnom)];
                    z = [0; 0; H-Tnom; H-Tnom];
                    if nargin > 4 && ~isempty(CA) && CA > 0
                        y = [y; -(A-Tnom)];
                        z = [z; H-Tnom - (CA-0.5*Tnom)];
                    
                    if nargin > 5 && ~isempty(CB) && CB > 0
                        y = [B-Tnom;y];
                        z = [CB-0.5*Tnom;z];
                    
                case 2
                case 3                    
                    rm = R+0.5*Tnom;
                    %  bottom flange edge stiffener
                    % (assume 90 degree angle)
                    if CB > 0                      
                        yBottomLip = B-Tnom;                      
                        y = [yBottomLip*[1;1;1];yBottomLip-0.7*rm;yBottomLip-rm];                      
                        z = [CB-0.5*Tnom;rm;0.7*rm;0;0];
                    else
                        y = B-0.5*Tnom;
                        z = 0;
                    
                    % bottom flange corner
                    y = [y;rm;0.7*rm;0;0];
                    z = [z;0;0;0.7*rm;rm];
                    % web;
                    hw = H-Tnom;
                    y = [y;0;0];
                    z = [z;hw-rm;hw-0.7*rm];
                    % top flange corner
                    y = [y;-0.7*rm;-rm];
                    z = [z;hw;hw];
                    % top flange                    
                    if CA > 0
                        yTopLip = A-Tnom;
                        y = [y;-(yTopLip-rm);-(yTopLip-0.7*rm);-yTopLip*[1;1;1]];
                        z = [z;hw;hw;hw-0.7*rm;hw-rm;hw-(CA-0.5*Tnom)];
                    else
                        y = [y;-(A-0.5*Tnom)];
                        z = [z;hw];                       
                    
                    %size(y)
                    %size(z)
            """        
            
        """
        def self = set.Tnom(self,val)
            self.Tnom = val;
            self.T = val-self.Tcoat;
            self.t(:) = self.T;
        
        
        def self = set.A(self,val)
            self.A = val;
            if self.CA > 0
                rm = self.rm;
                y0 = val-self.Tnom;
                self.y(-4:) = -[y0-rm,y0-0.7*rm,y0,y0,y0];
            else
                self.y() = -(val-0.5*self.Tnom);
        """
        
                
        def hWeb(self):
            """ Web heigth """
            return self.h-self.tnom
        
        
        def bTop(self):
            """ Width of the top flange (center line) """
            if self.ca == 0:
                r = self.a-0.5*self.tnom
            else:
                r = self.a-self.tnom
            return r
                    
        def bBottom(self):
            """ Width of the bottom flange (center line) """
            if self.cb == 0:
                r = self.b-0.5*self.tnom
            else:
                r = self.b-self.tnom
            return r
            
        def rounding_center(self):
            """ Corner rounding to the center line """
            return self.r+0.5*self.tnom
                
        def gr(self):
            """ Length gr from EN 1993-1-3
            Rm = self.rounding_center()
            a = 90
            gr = Rm*(tand(0.5*a)-sind(0.5*a))
            return gr
                
        def center_line_length(self):
            return self.ca+self.a+self.h+self.b+self.cb-8*(self.r+self.tnom)+2*math.pi()*(self.r+0.5*self.tnom)
                
        def average_yield(self):
            k = 7  # coefficient for roll forming
            n = 4  # number of 90 degree bs with r <= 5t
            fave = 0.5*(self.fyb+self.fub)
            fya = min(self.fyb + (self.fub-self.fyb)*k*n*self.t**2/self.Area,fave)
            
        
        
        def stCom = StressTopFlange(self)
            gc = self.Centroid;
            zgc = gc(2);
            zb1 = self.hWeb-zgc;
            r = min(zb1/(self.hWeb-zb1),1.0);
            stCom = self.fyb*r;
        
        
        % computes the notional width of a flange
        % input:
        %        C -- length of edge fold
        %  Bflange -- flange width
        def bp = NotionalWidthTopFlange(self)             
            bp = self.NotionalWidthFlange(self.A,self.CA,self.Tnom,self.gr);
        
        
        def bp = NotionalWidthBottomFlange(self)
            bp = self.NotionalWidthFlange(self.B,self.CB,self.T,self.gr);
        
        
        def bcp = NotionalWidthTopLip(self) 
            if self.CA > 1e-6
                bcp = self.CA-0.5*self.T-self.gr;
            else
                bcp = 0;
            
        
        
        def bcp = NotionalWidthBottomLip(self) 
            if self.CB > 1e-6
                bcp = self.CB-0.5*self.T-self.gr;
            else
                bcp = 0;
            
        
        
        def bp = NotionalWidthWeb(self)             
            bp = self.hWeb-2*self.gr;            
        
        
        
        def [beff,be,rho] = EffectiveFlange(self,flange,e)
            if nargin < 2 || isempty(flange)
                flange = 'top';
            
            
            switch flange
                case {'top',1}   
                    bp = self.NotionalWidthTopFlange;
                    C = self.CA;
%                     C = self.CA;
%                     Bflange = self.A;
                case {'bottom',0}
                    bp = self.NotionalWidthBottomFlange;
                    C = self.CB;
%                     C = self.CB;
%                     Bflange = self.B;
            
            
           % bp = self.NotionalWidthFlange(C,Bflange);
                        
            r = bp/self.T;    
            if nargin < 3
                e = Epsilon(self.fyb);
            
                                                    
            if abs(C) < 1e-6
                % no edge stiffener
                % assume constant stress in top flange:
                % this implies that stress ratio is 1.
                stressRatio = 1.0;
                ksigma = BucklingFactorOutstand(stressRatio);
                slred = PlateSlerness(r,e,ksigma);
                rho = OutstandReduction(slred);
                beff = rho*bp;
                be = beff;
            else
                % assume constant stress in top flange:
                % this implies that stress ratio is 1.
                stressRatio = 1.0;
                ksigma = BucklingFactorInternal(stressRatio);
                slred = PlateSlerness(r,e,ksigma);
                rho = InternalReduction(slred,stressRatio);
                beff = rho*bp;
                be = 0.5*beff*[1;1];
            
                        
            %self.DrawEffectiveFlange(flange,be(1),bp,stressRatio)
            %self.DrawEffectiveFlange(flange,be(1),bp,rho);            
        
        
        def [ceff,rho,bcp] = EffectiveLip(self,flange,e)
            verb = 0;
           
            if nargin < 2 || isempty(flange)
                flange = 'top';
            
            
            switch flange
                case {'top',1}               
                    bp = self.NotionalWidthTopFlange;
                    C = self.CA;
                    %Bflange = self.A;
                case {'bottom',0}
                    bp = self.NotionalWidthBottomFlange;
                    C = self.CB;
                    %Bflange = self.B;
            
            % unstiffened top flange:
            % plane element without stiffener
            %bp = self.NotionalWidthFlange(C,Bflange);
            bcp = self.NotionalWidthLip(C);
            r = bcp/self.T;                
            if nargin < 3
                e = Epsilon(self.fyb);
            
                                    
            rb = bcp/bp;
            ksigma = BucklingFactorLip(rb);            
            slred = PlateSlerness(r,e,ksigma);
            rho = OutstandReduction(slred);
            ceff = rho*bcp;
            
            if verb
                fprintf(1,'Tehollinen reunajäykiste\n');
                fprintf(1,'bpc = %6.3f\n',bcp);
                fprintf(1,'rc = %6.3f\n',r);
                fprintf(1,'bcp/bp = %6.3f\n',rb);
                fprintf(1,'ksigma = %6.3f\n',ksigma);
                fprintf(1,'slred = %6.3f\n',slred);
                fprintf(1,'rho = %6.3f\n',rho);
                fprintf(1,'ceff = %6.3f\n',ceff);
            
        
        
        def bcp = NotionalWidthLip(self,C)            
            bcp = C-0.5*self.Tnom-self.gr;
        
        
        def d = ReductionDelta(self)
            % Compute the delta to be used in reduction of 
            % cross-sectional properties taking into account
            % rounded corners
            bp(1) = self.NotionalWidthTopLip;
            bp(2) = self.NotionalWidthTopFlange;            
            bp(3) = self.NotionalWidthWeb;
            bp(4) = self.NotionalWidthBottomFlange;
            bp(5) = self.NotionalWidthBottomLip;
            d = 0.43*4*self.R/sum(bp);
        
        
        % Calculates the required shear stiffness to
        % stabilize the purlin
        def Sreq = RequiredShearStiffness(self,L)
            E = 210; % GPa
            G = 81;  % GPa
            Iw = self.WarpingConstant;
            It = self.TorsionConstant;
            I2 = self.SecondMomentArea;
            Iz = I2(2);
            h = self.H;
            Sreq = (E*Iw*pi^2/L^2 + G*It + 0.25*E*Iz*pi^2/L^2*h^2)*70/h^2;
        
        
        % Calculate effective cross section in bing
        def EffectiveSectionBing(self)
            
            % Start with effective top flange
            [beff,be,rhoF] = self.EffectiveFlange('top')
            
            % Get effective lip
            [ceff,rhoL,bcp] = self.EffectiveLip('top')
            
            if rhoF < 1.0
                % here, two points should added to mark
                % the non-effective part of the flange
            else
                % if the entire flange is effective,
                % one point should be added to the flange
            
        
        
    
    
    methods (Static = true)
        def bp = NotionalWidthFlange(Bflange,C,T,gr)
            if abs(C) < 1e-6
                b = Bflange-0.5*T;                
                bp = b-gr;
            else
                b = Bflange-T;
                bp = b-2*gr;
               
        
    
    

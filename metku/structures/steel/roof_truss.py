# Author(s): Kristo Mela
# This source code is licensed under the MIT license. See LICENSE in the repository root directory.
# Copyright 2022 Kristo Mela
# -*- coding: utf-8 -*-
""" Roof truss class
    Classes:
        RoofTruss -- subclass of Truss (see truss.py)

    Written by Kristo Mela
    18.4.2017
"""

import numpy as np
import math

from structures.truss import Truss, TrussMember, TrussNode
#from truss import TrussMember
#import truss as tr
from sections.steel.RHS import SHS

top_node_label = "TOP"
bottom_node_label = "BOTTOM"
    
top_member_label = "TC"
bottom_member_label = "BC"
brace_label = "BR"

class RoofTruss(Truss):
    """ Class for single-span steel roof trusses
        Subclass of 'Truss'
        Uses classes:
        'Cross-Section' (and its subclasses)
        'RHSJoint' (and its subclasses)
        'FrameFEM'
        'SteelMember'
        
        Attributes:
            span -- truss span
            h1   -- height from bottom chord to support
            slope -- top chord slope (e.g. slope = 1/20)
            braces -- list of brace members
            top_chord -- list of top chord members
            bottom_chord -- list of bottom chord members
            joints -- list of joints
            column -- column profile (optional)
            fyt -- top chord yield strength
            fyb -- bottom chord yield strength
            fy -- brace yield strength
        
        Methods:
            top_chord_angle -- determine top chord angle
            nbraces -- number of braces
            add_roof_member -- adds roof truss member
            roof_member_coord -- get roof member coordinates
            roof_truss_nodes -- create roof truss nodes
            top_chord_direction -- top chord direction vector
            brace_angle_top_chord -- angle between brace and top chord
            brace_angle_bottom_chord
    """
    
    top_node_label = "top"
    bottom_node_label = "bottom"
    
    top_member_label = "top"
    bottom_member_label = "bottom"
    brace_label = "brace"
    
    def __init__(self,span,h1,slope,topology="K",divX=4):
        """ Constructor
            
            NOTE: it might be better to leave creation of members and nodes
            out of the constructor in order to better allow free topologies
            
             span .. span of the truss (user defines units)
             h1 .. height from bottom chord to support
             slope .. slope of the top chord (e.g. slope = 1/20).
             members .. dictionary of truss members. The key:value pair is
                 member_id:member_object, where member_id is a string with
                 the unique identifier of the member and member_object is the
                 actual member object.
                 
                 the member_id identifies a member as a top chord member,
                 bottom chord member or a brace.
            
        """
        Truss.__init__(self)
        self.span = span
        self.h1 = h1
        # self.h2 = h2
        self.slope = slope
        
        
        self.members = {}
        
        self.braces = []
        self.top_chord = Chord()
        self.bottom_chord = Chord()
        self.joints = []
        self.column = 0.0
        self.fyt = 355.0
        self.fyb = 355.0
        self.fy = 355.0
        
        self.fem = []
    
    #THIS STUFF CAN BE LEFT OUT OF THE CONSTRUCTOR!
        # Create nodal coordinates
        self.roof_truss_nodes(divX,topology)
        
        # RETHINK these attributes!
        # Truss has 'nodes' now instead of 'x' and 'xmod'
        # 'nodes' might also include 'xmod' information!        
        self.roof_truss_members(topology)

        #self.symmetry_nodes()              
        #self.symmetry_pairs()        
        
        #self.determine_joints()
        
    
    def top_chord_angle(self):
        """ Angle of the top chord with respect to the horziontal axis. """    
        a = math.degrees(math.atan(self.slope))
        return a
   
    def nbraces(self):
        """ Number of brace members """
        NB = len(self.braces)
        return NB
    
    def add_node(self,x,y,chord):
        """ Create node and add it to the list of nodes
            
            chord -- "top" or "bottom"
        """
        newNode = RoofTrussNode(x,y,chord)
        self.nodes.append(newNode)
        
        if chord == "top":
            self.top_chord.add_node(newNode)
        else:
            self.bottom_chord.add_node(newNode)    
    
    def add_member(self,n1,n2,profile,member_type):
        """  Add new member 
        
            Calls the AddMember method of class 'Truss', with additional
            information on member type
        """        
        if member_type == "brace":
            mNew = Brace(n1,n2,profile)
            self.braces.append(mNew)
        elif member_type == "top":
            mNew = TrussMember([n1,n2],profile)
            self.top_chord.members.append(mNew)
        else:
            mNew = TrussMember([n1,n2],profile)
            self.bottom_chord.members.append(mNew)
    
        self.members.append(mNew)
    
    def roof_truss_nodes(self,n,topology):
        """ Create roof truss nodes for given topology and division
            n -- number of divisions on one half of the truss
            topology -- one of the standard topologies, 'K', 'N', or KT
        """
        L = self.span
        h = self.h1
        k = self.slope
        x = 0.0
        dx = 0.5*L/n
        
        """ First (left-most) supported node """
        self.add_node(0.0,0.0,top_node_label)
        
        if topology == "N":            
            for i in range(n):
                x = x+dx
                self.add_node(x,-h,bottom_node_label)
                self.add_node(x,k*x,top_node_label)
            
            for i in range(n-1):
                x = x+dx
                self.add_node(x,-h,self.bottom_node_label)
                self.add_node(x,-k*(x-L),self.top_node_label)

            self.add_node(L,0,top_node_label)
            
        elif topology == "K":
            for i in range(n):
                x = x+dx
                self.add_node(x-0.5*dx,-h,bottom_node_label)
                self.add_node(x,k*x,top_node_label)
            
            for i in range(n-1):
                x = x+dx
                self.add_node(x-0.5*dx,-h,bottom_node_label)
                self.add_node(x,-k*(x-L),top_node_label)

            self.add_node(L-0.5*dx,-h,bottom_node_label)
            self.add_node(L,0,top_node_label)
                
        elif topology == "KT":
            for i in range(n):
                x = x+dx
                self.add_node(x-0.5*dx,k*(x-0.5*dx),top_node_label)
                self.add_node(x-0.5*dx,-h,bottom_node_label)
                self.add_node(x,k*x,top_node_label)
            
            for i in range(n-1):
                x = x+dx
                self.add_node(x-0.5*dx,-k*(x-0.5*dx-L),top_node_label)
                self.add_node(x-0.5*dx,-h,bottom_node_label)
                self.add_node(x,-k*(x-L),top_node_label)

            self.add_node(L-0.5*dx,-k*(-0.5*dx),top_node_label)
            self.add_node(L-0.5*dx,-h,bottom_node_label)
            self.add_node(L,0,top_node_label)
    
    def roof_truss_members(self,topology):
        """ Creates members based on given nodal coordinates """        

        """ Default profiles """
        defTop = SHS(150,6)
        defBottom = SHS(140,5)
        defBrace = SHS(60,4)

        N = self.nnodes()

        if topology == "N":
            half_nodes = int(N/2)
            for i in range(half_nodes):
                if self.nodes[i].chord == top_node_label:
                    """ node is located at the top chord """
                    n1 = i+2
                    self.add_member(i,n1,defTop,top_member_label)
                    n2 = i+1
                    self.add_member(i,n2,defBrace,brace_label)
                else:
                    """ node is located at the bottom chord """
                    n1 = i+1
                    self.add_member(i,n1,defBrace,brace_label)
                    n2 = i+2
                    self.add_member(i,n2,defBottom,bottom_member_label)

            self.add_member(half_nodes-1,half_nodes+2,defBottom,bottom_member_label)

            for i in range(int(N/2),int(N-3)):
                if self.nodes[i].chord == top_node_label:
                    # node is located  at the top chord                                
                    n1 = i+2
                    self.add_member(i,n1,defTop,top_member_label)
                else:
                    # node is located at the bottom chord        
                    n1 = i+1
                    self.add_member(i,n1,defBrace,brace_label)
                    n2 = i+2
                    self.add_member(i,n2,defBottom,bottom_member_label)
                    n3 = i+3
                    self.add_member(i,n3,defBrace,brace_label)

            # Final two nodes:
            # 1) Last node of the bottom chord
            n1 = N-2
            n2 = N-1
            self.add_member(N-3,n1,defBrace,brace_label)
            self.add_member(N-3,n2,defBrace,brace_label)

            # 2) Last node of the top chord
            n1 = N-1
            self.add_member(N-2,n1,defTop,top_member_label)
        elif topology == "K":
            for i in range(N-2):
                if self.nodes[i].chord == top_node_label:
                    # node is located  at the top chord        
                    n1 = i+2
                    n2 = i+1
                    self.add_member(i,n1,defTop,top_member_label)
                    self.add_member(i,n2,defBrace,brace_label)
                else:
                    # node is located at the bottom chord
                    n1 = i+1
                    self.add_member(i,n1,defBrace,brace_label)
                    n2 = i+2
                    self.add_member(i,n2,defBottom,bottom_member_label)
        
            self.add_member(N-2,N-1,defBrace,bottom_member_label)
        
            
        elif topology == "KT":
            for i in range(N-1):
                if self.nodes[i].type == top_node_label:
                    self.add_member(i,i+1,defTop,top_member_label)
                    self.add_member(i,i+2,defTop,top_member_label)
                else:
                    """ bottom chord node """
                    self.add_member(i,i+1,defBrace,brace_label)
                    if i+3 < N-1:
                        self.add_member(i,i+3,defBottom,bottom_member_label)
                            

    def top_chord_direction(self,up_slope=True):
        """ Direction vector of the top chord:
            up_slope -- slope going upwards
        """
        a = self.top_chord_angle()
        
        ca = math.cos(math.radians(a))
        sa = math.sin(math.radians(a))
        
        if up_slope == True:
            v = np.array([ca,sa,0])
        else:
            v = np.array([ca,-sa,0])
        
        return v
    
    def brace_angle_top_chord(self,m):
        """ Calculate angle betweena brace and top chord """
        
        # Get the direction vector of member m
        #v = self.member_direction_vector(m)
        v = self.members[m].direction_vector()

        #X = self.member_coord(m)
        X = self.members[m].coord()
        aTop = self.top_chord_angle()
            
        """ if the brace is located on the left hand side
            of the ridge, then the direction vector of the
            top chord has a positive y-component.
        """
        if max(X[:,0]) <= 0.5*self.span:
            t = self.top_chord_direction(True)
        else:
            t = self.top_chord_direction(False)
                
        dvt = np.dot(v,t)
        if dvt < 0:
            a = 180-math.degrees(math.acos(dvt))
        else:
            a = math.degrees(math.acos(dvt))

        return a

    def brace_angle_bottom_chord(self,m):
        """ Determine the angle between brace and the bottom chord """

        #v = self.member_direction_vector(m)
        v = self.members[m].direction_vector()
        """ angle of the bottom chord with x-axis.
            it is assumed that the bottom chord is horizontal.
        """

        aBottom = 0
            
        """ direction vector of the bottom chord """
        t = np.array([math.cos(aBottom),math.sin(aBottom),0])
                
        dvt = np.dot(v,t)
        if dvt < 0:
            a = 180-math.degress(math.acos(dvt))
        else:
            a = math.degrees(math.acos(dvt))

        return a
"""
    %% Properties
    properties
        
        span    % span
        h1      % height from bottom chord to support
        h2      % height from bottom chord to top chord
        slope   % slope of the top chord
        UChord  % top chord members
        LChord  % bottom chord members
        Braces  % braces
        Joints  % joints, of class 'RHSJoint' (and its subclasses)
        Column = SHS(200,8) % profile of the column.
        fyt = 355 % top chord strength [MPa]
        fyb = 355 % bottom chord strength [MPa]
        fy = 355 % brace strength [MPa]
        
        % uniform loads on the chords
        TopChordLoad = 0;
        BottomChordLoad = 0;
        
        xmod % modified nodal coordinates.
        symmNode % pairs of symmetric nodes
        
        fem % finite element model of the truss
        jointFEMNodes % nodes of the fem model located at the joints
        
    end
    
    %% Methods
    
        % moves node 'n'.
        % 'value' is a 2x1 vector with new location
        % note: n may be a vector in which case 'value' must be a matrix
        def self = MoveNode(self,n,value)            
            self.xmod(:,n) = value;
            nsymm = [];
            
            for i = 1:length(n)
                [row,col] = find(self.symmNode == n(i));
                if ~isempty(row)
                    nnew = self.symmNode(row,setdiff(1:2,col));
                    self.xmod(:,nnew) =  value(:,i);
                    if value(1,i) < 0.5*self.span
                        self.xmod(1,nnew) = self.span-self.xmod(1,nnew);
                    else
                        self.xmod(1,nnew) = self.span-self.xmod(1,nnew);
                    end
                    nsymm = [nsymm,nnew];
                end
            end
            
            n = [n,nsymm];
            modNodes = length(n);
            for r = 1:modNodes
                % update the angles of the members
                % first, detect which braces are connected to the node:
                b = self.Joints{n(r)}.Braces;
                % then, compute the new angles of the braces:            
                for i = 1:length(b)
                    % find the nodes of member b(i)
                    nodes = self.members(b(i)).n;
                    for j = 1:2                    
                        chord = self.ChordNode(nodes(j));
                        % Determine the angle between the brace and
                        % the chord corresponding to node 'n'
                        switch chord
                            case 'top'
                                a(i,1) = self.TopAngle(b(i));
                            case 'bottom'    
                                a(i,1) = self.BottomAngle(b(i));
                        end
                        % determine the index of brace b(i) in joint n
                        nd = find(self.Joints{nodes(j)}.Braces == b(i));
                        self.Joints{nodes(j)}.Joint.theta(nd) = a(i);
                        B = self.Joints{nodes(j)}.Braces;
                        switch length(B)
                            case 2                            
                                % Update gap
                                %self.GapBetweenBraces(B(1),B(2),chord);
                                self.Joints{nodes(j)}.Joint.gap = self.GapBetweenBraces(B(1),B(2),chord);
                            case 3
                                % Update gaps
                                for k = 1:2                                
                                    %self.GapBetweenBraces(B(k),B(3),chord);
                                    self.Joints{nodes(j)}.Joint.gap(k) = self.GapBetweenBraces(B(k),B(3),chord);
                                end
                        end
                    end
                end
            end
%             if length(b) > 1
%                 self.Joints{n}.Joint.gap = self.GapBetweenBraces(b(1),b(2),chord);
%             end
            
        end
        
        
        def self = setTopChordProfile(self,p)            
            % set the corresponding profile to all 
            for i = 1:length(self.UChord)
                self.members(self.UChord(i)).profile = p;             
            end
            nTop = self.TopChordNodes;
            for i = 1:length(nTop)
                self.Joints{nTop(i)}.Joint.chord = p;
                switch length(self.Joints{nTop(i)}.Braces)
                    case 2
                    % For K- and N-joints determine the gaps again:
                        self.Joints{nTop(i)}.Joint.gap = self.GapBetweenBraces(self.Joints{nTop(i)}.Braces(1),self.Joints{nTop(i)}.Braces(2),'top');
                    case 3
                        for j = 1:2
                            self.Joints{nTop(i)}.Joint.gap(j) = self.GapBetweenBraces(self.Joints{nTop(i)}.Braces(j),self.Joints{nTop(i)}.Braces(3),'top');
                        end
                end
            end
        end
        
        % Set the bottom chord profile to 'p'
        def self = setBottomChordProfile(self,p)
            for i = 1:length(self.LChord)
                self.members(self.LChord(i)).profile = p;                
            end
            nBottom = self.BottomChordNodes;
            for i = 1:length(nBottom)
                self.Joints{nBottom(i)}.Joint.chord = p;                
                switch length(self.Joints{nBottom(i)}.Braces)
                    case 2
                        % For K- N- and KT-joints, determine the gaps again:
                        self.Joints{nBottom(i)}.Joint.gap = self.GapBetweenBraces(self.Joints{nBottom(i)}.Braces(1),self.Joints{nBottom(i)}.Braces(2),'bottom');
                    case 3
                        for j = 1:2
                            self.Joints{nBottom(i)}.Joint.gap(j) = self.GapBetweenBraces(self.Joints{nBottom(i)}.Braces(j),self.Joints{nBottom(i)}.Braces(3),'bottom');
                        end                
                end
            end
        end
        
        % get the profile of the chord to which node 'n' belongs to.
        def [p,chord] = ChordProfile(self,n)                
            chord = self.ChordNode(n);
            switch self.ChordNode(n)
                case 'top'
                    %p = self.members(self.UChord(1)).profile;
                    p = self.TopChordProfile;
                case 'bottom'
                    %p = self.members(self.LChord(1)).profile;
                    p = self.BottomChordProfile;
            end
        end
        
        % set profile 'p' to braces in vector 'b'
        % this leads to an update of the gaps to the corresponding
        % K and N joints.
        def self = BraceProfile(self,p,b)
            nb = length(b);
            for i = 1:nb
                % set for the member
                self.members(b(i)).profile = p;
                % the profile must also be set for the
                % joints that the member is a part of
                n = self.members(b(i)).n;
                for j = 1:2
                    %n(j)
                    %length(self.Joints{n(j)}.Braces)
                    nbraces = length(self.Joints{n(j)}.Braces);
                    if nbraces == 1                        
                            self.Joints{n(j)}.Joint.braces = p;
                    else                        
                        % determine which index is used for brace b(i)                        
                        nd = find(self.Joints{n(j)}.Braces == b(i));                           
                        % set profile of the corresponding brace                        
                        self.Joints{n(j)}.Joint.braces{nd} = p;                        
                        % determine chord of the node ('top' or 'bottom')                            
                        chord = self.ChordNode(n(j));
                        if nbraces == 2
                            self.Joints{n(j)}.Joint.gap = self.GapBetweenBraces(self.Joints{n(j)}.Braces(1),self.Joints{n(j)}.Braces(2),chord);
                        else
                            for k = 1:2
                                self.Joints{n(j)}.Joint.gap(k) = self.GapBetweenBraces(self.Joints{n(j)}.Braces(k),self.Joints{n(j)}.Braces(3),chord);
                            end
                        end
                    end
                end
            end
        end                               
        
        
        % Determines the symmetry pairs of roof truss members
        % Symmetry axis is located at midspan
        def self = SymmetryPairs(self)
            E = [self.members.n];
            S = zeros(self.nE,2);
            xS = 0.5*self.span;           
            for i = 1:self.nE               
                S(i,:) = sort(abs(self.x(1,E(:,i))-xS));
            end
            for i = 1:self.nE
                if isempty(self.members(i).pair)
                    for j = i+1:self.nE
                        if norm(S(i,:)-S(j,:),2) < 1e-3
                            self.members(i).pair = j;
                            self.members(j).pair = i;
                        end
                    end
                end
            end
        end
        
        % Find pairs of nodes that lie symmetric with respect to
        % the mid span of the truss
        def self = SymmetryNodes(self)
            S = [];
            xS = 0.5*self.span;
            n1 = find(self.x(1,:) < xS);
            n2 = find(self.x(1,:) > xS);
            S1 = abs(self.x(1,n1)-xS);
            S2 = abs(self.x(1,n2)-xS);
            for i = 1:length(n1)
                npair = find(abs(S2-S1(n1(i)))<1e-3);
                switch length(npair)
                    case 2
                    %n1(i)
                    %n2(npair(abs(self.x(2,n2(npair))-self.x(2,n1(i))) < 1e-3))
                        S = [S;n1(i),n2(npair(abs(self.x(2,n2(npair))-self.x(2,n1(i))) < 1e-3))];
                    case 1                
                        S = [S;[n1(i),n2(npair)]];                        
                end
            end
            self.symmNode = S;
        end
        
        % Determine eccentricity at node 'n'
        def e = Eccentricity(self,n)
            chord = self.ChordNode(n);
            switch chord
                case 'top'                    
                    v = self.TopDirection(n);                    
                case 'bottom'                    
                    v = [1;0];                    
            end
            
            e = self.Point2Line(self.xmod(:,n),self.x(:,n),v);

            switch chord
                case 'top'
                    if self.xmod(2,n) < self.x(2,n)
                        e = -e;
                    end
                case 'bottom'
                    if self.xmod(2,n) > self.x(2,n)
                        e = -e;
                    end
            end            
        end
        
 
        
    %% Methods in separate files
    methods
        % Calculate points where braces meet the surface of the chords
        [xChord,N,xProj] = ChordSurfacePoint(self,m,verb)
        % Determine joint types and generate joint selfects
        self = DetermineJoints(self);
        % Calculate gaps between braces
        g = GapBetweenBraces(self,m1,m2,chord);
        % Determine node locations such that eccentricity is minimized
        self = MinimizeEccentricity(self,x0);
        % for analytical expressions
        [self,P,xsol] = MinimizeEccentricity2(self,x0);
        % Creates a 'FrameFEM' model of the truss
        [fem,self,XndCenter] = CreateFEModel(self);
        % Writes a str file to be read in Robot Structural analysis
        WriteRobotFile(self,filename);
        % Perform member design according to EN 1993        
        [self,UtRatio] = MemberDesign(self,verb,allClass1or2)
    end
    
    %% Static methods
    methods (Static = true)
        def R = RotationMatrix(a)
        %
        % 2D rotation matrix
        % 
        % 28.12.2016 -k.m.
        %

        ca = cosd(a);
        sa = sind(a);

        R = [ca, -sa; sa, ca];
        end
        
        % distance of point r1 from line passing through point r0
        % with direction vector v
        def d = Point2Line(r1,r0,v)
           if length(r1) == 2
               r1 = [r1;0];
               r0 = [r0;0];
               v = [v;0];
           end
           d = norm(cross(r1-r0,v),2)/norm(v,2);
        end
    end
end

"""

class Chord:
    """ Class for storing chord data and manipulating things """
    def __init__(self):
        self.members = []
        self.nodes = []
        self.profile = None
        self.fy = 355
        self.slope = 0
    def angle(self,position):
        """ chord angle """
        
    def add_node(self,newNode):
        self.nodes.append(newNode)
        
    def add_member(self,newMember):
        """ Add member to chord """
        self.members.append(newMember)

    def set_profile(self,p):
        """ Sets a new profile for all members """
        self.profile = p

        for m in self.members:
            m.profile = p

    def set_fy(self,fy):
        """ Sets a new yield strength for all members """
        self.profile = fy
        
        for m in self.members:
            m.profile.fy = fy

class Brace(TrussMember):
    """ Class for braces 
        
        Subclass of TrussMember of 'truss' module
        The main difference is that the coordinates are taken from the
        variable 'xmod' instead of 'x'
    """

    def __init__(self,n1,n2,p):
        TrussMember.__init__(self,[n1,n2],p)

    def coord(self):
        """ Nodal coordinates of brace member
        
        """
        X = np.array([self.nodes[0].xmod,self.nodes[1].xmod])
        return X

    def is_vertical(self):
        res = False
        X = self.coord()
        if abs(X[0,0]-X[1,0]) < 1e-6:
            res = True

        return res

class RoofTrussNode(TrussNode):
    
    def __init__(self,x,y,chord):
        """ RoofTrussNode
            
            Subclass of TrussNode
            x, y -- coordinates
            chord -- "top" or "bottom", denotes the chord
        """
        TrussNode.__init__(self,x,y)
        
        self.xmod = np.array([x,y,0])
        self.chord = chord

""" For testing roof trusses """
if __name__ == "__main__":
    import sys
    topology = sys.argv[1]
    divX = int(sys.argv[2])
    t = RoofTruss(24000,2200,1/20,topology,divX)
    print("Number of nodes:",t.nnodes())
    print("Number of members:",t.nmembers())
    print("Slope (degrees):",t.top_chord_angle())
    
    nodes = [0,2,3,4]
    for n in nodes:
        print(t.nodes[n].x)
        print(t.nodes[n].xmod)
        print(t.nodes[n].chord)
    """
    braces = [2,5]
    for b in braces:
        print("Angle of member", b, "with top chord:",t.brace_angle_top_chord(b))
        print("Angle of member", b, "with bottom chord:",t.brace_angle_bottom_chord(b))
    """
    t.draw()
    
class Joint:
    """ Class for truss joints 
        properties:
            init_node .. initial node, where the center lines of the
                         braces meet at the center line of the chord
            eccentricity_node .. if the structural model contains a single eccentricity element,
                                the eccentricity node is stored here
            chord_face_nodes .. if the braces in the structural model end at
                                chord surfaces, the nodes are placed here
            eccentricity_elements .. element(s) from chord center line to eccentricity node
                                     or chord face
            gap_member .. member at the gap
            braces .. brace members connected to the joint
            chord_members .. chord members connected to the joint
    
    """
    
    def __init__(self,init_node):
        
        self.init_node = init_node
        self.eccentricity_node = None
        
        self.eccentricity_elements = None
        self.chord_face_nodes = None
        self.gap_member = None

        self.braces = []
        self.chord_members = []
    
    def add_brace(self,new_brace):
        self.braces.append(new_brace)
    
    def add_chord_member(self,new_chord_member):
        self.chord_members.append(new_chord_member)
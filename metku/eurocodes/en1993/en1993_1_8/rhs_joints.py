# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:11:28 2018

EN 1993-1-8 Rules for connections

@author: kmela
"""

from copy import deepcopy
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import numpy as np

try:
    from metku.eurocodes.en1993.constants import gammaM5, E
except:
    from eurocodes.en1993.constants import gammaM5, E

def eccentricity(h0,h1,h2,t1,t2,g):
    """ Eccentricity of a K gap joint """
    rt1 = math.radians(t1)
    rt2 = math.radians(t2)
    sin1 = math.sin(rt1)
    sin2 = math.sin(rt2)
    sin12 = math.sin(rt1+rt2)
    return (0.5*h1/sin1 + 0.5*h2/sin2 + g)*sin1*sin2/sin12 - 0.5*h0

def rotate(v,theta=0.0):
    """ Rotates vector 'v' by the angle 'theta'
        Works only in 2D
    """
    
    if len(v) != 2:
        print("Error: length of vector to be rotated must be 2.")
    else:
        trad = np.radians(theta)
        st = np.sin(trad)
        ct = np.cos(trad)
        R = np.array([[ct,-st],[st,ct]])
        
        v = np.asarray(v)
        vrot = R.dot(v)
    
    return vrot

class Line:
    """ Small class for lines in 2D """
    
    def __init__(self,**kwargs):
        """ Constructor
            Possible initialization calls:
                p1, p2: two points through which the line passes
                p1, v: point and direction vector
                
            v ..unit vector
                 
        """
        
        self.p1 = None
        self.p2 = None
        self.v = None
        
        for key, value in kwargs.items():
            if key == 'p1':
                self.p1 = np.asarray(value)
            elif key == 'p2':
                self.p2 = np.asarray(value)
            elif key == 'v':
                self.v = np.asarray(value)
    
        if self.p2 is None:
            self.p2 = self.p1 + self.v
        elif self.v is None:
            self.v = (self.p2-self.p1)/np.linalg.norm(self.p2-self.p1)
            
    def __call__(self,t):
        """ Calculates the coordinates for local coordinate t """
        
        return self.p1 + t*self.v
    
    def __repr__(self):
        
        return f"Line: p1 = ({self.p1[0]:.3f},{self.p1[1]:.3f}), v = ({self.v[0]:.3f},{self.v[1]:.3f})"
            
    def intersect(self,line):
        """ Finds intersection point between two lines
        """
        v = self.v
        p = self.p1
        u = line.v
        r = line.p1
        
        A = np.array([[v[0],-u[0]],[v[1],-u[1]]])
        b = r-p
                        
        t = np.linalg.solve(A,b)
        
        x_intersect = p + t[0]*v
        
        return t, x_intersect
    
    def distance(self,xR):
        """ Finds distance from a point with coordinates in xR
            Returns alsto the nearest point on the line
        """
        dx = self.p1-xR
        t_close = -self.v.dot(dx)
                
        
        x_closest = self.__call__(t_close)

        distance = np.linalg.norm(x_closest-xR)
        
        return distance, x_closest
    
class RHSJoint:
    """ Class for joints between RHS members """
    
    def __init__(self,chord_profile, N0=0.0, V0=0.0, M0=0.0):
        """ Constructor
            Input:
                
            :param: chord_profile .. RHS or SHS class member, sets the chord profile
            :param: N0 .. magnitude of the largest axial force in the chord around the joint.
            :param: V0 .. magnitude of the largest shear force in the chord around the joint.
            :param: M0 .. magnitude of the largest bending moment in the chord.
        
        
        """
        self.chord = chord_profile
        
        self.N0 = N0
        self.V0 = V0
        self.M0 = M0
        
        self.r = 0

        self.chord_face_plate = {'lp': 0, 'bp': 0,'tp': 0}
        self.chord_web_plate = {'lp': 0, 'bp': 0,'tp': 0}
             
    
    @property    
    def h0(self):
        return self.chord.H
        
    @property
    def b0(self):
        return self.chord.B
        
    @property
    def t0(self):
        return self.chord.T
    
    @property
    def r0(self):
        return self.chord.R
    
    @property
    def fy0(self):
        return self.chord.fy

    def gamma(self):
        return 0.5*self.b0/self.t0
        
    def b_straight(self):
        """ Length of the straight part of the chord profile """
        return self.b0-2*self.r0
        
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
            kn = min(1.3-0.4*n/b, 1.0)
        
        return kn
        
    def strength_reduction(self):
        """ Reduce the strength of joint for steel greater than S355
        
            NOTE! Strength reduction factor can be separately specified
            in attribute 'r'.
        """
        
        if self.r == 0:
            r = 1.0
            if self.fy0 <= 355:
                r = 1.0
            elif self.fy0 <= 460:            
                r = 0.9
            else:
                r = 0.8
        else:
            r = self.r
        return r
        
    def eval_N(self):
        """ Evaluate chord stress function """
        # s0Ed > 0 for compression
        s0Ed = self.chord_stress()
        n = s0Ed/self.fy0/gammaM5
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
        n = -s0Ed/self.fy0/gammaM5
        
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

    def info(self):
        """ Prints information on the joint """
        
        print("  Chord: {0:s}".format(self.chord.__repr__()))
        print("    N_0,Ed = {0:4.2f} kN".format(self.N0*1e-3))
        print("    M_0,Ed = {0:4.2f} kNm".format(self.M0*1e-6))
    
    def draw(self,fig=None,ax=None,chord_angle = 0, dx=500):
        """ Draw joint 
            For the main class, only the chord is plotted
            
            input:
                chord_angle .. angle of the chord and the horizontal axis (deg)
                dx .. length (mm) of the chord in each side of the origin
                
            origin is in the middle of the center line of the chord
        """
        
        grey = '0.8'
        if fig is None:
            fig, ax = plt.subplots(1)
        
        # y coordinate of the top face
        y = 0.5*self.h0
        yinner = y-self.t0
        
        ax.hlines([0],-dx,dx,'k','dashdot',zorder=2,linewidth=1)
        ax.hlines([-y,y],-dx,dx,'k','solid',zorder=2,linewidth=1.75)
        ax.hlines([-yinner,yinner],-dx,dx,grey,'solid',zorder=2,linewidth=1.25)
        
        
        ax.set_aspect('equal', 'datalim')
    
class RHSKGapJoint(RHSJoint):
    """ For K and N Gap joints between RHS members """
    
    def __init__(self,chord,braces,angles,gap=20, **kwargs):
        """ Constructor
            input:
                chord .. chord profile (RHS or SHS class object)
                braces .. list of brace profiles (RHS or SHS class objects)
                angles .. list of angles (degrees) between the braces and the chord
                gap .. gap between the braces
        """
        RHSJoint.__init__(self,chord, **kwargs)
                
        self.braces = braces
        self.gap = gap
        self.angles = angles
    
    @property
    def NEd(self):
        """ Numpy array of axial forces in braces """
        return np.array([self.braces[0].Ned,self.braces[1].Ned])
    
    @property
    def N1(self):
        """ Axial force in the compression member """
        if self.braces[0].Ned <= 0.0:
            return self.braces[0].Ned
        else:
            return self.braces[1].Ned
        
    @property
    def angle1(self):
        """ Angle of the compression member """
        if self.braces[0].Ned <= 0.0:
            return self.angles[0]
        else:
            return self.angles[1]
    
    @property
    def N2(self):
        """ Axial force in the tension member """
        if self.braces[0].Ned <= 0.0:
            return self.braces[1].Ned
        else:
            return self.braces[0].Ned
    
    @property
    def angle2(self):
        """ Angle of the tension member """
        if self.braces[0].Ned <= 0.0:
            return self.angles[1]
        else:
            return self.angles[0]
        
    @property
    def N0gap(self):
        """ Axial force in the gap 
            NOTE: The sign convention is that N0gap < 0 when there is compression
            in the gap and N0gap > 0 for tension.
            
            The expressions are obtained from horizontal equilibrium on the
            N0 side of the joint. When the chord is in compression, 
            N0 is on the same side as the tension brace. When the chord is in tension,
            N0 is on the same side as the compression brace.
        """
        if self.N0 < 0.0:
            # Chord in compression
            N0gap = self.N0+self.N2*np.cos(np.radians(self.angle2))
        else: 
            # Chord in tension            
            N0gap = self.N0+self.N1*np.cos(np.radians(self.angle1))
    
        return N0gap
    
    @property
    def V0gap(self):
        """ Shear force in the gap 
        
            The shear force is in the gap is obtained as the vertical component
            of either of the braces (see vertical equilibrium of the joint).
        """
        
        return abs(self.N1*np.sin(np.radians(self.angle1)))
    
    @property
    def h(self):
        """ Brace heights """
        return np.asarray([self.braces[0].H,self.braces[1].H])
    
    @property
    def b(self):
        """ Brace widths """
        return np.asarray([self.braces[0].B,self.braces[1].B])
    
    @property
    def t(self):
        """ Brace thickness """
        return np.asarray([self.braces[0].T,self.braces[1].T])
    
    @property
    def fy(self):
        """ Brace yield strength """
        return np.asarray([self.braces[0].fy,self.braces[1].fy])
    
    def beta(self):
        return (sum(self.b)+ sum(self.h)) / (4*self.b0)
    
    
    def eccentricity(self):
        t = self.angles
        h = self.h
        e0 = 0.5*h[0]/math.sin(math.radians(t[0])) + 0.5*h[1]/math.sin(math.radians(t[1])) + self.gap
        e = e0*math.sin(math.radians(t[0]))*math.sin(math.radians(t[1]))/math.sin(math.radians(sum(t)))-0.5*self.h0
        return e
    
    # Resistance formulas
    def chord_face_failure(self):
                        
        b = self.beta()
        g = self.gamma()                           
        kn = self.eval_KN()
        r = self.strength_reduction()
        fy0 = self.fy0
        if self.chord_face_plate['tp'] == 0:
            t0 = self.t0
        else:
            t0 = self.chord_face_plate['tp']
            
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
        fy0 = self.fy0
        # braces
        NiRd = fy0*Aeta/math.sqrt(3)/s/gammaM5
        # chord
        # Check first, if the chord can sustain the shear force
        if self.V0gap/self.chord.shear_force_resistance() > 1:
            # Chord profile is not sufficient, so set resistance to 0.0.
            #print("Chord is not strong enough for shear forces")
            N0Rd = 0.0

        else:
            N0Rd = fy0*((self.chord.A-Aeta)+Aeta*math.sqrt(1-(self.V0gap/self.chord.VRd)**2))/gammaM5
            
        r = self.strength_reduction()
        NiRd = r*NiRd
        N0Rd = r*N0Rd
        return NiRd, N0Rd
        
    def a_eta(self):
        a = math.sqrt(1/(1+4*self.gap**2/3/self.t0**2))
        return (2*self.h0+a*self.b0)*self.t0
        
    def brace_failure(self):
        r = self.strength_reduction()
        b = np.array(self.b)
        h = np.array(self.h)
        t = np.array(self.t)
        fy = np.array(self.fy)
        fy0 = self.fy0
        if self.chord_face_plate['tp'] == 0:
            t0 = self.t0
        else:
            t0 = self.chord_face_plate['tp']
        beff = np.minimum(10/(self.b0/t0)*fy0*t0/(fy*t),1)*b
        #beff_arr = np.array([min(x,beff[i]) for i, x in enumerate(b)])
            
        NiRd = r*self.fy*t*(2*h-4*t+b+beff)/gammaM5
        return NiRd
        
    def punching_shear(self):
        b = self.b
        fy0 = self.fy0
        if self.chord_face_plate['tp'] == 0:
            t0 = self.t0
        else:
            t0 = self.chord_face_plate['tp']
        bep = np.minimum(10/(self.b0/t0),1)*b 
        r = self.strength_reduction()
        
        s = np.sin(np.radians(np.array(self.angles)))
            
        NRd = r*fy0*t0/math.sqrt(3)/s*(2*self.h/s+b+bep)/gammaM5
        return NRd
    
    def add_chord_face_plate(self,lp=0,bp=0,tp=0):
        """
        Adds a reinforcement plate to the chord face

        Parameters
        ----------
        lp : float
            length of the plate.
        bp : float
            width of the plate.
        tp : float
            thickness of the plate.

        Returns
        -------
        None.

        """
        
        lp_min = 1.5*(self.braces[0].h/np.sin(np.radians(self.angles[0])) + self.gap + 
                                              self.braces[1].h/np.sin(np.radians(self.angles[1])))
        
        lp = max(lp,lp_min)
        
        bp = max(bp,self.b0-2*self.t0)

        tp = max(tp,2*max(self.t))
        
        self.chord_face_plate = {'lp':lp, 'bp':bp,'tp':tp}
        
        
    
    def validity(self,verb=False):
        """ Check the conditions of Table 7.8 regarding the geometry of the joint  
            TODO: i) Gap condition,
                 ii) h0/b0 and hi/bi condition
                iii) bi/ti and hi/ti condition
                    
        
        """
        
        biti_ratio = [True,True]
        hib0_ratio = [True,True]
        b0to_ratio = [True,True]
        gap_ratio = True
        
        b0t0 = self.b0/self.t0
        h0t0 = self.h0/self.t0
        bib0 = self.b/self.b0
        gap_b0 = self.gap/self.b0
        
        bib0_ratio = self.b/self.b0 > max(0.35,0.1+0.01*b0t0)
        
        
        b = self.b
        if self.braces[0].Ned < 0:
            b1 = b[0]
            b2 = b[1]
        else:
            b1 = b[1]
            b2 = b[0]
                            
        if verb:
            print("  Geometry conditions:")
            print("    a) h0/t0 = {0:4.2f} <= 35".format(h0t0))
            print("    b) b0/t0 = {0:4.2f} <= 35".format(b0t0))
            print("    c) bi/b0 > max(0.35, 0.1 + 0.01*b0/t0) = {0:4.2f}".format(max(0.35,0.1+0.01*b0t0)))
            for i in range(0,2):
                if bib0_ratio[i]:
                    OKNOTOK = "OK!"
                else:
                    OKNOTOK = "NOT OK!"
                print("  Brace {0:.0f}: bi/b0 = {1:4.2f}, {2:s}".format(i,bib0[i],OKNOTOK))
                
            print("  Further conditions for the use of Table 7.10 of EN 1993-1-8:")
            print(" 0.6 <= (b1+b2)/2b1 = {0:4.2f} <= 1.3".format((b1+b2)/2/b1))
            print(" b0/t0 = {0:4.2f} >= 15".format(b0t0))
    
    def info(self):
        """ Prints information on the joint """
        
        print("*** K gap Joint ***")
        
        super().info()
        
        for brace, angle in zip(self.braces,self.angles):
            print("  Brace: {0:s}".format(brace.__repr__()))
            if brace.Ned < 0:
                Nsense = "compression"
            elif brace.Ned > 0:
                Nsense = "tension"
            else:
                Nsense = "no axial force"
            print("  NEd: {0:4.2f} kN ({1:s})".format(brace.Ned*1e-3,Nsense))
            print("    Angle to chord: {0:4.2f} deg".format(angle))
        
        beta = self.beta()
        print("  Beta = {0:4.2f}".format(beta))
        print("  Gamma = {0:4.2f}".format(self.gamma()))        
        print("  Gap = {0:4.2f} mm".format(self.gap))
        print("    t1 + t2 = {0:4.2f} mm".format(sum([brace.t for brace in self.braces])))
        print("    0.5*b0*(1-beta) = {0:4.2f} mm".format(0.5*self.b0*(1-beta)))
        print("    1.5*b0*(1-beta) = {0:4.2f} mm".format(1.5*self.b0*(1-beta)))
        
        print("  Eccentricity = {0:4.2f} mm (Limit = {1:4.2f} mm)".format(self.eccentricity(),0.25*self.h0))
        
        self.validity(True)
        
        print(" -- RESISTANCE -- ")
        
        NRD_chord_face = self.chord_face_failure()
        
        print("    Chord face failure:")
        for i, NRd in enumerate(NRD_chord_face):
            print("      N_{0:.0f}_Rd = {1:4.2f} kN".format(i+1,NRd*1e-3))
            
        
        NiRd, N0Rd = self.chord_shear()
        
        print("    Chord shear failure:")
        for i, NRd in enumerate(NiRd):
            print("      N_{0:.0f}_Rd = {1:4.2f} kN".format(i+1,NRd*1e-3))

        print("      N_0_Rd = {0:4.2f} kN".format(N0Rd*1e-3))
        
        NiRd = self.brace_failure()
        
        print("    Brace failure:")
        for i, NRd in enumerate(NiRd):
            print("      N_{0:.0f}_Rd = {1:4.2f} kN".format(i+1,NRd*1e-3))
            
        NiRd = self.punching_shear()
        
        print("    Punching shear:")
        for i, NRd in enumerate(NiRd):
            print("      N_{0:.0f}_Rd = {1:4.2f} kN".format(i+1,NRd*1e-3))
            
    def draw(self,chord_type="bottom",chord_angle=0,dx=300,dy=250):
        """ Draw the joint 
            input:
                dy .. height of the braces above (below) the chord face
        """
            
        fig, ax = plt.subplots(1)
        
        # draw chord
        super().draw(fig,ax,chord_angle,dx)
        
        """ draw braces
            For the location of the braces, assume
            that the gap is divided equally between the two sides
            of the chord
            
            Assume brace 1 is on the left
        """
        WALL = 1.75
        INNER = 1.25
        CENTER = 1.5
                
        y0 = 0.5*self.h0
        dx0 = 0.5*self.gap
        
        xC = []
        vC = []
        
        for a, brace, K in zip(self.angles,self.braces,[-1,1]):
            arad = math.radians(a)
            
            # Draw the lines of the brace            
            x0 = np.array([K*dx0,y0])
            # Direction vector of the brace
            v = rotate([K,0],K*a)
            vC.append(deepcopy(-v))
            vortho = rotate(v,-K*90)
            
            # First line (inner edge)
            x1 = x0 + dy*v
            ax.plot([x0[0],x1[0]],[x0[1],x1[1]],linewidth=WALL,color='k')
            
            # Second pair of points
            sina = math.sin(arad)
            x0[0] += K * brace.t / sina
            x2 = x1 + brace.t * vortho
            ax.plot([x0[0],x2[0]],[x0[1],x2[1]],linewidth=INNER,color='0.8')
            
            # Center line
            dxin = 0.5*brace.h-brace.t
            x0[0] += K*dxin/sina
            x3 = x2 + dxin*vortho
            xC.append(deepcopy(x0))
            ax.plot([x0[0],x3[0]],[x0[1],x3[1]],linewidth=CENTER,color='k',linestyle='dashdot')
            
            # Second inner wall
            x0[0] += K*dxin/sina
            x4 = x3 + dxin*vortho
            ax.plot([x0[0],x4[0]],[x0[1],x4[1]],linewidth=INNER,color='0.8')
            
            # Second outer wall
            x0[0] += K * brace.t / sina
            x5 = x4 + brace.t / sina * vortho
            ax.plot([x0[0],x5[0]],[x0[1],x5[1]],linewidth=WALL,color='k')
            
            # Plot dashed line to the end of the brace
            ax.plot([x1[0],x5[0]],[x1[1],x5[1]],linewidth=INNER,color='k',linestyle='dashed')
        
        # Find the intersection of the center lines of the braces
        l1 = Line(p1=xC[0],v=vC[0])
        l2 = Line(p1=xC[1],v=vC[1])
        t, p = l1.intersect(l2)
        
        
        # Plot the point and the center lines
        ax.plot(p[0],p[1],'ko')
        for x in xC:
            ax.plot([x[0],p[0]],[x[1],p[1]],linewidth=CENTER,color='k',linestyle='dashdot')
        
        lchord = Line(p1=[0,0],v=[1,0])
        d, p0 = lchord.distance(p)
        
        ax.plot(p0[0],p0[1],'ko')
        ax.plot([p0[0],p[0]],[p0[1],p[1]],linewidth=INNER,color='k',linestyle='dashed')
    
    def design(self):
        """
        Designs the joint

        Returns
        -------
        Utilization ratios N1Ed/NjRd and N2Ed/NjRd

        """
        
        b = self.beta()
        g = self.gamma()
        
        NiEd = self.NEd
        
        N_cff_Rd = self.chord_face_failure()
        N_cs_Rd, N0Rd = self.chord_shear()
        N_bf_Rd = self.brace_failure()
        
        NjRd = np.minimum(N_cff_Rd,N_cs_Rd)
        NjRd = np.minimum(NjRd,N_bf_Rd)
        if b <= 1-1/g:
            N_ps_Rd = self.punching_shear()        
            NjRd = np.minimum(NjRd,N_ps_Rd)

        res = np.array(abs(self.NEd)/NjRd,abs(self.N0/N0Rd))

        return res

class RHSKTGapJoint(RHSJoint):
    """ For KT Gap joints between RHS members """
    
    def __init__(self,chord,braces,angles,gaps=[20,20], **kwargs):
        """ Constructor
            input:
                chord .. chord profile (RHS or SHS class object)
                braces .. list of brace profiles (RHS or SHS class objects)
                angles .. list of angles (degrees) between the braces and the chord
                gap .. gap between the braces
        """
        RHSJoint.__init__(self,chord, **kwargs)
                
        self.braces = braces
        self.gaps = gaps
        self.angles = angles

class RHSYJoint(RHSJoint):
    
    def __init__(self, chord_profile, brace, angle, **kwargs):
        super().__init__(chord_profile, **kwargs)
        self.brace = brace
        self.angle = angle
        
    @property
    def N1(self):
        """ Force in the brace """
        return self.brace.Ned
    
    @property
    def h(self):
        """ Brace heights """
        return self.brace.H
    
    @property
    def b(self):
        """ Brace widths """
        return self.brace.B
    
    @property
    def t(self):
        """ Brace thickness """
        return self.brace.T
    
    @property
    def fy(self):
        """ Brace yield strength """
        return self.brace.fy
    
    def beta(self):
        return self.b / self.b0
    
    @property
    def eta(self):
        return self.h/self.b0
    
    def fb(self):
        """ Yield strength to be used in chord side wall buckling. """
        
        if self.N1 < 0.0:
            alpha = self.chord.imp_factor[0]
            slend = self.slend()
            phi = 0.5*(1 + alpha*(slend - 0.2) + slend **2)
            chi = 1 / (phi + np.sqrt(phi**2 - slend**2))
            
            fb = chi * self.fy0
            #print(chi,fb)
        else:
            fb = self.fy0
        
        return fb
            
    def chord_face_failure(self):
                        
        b = self.beta()
        if b >= 1:
            # If beta is grater than 1, artificially reduce it to "nearly 1".
            b = 0.99
        kn = self.eval_KN()
        
        r = self.strength_reduction()
        fy0 = self.fy0
        t0 = self.t0
        NRd = r**kn*fy0*t0**2 / ((1-b)*np.sin(np.radians(self.angle)))
        NRd = NRd * (2*self.eta / np.sin(np.radians(self.angle)) + 4*np.sqrt(1-b))
        NRd = NRd/ gammaM5

        return NRd
    
    
    def slend(self):

        sl = 3.46*(self.h0/self.t0 -2)*np.sqrt(1/np.sin(np.radians(self.angle)))
        sl = sl / (np.pi * np.sqrt(E / self.fy0))
        return sl
   
    def chord_web_buckling(self):
        
        fb = self.fb()
        
        kn = self.eval_KN()
        
        r = self.strength_reduction()
        t0 = self.t0
        h = self.h
        NRd = r**kn*fb*t0 / (np.sin(np.radians(self.angle)))
        NRd = NRd * (2*h / np.sin(np.radians(self.angle)) + 10*t0)
        NRd = NRd / gammaM5
        return NRd
    
    def brace_failure(self):
        r = self.strength_reduction()
        b = self.b
        h = self.h
        t = self.t
        fy = self.fy
        fy0 = self.fy0
        t0 = self.t0
        beff = min(10/(self.b0/t0)*fy0*t0/(fy*t),1)*b
            
        NiRd = r*self.fy*t*(2*h-4*t+2*beff)/gammaM5
        return NiRd
    
    def punching_shear(self):
        b = self.b
        fy0 = self.fy0
        t0 = self.t0
        bep = min(10/(self.b0/t0),1.0)*b 
        r = self.strength_reduction()
        
        s = np.sin(np.radians(self.angle))
            
        NRd = r*fy0*t0/math.sqrt(3)/s*(2*self.h/s+2*bep)/gammaM5
        return NRd
    
    def validity(self,verb=False):
        """ Check the conditions of Table 7.8 regarding the geometry of the joint  
            
        """
        
        biti_ratio = True
        hib0_ratio = True
        b0to_ratio = True        
        
        b0t0 = self.b0/self.t0
        h0t0 = self.h0/self.t0
        bib0 = self.b/self.b0        
        
        bib0_ratio = self.b/self.b0 > 0.25
                
        b = self.b                
        if verb:
            
            if bib0_ratio:
                OKNOTOK = "OK!"
            else:
                OKNOTOK = "NOT OK!"
            
            print("  Geometry conditions:")
            print("    a) h0/t0 = {0:4.2f} <= 35".format(h0t0))
            print("    b) b0/t0 = {0:4.2f} <= 35".format(b0t0))
            print("    c) bi/b0 = {0:4.2f} > 0.25, {1:s}".format(bib0,OKNOTOK))            
                        
                
            print("  Further conditions for the use of Table 7.10 of EN 1993-1-8:")
            print(" bi/b0 = {0:4.2f} <= 0.85".format(bib0))
            print(" b0/t0 = {0:4.2f} >= 10".format(b0t0))
    
    def info(self):
        """ Prints information on the joint """
        
        print("*** Y Joint ***")
        
        super().info()
        
        print("  Brace: {0:s}".format(self.brace.__repr__()))
        if self.brace.Ned < 0:
            Nsense = "compression"
        elif self.brace.Ned > 0:
            Nsense = "tension"
        else:
            Nsense = "no axial force"
        print("  NEd: {0:4.2f} kN ({1:s})".format(self.brace.Ned*1e-3,Nsense))
        print("    Angle to chord: {0:4.2f} deg".format(self.angle))
        
        beta = self.beta()
        print("  Beta = {0:4.2f}".format(beta))
        print("  Gamma = {0:4.2f}".format(self.gamma()))                        
        
        self.validity(True)
        
        print(" -- RESISTANCE -- ")
        
        NRD_chord_face = self.chord_face_failure()
        
        print("    Chord face failure:")
        print("      N_1_Rd = {0:4.2f} kN".format(NRD_chord_face*1e-3))
        
        NiRd = self.chord_web_buckling()
        
        print("    Chord side wall buckling:")        
        print("      N_1_Rd = {0:4.2f} kN".format(NiRd*1e-3))        
        
        NiRd = self.brace_failure()
        
        print("    Brace failure:")        
        print("      N_1_Rd = {0:4.2f} kN".format(NiRd*1e-3))
            
        NiRd = self.punching_shear()
        
        print("    Punching shear:")        
        print("      N_1_Rd = {0:4.2f} kN".format(NiRd*1e-3))
    
    def design(self):
        """
        Designs the joint

        Returns
        -------
        Utilization ratio N1Ed/NjRd.

        """
        
        b = self.beta()
        g = self.gamma()
        
        
        if b < 0.85:
            NjRd = self.chord_face_failure()
        
        else:
            N_cff_Rd = self.chord_face_failure()
            N_csw_Rd = self.chord_web_buckling()
            N_bf_Rd = self.brace_failure()
            
            # Linear interpolation between chord face failure at b = 0.85
            # and chord side wall failure at b = 1.0.
            # The denominator is 1-0.85 = 0.15.            
            N_interp_Rd = N_cff_Rd + (b-0.85)/0.15*(N_csw_Rd-N_cff_Rd)
            
            NjRd = min(N_interp_Rd,N_bf_Rd)
            
            if b <= 1-1/g:
                N_ps_Rd = self.punching_shear()                
                NjRd = min(NjRd,N_ps_Rd)

        res = abs(self.N1)/NjRd        

        return res
            
            
    
class RHSXJoint(RHSYJoint):
    
    def __init__(self, chord_profile, brace, angle, **kwargs):
        super().__init__(chord_profile, brace, angle, **kwargs)
    
    def fb(self):
        """ Yield strength to be used in chord side wall buckling. """        
                
                
        if self.N1 < 0.0:            
            alpha = self.chord.imp_factor[0]
            slend = self.slend()
            phi = 0.5*(1 + alpha*(slend - 0.2) + slend **2)
            chi = 1 / (phi + np.sqrt(phi**2 - slend**2))
            
            res = 0.8*chi * self.fy0*np.sin(np.radians(self.angle))            
        else:
            res = self.fy0
        
        return res
    
    def chord_web_buckling(self):
        
        # Using the chord web buckling of the Y joint class, but
        # with fb calculated for X joints.
        NRd = super().chord_web_buckling()
        
        if np.cos(np.radians(self.angle)) > self.h/self.h0:
            # Take the minimium of chord web buckling for X joints
            # and chord shear for K and N joints
            # 
            Nshear = self.fy0*self.chord.Ashear/np.sqrt(3)/np.sin(np.radians(self.angle))
            # The strength reduction factor is already taken into account
            # for the chord web buckling of the Y joint class method.
            # It must be applied for the chord shear.
            r = self.strength_reduction()
            NRd = min(NRd,r*Nshear)
        
        NRd = NRd / gammaM5
        return NRd
        
    def info(self):
        """ Prints information on the joint """
        
        print("*** X Joint ***")
        
        RHSJoint.info(self)
        
        print("  Brace: {0:s}".format(self.brace.__repr__()))
        if self.brace.Ned < 0:
            Nsense = "compression"
        elif self.brace.Ned > 0:
            Nsense = "tension"
        else:
            Nsense = "no axial force"
        print("  NEd: {0:4.2f} kN ({1:s})".format(self.brace.Ned*1e-3,Nsense))
        print("    Angle to chord: {0:4.2f} deg".format(self.angle))
        
        beta = self.beta()
        print("  Beta = {0:4.2f}".format(beta))
        print("  Gamma = {0:4.2f}".format(self.gamma()))                        
        
        self.validity(True)
        
        print(" -- RESISTANCE -- ")
        
        NRD_chord_face = self.chord_face_failure()
        
        print("    Chord face failure:")
        print("      N_1_Rd = {0:4.2f} kN".format(NRD_chord_face*1e-3))
        
        NiRd = self.chord_web_buckling()
        
        print("    Chord side wall buckling:")        
        print("      N_1_Rd = {0:4.2f} kN".format(NiRd*1e-3))        
        
        NiRd = self.brace_failure()
        
        print("    Brace failure:")        
        print("      N_1_Rd = {0:4.2f} kN".format(NiRd*1e-3))
            
        NiRd = self.punching_shear()
        
        print("    Punching shear:")        
        print("      N_1_Rd = {0:4.2f} kN".format(NiRd*1e-3))

def Y_example():
    
    
    chord = SHS(140,6,fy=700)

    brace = SHS(140,8,fy=700)
    brace.Ned = -50e3
    
    
    YJoint = RHSYJoint(chord,brace,45,N0=-887.5e3)

    YJoint.r = 1.0
    
    YJoint.info()
    
    return YJoint

def K_example():
    
    chord = RHS(80,140,6,fy=700)
    b1 = RHS(70,70,3,fy=700)
    b2 = RHS(40,60,4,fy=700)
    
    b1.Ned = 123.743e3
    b2.Ned = -123.743e3
    
    Kjoint = RHSKGapJoint(chord, [b1,b2], [45,45],40,N0=-362500.0)
    
    print(Kjoint.eval_KN())
    
    return Kjoint
    
if __name__ == '__main__':
    
    from sections.steel.RHS import RHS, SHS
    #YJ = Y_example()
    
    #KJ = K_example()
    
    """
    C = RHS(60,40,3,fy=700)
    C.material.fy = 770
    B1 = RHS(40,40,4,fy=700)
    B1.material.fy = 770
    B2 = RHS(40,40,4,fy=700)
    B2.material.fy = 770
    """
    
    C = RHS(60,40,3,fy=700)
    C.material.fy = 770
    B1 = RHS(40,40,4,fy=700)
    B1.material.fy = 770
    B2 = RHS(40,40,4,fy=700)
    B2.material.fy = 770
    
    C.Ned = -500e3
    B1.Ned = 200e3
    B2.Ned = -200e3
    
    
    
    K = RHSKGapJoint(C, [B1,B2], [45,45],15,N0=C.Ned)
    
    K.info()
    
    """
    
    
    
    
    
    chord = SHS(200,8,fy=420)
    
    chord.Ned = -1364e3
    b1 = SHS(150,6,fy=420)
    b2 = SHS(150,6,fy=420)
        
    
    J = RHSKGapJoint(chord,[b1,b2], [45,45], 28)
    
    J.N0 = -1364e3
    J.braces[0].Ned = 600e3
    J.braces[1].Ned = -600e3
    
    #print(J.N0gap*1e-3)
    #print(J.V0gap*1e-3)
    
    #J.info()
    #J.draw()
    
    b1.Ned = -600e3
    
    #X = RHSXJoint(chord,b1,45)
    
    #X.chord_web_buckling()
    
    #X.info()
    
    #l1 = Line(p1=[0,0],v=[1,0])
    #l2 = Line(p1=[1,1],v=[0,1])
    
    #t, x0 = l1.intersect(l2)
    
    # Y Joint example from SSAB's book  (Ex. 3.1) 
    chord = SHS(200,8,fy=420)
    brace = SHS(100,5,fy=420)
    
    brace.Ned = 590e3
    
    Y = RHSYJoint(chord,brace,45,N0=-936e3)
    
    # Y Joint example from SSAB's book  (Chapter 7) 
    chord = SHS(200,8,fy=420)
    brace = SHS(120,5,fy=420)
    
    brace.Ned = 750e3
    
    Y = RHSYJoint(chord,brace,54,N0=-527e3,M0=-653e3*200)
    
    Y.info()
    """
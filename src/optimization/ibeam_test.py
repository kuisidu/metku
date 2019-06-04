# -*- coding: utf-8 -*-
"""
Created on Tue May 21 19:24:45 2019

@author: kmela
"""

import math
from functools import partial
from sections.steel.wi import WISection
#from structopt import Variable, LinearConstraint, Constraint, OptimizationProblem
import structopt as sop
import numpy as np



def IBeamWeight(p,x):
    """ Weight of I-beam
        
        Parameters:
        -----------
            param:x: [h,tw,bf,tf]
            param:beam: SteelMember object
    
        Output:
        ------
            f: weight of the beam
            
    """
    L = 6000
    p.h = x[0]
    p.tw = x[1]
    
    p.b[0] = x[2]
    p.b[1] = x[2]
    
    p.tf[0] = x[3]
    p.tf[1] = x[3]
    
    return p.weight()*L

def IBeamWeightGrad(p,x):
    """ Gradient of the weight of I-beam
        
        Parameters:
        -----------
            param:x: [h,tw,bf,tf]
            param:beam: SteelMember object
    
        Output:
        ------
            f: weight of the beam
            
    """
    L = 6000
    """
    p.h = x[0]
    p.tw = x[1]
    
    p.b[0] = x[2]
    p.b[1] = x[2]
    
    p.tf[0] = x[3]
    p.tf[1] = x[3]
    """
    
    df = p.density*L*np.array([x[1],x[0]-2*x[3],2*x[3],2*x[2]-2*x[1]])
    
    return df

def IBeamWeightHessian(p,x):
    """ Hessian of the weight of I-beam
        
        Parameters:
        -----------
            param:x: [h,tw,bf,tf]
            param:beam: SteelMember object
    
        Output:
        ------
            f: weight of the beam
            
    """
    L = 6000
    """
    p.h = x[0]
    p.tw = x[1]
    
    p.b[0] = x[2]
    p.b[1] = x[2]
    
    p.tf[0] = x[3]
    p.tf[1] = x[3]
    """
    
    H = np.zeros([4,4])
    H[0,1] = 1
    H[1,0] = 1
    H[1,3] = -2
    H[2,3] = 2
    H[3,1] = -2
    H[3,2] = 2

    return p.density*L*H


def IBeamFlangeClassCon(p,flange_class=2):
    """ Builds a constraint for cross-section class of an I beam """
    
    a = [0,0,0,0]
    b = 0
    con_type = '<'
    
    e = p.eps
    
    if flange_class == 1:
        C0 = 9*e
    elif flange_class == 2:
        C0 = 10*e
    elif flange_class > 2:
        C0 = 14*e
    
    b = math.sqrt(2)*p.weld_throat
    a[1] = -C0*e
    a[2] = 0.5
    a[3] = -0.5
    
    if flange_class > 3:
        """ if class 4 is required, the cf/tf ratio needs to
            be greater than the class 3 limit. This  changes the
            direction of the constraint from < to >.
        """
        con_type = '>'
    
    con_name = "Flange in class " + str(flange_class)
    con = sop.LinearConstraint(a,b,con_type,name=con_name)
    
    return con

def IBeamWebClassCon(p,web_class=2):
    """ Builds a constraint for cross-section class of an I beam """
    
    a = [0,0,0,0]
    b = 0
    con_type = '<'
    
    e = p.eps
     
    """ web is assumed to be in bending """
    if web_class == 1:
        C1 = 72*e
    elif web_class == 2:
        C1 = 83*e
    elif web_class > 2:
        C1 = 124*e
    
    b = 2*math.sqrt(2)*p.weld_throat
    a[0] = 1
    a[1] = -C1*e
    a[3] = -2
    
    if web_class > 3:
        """ if class 4 is required, the cw/tw ratio needs to
            be greater than the class 3 limit. This  changes the
            direction of the constraint from < to >.
        """
        con_type[1] = '>'
    
    con_name = "Web in class " + str(web_class)
    con = sop.LinearConstraint(a,b,con_type,name=con_name)
    
    return con

def IBeamBendingCon(p,section_class,x):
    """ Constraint for the bending moment resistance """
    p.h = x[0]
    p.tw = x[1]
    
    p.b[0] = x[2]
    p.b[1] = x[2]
    
    p.tf[0] = x[3]
    p.tf[1] = x[3]

    if section_class < 3:
        MRd = p.plastic_bending_resistance()
    elif section_class == 3:
        MRd = p.elastic_bending_resistance()
    
    return p.Med - MRd
        

if __name__ == '__main__':
    
    p = WISection(h=300,tw=10,b=[100,100],tf=[15,15])
    
    """ Design variables:
        x[0] .. height of the beam
        x[1] .. thickness of web
        x[2] .. width of the flanges
        x[3] .. thickness of flanges
    """
    dvars = []
    dvars.append(sop.Variable("Height",150,1000))
    dvars.append(sop.Variable("Web thickness",5,40))
    dvars.append(sop.Variable("Flange width",150,500))
    dvars.append(sop.Variable("Flange thickness",5,40))
    #dvars.append(sop.Variable("Web thickness",5,40))
    
    obj = partial(IBeamWeight,p)
    obj_grad = partial(IBeamWeightGrad,p)
    obj_hess = partial(IBeamWeightHessian,p)
    
    prob = sop.OptimizationProblem(name="I-Beam Weight Minimization",\
                            variables=dvars,objective=obj,gradient=obj_grad)
    
    flange_class = 2
    web_class = 2
    sec_class = max(flange_class,web_class)
    
    flange_con = IBeamFlangeClassCon(p,flange_class)
    web_con = IBeamWebClassCon(p,web_class)
    
    prob.add_constraints(flange_con)
    prob.add_constraints(web_con)
    
    bending_con_fun = partial(IBeamBendingCon,p,sec_class)
    bending_con = sop.NonLinearConstraint(bending_con_fun,'<',\
                                      "Bending resistance (mid-span)")
        

    prob.add_constraints(bending_con)
    
    x = [280,8,110,12]
    #prob(x)
    #print(prob.grad(x))
    #x0 = np.array(x)
    #print(sop.NumGrad(prob.obj,0.01*x0,x0))
    
    r0 = prob.solve("slsqp",x0=x)
    r1 = prob.solve("trust-constr",x0=x)
    
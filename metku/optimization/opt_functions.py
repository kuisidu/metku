# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:33:14 2022

Various functions to be used in optimization

@author: kmela
"""

import numpy as np
from copy import copy

from metku.optimization.variables import Variable, DiscreteVariable, IndexVariable, BinaryVariable
from metku.optimization.variables import CrossSectionVariable
#from metku.optimization.constraints import LinearConstraint

import metku.optimization.constraints as constr

def NumGrad(fun, x):
    """ Gradient by finite difference """
    fx = fun(x)
    n = len(x)

    """ if only a single value for step length 'h' is given,
        transform it to a list
    """

    # if isinstance(h, float):
    #    h = h * np.ones(n)

    df = np.zeros(n, dtype=np.float)
    #    
    for i in range(len(x)):
        # print("Täällä")
        # xh = deepcopy(x)
        xh = np.asarray([val for val in x], dtype=np.float)
        h = max(1e-6 * abs(xh[i]), 1e-8)
        xh[i] += h

        fh = fun(xh)

        df[i] = (fh - fx) / h

    return df


def conlin_fun(xk, fk, dfk, pos_vars=True):
    """ Generates convex-linear approximation for a function
        at point 'xk'
        
        :param xk: approximation point
        :param fk: value of the function at 'xk'
        :parm dfk: value of the gradient at 'xk'
    """
    if isinstance(xk, float) or isinstance(xk, int):
        xk = np.array([xk])
        dfk = np.array([dfk])

    n = len(xk)

    """ Initialize 'alpha' as zeros. """
    a = np.zeros(n)
    for i, Xk in enumerate(xk):
        """ If the derivative of f with respect to xk is non-negative,
            set a[i] = 1
        """
        if dfk[i] >= 0:
            a[i] = 1

    def generate_a(x):

        a = np.ones(n)
        for i, Xk in enumerate(xk):
            if x[i] * dfk[i] < 0:
                a[i] = Xk / x[i]

    def eval_conlin(x):
        fC = fk
        if isinstance(x, float) or isinstance(x, int):
            x = np.array([x])

        for i in range(n):
            if dfk[i] >= 0:
                a[i] = 1.0
            else:
                a[i] = xk[i] / x[i]

        fC += (a * dfk).dot(x - xk)

        return fC

    return eval_conlin


def mma_fun(xk, fk, dfk, L, U):
    """ Generates the MMA approximation of a function """

    if isinstance(xk, float) or isinstance(xk, int):
        xk = np.array([xk])
        dfk = np.array([dfk])

        L = np.array([L])
        U = np.array([U])

    n = len(xk)

    """ Initialize 'p' and 'q' as zeros. """
    p = np.zeros(n, dtype=np.float)
    q = np.zeros(n, dtype=np.float)
    for i, (Xk, Lk, Uk) in enumerate(zip(xk, L, U)):
        # Determine the constants p and q
        if dfk[i] >= 0:
            p[i] = (Uk - Xk) ** 2 * dfk[i]
        else:
            q[i] = -(Xk - Lk) ** 2 * dfk[i]

    # Evaluate r(xk)
    rk = fk - sum(p / (U - xk) + q / (xk - L))

    def eval_mma(x):
        # Evaluate MMA approximation
        return rk + sum(p / (U - x) + q / (x - L))

    return eval_mma


class MMAFunction:
    """ Class for creating and evaluating MMA approximations """

    def __init__(self, xk, fk, dfk, L, U):
        """

        Parameters
        ----------
        xk : numpy array
            Point, where the MMA approximation is made.
        fk : float
            Value of the approximated function at xk, i.e. fk = f(xk).
        dfk : numpy array
            Gradient of the approximated function at xk.
        L : numpy array
            lower asymptotes.
        U : numpy array
            upper asymptotes.

        Returns
        -------
        None.

        """

        self.n = len(xk)
        self.xk = np.array(xk)
        self.fk = np.array(fk)
        self.dfk = np.array(dfk)
        self.L = np.array(L)
        self.U = np.array(U)

        self.p = np.zeros(self.n)
        self.q = np.zeros(self.n)
        for i, (Xk, Lk, Uk) in enumerate(zip(self.xk, self.L, self.U)):
            # Determine the constants p and q
            if self.dfk[i] >= 0:
                self.p[i] = (Uk - Xk) ** 2 * self.dfk[i]
            else:
                self.q[i] = -(Xk - Lk) ** 2 * self.dfk[i]

        # Evaluate r(xk)
        self.rk = self.fk - sum(self.p / (self.U - self.xk) + self.q / (self.xk - self.L))

    def __call__(self, x):
        # Evaluate MMA approximation
        return self.rk + sum(self.p / (self.U - x) + self.q / (x - self.L))

    def grad(self, x):
        return self.p / (self.U - x) ** 2 - self.q / (x - self.L) ** 2


def Linearize(fun, x, grad=None):
    """ Linearization of a function

        input:
            fun .. function that takes x as input
            x .. point of linearization (array)
            grad .. returns the gradient of fun at x

        output:
            a .. grad(x)
            b .. fun(x)
            c .. grad(x)*x

    """

    b = fun(x)

    if grad == None:
        h = 0.01 * x  # check if this is OK
        # h = np.ones_like(x) * 1e-2
        a = NumGrad(fun, x, h)
    else:
        a = grad(x)

    c = a.dot(x)

    return a, b, c

def index2binary(P,disc_vars,var_type='continuous'):
    """ Changes problem formulation from index variables to binary variables
        to be used in selecting the cross-section.
        
        disc_vars .. list with names of attributes to be added as discrete 
                     variables. For example, disc_vars = ['A', 'Iy', 'H', 'Wply']
    """
    
    Pnew = P
    
    var_pack = []
    
    
    for i, var in enumerate(P.vars):
        if isinstance(var,IndexVariable):
            # Operaatiot
            # 1. Luo muuttujat disc_varista. Niillä on sama target[object] kuin
            # 'var'-muuttujalla, ja property otetaan disc_varista
            # 2. Luo binäärimuuttujat. Binäärimuuttujia luodaan yhtä monta kuin on
            # profiileja, eli 'values'-attribuutissa alkioita.
            # 3. Luo lineaariset rajoitusehdot poikkileikkaus- ja binäärimuuttujien välille.
            # Tässä haasteeksi tulee muuttujien indeksin kanssa pelaaminen.
            ndx_set = {'props':[],'binary':[]}
            
            
            for sec_prop in disc_vars:
                # props sisältää poikkileikkaussuureen 'sec_prop' arvot
                # jotka vastaavat indeksimuuttujan katalogia
                # var.values = indeksimuuttujan profiililista
                props = [round(s.__getattribute__(sec_prop),3) for s in var.values]
                var_name = sec_prop + str(i)
                # ONGELMA:
                # Indeksimuuttujalla on target['objects'] lista sauvoja, joiden poikkileikkausta muutetaan
                # 'property' on 'cross_section'. Nyt pitäisi saada target['objects'] edelleen lista sauvoja, mutta
                # niiden tiettyä poikkileikkaussuuretta muutetaan.
                for mem in var.target['objects']:                    
                    mem.cross_section = copy(mem.cross_section)
                    
                #target = {'objects': [mem.cross_section for mem in var.target['objects']],
                #          'property': sec_prop}
                target = {'objects': [mem for mem in var.target['objects']],
                          'property': sec_prop}
                
                #for mem in var.target['objects']:
                #    mem.cross_section = target['objects'][0]
                
                if var_type == 'continuous':
                    lb = min(props)
                    ub = max(props)
                
                    new_var = CrossSectionVariable(var_name,lb,ub,target=target,value=props[var.value])                    
                    
                else:
                    new_var = DiscreteVariable(var_name,props,target,value=props[var.value])
                
                Pnew.add(new_var)
                ndx_set['props'].append(new_var)
            
            # Tehdään binäärimuuttujat:
            # Jokaiselle valittavissa olevalle profiilille luodaan omansa
            for j, sec in enumerate(var.values):
                var_name = var.name + ' profile ' + str(j)
                if var.value == j:
                    val = 1
                else:
                    val = 0
                #new_var = BinaryVariable(var_name,section=sec,target=var.target,value=val)
                new_var = BinaryVariable(var_name,section=sec,objects=var.target['objects'],target=None,value=val)
                Pnew.add(new_var)
                ndx_set['binary'].append(new_var)
            
            var_pack.append(ndx_set)
            
            Pnew.del_var(var)
    
    
    # Luodaan rajoitusehdot
    nx = Pnew.nvars()
    for pack in var_pack:
        # Jokainen poikkileikkausmuuttuja sidotaan ko. sauvan binäärimuuttujiin
        for prop_var in pack['props']:               
            a = np.zeros(nx)
            a[Pnew.vars.index(prop_var)] = 1
            for bin_var in pack['binary']:
                a[Pnew.vars.index(bin_var)] = -round(bin_var.section.__getattribute__(prop_var.target['property']),3)
            
            con_name = prop_var.name + " unique profile"
            Pnew.add(constr.LinearConstraint(a,0,'=',con_name,Pnew))
        
        # Kullekin sauvalle/ryhmälle binäärimuuttujien summa on tasan 1
        a = np.zeros(nx)
        for bin_var in pack['binary']:
            a[Pnew.vars.index(bin_var)] = 1
        
        con_name = "Unique profile"
        Pnew.add(constr.LinearConstraint(a,1,'=',con_name,Pnew))
           
    return Pnew, var_pack

def potentially_active_constraints(prob,x):
    
    ACT_EPS = 0.25
    
    act_cons = []
    
    for con in prob.cons:
        if con.type == '=':
            act_cons.append(con)
        elif con.type == '<':
            if con(x) > -ACT_EPS:
                act_cons.append(con)
        else:
            if con(x) < ACT_EPS:
                act_cons.append(con)
    
    return act_cons

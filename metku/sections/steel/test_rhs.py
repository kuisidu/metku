# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:43:48 2018

@author: huuskoj
"""

import unittest
from RHS import RHS
# 0  1   2   3    4     5    6   7    8     9    10
# A, It, Wt, Iy, Wely, Wply, iy, Iz, Welz, Wplz, iz
ARGS1 = [9.59*1e2, 215.8*1e4, 43.23*1e3, 195.8*1e4, 32.63*1e3,
         39.07*1e3, 45.2, 105.2*1e4, 26.30*1e3, 29.65*1e3, 33.1]
ARGS2 = [2.14*1e2,3.45*1e4, 2.36*1e3, 4.05*1e4, 2.02*1e3,
         2.61*1e3, 13.8, 1.34*1e4, 1.34*1e3, 1.60*1e3, 7.9]
# NcRd, McyRd, MczRd, VplyRd, VplzRd [kN /kNm]
RD1 = [299.5,13.87,8.02,78.61,117.9]
RD2 = [75.86,0.93,0.57,14.60,29.20]
"""
!!! NOTE !!!
rhs1 is cross-section class 4 and class 4 
calculations haven't been implemented yet
!!! NOTE !!!
"""
class TestRHS(unittest.TestCase):
    
    def test_area(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.A, ARGS1[0], delta=0.01*ARGS1[0])
        self.assertAlmostEqual(rhs2.A, ARGS2[0], delta=0.01*ARGS2[0])
    """
    def test_Wt(self):
        rhs1 = RHS(120,80,3)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.It, ARGS1[1], delta=1)
        self.assertAlmostEqual(rhs2.It, ARGS2[1], delta=1)
    """
    def test_It(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.It, ARGS1[2], delta=0.01*ARGS1[2])
        self.assertAlmostEqual(rhs2.It, ARGS2[2], delta=0.01*ARGS2[2])
    
    def test_Iy(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.I[0], ARGS1[3], delta=0.01*ARGS1[3])
        self.assertAlmostEqual(rhs2.I[0], ARGS2[3], delta=0.01*ARGS2[3])
    
    def test_Wely(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.Wel[0], ARGS1[4], delta=0.01*ARGS1[4])
        self.assertAlmostEqual(rhs2.Wel[0], ARGS2[4], delta=0.01*ARGS2[4])
        
    def test_Wply(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.Wpl[0], ARGS1[5], delta=0.01*ARGS1[5])
        self.assertAlmostEqual(rhs2.Wpl[0], ARGS2[5], delta=0.01*ARGS2[5])
    """   
    def test_iy(self):
        rhs1 = RHS(120,80,3)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.Wpl[1], ARGS1[6], delta=1)
        self.assertAlmostEqual(rhs2.Wpl[1], ARGS2[6], delta=1)
    """    
    def test_Iz(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.I[1], ARGS1[7], delta=0.01*ARGS1[7])
        self.assertAlmostEqual(rhs2.I[1], ARGS2[7], delta=0.01*ARGS2[7])
        
    def test_Welz(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.Wel[1], ARGS1[8], delta=0.01*ARGS1[8])
        self.assertAlmostEqual(rhs2.Wel[1], ARGS2[8], delta=0.01*ARGS2[8])
        
    def test_Wplz(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.Wpl[1], ARGS1[9], delta=0.01*ARGS1[9])
        self.assertAlmostEqual(rhs2.Wpl[1], ARGS2[9], delta=0.01*ARGS2[9])
    
    """
    def test_iz(self):
        rhs1 = RHS(120,80,3)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.Wpl[1], ARGS1[10], delta=1)
        self.assertAlmostEqual(rhs2.Wpl[1], ARGS2[10], delta=1)
    """   
       
    def test_NcRd(self):
        rhs1 = RHS(120,80,2.5)
        rhs1.Ned = -1
        rhs2 = RHS(40,20,2)
        rhs2.Ned = -1
        #self.assertAlmostEqual(rhs1.axial_force_resistance()/1e3, RD1[0], delta=1)
        self.assertAlmostEqual(rhs2.axial_force_resistance()/1e3, RD2[0], delta=1)
        
    def test_McyRd(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        #self.assertAlmostEqual(rhs1.bending_resistance()[0]/1e6, RD1[1], delta=1)
        self.assertAlmostEqual(rhs2.bending_resistance()[0]/1e6, RD2[1], delta=1)
        
    def test_MczRd(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        #self.assertAlmostEqual(rhs1.bending_resistance(axis='z')[0]/1e6, RD1[2], delta=1)
        self.assertAlmostEqual(rhs2.bending_resistance(axis='z')[0]/1e6, RD2[2], delta=1)
        
    def test_VcplyRd(self):
        rhs1 = RHS(120,80,2.5)
        rhs2 = RHS(40,20,2)
        self.assertAlmostEqual(rhs1.shear_force_resistance()/1e3, RD1[3], delta=1)
        self.assertAlmostEqual(rhs2.shear_force_resistance()/1e3, RD2[3], delta=1)       
        
        
if __name__ == '__main__':
    unittest.main()
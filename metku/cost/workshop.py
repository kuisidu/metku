# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 13:12:18 2019

@author: kmela

Contents

Class Workshop

Class CostCentre
Class BlastingCost(CostCentre)
Class CuttingCost(CostCentre)
Class PaintingCost(CostCentre)
Class AssemblingCostWeld(CostCentre)

TODO:
    Beam Welding .. Haapio 4.5.4, for welded beams
    Sawing .. Haapio 4.5.5
    Drilling .. Haapio 4.5.6
    Coping .. Haapio 4.5.7
    Fabrication of parts .. Haapio 4.5.8 (e.g. angle profiles, plates)
    Assembling (Bolts) .. Haapio 4.5.9
    Post-treatment and inspection .. Haapio 4.5.10
    Drying
    Erecting

"""

import copy
from colorama import Fore
from cost.cost_data import BASIC_STEEL_PRICE

NOT_IMPLEMENTED = "NOT IMPLEMENTED YET! Add this feature: "

work_days_per_year = 252
work_hours_per_day = 8
work_minutes_per_year = work_days_per_year*work_hours_per_day*60

wages = {"plate_and_construction": 18.37,
         "welding": 18.13,
         "machine": 19.76,
         "quality_control":19.16,
         "erecting_work":18.76
         }

overhead = 0.7

wage_per_minute = copy.deepcopy(wages)

# change hourly rates to euros per minute, including overheads
for key in wage_per_minute:
    wage_per_minute[key] *= (1+overhead)/60

#total_wages_per_min = 

""" Construction of the workshop, based on
    Haahtela & Kiiras, Talonrakennuksen kustannustieto 2015
    Alue 4, Kasvukeskukset, hintataso 78
""" 
construction_cost = 1000 # €/m2

""" Cost of plot
    It is assumed that the workshop is located in Seinäjoki, Finland
"""
plot_efficiency = 0.5
raw_land_cost = 6 # €/m2
land_cost = raw_land_cost/plot_efficiency

""" Maintenance cost of real estate
    KH-cards KH X1-00244, KH X1-00291 and KH X1-00379
    Teollisuushallin ylläpito 2000/4, KH X1-00291 taulukko 6
"""
industrial_hall_maintenance = 6.385 # €/(m3*a)

""" Kiinteistön ylläpidon kustannusindeksi 
    rakennetyypeittäin (Tilastokeskus), 2000=100, Teollisuushalli
"""
maintenance_index = 181 # Q4/2014

building_height = 8.0 # m
real_estate_maintenance_cost = building_height*industrial_hall_maintenance*maintenance_index/100

#print(f'{real_estate_maintenance_cost:.2f} €')

""" Cost of electricity """
cost_energy = 0.1 # €/kWh
cost_energy_min = cost_energy/60

""" Investements 
    IEq .. interest rate of equipment
    nEq .. life of equipment (years)
    IRealS .. interest reate of real estate
    nRealS .. life of real estate (years)
"""

IEq = 0.05
nEq = 20
IRealS = 0.05
nRealS = 50

""" Paint data
    vol_solid .. volume of solids
    DFT .. dry film thickness
    price .. €/l
"""
tip_flow = 900000 # [mm3/min]

paints = {"alkyd": {"vol_solid":[0.48,0.45,0.45], "DFT":[80,40,40],"price":[2.8,3.6,3.6],"time":3,"TPP":0.0,"CP":0.0},
          "epoxy": {"vol_solid":[0.50,0.50,0.48], "DFT":[60,60,40],"price":[4.2,4.2,4.9],"time":14,"TPP":0.0,"CP":0.0},
          "polyurethane": {"vol_solid":[0.50,0.50], "DFT":[60,60],"price":[4.2,5.9],"time":27,"TPP":0.0,"CP":0.0},
          "acryl": {"vol_solid":[0.48,0.48,0.40], "DFT":[60,60,40],"price":[6.4,6.4,6.4],"time":7,"TPP":0.0,"CP":0.0},
          }

# Calculate basic values for productive time and cost of consumables
for paint, films in paints.items():    
    nfilms = len(films["vol_solid"])
    TPP = 0.0
    CP = 0.0
    for i in range(nfilms):
        tloss = 0.01+0.3*films["DFT"][i]*1e-3
        P0 = 1/films["vol_solid"][i]*(films["DFT"][i]*1e-3+tloss)
        TPP += P0
        CP += films["price"][i]*1e-6*P0

    TPP = TPP/tip_flow
    paints[paint]["TPP"] = TPP
    paints[paint]["CP"] = CP


""" Welding consumables, prices """
""" Welding process 135 (MAG with solid wire) """
wire_prices = {"S355": 1.66,
               "S420": 1.66,
               "S550": 3.5,
               "S690": 2.34}

sawing_times = {"setup":2,
                "blade_descent":0.25,
                "blade_raising":0.5,
                "profile_fixing":0.5,
                }

sawing_base_time = sawing_times["setup"]+\
            2*(sawing_times["blade_descent"]+sawing_times["blade_raising"]+\
               sawing_times["profile_fixing"])

saw_bevel_time = 1.0 # min

sawing_parameters = {"S235":{"Sm":2.0, "Q":1800},
                     "S355":{"Sm":2.0, "Q":1800},
                     "S420":{"Sm":2.0, "Q":1800},
                     "S500":{"Sm":2.0, "Q":1800},
                     "S550":{"Sm":2.0, "Q":1800},
                     "S690":{"Sm":2.0, "Q":1800},
                     "S700":{"Sm":2.0, "Q":1800}}

""" Cost of one saw blade [€] """
#saw_blade_cost = 100.0
saw_blade_cost = 170.0

def print_initial_cost_data():
    """ Prints a set of initial data used in cost calculations """
    
    header = "-----"
    
    print("Cost analysis following Haapio (2012)")
    print(header + " Working time " + header)
    print(f"Working days per year: {work_days_per_year} d")
    print(f"Working hours per shift: {work_hours_per_day} h")
    print(f"Working minutes per year: {work_minutes_per_year}")
    
    print('\n' + header + " Wages " + header)
    
    for key, value in wages.items():
        print(f'{key}: {value:.2f} €/h')
    
    print(f'Overhead: {overhead}')

    for key, value in wage_per_minute.items():
        print(f'{key}: {value:.3f} €/min')
    
    print('\n' + header + " Real estate " + header)
    
    print(f"Construction cost: {construction_cost} €/m2") 
    print(f"Plot efficiency: {plot_efficiency} ") 
    print(f"Raw land cost: {raw_land_cost} €/m2") 
    print(f"Land cost (Seinäjoki, Finland): {land_cost} €/m2") 
    print(f"Interest rate of real estate: {IRealS*100:.1f}%")
    print(f"Life time of real estate: {nRealS} a")
   
    print(f"Maintenance cost, industrial hall: {industrial_hall_maintenance} €/(m3*a)") 
    print(f"Maintenance index (Q4/2014): {maintenance_index} ") 
    print(f"Building height: {building_height} m")
    print(f"Real estate maintenance cost: {real_estate_maintenance_cost:.2f} €/(m2*a)")
    
    print('\n' + header + " Cost of electricity " + header)
    print(f'Cost of energy: {cost_energy} €/kWh')
    print(f'Cost of energy per minute: {cost_energy_min:.2f} €/kWh/min')

    print('\n' + header + " Investment of equipment " + header)
    print(f'Interest rate of equipment: {IEq*100} %')
    print(f'Life time of equipment: {nEq} a')
    
    print('\n' + header + " Painting " + header)
    for paint, values in paints.items():
        print(f"{paint}: TPP = {values['TPP']:.2e}, CP = {values['CP']:.2e} ")
        
    
    
    """
    paints = {"alkyd": {"vol_solid":[0.48,0.45,0.45], "DFT":[80,40,40],"price":[2.8,3.6,3.6],"time":3,"TPP":0.0,"CP":0.0},
              "epoxy": {"vol_solid":[0.50,0.50,0.48], "DFT":[60,60,40],"price":[4.2,4.2,4.9],"time":14,"TPP":0.0,"CP":0.0},
              "polyurethane": {"vol_solid":[0.50,0.50], "DFT":[60,60],"price":[4.2,5.9],"time":27,"TPP":0.0,"CP":0.0},
              "acryl": {"vol_solid":[0.48,0.48,0.40], "DFT":[60,60,40],"price":[6.4,6.4,6.4],"time":7,"TPP":0.0,"CP":0.0},
              }
    """
# Annual cost, Haapio, eq (2)
# Uniform series end-of-period installment [e/a]
def annual_cost(P, I, n):
    # P = a sum of money invested in the initial year less resale value [e]
    # I = interest rate (cost of capital) [%]
    # n = time, the number of units (years) over which interest accumulates

    return P*((I*(1.0 + I)**n)/((1.0 + I)**n - 1.0))       # [e/a]

class CostCentre:
    """
        Class for cost centres
    """
    def __init__(self,length=20,width=10,workers=1,
                 eq_invest=1000,eq_resale=0,eq_maint=1000,eq_year=nEq):
        
        self.name = "Cost Centre"
        self.length = length
        self.width = width
        self.nw = workers        
        
        self.times = {"productive":0.0,"non-productive":0.0}
        #self.productive_time = 0.0
        #self.non_productive_time = 0.0
        
        self.costs = {"labour":0.0,
                      "equipment":0.0,
                      "equip_maintenance":0.0,
                      "real_estate":0.0,
                      "real_maintenance":0.0,                     
                      }
        
        self.productive_costs = {"consumables":0.0,
                                 "energy":0.0
                                 }
        
        self.consumables_cost = 0.0
        self.utility = 1.0             
        
        """ Unit costs of real estate can be calculated similarly to 
            all cost centres.
            
            Other unit costs depend on the specifics.
        """
        # Investment cost of workspace area [e]
        P_invest_RealS = self.Area*construction_cost
        # Resale value of workspace [e], equals value of land
        P_resale_RealS = self.Area*land_cost          
        self.costs["real_estate"] = annual_cost(P_invest_RealS - P_resale_RealS, IRealS, nRealS)/work_minutes_per_year
        
        # Real estate maintenance cost
        c_maintenance_RealS = self.Area*real_estate_maintenance_cost
        # Real estate maintenance cost[e/min]
        self.costs["real_maintenance"] = c_maintenance_RealS/work_minutes_per_year
        
        # Equipment installment cost [e/min]
        # Investment of equipment [e]
        
        self.costs["equipment"] = annual_cost(eq_invest - eq_resale, IEq, eq_year)/work_minutes_per_year

        # Maintenance cost of equipment        
        self.costs["equip_maintenance"] = eq_maint/work_minutes_per_year
        
    @property
    def Area(self):
        return self.length*self.width
        
    # Total cost of cost center [e]
    def cost(self):
        return sum(self.times.values())*sum(self.costs.values())/self.utility + \
                self.times["productive"]*sum(self.productive_costs.values()) + \
                self.consumables_cost
    
    def info(self):
        
        print(Fore.GREEN + f"Cost centre: {self.name}")
        print(Fore.RESET + f"Area: {self.Area:.2f} m2")
        
        print(f"Labour costs: {self.costs['labour']:.2f} €/min")
        print(f"Equipment costs: {self.costs['equipment']:.2f} €/min")
        print(f"Equipment maintenance costs: {self.costs['equip_maintenance']:.2f} €/min")
        print(f"Real estate costs: {self.costs['real_estate']:.2f} €/min")
        print(f"Real estate maintenance costs: {self.costs['real_maintenance']:.2f} €/min")
# Cost center for beam welding
class BeamWeldingCost(CostCentre):

    def __init__(self,length=20,width=15,workers=2,
                 eq_invest=200000.0,eq_resale=0.0,eq_maint=1000.0):
    
        CostCentre.__init__(self, length,width,workers,
                            eq_invest,eq_resale,eq_maint)
        
        self.name = "Beam Welding"
        
        self.times["non-productive"] = 6.25 # minutes
        
        self.nwelds = 2 # number of simultaneous welds
        self.deposit_rate = 0.14 # kg/min
        self.wire_price = 1.84 # €/kg
        self.flux_price = 1.46 # €/kg
        self.flux_wire_ratio = 1.0
        

        # Labour cost
        self.costs["labour"] = self.nw*wage_per_minute["machine"] # Hourly cost [e/h]
         
                     
        # Cost of consumables [e/min] TODO!
        self.productive_costs["consumables"] = 1.0
        # Energy
        # Total power consumption of one weld head [kW]
        self.power = 24                       
        # Cost of energy [e/min]
        self.productive_costs["energy"] = self.power*cost_energy_min                

    def ProductiveTime(self,weld_length,weld_size,weld_type="fillet"):
        """ Calculate the productive time for welding
            input: Lb .. length of steel piece (mm)
            output: productive time (min)
        """
        
        
        TP = Lb/self.vc
        self.times["productive"] = TP
        return TP  


    def welding_time(self):
        """ Time required for welding [min]
            input: weld_size .. throat thickness of fillet weld or
                                depth of single-bevel butt weld
                   weld_length [mm]
                   weld_type .. "fillet" or "butt"
        """
        if self.weld_type == "fillet":
            a = self.weld_size            
            tweld = self.weld_length*1e-3*(0.4988*a**2-0.0005*a+0.0021)            
        elif self.weld_type == "butt":
            b = self.weld_size
            tweld = self.weld_length*1e-3*(0.249*b**2+0.0096*b-0.0506)
        
        return tweld

    def productive_time(self,weld_length,weld_size,weld_type="fillet"):
        """ Productive time of welding
            
        """
        # Productive time [min], painting        
        if self.weld_type == "fillet":
            a = self.weld_size
            weld_weight = self.weld_length*a**2*7850e-9
        elif self.weld_type == "butt":
            b = self.weld_size
            weld_weight = self.weld_length*b**2/2*7850e-9
        
        self.times["productive"] = weld_weight/self.deposit_rate
        

    def consumables(self):
        """ Cost of welding consumables            
        """
        
        if self.weld_type == "fillet":
            a = self.weld_size
            weld_weight = self.weld_length*a**2*7850e-9
        elif self.weld_type == "butt":
            b = self.weld_size
            weld_weight = self.weld_length*b**2/2*7850e-9
                
        self.consumables_cost = weld_weight*(self.wire_price + self.gas_consumption*self.gas_price)
                    
        return self.consumables_cost      


    
# Cost center for blasting
class BlastingCost(CostCentre):

    def __init__(self,length=40,width=10,workers=1,speed=3000,
                 eq_invest=200000.0,eq_resale=0.0,eq_maint=1000.0):
    
        CostCentre.__init__(self, length,width,workers,
                            eq_invest,eq_resale,eq_maint)
        
        self.name = "Blasting"
        # Productive time
        self.vc = speed                         # Conveyor speer [mm/m]        

        # Labour cost
        self.costs["labour"] = self.nw*wage_per_minute["machine"] # Hourly cost [e/h]
         
        # Consumables: grains
        # Grain consumption on blasting [kg/min], equipment utilization ratio 1
        # For utilization ratio of 0.5, grain consumption is 10kg/8h
        grain_consump = 0.042                   
        # Grain weight price [e/kg]
        #grain_unit_cost = 0.50                  
        grain_unit_cost = 0.75 

        # Cost of consumables [e/min]
        self.productive_costs["consumables"] = grain_consump*grain_unit_cost
        # Energy
        # Total power consumption using four blasting turbines [kW]
        P_total = 40                                    
        # Cost of energy [e/min]
        self.productive_costs["energy"] = P_total*cost_energy_min                

    def ProductiveTime(self,Lb):
        """ Calculate the productive time for steel piece of length Lb
            input: Lb .. length of steel piece (mm)
            output: productive time (min)
        """
        TP = Lb/self.vc
        self.times["productive"] = TP
        return TP

    def cost(self,Lb):
        
        self.ProductiveTime(Lb)
        
        return super().cost()

# Cost center for cutting
class CuttingCost(CostCentre):

    def __init__(self, length=20.0,width=10.0, workers=2.0,
                 eq_invest=280000.0,eq_resale=0.0,eq_maint=1000.0):

        CostCentre.__init__(self, length,width,workers,
                            eq_invest,eq_resale,eq_maint)
        
        self.name = "Cutting"
        
        self.times["non-productive"] = 3.0          
        
        # Labour cost
        self.costs["labour"] = self.nw*wage_per_minute["machine"] # Hourly cost [e/h]
            
        

    def ProductiveTime(self,t_p,cut_length):
        """ Calculate the productive time
            input: Lb .. length of steel piece (mm)
            output: productive time (min)
        """
        # Cutting speed according to thickness [mm] of the plate
        if t_p <= 30.0:   # Plasma cutting
            v_cut = (8.9212*t_p**2.0 - 486.87*t_p + 8155.8)     # [mm/min]
        else:   # Flame cutting
            v_cut = -4.1939*t_p + 658.67                        # [mm/min]
            
        TP = cut_length/v_cut                                   # [min]
        self.times["productive"] = TP
        return TP
    
    def Consumables(self,t_p):
        """ Calculate the cost of consumables
            input: t_p .. plate thickness
            
        """
        # Consumables: Gases, plasma electrodes and nozzles
        if t_p <= 30.0:   # Plasma cutting
            # Total consumables cost [e/min]
            cutting_oxygen = 0.071                                       # [m^3/min]
            oxygen_cost = 0.3                                            # [e/min]
            nozzle_unit_cost = 40                                        # [e]
            nozzle_life = 480                                            # [min]
            nozzle_cost = nozzle_unit_cost/nozzle_life                   # [e/min]
            self.productive_costs["consumables"] = cutting_oxygen*oxygen_cost + nozzle_cost
            
        else:   # Flame cutting
            # Total consumables cost [e/min]
            propane_consump = 0.0063                                        # [m^3/min]
            propane_cost = 18.40                                            # [e/m^3]
            preheat_oxygen_consump = 0.025                                  # [m^3/min]
            cutting_oxygen_consump = 1.0e-5*t_p**2.0 + 0.001*t_p + 0.0224   # [m^3/min]
            oxygen_cost = 4.18                                              # [e/m^3]
            self.productive_costs["consumables"] = propane_consump*propane_cost + (preheat_oxygen_consump + cutting_oxygen_consump)*oxygen_cost        
    
    def Energy(self,t_p):
        """ Cost of energy used by the machine """
        
        if t_p <= 30:
            # Plasma cutting
            P_total = 72.0
        else:
            # Flame cutting
            P_total = 0.0
            
        # Cost of energy [e/min]
        self.productive_costs["energy"] = P_total*cost_energy_min

    def cost(self,t_p,cut_length):
        """ Evaluate different costs and run the total cost
            calculation
        """
        self.ProductiveTime(t_p, cut_length)
        self.Consumables(t_p)
        self.Energy(t_p)
        return super().cost()

class PaintingCost(CostCentre):    

    def __init__(self, length=15.0,width=5.0, workers=1.0,
                 paint = "alkyd"):
        """ NOTE: Drying costs are missing, 26.1.2022 """

        CostCentre.__init__(self, length,width,workers,
                            eq_invest=0.0,eq_resale=0.0,eq_maint=0.0)        
        
        self.name = "Painting"
        
        self.paint = paint
        # flow of paint through the tip of the paint gun
        # (mm3/min)
        self.tip_flow = tip_flow
            
        self.costs["labour"] = self.nw*wage_per_minute["machine"]
            

    def productive_time(self,Ap):
        """ Productive time of painting
            input: Ap .. area to be painted (mm2)
        """
        # Productive time [min], painting
        T_P = paints[self.paint]["TPP"]*Ap       
        
        self.times["productive"] = T_P
        return T_P

    def consumables(self,Ap):
        """ Cost of painting consumables
            input: Ap .. area to be painted (mm2)
        """
        # Consumables: Paint
        # Total cost of non-time-related consumables [e]
        C_C = paints[self.paint]["CP"]*Ap
       
        self.consumables_cost = C_C
        
    def cost(self,Ap):
        """ Evaluate total painting cost """
        
        self.productive_time(Ap)
        self.consumables(Ap)
        
        return super().cost()
        
class AssemblingCostWeld(CostCentre):    
    """ Cost of welding an assembly """
    
    def __init__(self, length=15.0,width=5.0, workers=1.0,
                 eq_invest=5000.0,eq_resale=0.0,eq_maint=0.0,
                            eq_year=10,
                            weld_type="fillet",
                            weld_length=250.0,
                            weld_size=4.0,
                            steel_grade="S355"):
        """
            input:  length, width .. dimensions of cost centre
                    workers .. number of workers
                    eq_invest .. equipment investment cost
                    eq_resale .. equipment resale value
                    eq_maint .. equipment maintenance cost (€/year)
                    eq_year .. life of equipment (years)
                    weld_type .. "fillet" or "butt"
                    weld_length .. (mm)
                    weld_size .. (mm)
                    steel_grade .. of steel piece to be welded. This determines
                                    the cost of the wire
        """
        CostCentre.__init__(self, length,width,workers,
                            eq_invest,eq_resale,eq_maint,eq_year)                            
        
        self.name = "Assembly welding"
        
        self.weld_type = weld_type
        self.weld_length = weld_length
        self.weld_size = weld_size
        
        number_of_tacks = 4
        
        self.tacking_times = {"debur":0.5,
                         "insertion":0.5,
                         "clamp":0.27,
                         "tack":0.09*number_of_tacks+0.05}
        
        self.power = 200*30*1e-3 # current * voltage (kW)
        self.productive_costs["energy"] = self.power*cost_energy_min
        
        self.gas_consumption = 0.4  # m3/kg wire
        self.gas_price = 3.79      # euro/m3 gas        
        self.wire_price = wire_prices[steel_grade] # euro/kg wire            
        
        self.tacking_time = sum(self.tacking_times.values())
        
        self.costs["labour"] = self.nw*wage_per_minute["machine"]
        self.times["productive"] = self.welding_time()+self.tacking_time
        self.consumables()
            
    def welding_time(self):
        """ Time required for welding [min]
            input: weld_size .. throat thickness of fillet weld or
                                depth of single-bevel butt weld
                   weld_length [mm]
                   weld_type .. "fillet" or "butt"
        """
        if self.weld_type == "fillet":
            a = self.weld_size            
            tweld = self.weld_length*1e-3*(0.4988*a**2-0.0005*a+0.0021)            
        elif self.weld_type == "butt":
            b = self.weld_size
            tweld = self.weld_length*1e-3*(0.249*b**2+0.0096*b-0.0506)
        
        return tweld

    def productive_time(self):
        """ Productive time of welding
            input: Ap .. area to be painted (mm2)
        """
        # Productive time [min], painting        
        self.times["productive"] = self.welding_time()+self.tacking_time
        

    def consumables(self):
        """ Cost of welding consumables            
        """
        
        if self.weld_type == "fillet":
            a = self.weld_size
            weld_weight = self.weld_length*a**2*7850e-9
        elif self.weld_type == "butt":
            b = self.weld_size
            weld_weight = self.weld_length*b**2/2*7850e-9
                
        self.consumables_cost = weld_weight*(self.wire_price + self.gas_consumption*self.gas_price)
                    
        return self.consumables_cost
    
    def cost(self,weld_size,weld_length,weld_type='fillet'):
        """ Evaluate total painting cost """
        
        self.weld_size = weld_size
        self.weld_length = weld_length
        self.weld_type = weld_type
        
        self.consumables()
        self.productive_time()
        
        return super().cost()
        
class SawingCost(CostCentre):    
    """ Cost of sawing """
    
    def __init__(self, length=35.0,width=15.0, workers=1.0,
                 eq_invest=310e3,eq_resale=0.0,eq_maint=1e3,eq_year=20):
        """
            input:  length, width .. dimensions of cost centre
                    workers .. number of workers
                    eq_invest .. equipment investment cost
                    eq_resale .. equipment resale value
                    eq_maint .. equipment maintenance cost (€/year)
                    eq_year .. life of equipment (years)                    
        """
        CostCentre.__init__(self, length,width,workers,
                            eq_invest,eq_resale,eq_maint,eq_year)                            
        
        self.name = "Sawing"
        
        # Here, a basic value for the non-productive time is calculated.
        # It consists of setting up the machine, fixing the profile,
        # descending and ascending the blade for both ends
        self.times["non-productive"] = sawing_times["setup"]+\
            2*(sawing_times["blade_descent"]+sawing_times["blade_raising"]+\
               sawing_times["profile_fixing"])
          
        self.conveyor_speed = 20e3 # mm/min
        self.blade_cost = 170 # euro
        self.Fs = 1350
        
        self.power = 10 # (kW)
        self.productive_costs["energy"] = self.power*cost_energy_min
                
        self.costs["labour"] = self.nw*wage_per_minute["machine"]
        
        #self.consumables()
        
    def non_productive_time(self,L,bevels=0):
        """ Calculate non-productive time 
            input: bevels .. number of bevels (0,1 or 2)
                    L .. length of the member [mm]
        """
        
        self.times["non-productive"] = sawing_base_time + saw_bevel_time*bevels + L/self.conveyor_speed
        
        return self.times["non-productive"]    
        
    def productive_time(self,steel_grade,h,t,Ah):
        """ Productive time of sawing
            input: h .. vertical part of the material, where thickness of sawed
                        material is comparable to plate thickness [mm]
                   t ..plate thickness corresponding to h
                Ah .. area of the horizontal part of the profile [mm2]
        """
        #self.times["productive"] = 0.0
        Sm = sawing_parameters[steel_grade]["Sm"]
        Q = sawing_parameters[steel_grade]["Q"]
        # calculate productive time for obth ends
        #for i, T in enumerate(t):
        #S = self.feeding_speed(T)
            #TPadd = h[i]/Sm/S + Ah[i]/Q
            #self.times["productive"] += TPadd            
        S = self.feeding_speed(t)
        TP = h/Sm/S + Ah/Q
        return TP
        
        
    def consumables(self,At,t,TP,steel_grade='S355',verb=False):
        """ Cost of sawing consumables 
            Input:
                At .. cross-sectional area of the sawed profile
                t .. wall thickness of the sawed profile
        """
        
        St = self.blade_durability(t,steel_grade)
        
        if verb:
            print(f'St = {St:.2f}')
            print(f'At = {At:.2f}')
            print(f'TP = {TP:.2f}')
            print(f'Saw blade = {self.blade_cost:.2f}')
        
        return At/St/TP*self.blade_cost
        
        #return self.consumables_cost
    
    def feeding_speed(self,t):
        """ Speed of the blade
            input: t .. thickness through which the blade travels [mm]
            output: speed [mm/min]
        """
        
        if t <= 5.0:
            S = 120
        elif t <= 10:
            S = 100
        elif t <= 15:
            S = 90
        elif t <= 20:
            S = 80
        elif t <= 25:
            S = 70
        elif t <= 30:
            S = 60
        elif t <= 35:
            S = 50
        else:
            S = 40
        
        return S
    
    def Fsp(self,t):
        """ Sawing parameter Fsp
            input: t .. thickness through which the blade travels [mm]
            output: parameter Fsp
        """
        
        if t <= 5.0:
            Fsp = 0.4
        elif t <= 10:
            Fsp = 0.45
        elif t <= 15:
            Fsp = 0.5
        elif t <= 20:
            Fsp = 0.55
        elif t <= 25:
            Fsp = 0.6
        elif t <= 30:
            Fsp = 0.65
        elif t <= 35:
            Fsp = 0.7
        else:
            Fsp = 0.8
        
        return Fsp
    
    def blade_durability(self,t,steel_grade="S355"):
        """ Durability of the blade 
            t .. thickness of the plate to be sawn
        """
        
        Q = sawing_parameters[steel_grade]["Q"]
        Fsp = self.Fsp(t)
        St = Q*self.Fs*Fsp
        
        return St
    
    def cost(self,L,bevels, steel_grade, H, T, AH, AT):
        """
        

        Parameters
        ----------
        L : float
            Length of member.
        bevels : float/int
            Number of bevels (0, 1 or 2).
        steel_grade : string
            Steel grade, one of the values in sawing_parameters dict.
        h : float
            height of the sawn portion, where the thickness of the plate
            is close to wall thickness.
        T : list
            wall thicknesses of the "thin" portion of the sawn part at each end of the member
        AH : list
            area of the "solid" part of the section at each end.
        AT : list
            total cross_sectional area of the profile at each end.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        
        
        TNP = self.non_productive_time(L,bevels)
        #TP  = self.productive_time(steel_grade, H, T, AH)
        TP = 0.0
        cCS = 0.0
        # calculate productive time and cost of consumables for each end.   
        for t, h, Ah, At in zip(T,H,AH,AT):     
            TPnew = self.productive_time(steel_grade, h, t, Ah)
            cCS += self.consumables(At, t, TPnew, steel_grade)
            TP += TPnew            
                    
        self.times["productive"] = TP
        self.consumables_cost = cCS
        return super().cost()        

class Workshop:
    """
        Class for storing workshop data and doing cost calculations
        
        Data to be contained:
            Dict of cost centres.
        
    """
    
    def __init__(self,steel_price={"plates":BASIC_STEEL_PRICE},cost_centres=None):
        """ Constructor """        
    
        self.name = "Steel Workshop"
        
        self.steel_price = steel_price
    
        if cost_centres is None:
            self.cost_centres = {'blasting':BlastingCost(),
                                 'cutting':CuttingCost(),
                                 'assembly_welding':AssemblingCostWeld(),
                                 'sawing':SawingCost(),
                                 'painting':PaintingCost(),
                                 'beam_welding':BeamWeldingCost()}
        else:
            self.cost_centres = {}
            for key, value in cost_centres.items():
                self.cost_centres[key] = value
        
    
    def info(self):
        """ Prints a bunch of information regarding the workshop """
        
        for centre in self.cost_centres.values():
            centre.info()
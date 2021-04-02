# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 08:20:46 2020

Data for loads and load combinations

EN 1990 and EN 1991

@author: kmela
"""

import numpy as np
import copy
from abc import ABC, abstractmethod

# RIL 201-1-2008 Osa 0 Taulukko A1.1
comb_factors = {
    'A': [0.7, 0.5, 0.3],
    'B': [0.7, 0.5, 0.3],
    'C': [0.7, 0.7, 0.3],
    'D': [0.7, 0.7, 0.6],
    'E': [1.0, 0.9, 0.8],
}

# RIL 205-2-2009 Taulukko 2.2
comb_factors_fire = {
    'A': [0.3, 0.3],
    'B': [0.3, 0.3],
    'C': [0.3, 0.3],
    'D': [0.6, 0.6],
    'E': [0.8, 0.8],
    'F': [0.6, 0.6],
}


class Load:
    """ General class for loads """
    
    def __init__(self,value=1.0,comb_factors=[1.0,1.0,1.0]):
        """ Constructor """
        
        self.value = value
        self.comb_factors = comb_factors
        self.final_value = None
        
    def combined(self,combi=0):
        """ Returns the combined value
            :param combi: integer (0, 1, 2) for choosing the combination factor
        """
        return self.value*self.comb_factors[combi]

    def __copy__(self):
        return Load(self.value, self.comb_factors)

    def __repr__(self):
        return f"Load"

class AreaLoad(Load):
    def __init__(self, value=1.0, comb_factors=[1.0, 1.0, 1.0]):
        super().__init__(value, comb_factors)

    def __copy__(self):
        return AreaLoad(self.value, self.comb_factors)

    def __repr__(self):
        return "Area Load"

class ImposedLoad(Load):
    """ Class for imposed (live) loads """
    
    def __init__(self,value_ground=2.0,comb_factors=[0.7,0.5,0.2]):
        """ Constructor
        
            :param value_ground: value of snow load on ground [kN/m2]
        """
        
        super().__init__(value_ground,comb_factors)

    def __copy__(self):
        return ImposedLoad(self.value, self.comb_factors)

    def __repr__(self):
        return f"Imposed Load"


class SnowLoad(Load):
    """ Class for snow loads """
    
    def __init__(self,value_ground=2.0,comb_factors=[0.7,0.5,0.2]):
        """ Constructor
        
            :param value_ground: value of snow load on ground [kN/m2]
        """
        
        if value_ground < 2.75:
            comb_factors = [0.7,0.4,0.2]
        else:
            comb_factors = [0.7,0.5,0.2]
            
        super().__init__(value_ground,comb_factors)

    def __copy__(self):
        return SnowLoad(self.value, self.comb_factors)

    def __repr__(self):
        
        return "Snow load"


class WindLoad(Load):
    """ Class for wind loads """
    
    def __init__(self,value=0.6,comb_factors=[0.6,0.2,0.0]):
        """ Constructor
        
            :param value: basic value of wind load [kN/m2]
        """
            
        super().__init__(value,comb_factors)

    def __copy__(self):
        return WindLoad(self.value, self.comb_factors)

    def __repr__(self):
        
        return "Wind load"


class Combiner(ABC):
    def __init__(self, structure_class='E', CC=2, kk_jako=1):
        self.structure_class = structure_class
        self.CC = CC
        if CC == 3:
            self.kfi = 1.1
        elif CC == 2:
            self.kfi = 1.0
        elif CC == 1:
            self.kfi = 0.9
        self.kk_jako = kk_jako

        self.G_k = []
        self.Q_k = []

    def get_comb_factors(self):
        return comb_factors[self.structure_class]

    def add(self, *args):
        for l in args:
            load = [copy.copy(l[0]), l[1], l[2]]
            if isinstance(load[0], (SnowLoad, WindLoad, ImposedLoad)):
                if isinstance(load[0], ImposedLoad):
                    load[0].comb_factors = self.get_comb_factors()
                load[0].value *= self.kk_jako
                self.Q_k.append(load)
            elif isinstance(load[0], (Load, AreaLoad)):
                if isinstance(load[0], AreaLoad):
                    load[0].value *= self.kk_jako
                load[0].comb_factors = self.get_comb_factors()
                self.G_k.append(load)

    def sum_G_k(self):
        s = sum([G[0].value for G in self.G_k])
        return s

    def sum_Q_k(self, kahveliID):
        s = sum([Q[0].combined(kahveliID) for Q in self.Q_k])
        return s

    def sum_Q_k_excluding(self, exclude_load, kahveliID):
        mod = self.Q_k.copy()
        filtered = filter(lambda q: q[0] != exclude_load, mod)
        s = sum([Q[0].combined(kahveliID) for Q in filtered])
        return s

    def use_reduction_factor(self, length):
        """
        Reduce Imposed load by a factor alpha_A

        @param length: In meters
        """
        if self.structure_class in ('A', 'B', 'C', 'D'):
            for load in self.Q_k:
                if isinstance(load[0], ImposedLoad):
                    alpha = (5/7) * load[0].comb_factors[0] + 10 / (self.kk_jako * length)
                    alpha = np.clip(alpha, 0.7, 1.0)
                    load[0].value *= alpha

    @staticmethod
    def get_single_result(res):
        result = []
        if 'G' in res['combination']:
            for load in res['combination']['G']:
                load['load'][0].final_value = load['load'][0].value * load['factor']
            result.extend(list(map(lambda x: x['load'], res['combination']['G'])))
        if 'Q1' in res['combination']:
            for load in res['combination']['Q1']:
                load['load'][0].final_value = load['load'][0].value * load['factor']
            result.extend(list(map(lambda x: x['load'], res['combination']['Q1'])))
        if 'Qi' in res['combination']:
            for load in res['combination']['Qi']:
                load['load'][0].final_value = load['load'][0].value * load['factor']
            result.extend(list(map(lambda x: x['load'], res['combination']['Qi'])))


        max_load_id = None
        for load in result:
            if max_load_id is None:
                max_load_id = max(load[2])
            if max(load[2]) > max_load_id:
                max_load_id = max(load[2])

        loads_by_id = [None] * (max_load_id + 1)
        qval_y = [None] * (max_load_id + 1)
        qval_x = [None] * (max_load_id + 1)
        for load in result:
            for id in load[2]:
                if loads_by_id[id] is None:
                    loads_by_id[id] = [load]
                else:
                    loads_by_id[id].append(load)
        for i, loads in enumerate(loads_by_id):
            if loads is None:
                continue
            else:
                qx = 0
                qy = 0
                for load in loads:
                    if load[1] == 'x':
                        qx += round(load[0].final_value, 2)
                    elif load[1] == 'y':
                        qy += round(load[0].final_value, 2)
                qval_y[i] = qy
                qval_x[i] = qx
        return [qval_x, qval_y]

    @staticmethod
    def get_single_split_result(self, res):
        def partial_result(part):
            result = []
            for load in res['combination'][part]:
                load['load'][0].final_value = load['load'][0].value * load['factor']
            result.extend(list(map(lambda x: x['load'], res['combination'][part])))

            max_load_id = None
            for load in result:
                if max_load_id is None:
                    max_load_id = max(load[2])
                if max(load[2]) > max_load_id:
                    max_load_id = max(load[2])

            loads_by_id = [None] * (max_load_id + 1)
            qval_y = [None] * (max_load_id + 1)
            qval_x = [None] * (max_load_id + 1)

            for load in result:
                for id in load[2]:
                    if loads_by_id[id] is None:
                        loads_by_id[id] = [load]
                    else:
                        loads_by_id[id].append(load)
            for i, loads in enumerate(loads_by_id):
                if loads is None:
                    continue
                else:
                    qx = 0
                    qy = 0
                    for load in loads:
                        if load[1] == 'x':
                            qx += round(load[0].final_value, 2)
                        elif load[1] == 'y':
                            qy += round(load[0].final_value, 2)
                    qval_y[i] = qy
                    qval_x[i] = qx
            return [qval_x, qval_y]

        split_result = {}
        if 'G' in res['combination']:
            split_result['G'] = partial_result('G')
        if 'Q1' in res['combination']:
            split_result['Q1'] = partial_result('Q1')
        if 'Qi' in res['combination']:
            split_result['Qi'] = partial_result('Qi')
        return split_result

    @abstractmethod
    def combine(self):
        pass

    @abstractmethod
    def get_result(self):
        pass


class ULSCombiner(Combiner):
    def combine(self):
        vals = []
        case1 = []
        case3 = []
        for l in self.G_k:
            case1.append({'load': l, 'factor': 1.35 * self.kfi})
            case3.append({'load': l, 'factor': 1.15 * self.kfi})

        vals.append({'value': 1.35 * self.kfi * self.sum_G_k(), 'combination': {'G': case1}})

        for l in self.Q_k:
            case2q1 = []
            case2qi = []
            case2q1.append({'load': l, 'factor': 1.5 * self.kfi})
            for lo in self.Q_k:
                if l[0] is not lo[0]:
                    case2qi.append({'load': lo, 'factor': 1.5 * self.kfi * lo[0].comb_factors[0]})

            vals.append({'value': 1.15 * self.kfi * self.sum_G_k() + 1.5 * self.kfi * l[0].value +
                                  1.5 * self.kfi * self.sum_Q_k_excluding(l[0], 0), 'combination': {'G': case3, 'Q1': case2q1, 'Qi': case2qi}})

        return max(vals, key=lambda x: x['value'])

    def get_result(self):
        """

        @return: returns the combination with the highest value
        """
        if len(self.Q_k) + len(self.G_k) <= 0:
            raise Exception('Give at least one Load with add method')
        res = self.combine()
        return self.get_single_result(res)


class SLSCombiner(Combiner):
    def combine(self):
        characteristic_combination = []
        frequent_combination = []
        quasi_permanent_combination = []
        case1 = []
        for l in self.G_k:
            case1.append({'load': l, 'factor': 1})

        if len(self.Q_k) == 0:
            characteristic_combination.append({'value': self.sum_G_k(),
                                               'combination': {'G': case1}})
            frequent_combination.append({'value': self.sum_G_k(),
                                         'combination': {'G': case1}})
        case4 = []
        for l in self.Q_k:
            case2q1 = []
            case2qi = []
            case3q1 = []
            case3qi = []
            case2q1.append({'load': l, 'factor': 1})
            case3q1.append({'load': l, 'factor': l[0].comb_factors[1]})
            case4.append({'load': l, 'factor': l[0].comb_factors[2]})

            for lo in self.Q_k:
                if l is not lo:
                    case2qi.append({'load': lo, 'factor': lo[0].comb_factors[0]})
                    case3qi.append({'load': lo, 'factor': lo[0].comb_factors[2]})
            characteristic_combination.append({'value': self.sum_G_k() + l[0].value + self.sum_Q_k_excluding(l[0], 2),
                                               'combination': {'G': case1, 'Q1': case2q1, 'Qi': case2qi}})
            frequent_combination.append({'value': self.sum_G_k() + l[0].combined(1) + self.sum_Q_k_excluding(l[0], 2),
                                         'combination': {'G': case1, 'Q1': case3q1, 'Qi': case3qi}})

        quasi_permanent_combination.append({'value': self.sum_G_k() + self.sum_Q_k(2),
                                            'combination': {'G': case1, 'Qi': case4}})

        for cc in characteristic_combination:
            if 'G' in cc['combination']:
                if len(cc['combination']['G']) == 0:
                    cc['combination'].pop('G')
            if 'Q1' in cc['combination']:
                if len(cc['combination']['Q1']) == 0:
                    cc['combination'].pop('Q1')
            if 'Qi' in cc['combination']:
                if len(cc['combination']['Qi']) == 0:
                    cc['combination'].pop('Qi')

        for fc in frequent_combination:
            if 'G' in fc['combination']:
                if len(fc['combination']['G']) == 0:
                    fc['combination'].pop('G')
            if 'Q1' in fc['combination']:
                if len(fc['combination']['Q1']) == 0:
                    fc['combination'].pop('Q1')
            if 'Qi' in fc['combination']:
                if len(fc['combination']['Qi']) == 0:
                    fc['combination'].pop('Qi')

        for qc in quasi_permanent_combination:
            if 'G' in qc['combination']:
                if len(qc['combination']['G']) == 0:
                    qc['combination'].pop('G')
            if 'Qi' in qc['combination']:
                if len(qc['combination']['Qi']) == 0:
                    qc['combination'].pop('Qi')

        return [max(characteristic_combination, key=lambda x: x['value']),
                max(frequent_combination, key=lambda x: x['value']),
                max(quasi_permanent_combination, key=lambda x: x['value'])]

    def get_result(self, split_result=False):
        """
        EN 1990 (6.14) characteristic (6.15) frequent (6.16) Quasi-permanent

        Combines loads for Service Limit State and returns combination with highest value

        @return: returns three results in a list.
                 First result is for characteristic combination,
                 second is frequent combination,
                 third is quasi permanent combination.
                 Each result consists of two parts: x and y values. Those values are for given ids.
        """
        if len(self.Q_k) + len(self.G_k) <= 0:
            raise Exception('Give at least one Load with add method')
        if split_result:
            return [self.get_single_split_result(self, res) for res in self.combine()]
        final_result = [self.get_single_result(res) for res in self.combine()]
        return final_result


class ACCCombiner(Combiner):
    def combine(self):
        vals = []
        case1 = []
        for l in self.G_k:
            case1.append({'load': l, 'factor': 1})

        for l in self.Q_k:
            case2q1 = []
            case2qi = []
            if isinstance(l[0], (WindLoad, SnowLoad)):
                case2q1.append({'load': l, 'factor': l[0].comb_factors[1]})
            else:
                case2q1.append({'load': l, 'factor': l[0].comb_factors[2]})

            for lo in self.Q_k:
                if l is not lo:
                    case2qi.append({'load': lo, 'factor': lo[0].comb_factors[2]})
            if isinstance(l[0], (WindLoad, SnowLoad)):
                vals.append({'value': self.sum_G_k() + l[0].combined(1) + self.sum_Q_k_excluding(l[0], 2),
                            'combination': {'G': case1, 'Q1': case2q1, 'Qi': case2qi}})
            else:
                vals.append({'value': self.sum_G_k() + l[0].combined(2) + self.sum_Q_k_excluding(l[0], 2),
                             'combination': {'G': case1, 'Q1': case2q1, 'Qi': case2qi}})
        else:
            vals.append({'value': self.sum_G_k(), 'combination': {'G': case1}})

        for cc in vals:
            if 'G' in cc['combination']:
                if len(cc['combination']['G']) == 0:
                    cc['combination'].pop('G')
            if 'Q1' in cc['combination']:
                if len(cc['combination']['Q1']) == 0:
                    cc['combination'].pop('Q1')
            if 'Qi' in cc['combination']:
                if len(cc['combination']['Qi']) == 0:
                    cc['combination'].pop('Qi')

        return max(vals, key=lambda x: x['value'])

    def get_result(self):
        """

        @return: returns the combination with the highest value
        """
        if len(self.Q_k) + len(self.G_k) <= 0:
            raise Exception('Give at least one Load with add method')
        res = self.combine()
        return self.get_single_result(res)


class MultiCombiner:
    def __init__(self, structure_class='E', CC=2, kk_jako=1):
        self.ULS = ULSCombiner(structure_class, CC, kk_jako)
        self.SLS = SLSCombiner(structure_class, CC, kk_jako)
        self.Acc = ACCCombiner(structure_class, CC, kk_jako)

    def add(self, *args):
        for arg in args:
            self.ULS.add(arg)
            self.SLS.add(arg)
            self.Acc.add(arg)

    def use_reduction_factor(self, val):
        self.ULS.use_reduction_factor(val)
        self.SLS.use_reduction_factor(val)
        self.Acc.use_reduction_factor(val)

    def get_result(self):
        uls_res = self.ULS.get_result()
        sls_res = self.SLS.get_result()
        acc_res = self.Acc.get_result()
        return [uls_res, sls_res, acc_res]


# OP = 1      kN/m    1   G_k
# KR = 0.5    kN/m2   5   G_k
# LK = 2.5    kN/m2   25  Q_k
# TK = 1.5    kN/m2   15  Q_k
# SP = 0.2    kN/m2   2   Q_k
#
# kk-jako = 10 m
#
# CC2
#
# Class='A'
#
#
# CC3 = K_fi = 1.1
# CC2 = K_fi = 1.0
# CC1 = K_fi = 0.9
#
# Rakenneluokka(Class='A', CC3)
#
# Palkki = 1 kN/m
# Välipohja = 0.5 kN/m¨2
#
# Omatpainot(palkki, välipohja)
#
# Snow_load(<2,75kN/m2) = [0.7, 0.4, 0.2]
# Snow_load(>=2,75kN/m2) = [0.7, 0.5, 0.2]
#
# Wind_load = [0.6, 0.2, 0.0]
#
# ULS:            # Murtorajatila
# 1.15 * K_fi * summa(G_k) + 1.5 * K_fi * Q_k1 + 1.5 * K_fi * summa(kahveli_0 * Q_kj)
# 1.35 * K_fi * summa(G_k)            vähintään tämän verran
#
# SLS:            # Käyttörajatila
# a) ominaisyhdistelmä        summa(G_k) + Q_k1 + summa(kahveli_0 * Q_ki)
# b) tavallinen yhdistelmä    summa(G_k) + kahveli_1 * Q_k1 + summa(kahveli_2 * Q_ki)
# c) pitkäaikaisyhdistelmä    summa(G_k) + summa(kahveli_2 * Q_ki)
#
# ACC:            # Onnettomuusrajatila
# summa(G_k) + kahveli_1 * Q_k1 + summa(kahveli_2 * Q_ki)    jos Q_k1 = lumi tai tuuli
# summa(G_k) + kahveli_2 * Q_k1 + summa(kahveli_2 * Q_ki)     muussa tapauksessa

if __name__ == "__main__":
    import sections.timber.timber_section as ts

    SW = Load(9.31)

    ont = AreaLoad(2.6)
    vesi = AreaLoad(0.5)
    Qs = SnowLoad(2.0)
    Qw = WindLoad(0.5)
    Hk = ImposedLoad(3)

    M = MultiCombiner('A', kk_jako=8)
    M.add([SW, 'y', [0]])
    M.add([ont, 'y', [0]])
    M.add([Qs, 'y', [0]])
    M.add([vesi, 'y', [0]])
    M.add([Hk, 'y', [0]])
    M.use_reduction_factor(16)
    result = M.get_result()

    # result[ULS/SLS/Acc][x/y][index]
    #print(f'ULS result: {result[0]}')
    print(f'ULS result: x: {result[0][0][0]}, y: {result[0][1][0]}')
    #print(f'SLS result: {result[1]}')

    # result[SLS = 1][ominais/tavallinen/pitkäaikais][x/y][index]
    print(f'SLS result: ominais: x: {result[1][0][0][0]}, y: {result[1][0][1][0]}')
    print(f'SLS result: tavallinen: x: {result[1][1][0][0]}, y: {result[1][1][1][0]}')
    print(f'SLS result: pitkäaikais: x: {result[1][2][0][0]}, y: {result[1][2][1][0]}')
    #print(f'Acc result: {result[2]}')
    print(f'Acc result: x: {result[2][0][0]}, y: {result[2][1][0]}')

    # ULS = ULSCombiner('E', kk_jako=8)
    # ULS.add([SW, 'y', [0]])
    # #ULS.add([Hk, 'y', [0]])
    # ULS.add([Qs, 'y', [0]])
    # ULS.add([ont, 'y', [0]])
    # #ULS.add([Qw, 'x', [0, 4]])
    # #ULS.add([ont, 'y', [0]])
    # #ULS.use_reduction_factor(16)
    # res = ULS.get_result()
    # print(res)
    # #
    # SLS = SLSCombiner('A', kk_jako=8)
    # SLS.add([SW, 'y', [0]])
    # #SLS.add([Qw, 'x', [0,4]])
    # #SLS.add([Hk, 'y', [0]])
    # SLS.add([Qs, 'y', [0]])
    # SLS.add([ont, 'y', [0]])
    # #SLS.add([ont, 'y', [0]])
    # #SLS.use_reduction_factor(16)
    # ress = SLS.get_result()
    # print(ress)
    #
    # Acc = ACCCombiner('A', kk_jako=10)
    # Acc.add([SW, 'y', [0]])
    # #Acc.add([Qw, 'x', [0]])
    # Acc.add([Qs, 'y', [0]])
    # Acc.add([Hk, 'y', [0]])
    # Acc.add([ont, 'y', [0]])
    # #Acc.add([Hk, 'y', [0]])
    # #Acc.use_reduction_factor(16)
    # resss = Acc.get_result()
    # print(resss)
    # # #print(Acc.combine())


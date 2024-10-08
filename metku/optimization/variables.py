# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:56:32 2022

Optimization variable classes

@author: kmela
"""

import math

from typing import Optional
from metku.optimization.constants import XLB, XUB, INT_TOL, DISC_TOL


class Variable:
    """ Class for optimization variables
    """

    def __init__(self,
                 name: str,
                 lb: float = XLB,
                 ub: float = XUB,
                 target: Optional[dict[str, object | list]] = None,
                 target_property: Optional[str] = None,
                 target_properties: Optional[list[str]] = None,
                 target_object: Optional[object] = None,
                 target_objects: Optional[list[object]] = None,
                 profiles: Optional[list] = None,
                 value: Optional[float] = None,
                 target_fun: Optional[callable] = None):
        """
        Parameters:
            -----------
                :param name: string stating the name of the variable
                :param value: current value of the variable
                :param lb: lower bound
                :param ub: upper bound
                :param target: dict that provides information about the part of
                the structure that the variable is changing
                :param profiles: list of available profiles (optional)

                target = {"property": string, "objects": list}
                    property .. unique identifier that shows which property of
                                the structure is affected by the variable
                                
                                this can also be a list of strings
                    objects .. list of objects (e.g. frame members) that are
                                affected by the variable

                if target is 'None', then the variable does not change anything
                in the structure
                :param target_fun: (optional) dict that provides functions for evaluating a functional relationship
                                    between the variable and its target. For cases, where the variable substitution
                                    is not a direct attribute insertion. For example, when the cross-sectional area
                                    of a member is used as a variable, and it also affects the second moment of area
                                    by some function I = I(A)

                Allowed values for 'property':
                    'AREA' .. cross-sectional area of a member
                    'IY' .. second moment of area of a member (major axis)
                    'IZ' .. second moment of area of a member (minor axis)
                    'WELY' .. elastic section modulus of a member (major axis)
                    'WELZ' .. elastic section modulus of a member (minor axis)
                    'WPLY' .. plastic section modulus of a member (major axis)
                    'WPLZ' .. plastic section modulus of a member (minor axis)
                    'PROFILE' .. choosing a profile from a given catalogue
                    'BF' .. width of flange
                    'TF' .. thickness of flange
                    'H' .. height of the profile
                    'TW' .. thickness of the web
        """

        self.name = name
        self.value = value
        self.lb = lb
        self.ub = ub

        if any([target_property,target_properties, target_objects, target_object]):

                if not any([target_properties, target_property]):
                    raise ValueError("target_property/properties cannot be None "
                                     "if target_object/objects is defined!")

                if not any([target_objects, target_object]):
                    raise ValueError("target_object/objects cannot be None "
                                     "if target_property/properties is defined!")

                if target_properties is None and target_property is not None:
                    target_properties = [target_property]
                if target_objects is None:
                    target_objects = [target_object]

                target = {"property": target_properties,
                          "objects": target_objects}

        self.target = target
        self.profiles = profiles
        self.locked = False
        self.branch_priority = 0
        self.scaling = 1.0
        self.target_fun = target_fun

        # Lagrange multiplier for lower bound
        self.lb_mult = 1.0
        # Lagrange multiplier for upper bound
        self.ub_mult = 1.0

    def __repr__(self):

        # return self.name + ": [" + str(self.lb) + "," + str(self.ub) + "]"
        if self.locked:
            fixed = " (fixed to {0:g})".format(self.value)
        else:
            fixed = ""

        return str(self.lb) + " <= " + self.name + " <= " + str(self.ub) + fixed

    def lock(self, val=None):
        """
        Locks the variable, so it's value can't be changed

        Parameters:
        -----------
        :param val: value to give to variable before locking
        """
        if val is not None:
            self.substitute(val)

        self.locked = True

    def unlock(self):
        """
        Unlocks the variable so it's value can be changed
        """
        self.locked = False

    def scale(self, scaling_factor):
        """ Scale variable by the factor 'scaling_factor' 
            This alters the variable bounds. In 'substitute' method, scaling
            is reflected such that the variable value is multiplied by
            'scaling_factor'
        """

        self.scaling = scaling_factor
        self.lb = self.lb / scaling_factor
        self.ub = self.ub / scaling_factor

    def substitute(self, new_value):
        """
        Changes variable value and modifies target
        :param new_value:
        :return:
        """
        if not self.locked:
            """ Substitute a new value for the variable """
            self.value = new_value
            """ Modify target object(s) if any """
            if self.target is not None:
                for obj in self.target['objects']:
                    prop = self.target['property']
                    """ If the variable affects several properties at a time
                        'prop' is a list of property names
                    """
                    if isinstance(prop, list):
                        for p in prop:
                            """ If the variable has a 'target_fun',
                                evaluate it for right properties
                            """
                            if self.target_fun is not None:
                                """ If 'p' is not a property that has a functional dependency,
                                    then do regular substitution
                                """
                                try:
                                    new_val = self.target_fun[p](new_value * self.scaling)
                                except KeyError:
                                    new_val = new_value * self.scaling

                                obj.__setattr__(p, new_val)
                            else:
                                obj.__setattr__(p, new_value)

                    elif isinstance(self, IndexVariable):
                        obj.__setattr__(prop, new_value)
                    elif isinstance(self, BinaryVariable):
                        """ If a binary variable has a target, assume that target
                            is the section to be substituted. The section is
                            given in the 'section' attribute of the variable.
                        """
                        if new_value == 1:
                            if prop == 'profile' or prop == 'cross_section':
                                obj.__setattr__(prop, self.section)
                        else:
                            """ If the target of a binary variable is a member,
                                the binary variable is a member existence variable
                                used in topology optimization. In that case,
                                mark the member as passive.
                            """
                            if prop == 'member':
                                obj.__setattr__('active', False)
                    else:
                        obj.__setattr__(prop, new_value * self.scaling)

    def discrete_violation(self, x):
        """ Identify the violation of value 'x' from allowable (Discrete)
            value. Returns 0 for continuous variables
        :return:

        """
        return 0


class IntegerVariable(Variable):
    """ Class for integer variables

        Integer variables can take integer values from the interval lb, ub
    """

    def __init__(self,
                 name: str = "",
                 lb: int = 0,
                 ub: int = 1e5,
                 **kwargs):
        """ Constructor

            Arguments:
                name .. string stating the name of the variable
                lb .. lower bound
                ub .. upper bound
        """

        super().__init__(name, lb, ub, **kwargs)

    def is_allowed_value(self, x):
        """ Check if 'x' is an integer within the range of the variable """

        if self.lb <= x <= self.ub:
            """ If x is within INT_TOL of its rounded integer value,
                it is deemed integer
            """
            if math.isclose(x, round(x, 0), abs_tol=INT_TOL):
                return True
            else:
                return False
        else:
            return False

    def discrete_violation(self, x):
        """ Identify the violation of value 'x' from allowable (Discrete)
            value. Returns 0 for continuous variables
        """
        violation = abs(x - round(x, 0))
        if violation <= INT_TOL:
            violation = 0

        return violation


class BooleanVariable(IntegerVariable):
    """ Class for boolean variables

        Boolean variables can only have values 0 or 1

    """

    def __init__(self,
                 name="",
                 **kwargs):
        """ Constructor

            Arguments:
                name .. string stating the name of the variable
        """
        super(IntegerVariable, self).__init__(name, 0, 1, **kwargs)


class BinaryVariable(IntegerVariable):
    """ Class for binary variables

        Binary variables can only have values 0 or 1

    """

    def __init__(self, name="", section=None, objects=None, **kwargs):
        """ Constructor

            Arguments:
                name .. string stating the name of the variable
        """
        super(IntegerVariable, self).__init__(name, 0, 1, **kwargs)

        self.section = section
        self.objects = objects


class IndexVariable(IntegerVariable):
    """
    Class for index variables.
    
    An index variable is usually used to select a profile for
    a member.
    
    """

    def __init__(self, name="", values=None, target=None, **kwargs):
        """ Constructor:
            :param values: list of possible 
        """
        self.values = values

        if values is not None:
            lb = 0
            ub = len(values) - 1
        else:
            lb = None
            ub = None

        super().__init__(name, lb, ub, target=target, **kwargs)

    def substitute(self, new_value):
        """ Substitute the profile corresponding to the integer new_value
        """
        try:
            # Take into account that new_value can be a float
            new_value = int(new_value)
            super().substitute(self.values[new_value])
        except:
            raise ValueError(
                f"Input {new_value} is erroneous "
                "IndexVariable's value must be either"
                " index or a value from the given list!\n"
                f"Values: {self.values}")

        self.value = new_value

    # @property
    # def idx(self):
    #     if self.value in self.values:
    #         return self.values.index(self.value)
    #     return self.value


class DiscreteVariable(Variable):
    """ Class for general discrete variables """

    def __init__(self, name="", values=None, target=None, **kwargs):
        """ Constructor

            Parameters:
            -----------
                :param name: string stating the name of the variable
                :param values: list of allowable discrete values

                Variables:
                ----------
                :ivar values: list of allowable discrete values
                :ivar name: name of the variable (string)
        """

        self.values = values

        if values is not None:
            lb = min(values)
            ub = max(values)

        else:
            lb = None
            ub = None

        super().__init__(name, lb, ub, target=target, **kwargs)

    def is_allowed_value(self, x):
        """ Check if 'x' is one of the values listed in
            self.values
        """

        if self.lb <= x <= self.ub:
            """ If x is within INT_TOL of any of the allowed values,
                it is deemed allowed
            """
            if min([abs(val - x) for val in self.values]) < INT_TOL:
                return True

            """
            for val in self.values:
                if abs(x-val) < INT_TOL:            
                    return True        
            """
            return False
        else:
            return False

    def discrete_violation(self, x):
        """ Identify the violation of value 'x' from allowable (Discrete)
            value.
        """

        violation = min([abs(val - x) for val in self.values])

        if violation <= DISC_TOL:
            violation = 0

        return violation

    def smaller_discrete_value(self, x):
        """ Returns discrete value smaller than 'x' and closest to it  """

        return max(list(filter(lambda val: (val - x < 0), self.values)))

        # diff = [val-x for val in self.values]

    def larger_discrete_value(self, x):
        """ Returns discrete value smaller than 'x' and closest to it  """

        return min(list(filter(lambda val: (val - x > 0), self.values)))


class CrossSectionVariable(Variable):
    """ Class for cross-sectional variables (continuous) """

    def __init__(self, name, lb=XLB, ub=XUB, target=None, profiles=None, value=None, target_fun=None):

        super().__init__(name, lb, ub, target, profiles, value, target_fun=target_fun)

    def substitute(self, new_value):
        """
        Changes variable value and modifies target
        :param new_value:
        :return:
        """
        if not self.locked:
            """ Substitute a new value for the variable """
            self.value = new_value

            for obj in self.target['objects']:
                prop = self.target['property']
                """ If the variable affects several properties at a time
                    'prop' is a list of property names
                """
                if isinstance(prop, list):
                    for p in prop:
                        """ If the variable has a 'target_fun',
                            evaluate it for right properties
                        """
                        if self.target_fun is not None:
                            """ If 'p' is not a property that has a functional dependency,
                                then do regular substitution
                            """
                            try:
                                new_val = self.target_fun[p](new_value * self.scaling)
                            except KeyError:
                                new_val = new_value * self.scaling

                            obj.__setattr__(p, new_val)
                        else:
                            obj.__setattr__(p, new_value)
                else:
                    obj.cross_section.__setattr__(prop, new_value * self.scaling)


if __name__ == "__main__":
    cont_var = Variable(
        name="Jatkuva muuttuja",
        lb=0,
        ub=100,
        value=50
    )
    print(cont_var)
    print(cont_var.value)
    cont_var.substitute(75.123456)
    print(cont_var.value)


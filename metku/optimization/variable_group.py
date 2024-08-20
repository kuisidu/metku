from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional

from metku.optimization.variables import Variable, DiscreteVariable, IndexVariable, BinaryVariable, BooleanVariable


class VariableTypeEnum(Enum):
    CONTINUOUS = auto()
    DISCRETE = auto()
    INDEX = auto()
    BINARY = auto()
    BOOLEAN = auto()


def has_attribute(obj: object, attr: list | str) -> bool:
    if isinstance(attr, list):
        for a in attr:
            if not hasattr(obj, a):
                return False
    else:
        return hasattr(obj, attr)
    return True


@dataclass
class VariableGroup:
    name: str
    var_type: VariableTypeEnum
    object: object = None
    objects: list[object] = None
    attribute: Optional[str] = None
    attributes: Optional[list[str]] = None
    value: Optional[list | float] = None
    values: Optional[list] = None
    lower_bounds: Optional[list] = None
    upper_bounds: Optional[list] = None
    variables: list = field(init=False)

    def __post_init__(self):

        if not self.objects:
            self.objects = [self.object]
        if not self.attributes:
            self.attributes = [self.attribute]
        if self.value is None:
            self.value = [0] * len(self.attributes)

        # Attribute check
        for attr in self.attributes:
            for obj in self.objects:
                if not has_attribute(obj, attr):
                    raise TypeError(f"Object {type(obj)} has no attribute: {attr}")

        variables = []
        match self.var_type:

            case VariableTypeEnum.CONTINUOUS:
                for attr, lb, ub, value in zip(self.attributes,
                                               self.lower_bounds,
                                               self.upper_bounds,
                                               self.value):
                    var = Variable(name=f"{self.name} {attr}",
                                   value=value,
                                   lb=lb,
                                   ub=ub,
                                   target={
                                       'property': attr,
                                       'objects': self.objects
                                   })
                    variables.append(var)
            case VariableTypeEnum.DISCRETE:
                for attr, value, values in zip(self.attributes,
                                               self.value,
                                               self.values):
                    var = DiscreteVariable(name=f"{self.name} {attr}",
                                           value=value,
                                           values=values,
                                           target={
                                               'property': attr,
                                               'objects': self.objects
                                           })
                    variables.append(var)
            case VariableTypeEnum.INDEX:
                var = IndexVariable(name=f"{self.name} {self.attribute}",
                                    value=self.value,
                                    values=self.values,
                                    target={
                                        'property': self.attribute,
                                        'objects': self.objects
                                    })
                variables.append(var)
            case VariableTypeEnum.BINARY:
                raise NotImplementedError(f"{VariableTypeEnum.BINARY} Not Implemented")
                # var = BinaryVariable()
            case VariableTypeEnum.BOOLEAN:
                raise NotImplementedError(f"{VariableTypeEnum.BOOLEAN} Not Implemented")
            case _:
                raise ValueError(f"var_type must be VariableTypeEnum, not {self.var_type}")

        self.variables = variables


if __name__ == '__main__':
    from metku.frame2d import SteelBeam
    from metku.sections.steel.catalogue import shs_profiles

    sb = SteelBeam([[0, 0], [5000, 0]])

    vg = VariableGroup(name="ContinuousTest",
                       object=sb.cross_section,
                       attributes=[['Iy', 'Iz'], 'A'],
                       var_type=VariableTypeEnum.CONTINUOUS,
                       # value=[10, 100],
                       lower_bounds=[0, 123],
                       upper_bounds=[1e6, 456])

    vg2 = VariableGroup(name="IndexTest",
                        object=sb,
                        attribute='profile',
                        var_type=VariableTypeEnum.INDEX,
                        values=shs_profiles.keys())

    print(vg.variables)

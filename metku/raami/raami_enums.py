from enum import auto, StrEnum
from metku.helpers.upper_str_enum import UpperStrEnum


class MemberTypeEnum(UpperStrEnum):
    COLUMN = auto()
    BEAM = auto()
    BRACE = auto()
class SupportMethodEnum(UpperStrEnum):
    REM = auto()
    ZERO = auto()

class LoadCaseEnum(UpperStrEnum):
    PERMANENT = auto()
    IMPOSED_A = auto()
    IMPOSED_B = auto()
    IMPOSED_C = auto()
    IMPOSED_D = auto()
    IMPOSED_E = auto()
    SNOW = auto()
    WIND = auto()
class LoadCombinationEnum(UpperStrEnum):
    ULS = auto()
    SLS_CH = auto()
    SLS_FR = auto()
    SLS_QP = auto()
    FIRE = auto()
    ACC = auto()


if __name__ == "__main__":
    r = LoadCaseEnum.SNOW
    print(r)
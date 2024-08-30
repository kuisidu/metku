from enum import Enum, auto, StrEnum

class MemberTypeEnum(StrEnum):
    COLUMN = auto()
    BEAM = auto()
    BRACE = auto()
class SupportMethodEnum(StrEnum):
    REM = auto()
    ZERO = auto()

class LoadCaseEnum(StrEnum):

    PERMANENT = auto()
    IMPOSED_A = auto()
    IMPOSED_B = auto()
    IMPOSED_C = auto()
    IMPOSED_D = auto()
    IMPOSED_E = auto()
    SNOW = auto()
    WIND = auto()
class LoadCombinationEnum(StrEnum):
    ULS = auto()
    SLS_CH = auto()
    SLS_FR = auto()
    SLS_QP = auto()
    FIRE = auto()
    ACC = auto()



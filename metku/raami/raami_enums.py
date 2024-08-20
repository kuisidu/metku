from enum import Enum, auto

class MemberTypeEnum(Enum):
    COLUMN = "COLUMN"
    BEAM = "BEAM"
class SupportMethodEnum(Enum):
    REM = "REM"
    ZERO = "ZERO"

class LoadTypeEnum(Enum):

    PERMANENT = auto()
    IMPOSED_A = auto()
    IMPOSED_B = auto()
    IMPOSED_C = auto()
    IMPOSED_D = auto()
    IMPOSED_E = auto()
    SNOW = auto()
    WIND = auto()
class LoadCombinationTypeEnum(Enum):
    ULS = auto()
    SLS_CH = auto()
    SLS_FR = auto()
    SLS_QP = auto()
    ACC = auto()



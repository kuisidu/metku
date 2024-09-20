from enum import StrEnum
class UpperStrEnum(StrEnum):
    """
    Ensures every enum is in uppercase letters
    """
    def __new__(cls, value):
        # Convert the string value to uppercase
        obj = str.__new__(cls, value.upper())
        obj._value_ = obj
        return obj
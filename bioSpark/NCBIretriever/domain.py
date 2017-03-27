import re
from bioSpark.common.util import Constants

class NucleotidesFromNCBI(object):
    """NucleotidesFromNCBI model"""
    def __init__(self, id = '0', sequence = '_'):
        self.id = id
        self.sequence = sequence

    def get_id(self) -> str:
        return self._id

    def set_id(self, id) -> None:
        self._id = id

    def get_sequence(self) -> str:
        return self._sequence

    def set_sequence(self, sequence) -> None:
        # validate DNA sequence
        pattern = re.compile(Constants.REG_EXP_DNA_SEQUENCE)
        if not pattern.fullmatch(sequence):
            raise ValueError(Constants.MSG_ERROR_DNA_SEQUENCE)
        self._sequence = sequence

    id = property(fget= get_id, fset= set_id)
    sequence = property(fget= get_sequence, fset= set_sequence)

# test = NucleotidesFromNCBI(id= "01", sequence="GCATCA")
#
# print(test.sequence)






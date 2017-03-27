import re
from bioSpark.common.util import Constants

import mongoengine

# option 1 --> using ObjectID mongo

class NucleotidesFromNCBI(object):
    """NucleotidesFromNCBI model"""
    def __init__(self, idNcbi ='0', sequence ='_'):
        self.idNcbi = idNcbi
        self.sequence = sequence

    def get_idNcbi(self) -> str:
        return self._idNcbi

    def set_idNcbi(self, id) -> None:
        self._idNcbi = id

    def get_sequence(self) -> str:
        return self._sequence

    def set_sequence(self, sequence) -> None:
        # validate DNA sequence
        pattern = re.compile(Constants.REG_EXP_DNA_SEQUENCE)
        if not pattern.fullmatch(sequence):
            raise ValueError(Constants.MSG_ERROR_DNA_SEQUENCE)
        self._sequence = sequence

    idNcbi = property(fget= get_idNcbi, fset= set_idNcbi)
    sequence = property(fget= get_sequence, fset= set_sequence)

# option 2 --> using id from idNcbi for collections

class NucleotidesFromNCBIusingIdOfSequenceInMongo(object):
    # """NucleotidesFromNCBI model"""
    def __init__(self, id, sequence ='_'):
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





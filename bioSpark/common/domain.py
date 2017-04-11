import re
from bioSpark.common.util import Constants


class NucleotidesFromNCBI(object):
    """NucleotidesFromNCBI model"""

    def __init__(self):
        self._idNcbi = None
        self._sequence = None

    def _get_idNcbi(self) -> str:
        return self._idNcbi

    def _set_idNcbi(self, id) -> None:
        self._idNcbi = id

    def _get_sequence(self) -> str:
        return self._sequence

    def _set_sequence(self, sequence) -> None:
        # validate DNA sequence
        # pattern = re.compile(Constants.REG_EXP_DNA_SEQUENCE)
        # if not pattern.fullmatch(sequence):
        #     raise ValueError(Constants.MSG_ERROR_DNA_SEQUENCE)
        self._sequence = sequence

    idNcbi = property(fget=_get_idNcbi, fset=_set_idNcbi)
    sequence = property(fget=_get_sequence, fset=_set_sequence)


class NCBIsearcher(object):
    def __init__(self):
        self._search_by_NCBI_id_criteria = None

    def _get_NCBIIdSearchCriteria(self) -> str:
        return self._search_by_NCBI_id_criteria

    def _set_NCBIIdSearchCriteria(self, NCBI_id_criteria: str) -> None:
        self._search_by_NCBI_id_criteria = NCBI_id_criteria

    search_by_NCBI_id_criteria = property(fget=_get_NCBIIdSearchCriteria, fset=_set_NCBIIdSearchCriteria)

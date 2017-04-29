class GeneDTO(object):
    """GeneDTO model"""

    def __init__(self):
        self._gene_id = None
        self._sequence = None

    def _get_gene_id(self) -> str:
        return self._gene_id

    def _set_gene_id(self, id) -> None:
        self._gene_id = id

    def _get_sequence(self) -> str:
        return self._sequence

    def _set_sequence(self, sequence) -> None:
        self._sequence = sequence

    gene_id = property(fget=_get_gene_id, fset=_set_gene_id)
    sequence = property(fget=_get_sequence, fset=_set_sequence)


class GeneSearcher(object):
    def __init__(self):
        self._search_by_gene_id_criteria = None

    def _get_gene_id_search_criteria(self) -> str:
        return self._search_by_gene_id_criteria

    def _set_gene_id_search_criteria(self, NCBI_id_criteria: str) -> None:
        self._search_by_gene_id_criteria = NCBI_id_criteria

    search_by_gene_id_criteria = property(fget=_get_gene_id_search_criteria, fset=_set_gene_id_search_criteria)

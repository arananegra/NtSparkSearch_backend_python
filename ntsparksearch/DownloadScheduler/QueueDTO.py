class QueueDTO(object):
    """GeneDTO model"""

    def __init__(self):
        self._gene_id_list_to_download = None

    def _get_gene_id_list_to_download(self) -> str:
        return self._gene_id_list_to_download

    def _set_gene_id_list_to_download(self, list_of_genes_to_download: list) -> None:
        self._gene_id_list_to_download = list_of_genes_to_download

    gene_id_list_to_download = property(fget=_get_gene_id_list_to_download, fset=_set_gene_id_list_to_download)
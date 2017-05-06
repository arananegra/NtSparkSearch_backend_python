from datetime import datetime


class QueueDTO(object):
    """GeneDTO model"""

    def __init__(self):
        self._gene_id_list_to_download = None
        self._entry_date_in_collection = None

    def _get_gene_id_list_to_download(self) -> list:
        return self._gene_id_list_to_download

    def _set_gene_id_list_to_download(self, list_of_genes_to_download: list) -> None:
        self._gene_id_list_to_download = list_of_genes_to_download

    def _get_entry_date_in_collection(self) -> datetime:
        return self._entry_date_in_collection

    def _set_entry_date_in_collection(self, entry_date_in_collection: list) -> None:
        self._entry_date_in_collection = entry_date_in_collection

    gene_id_list_to_download = property(fget=_get_gene_id_list_to_download, fset=_set_gene_id_list_to_download)
    entry_date_in_collection = property(fget=_get_entry_date_in_collection, fset=_set_entry_date_in_collection)

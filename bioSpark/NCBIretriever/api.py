from abc import abstractmethod, ABCMeta
from bioSpark.common.domain import NucleotidesFromNCBI


class INCBIretriever(metaclass=ABCMeta):
    @abstractmethod
    def insert_in_collection_from_excel(self, sheet: int, column_name: str):
        """
        From an excel file, recover the genes of a certain column in a certain sheet and
        insert them in the ncbiunfiltered collection
        :param list_of_genes: list with id's of genes
        :return: A dictionary with the id:sequence from the list_of_genes target id's
        """

    @abstractmethod
    def download_sequences_from_list_to_dict(self, list_of_genes: list):
        """
        From a list of genes (strings representing id's from NCBI),this method
        retrieves their sequences and creates a dictionary.
        :param list_of_genes: list with id's of genes
        :return: A dictionary with the id:sequence from the list_of_genes target id's
        """

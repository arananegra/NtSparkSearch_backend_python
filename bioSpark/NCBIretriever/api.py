from abc import abstractmethod, ABCMeta


class INCBIretriever(metaclass=ABCMeta):
    @abstractmethod
    def insert_in_collection_from_excel(self, sheet: int, column_name: str) -> None:
        """
        From an excel file, recover the genes of a certain column in a certain sheet and
        insert them in the ncbiunfiltered collection
        :param sheet: sheet from which extract the genes
        :param column_name: name of the column with the genes
        """

    @abstractmethod
    def obtain_list_of_ids_from_mongo(self) -> list:
        """
        From a the unfiltered collection of mongo, extract the ids and returns a list
        with them
        :return: A list with the ids of the unfiltered genes
        """

    @abstractmethod
    def update_genes_from_dict(self, dict_of_genes: dict) -> None:
        """
        From a dict of genes (strings representing id's from NCBI and sequences), this method
        updates the mongo collection with the new information coming from the dict.
        :param dict_of_genes: dictionary with id's and sequences to update in the collection
        unfiltered
        """

    @abstractmethod
    def download_sequences_from_list_as_dict(self, list_of_genes: list) -> dict:
        """
        From a list of genes (strings representing id's from NCBI),this method
        retrieves their sequences and creates a dictionary (ids:sequences)
        :param list_of_genes: list with id's of genes
        :return: A dictionary with all the gene ids and their sequences
        """
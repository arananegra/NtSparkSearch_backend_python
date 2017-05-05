from abc import abstractmethod, ABCMeta


class IGeneRetriever(metaclass=ABCMeta):

    @abstractmethod
    def insert_in_collection_from_excel(self, file_path: str, sheet: str, column_name: str) -> None:
        """
        From an excel file, recover the genes of a certain column in a certain sheet and
        insert them in the ncbiunfiltered collection
        :param file_path: path of the excel file
        :param sheet: sheet from which extract the genes
        :param column_name: name of the column with the genes
        """

    @abstractmethod
    def insert_in_collection_from_fasta(self, file_path: str) -> None:
        """
        From an fasta file, recover the ids and sequences and
        insert them in the ncbiunfiltered collection
        :param file_path: path of the fasta file
        """
    @abstractmethod
    def insert_in_collection_from_list_of_ids(self, list_of_gene_ids: list) -> None:
        """
        From a list of gene ids, insert them into the unfiltered collection only if
        that record did not previously exists
        :param list_of_gene_ids: list of gene ids
        """

    @abstractmethod
    def export_unfiltered_genes_collection_to_file_with_just_ids(self, file_name: str) -> None:
         """
        Creates a text file with the whole unfiltered collection of genes (just ids)
        :param file_name: name of the text file that will be created
        """

    @abstractmethod
    def export_unfiltered_genes_collection_to_fasta(self, fasta_name: str) -> None:
        """
        Creates a fasta file with the whole unfiltered collection of genes
        :param fasta_name: name of the fasta file that will be created
        """

    @abstractmethod
    def get_list_of_ids_from_mongo(self) -> list:
        """
        From a the unfiltered collection of mongo, extract the ids and returns a list
        with them
        :return: A list with the ids of the unfiltered genes
        """

    @abstractmethod
    def get_list_of_ids_from_mongo_without_sequence(self) -> list:
        """
        From a the unfiltered collection of mongo, extract the ids and returns a list
        with them (only if their sequence is None)
        :return: A list with the ids of the unfiltered genes whose sequence is None
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
    def delete_unfiltered_collection(self) -> None:
        """
        Removes the whole unfiltered collection of genes
        """

    @abstractmethod
    def download_sequences_from_list_as_dict_from_NCBI(self, list_of_genes: list) -> dict:
        """
        From a list of genes (strings representing id's from NCBI),this method
        retrieves their sequences and creates a dictionary (ids:sequences)
        :param list_of_genes: list with id's of genes
        :return: A dictionary with all the gene ids and their sequences
        """
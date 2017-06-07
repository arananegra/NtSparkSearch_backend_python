from abc import abstractmethod, ABCMeta


class IGeneHandler(metaclass=ABCMeta):
    @abstractmethod
    def get_dict_from_filtered_with_sequences(self) -> dict:
        """
        From the filtered collection, obtain all data in dict format (only if
        the sequence field is not empty).
        :return dictionary with id:sequence
        """

    @abstractmethod
    def get_dict_from_unfiltered_with_sequences(self) -> dict:
        """
        From the unfiltered collection, obtain all data in dict format (only if
        the sequence field is not empty).
        :return dictionary with id:sequence
        """

    @abstractmethod
    def get_list_of_ids_from_mongo_unfiltered(self) -> list:
        """
        From a the unfiltered collection of mongo, extract the ids and returns a list
        with them
        :return: A list with the ids of the unfiltered genes
        """

    @abstractmethod
    def get_list_of_ids_from_mongo_filtered(self) -> list:
        """
        From a the filtered collection of mongo, extract the ids and returns a list
        with them
        :return: A list with the ids of the unfiltered genes
        """

    @abstractmethod
    def get_list_of_ids_from_mongo_unfiltered_without_sequence(self) -> list:
        """
        From a the unfiltered collection of mongo, extract the ids and returns a list
        with them (only if their sequence is None)
        :return: A list with the ids of the unfiltered genes whose sequence is None
        """

    @abstractmethod
    def get_list_of_ids_from_mongo_unfiltered_with_sequence(self) -> list:
        """
        From a the unfiltered collection of mongo, extract the ids and returns a list
        with them (only if their sequence is NOT None)
        :return: A list with the ids of the unfiltered genes whose sequence is NOT None
        """

    @abstractmethod
    def insert_in_unfiltered_collection_from_excel(self, file_path: str, sheet: str, column_name: str) -> None:
        """
        From an excel file, recover the genes of a certain column in a certain sheet and
        insert them in the ncbiunfiltered collection
        :param file_path: path of the excel file
        :param sheet: sheet from which extract the genes
        :param column_name: name of the column with the genes
        """

    @abstractmethod
    def insert_in_unfiltered_collection_from_fasta(self, file_path: str) -> None:
        """
        From an fasta file, recover the ids and sequences and
        insert them in the ncbiunfiltered collection
        :param file_path: path of the fasta file
        """

    @abstractmethod
    def insert_in_unfiltered_collection_from_list_of_ids(self, list_of_gene_ids: list) -> None:
        """
        From a list of gene ids, insert them into the unfiltered collection only if
        that record did not previously exists
        :param list_of_gene_ids: list of gene ids
        """

    @abstractmethod
    def insert_filtered_dict_in_filtered_collection(self, filtered_dict: dict) -> None:
        """
        Insert a filtered dictionary in the filtered collection
        :param filtered_dict: dictionary with id:sequence elements which has passed out the prune
        """

    @abstractmethod
    def delete_gene_document_by_gene_id(self, gene_id: str) -> None:
        """
        From gene ids, finds and removes that document in the database
        :param gene_id: gene_id to remove
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
    def export_filtered_genes_collection_to_file_with_just_ids(self, file_name: str) -> None:
        """
        Creates a text file with the whole filtered collection of genes (just ids)
        :param file_name: name of the text file that will be created
        """

    @abstractmethod
    def export_filtered_genes_collection_to_fasta(self, fasta_name: str) -> None:
        """
        Creates a fasta file with the whole filtered collection of genes
        :param fasta_name: name of the fasta file that will be created
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
    def delete_filtered_collection(self) -> None:
        """
        Removes the whole unfiltered collection of genes
        """

    @abstractmethod
    def check_gene_id_list_existance_on_ncbi_from_list(self, list_of_genes: list) -> list:
        """
        From a list of genes (strings representing id's from NCBI),this method
        checks if the gene exists on the NCBI
        :param list_of_genes: list with id's of genes
        :return: A list with genes which exists from the original param list
        """

    @abstractmethod
    def download_sequences_from_list_as_dict_from_NCBI(self, list_of_genes: list) -> dict:
        """
        From a list of genes (strings representing id's from NCBI),this method
        retrieves their sequences and creates a dictionary (ids:sequences)
        :param list_of_genes: list with id's of genes
        :return: A dictionary with all the gene ids and their sequences
        """
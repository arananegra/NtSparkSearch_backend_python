from abc import abstractmethod, ABCMeta


class ISubSequenceSparkMatcher(metaclass=ABCMeta):
    @abstractmethod
    def get_dict_from_unfiltered_with_sequences(self) -> dict:
        """
        From the unfiltered collection, obtain all data in dict format (only if
        the sequence field is not empty).
        :return dictionary with id:sequence
        """

    @abstractmethod
    def filter_sequences_by_sequence_string_to_dict(self, sequence_to_filter: str,
                                                    remove_previous_result: str) -> dict:
        """
        From the get_dict_from_unfiltered_with_sequences method, this method prune the
        dictionary following the parameter pattern
        :param sequence_to_filter: nucleotide sequence used to prune the dictionary
        :param remove_previous_result: value to decide if the filtered collection is removed or not
        :return dictionary with id:sequence after the prune
        """

    @abstractmethod
    def insert_filtered_dict_in_filtered_collection(self, filtered_dict: dict) -> None:
        """
        Insert a filtered dictionary in the filtered collection
        :param filtered_dict: dictionary with id:sequence elements which has passed out the prune
        """

    @abstractmethod
    def get_list_of_ids_from_mongo_filtered(self) -> list:
        """
        Obtain the whole filtered collection (just the ids)
        :return list with the ids whose sequence is an exact match with the pattern provided.
        """

    @abstractmethod
    def delete_filtered_collection(self) -> None:
        """
        Removes the entire filtered collection
        """

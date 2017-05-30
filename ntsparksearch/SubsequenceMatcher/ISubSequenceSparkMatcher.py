from abc import abstractmethod, ABCMeta


class ISubSequenceSparkMatcher(metaclass=ABCMeta):

    @abstractmethod
    def filter_sequences_by_sequence_string_to_dict(self, sequence_to_filter: str) -> dict:
        """
        From the get_dict_from_unfiltered_with_sequences method, this method prune the
        dictionary following the parameter pattern
        :param sequence_to_filter: nucleotide sequence used to prune the dictionary
        :param remove_previous_result: value to decide if the filtered collection is removed or not
        :return dictionary with id:sequence after the prune
        """
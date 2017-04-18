from abc import abstractmethod, ABCMeta
from pyspark.sql import SparkSession


class ISubSequenceSparkMatcher(metaclass=ABCMeta):
    @abstractmethod
    def get_dict_from_unfiltered_with_sequences(self) -> dict:
        """
        From the unfiltered collection, obtain all data in dict format (only if
        the sequence field is not empty).
        :return dictionary with id:sequence
        """
    @abstractmethod
    def filter_sequences_by_sequence_string_to_dict(self, sequence_to_filter: str, spark_session: SparkSession) -> dict:
        """
        From the get_dict_from_unfiltered_with_sequences method, this method prune the
        dictionary following the parameter pattern
        :param sequence_to_filter: nucleotide sequence used to prune the dictionary
        :param spark_session: spark_session in which execute the method
        :return dictionary with id:sequence after the prune
        """
    @abstractmethod
    def insert_filtered_dict_in_filtered_collection(self, filtered_dict: dict) -> None:
        """
        Insert a filtered dictionary in the filtered collection
        :param filtered_dict: dictionary with id:sequence elements which has passed out the prune
        """
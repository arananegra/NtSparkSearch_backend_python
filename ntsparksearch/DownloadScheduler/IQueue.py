from abc import abstractmethod, ABCMeta


class IQueue(metaclass=ABCMeta):

    @abstractmethod
    def insert_queue_document_from_object(self, list_of_genes: list) -> None:
        """
        From a list of gene ids, insert them in the queue collection
        :param list_of_genes: path of the excel file
        """

    @abstractmethod
    def delete_collection(self):
        """
        Removes the whole queue collection
        """
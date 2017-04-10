from abc import abstractmethod, ABCMeta
from bioSpark.common.domain import NucleotidesFromNCBI

class INCBIretriever(metaclass=ABCMeta):

    @abstractmethod
    def download_sequences_to_dict(self):
       ""


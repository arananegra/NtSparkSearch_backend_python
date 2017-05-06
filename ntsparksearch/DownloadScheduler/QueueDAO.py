from ntsparksearch.DownloadScheduler.QueueDTO import QueueDTO
from pymongo import MongoClient
from pymongo import TEXT as index_text
from pymongo.collection import Collection
from ntsparksearch.Common.Constants import Constants
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet


class QueueDAO(object):
    def __init__(self, client_reference: MongoClient, database_name: str, collection_name: str):
        self._client_reference = client_reference
        self._database_name = database_name
        self._collection_name = collection_name

    def get_collection(self) -> Collection:
        try:
            collection = self._client_reference[
                self._database_name][self._collection_name]

            return collection
        except Exception as error:
            print('Caught exception obtaining the collection: ' + repr(error))

# insertar desde objeto luego llamada a actualizar el otro

    def insert_queue_document_from_object(self, queue_object: QueueDTO) -> None:
        collection_from_client_reference = None
        try:

            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.insert_one(queue_object.__dict__)

        except Exception as error:
            print('Caught exception at insert from object operation at queue collection: ' + repr(error))

    def delete_collection(self) -> None:

        try:
            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.delete_many({})

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))
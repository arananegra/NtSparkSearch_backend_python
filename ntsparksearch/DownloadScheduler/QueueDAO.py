from ntsparksearch.DownloadScheduler.QueueDTO import QueueDTO
from pymongo import MongoClient
from pymongo.collection import Collection
from datetime import datetime
from ntsparksearch.Common.Constants import Constants


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

    def insert_queue_document_from_object(self, queue_object: QueueDTO) -> None:
        collection_from_client_reference = None
        try:

            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.insert_one(queue_object.__dict__)

        except Exception as error:
            print('Caught exception at insert from object operation at queue collection: ' + repr(error))

    def get_most_recent_document_by_date(self) -> QueueDTO:
        collection_from_client_reference = None
        try:

            collection_from_client_reference = self.get_collection()
            collection_most_recent_document = collection_from_client_reference.find() \
                .sort(Constants.QUEUE_ENTRY_DATE,
                      Constants.MONGODB_COLLECTION_ASCENDING_ORDER) \
                .limit(1)

            most_recent_queue = QueueDTO()
            for document in collection_most_recent_document:
                list_of_genes = document[Constants.QUEUE_LIST_OF_GENES]
                entry_date = document[Constants.QUEUE_ENTRY_DATE]

                most_recent_queue.gene_id_list_to_download = list_of_genes
                most_recent_queue.entry_date_in_collection = entry_date

            return most_recent_queue
        except Exception as error:
            print('Caught exception at insert from object operation at queue collection: ' + repr(error))

    def delete_collection(self) -> None:

        try:
            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.delete_many({})

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))


# test = QueueDAO(MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
#                 Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_QUEUE)
#
# q = QueueDTO()
# l = ["3", "3", "3"]
# currentD = datetime.now()
#
# q.gene_id_list_to_download = l
# q.entry_date_in_collection = currentD
#
# queue = test.get_most_recent_document_by_date()
# print(queue)
# # test.insert_queue_document_from_object(q)

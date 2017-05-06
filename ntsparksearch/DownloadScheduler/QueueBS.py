from ntsparksearch.DownloadScheduler.IQueue import IQueue
from ntsparksearch.DownloadScheduler.QueueDTO import QueueDTO
from ntsparksearch.DownloadScheduler.QueueDAO import QueueDAO
from pymongo import MongoClient
from ntsparksearch.Common.Constants import Constants


class QueueBS(IQueue):
    def insert_queue_document_from_object(self, list_of_genes: list) -> None:
        mongo_queue_dao = QueueDAO(
            MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
            Constants.MONGODB_DB_NAME,
            Constants.MONGODB_COLLECTION_QUEUE)

        queueDTO = QueueDTO()

        try:
            if list_of_genes is not None:
                queueDTO.gene_id_list_to_download = list_of_genes
                mongo_queue_dao.insert_queue_document_from_object(queueDTO)

        except Exception as error:
            print('Caught exception when inserting at queue the list: ' + repr(error))

    def delete_collection(self):
        mongo_queue_dao = QueueDAO(
            MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
            Constants.MONGODB_DB_NAME,
            Constants.MONGODB_COLLECTION_QUEUE)

        try:
            mongo_queue_dao.delete_collection()

        except Exception as error:
            print('Caught exception when deleting queue collection: ' + repr(error))
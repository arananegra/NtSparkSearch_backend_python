from ntsparksearch.Common.GeneDTO import GeneDTO, GeneSearcher
from pymongo import MongoClient
from pymongo import TEXT as indexText
from pymongo.collection import Collection
from ntsparksearch.Common.Constants import Constants


class GeneDAO(object):
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

    @staticmethod
    def create_text_index_in_collection(collection_from_client_reference: Collection, attribute_of_mongo: str,
                                        name_of_index: str) -> None:
        try:
            collection_from_client_reference.create_index([(attribute_of_mongo, indexText)], name=name_of_index,
                                                          default_language='english')
        except Exception as error:
            print('Caught exception creating the index at collection: ' + repr(error))

    def search_gene_objects_and_return_as_list(self, search_criteria: GeneSearcher) -> list:
        collection_from_client_reference = None
        collection_cursor = None

        single_gene_record = None
        list_gene_records = None

        mongodb_find = None

        collection_from_client_reference = self.get_collection()

        try:

            list_gene_records = []

            if (search_criteria.search_by_gene_id_criteria is not None):
                mongodb_find = search_criteria.search_by_gene_id_criteria

            self.create_text_index_in_collection(collection_from_client_reference, Constants.gene_id, '_gene_id_text')

            collection_cursor = collection_from_client_reference. \
                find({"$text": {"$search": "\"%s\"" % mongodb_find}})

            for document in collection_cursor:
                single_gene_record = GeneDTO()

                single_gene_record.gene_id = document[Constants.gene_id]
                single_gene_record.sequence = document[Constants.sequence]

                list_gene_records.append(single_gene_record)

            if (len(list_gene_records) == 0):
                list_gene_records = None

            return list_gene_records

        except Exception as error:
            print('Caught exception searching by id: ' + repr(error))

    def get_all_gene_objects_as_list(self) -> list:
        collection_from_client_reference = None
        list_gene_records = None
        single_gene_record = None

        collection_from_client_reference = self.get_collection()

        try:
            collection_cursor = collection_from_client_reference. \
                find({})

            list_gene_records = []

            for document in collection_cursor:
                single_gene_record = GeneDTO()

                single_gene_record.gene_id = document[Constants.gene_id]
                single_gene_record.sequence = document[Constants.sequence]

                list_gene_records.append(single_gene_record)

            if len(list_gene_records) == 0:
                list_gene_records = None

            return list_gene_records

        except Exception as error:
            print('Caught exception getting all elements as list: ' + repr(error))

    def get_all_gene_objects_as_dict(self) -> dict:
        collection_from_client_reference = None
        dict_gene_records = None
        single_gene_record = None

        collection_from_client_reference = self.get_collection()

        try:
            collection_cursor = collection_from_client_reference. \
                find({})

            dict_gene_records = {}

            for document in collection_cursor:
                single_gene_record = GeneDTO()

                single_gene_record.gene_id = document[Constants.gene_id]
                single_gene_record.sequence = document[Constants.sequence]

                dict_gene_records[single_gene_record.gene_id] = single_gene_record.sequence

            if len(dict_gene_records) == 0:
                dict_gene_records = None

            return dict_gene_records

        except Exception as error:
            print('Caught exception getting all elements as dict: ' + repr(error))

    def delete_gene_document_by_gene_id(self, gene_id: str) -> None:
        collection_from_client_reference = None
        try:
            collection_from_client_reference = self.get_collection()
            criteria = GeneSearcher()
            criteria.search_by_gene_id_criteria = gene_id

            # !!!!!!!
            # llevar esta condición a la bs --> es esa capa la que se tiene que encargar
            # de hacer esas comprobaciones
            #if self.search_gene_objects_and_return_as_list(criteria) is None:
            collection_from_client_reference.delete_one({Constants.gene_id: gene_id})

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))

    def delete_collection(self) -> None:

        try:
            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.delete_many({})

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))

    def update_gene_element_from_object(self, gene_record: GeneDTO, upsert: bool) -> None:
        collection_from_client_reference = None
        try:
            collection_from_client_reference = self.get_collection()

            self.create_text_index_in_collection(collection_from_client_reference, Constants.gene_id, Constants.gene_id_index)

            collection_from_client_reference.update_one({"$text": {"$search": "\"%s\"" % gene_record.gene_id}}, {
                '$set': {
                    Constants.gene_id: gene_record.gene_id,
                    Constants.sequence: gene_record.sequence
                }
            }, upsert=upsert)

        except Exception as error:
            print('Caught exception at update operation: ' + repr(error))

    def insert_gene_document_from_object(self, gene_object: GeneDTO) -> None:
        collection_from_client_reference = None
        try:

            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.insert_one(gene_object.__dict__)

        except Exception as error:
            print('Caught exception at insert from object operation: ' + repr(error))

    def insert_gene_document_from_list_of_objects(self, list_of_gene_objects: list) -> None:
        try:

            for gene_object in list_of_gene_objects:
                self.insert_gene_document_from_object(gene_object)

        except Exception as error:
            print('Caught exception at insert from list operation: ' + repr(error))

    def insert_gene_document_from_non_object_dict(self, dict_with_gene_raw_data: dict) -> None:

        try:

            for k, v in dict_with_gene_raw_data.items():
                gene_object = GeneDTO()

                gene_object.gene_id = k
                gene_object.sequence = v

                self.insert_gene_document_from_object(gene_object)

        except Exception as error:
            print('Caught exception at insert from dict operation: ' + repr(error))
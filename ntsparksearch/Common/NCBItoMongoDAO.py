from ntsparksearch.Common.NucleotidesFromNCBIDTO import NucleotidesFromNCBIDTO, NCBIsearcher
from pymongo import MongoClient
from pymongo import TEXT as indexText
from pymongo.collection import Collection


class NCBItoMongoDAO(object):
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

    def search_ncbi_objects_and_return_as_list(self, search_criteria: NCBIsearcher) -> list:
        collection_from_client_reference = None
        collection_cursor = None

        single_ncbi_record = None
        list_ncbi_records = None

        mongodb_find = None

        collection_from_client_reference = self.get_collection()

        try:

            list_ncbi_records = []

            if (search_criteria.search_by_NCBI_id_criteria is not None):
                mongodb_find = search_criteria.search_by_NCBI_id_criteria

            self.create_text_index_in_collection(collection_from_client_reference, '_idNcbi', '_idNcbi_text')

            collection_cursor = collection_from_client_reference. \
                find({"$text": {"$search": "\"%s\"" % mongodb_find}})

            for document in collection_cursor:
                single_ncbi_record = NucleotidesFromNCBIDTO()

                single_ncbi_record.idNcbi = document["_idNcbi"]
                single_ncbi_record.sequence = document["_sequence"]

                list_ncbi_records.append(single_ncbi_record)

            if (len(list_ncbi_records) == 0):
                list_ncbi_records = None

            return list_ncbi_records

        except Exception as error:
            print('Caught exception searching by id: ' + repr(error))

    def get_all_ncbi_objects_as_list(self) -> list:
        collection_from_client_reference = None
        list_ncbi_records = None
        single_ncbi_record = None

        collection_from_client_reference = self.get_collection()

        try:
            collection_cursor = collection_from_client_reference. \
                find({})

            list_ncbi_records = []

            for document in collection_cursor:
                single_ncbi_record = NucleotidesFromNCBIDTO()

                single_ncbi_record.idNcbi = document["_idNcbi"]
                single_ncbi_record.sequence = document["_sequence"]

                list_ncbi_records.append(single_ncbi_record)

            if len(list_ncbi_records) == 0:
                list_ncbi_records = None

            return list_ncbi_records

        except Exception as error:
            print('Caught exception getting all elements as list: ' + repr(error))

    def get_all_ncbi_objects_as_dict(self) -> dict:
        collection_from_client_reference = None
        dict_ncbi_records = None
        single_ncbi_record = None

        collection_from_client_reference = self.get_collection()

        try:
            collection_cursor = collection_from_client_reference. \
                find({})

            dict_ncbi_records = {}

            for document in collection_cursor:
                single_ncbi_record = NucleotidesFromNCBIDTO()

                single_ncbi_record.idNcbi = document["_idNcbi"]
                single_ncbi_record.sequence = document["_sequence"]

                dict_ncbi_records[single_ncbi_record.idNcbi] = single_ncbi_record.sequence

            if len(dict_ncbi_records) == 0:
                dict_ncbi_records = None

            return dict_ncbi_records

        except Exception as error:
            print('Caught exception getting all elements as dict: ' + repr(error))

    def delete_ncbi_document_by_idNcbi(self, idNcbi: str) -> None:
        collection_from_client_reference = None
        try:
            collection_from_client_reference = self.get_collection()
            criteria = NCBIsearcher()
            criteria.searchByNCBIidCriteria = idNcbi

            # !!!!!!!
            # llevar esta condición a la bs --> es esa capa la que se tiene que encargar
            # de hacer esas comprobaciones
            #if self.search_ncbi_objects_and_return_as_list(criteria) is None:
            collection_from_client_reference.delete_one({'_idNcbi': idNcbi})

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))

    def delete_collection(self) -> None:

        try:
            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.delete_many({})

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))

    def update_ncbi_element_from_object(self, ncbi_record: NucleotidesFromNCBIDTO, upsert: bool) -> None:
        collection_from_client_reference = None
        try:
            collection_from_client_reference = self.get_collection()

            self.create_text_index_in_collection(collection_from_client_reference, '_idNcbi', '_idNcbi_text')

            collection_from_client_reference.update_one({"$text": {"$search": "\"%s\"" % ncbi_record.idNcbi}}, {
                '$set': {
                    '_idNcbi': ncbi_record.idNcbi,
                    '_sequence': ncbi_record.sequence
                }
            }, upsert=upsert)

        except Exception as error:
            print('Caught exception at update operation: ' + repr(error))

    def insert_ncbi_document_from_object(self, ncbi_object: NucleotidesFromNCBIDTO) -> None:
        collection_from_client_reference = None
        try:

            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.insert_one(ncbi_object.__dict__)

        except Exception as error:
            print('Caught exception at insert from object operation: ' + repr(error))

    def insert_ncbi_document_from_list_of_objects(self, list_of_ncbi_objects: list) -> None:
        try:

            for ncbi_object in list_of_ncbi_objects:
                self.insert_ncbi_document_from_object(ncbi_object)

        except Exception as error:
            print('Caught exception at insert from list operation: ' + repr(error))

    def insert_ncbi_document_from_non_object_dict(self, dict_with_ncbi_raw_data: dict) -> None:

        try:

            for k, v in dict_with_ncbi_raw_data.items():
                ncbi_object = NucleotidesFromNCBIDTO()

                ncbi_object.idNcbi = k
                ncbi_object.sequence = v

                self.insert_ncbi_document_from_object(ncbi_object)

        except Exception as error:
            print('Caught exception at insert from dict operation: ' + repr(error))
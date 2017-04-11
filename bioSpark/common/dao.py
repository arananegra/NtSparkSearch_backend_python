from bioSpark.common.domain import NucleotidesFromNCBI, NCBIsearcher
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

    def create_text_index_in_collection(self, collection_from_client_reference: Collection, attribute_of_mongo: str,
                                        name_of_index: str) -> None:
        collection_from_client_reference.create_index([(attribute_of_mongo, indexText)], name=name_of_index,
                                                      default_language='english')

    def search_ncbi_objects_and_return_as_list(self, search_criteria: NCBIsearcher) -> list:
        collection_from_client_reference = None
        collection_cursor = None

        single_ncbi_record = None
        list_ncbi_records = None

        mongodbFind = None

        collection_from_client_reference = self.get_collection()

        try:

            list_ncbi_records = []

            if (search_criteria.search_by_NCBI_id_criteria is not None):
                mongodbFind = search_criteria.search_by_NCBI_id_criteria

            self.create_text_index_in_collection(collection_from_client_reference, '_idNcbi', '_idNcbi_text')

            collection_cursor = collection_from_client_reference. \
                find({"$text": {"$search": "\"%s\"" % mongodbFind}})

            for document in collection_cursor:
                single_ncbi_record = NucleotidesFromNCBI()

                single_ncbi_record.idNcbi = document["_idNcbi"]
                single_ncbi_record.sequence = document["_sequence"]

                list_ncbi_records.append(single_ncbi_record)

            if (len(list_ncbi_records) == 0):
                list_ncbi_records = None

            return list_ncbi_records

        except Exception as error:
            print('Caught exception when searching by id: ' + repr(error))

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
                single_ncbi_record = NucleotidesFromNCBI()

                single_ncbi_record.idNcbi = document["_idNcbi"]
                single_ncbi_record.sequence = document["_sequence"]

                list_ncbi_records.append(single_ncbi_record)

            if len(list_ncbi_records) == 0:
                list_ncbi_records = None

            return list_ncbi_records

        except Exception as error:
            print('Caught exception when getting all elements as list: ' + repr(error))

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
                single_ncbi_record = NucleotidesFromNCBI()

                single_ncbi_record.idNcbi = document["_idNcbi"]
                single_ncbi_record.sequence = document["_sequence"]

                dict_ncbi_records[single_ncbi_record.idNcbi] = single_ncbi_record.sequence

            if len(dict_ncbi_records) == 0:
                dict_ncbi_records = None

            return dict_ncbi_records

        except Exception as error:
            print('Caught exception when getting all elements as dict: ' + repr(error))

    def delete_ncbi_document_by_idNcbi(self, idNcbi: str) -> None:
        collection_from_client_reference = None
        try:
            collection_from_client_reference = self.get_collection()
            criteria = NCBIsearcher()
            criteria.searchByNCBIidCriteria = idNcbi

            # !!!!!!!
            # llevar esta condición a la bs --> es esa capa la que se tiene que encargar
            # de hacer esas comprobaciones
            if self.search_ncbi_objects_and_return_as_list(criteria) is None:
                collection_from_client_reference.delete_one({'_idNcbi': idNcbi})

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))

    def update_ncbi_element_from_object(self, ncbi_record: NucleotidesFromNCBI, upsert: bool) -> None:
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
            print('Caught exception at delete operation: ' + repr(error))

    def insert_ncbi_document_from_object(self, ncbi_object: NucleotidesFromNCBI) -> None:
        collection_from_client_reference = None
        try:

            collection_from_client_reference = self.get_collection()
            collection_from_client_reference.insert_one(ncbi_object.__dict__)

        except Exception as error:
            print('Caught exception at insert operation: ' + repr(error))

    def insert_ncbi_document_from_list_of_objects(self, list_of_ncbi_objects: list) -> None:
        try:

            for ncbi_object in list_of_ncbi_objects:
                self.insert_ncbi_document_from_object(ncbi_object)

        except Exception as error:
            print('Caught exception at insert operation: ' + repr(error))

    def insert_ncbi_document_from_non_object_dict(self, dict_with_ncbi_raw_data: dict) -> None:

        try:

            for k, v in dict_with_ncbi_raw_data.items():
                ncbi_object = NucleotidesFromNCBI()

                ncbi_object.idNcbi = k
                ncbi_object.sequence = v

                self.insert_ncbi_document_from_object(ncbi_object)

        except Exception as error:
            print('Caught exception at insert operation: ' + repr(error))

# client = MongoClient('localhost', 27017)
#
# testDao = NCBItoMongoDAO(client, "SparkTest", "test01")
#
# criteria = NCBIsearcher()
# criteria.search_by_NCBI_id_criteria = "002"
#
# # Example search_ncbi_objects_and_return_as_list
#
# listResult = testDao.search_ncbi_objects_and_return_as_list(criteria)
# print(listResult[0].idNcbi)
#
# # Example get_all_ncbi_objects_as_dict
#
# dictResult = testDao.get_all_ncbi_objects_as_dict()
# print(dictResult)
#
# # Example delete_ncbi_document_by_idNcbi
#
# # testDao.delete_ncbi_document_by_idNcbi('003')
#
# # Example insert_ncbi_document
#
# testNcbi1 = NucleotidesFromNCBI()
#
# testNcbi1.idNcbi = "004"
# testNcbi1.sequence = "GCTACG"
#
# testNcbi2 = NucleotidesFromNCBI()
#
# testNcbi2.idNcbi = "005"
# testNcbi2.sequence = "GCTACG"

# testDao.insert_ncbi_document_from_object(testNcbi1)

# Example of insert_ncbi_document_from_list_of_objects

# list_of_ncbi_objects = [testNcbi1, testNcbi2]
# testDao.insert_ncbi_document_from_list_of_objects(list_of_ncbi_objects)

# Example of insert_ncbi_document_from_non_object_dict


# dict_test = {"006": "GCTAGCA", "007": "ATGCTAG"}
#
# # testDao.insert_ncbi_document_from_non_object_dict(dict_test)
#
# testNcbi3 = NucleotidesFromNCBI()
# testNcbi3.idNcbi = "013"
# testNcbi3.sequence = None
#
# testDao.update_ncbi_element_from_object(testNcbi3, upsert=True)

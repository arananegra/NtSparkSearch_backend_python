from bioSpark.common.domain import NucleotidesFromNCBI, NCBIsearcher
from pymongo import MongoClient
from bson.objectid import ObjectId


class NCBIDao(object):
    def __init__(self, client_reference: MongoClient, database_name: str, collection_name: str):
        self._client_reference = client_reference
        self._database_name = database_name
        self._collection_name = collection_name

    def get_collection(self):
        return self._client_reference[
            self._database_name][self._collection_name]

    def search_ncbi_objects_and_return_as_list(self, searchCriteria: NCBIsearcher):
        collection_from_client_reference = None
        collection_cursor = None

        single_ncbi_record = None
        list_ncbi_records = None

        mongodbFind = None
        list_ncbi_records = []

        if (searchCriteria.searchByNCBIidCriteria != None):
            mongodbFind = searchCriteria.searchByNCBIidCriteria

        collection_from_client_reference = self.get_collection()

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

    def get_all_ncbi_objects_as_dict(self):
        collection_from_client_reference = None
        dict_ncbi_records = None
        single_ncbi_record = None

        collection_from_client_reference = self.get_collection()
        collection_cursor = collection_from_client_reference. \
            find({})

        dict_ncbi_records = {}

        for document in collection_cursor:
            single_ncbi_record = NucleotidesFromNCBI()

            single_ncbi_record.idNcbi = document["_idNcbi"]
            single_ncbi_record.sequence = document["_sequence"]

            dict_ncbi_records[single_ncbi_record.idNcbi] = single_ncbi_record.sequence

        if (len(dict_ncbi_records) == 0):
            dict_ncbi_records = None

        return dict_ncbi_records

    def delete_ncbi_document_by_idNcbi(self, idNcib: str):
        collection_from_client_reference = None
        try:
            collection_from_client_reference = self.get_collection()
            criteria = NCBIsearcher()
            criteria.searchByNCBIidCriteria = idNcib
            if self.search_ncbi_objects_and_return_as_list(criteria) != None:
                collection_from_client_reference.delete_one({'_idNcbi': idNcib })

            # cómo tratar el caso de que se llame a borrar un elemento que no existe
        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))

    def insert_ncbi_document(self, ncbi_object: NucleotidesFromNCBI):
        collection_from_client_reference = None
        try:

            collection_from_client_reference = self.get_collection()
            criteria = NCBIsearcher()
            criteria.searchByNCBIidCriteria = ncbi_object.idNcbi

            if self.search_ncbi_objects_and_return_as_list(criteria) == None:
                collection_from_client_reference.insert(ncbi_object.__dict__)
                print("insertado")
        except Exception as error:
            print('Caught exception at insert operation: ' + repr(error))


client = MongoClient('localhost', 27017)

testDao = NCBIDao(client, "SparkTest", "test01")

criteria = NCBIsearcher()
criteria.searchByNCBIidCriteria = "002"

# Example search_ncbi_objects_and_return_as_list

listResult = testDao.search_ncbi_objects_and_return_as_list(criteria)
print(listResult[0].idNcbi)

# Example get_all_ncbi_objects_as_dict

dictResult = testDao.get_all_ncbi_objects_as_dict()
print(dictResult)


# Example delete_ncbi_document_by_idNcbi

#testDao.delete_ncbi_document_by_idNcbi('003')

# Example insert_ncbi_document

testNcbi = NucleotidesFromNCBI()

testNcbi.idNcbi = "004"
testNcbi.sequence = "GCTACG"

testDao.insert_ncbi_document(testNcbi)
from bioSpark.common.domain import NucleotidesFromNCBI, NCBIsearcher
from pymongo import MongoClient
import json


class NCBIDao(object):
    def __init__(self, connectionReference: MongoClient, collectionName: str):
        self._connectionReference = connectionReference
        self._collectionName = collectionName

    def searchNCBISequence(self, searchCriteria: NCBIsearcher):
        collection_from_connectionInstance = None
        collection_cursor = None

        single_ncbi_record = None
        list_ncbi_records = None

        list_ncbi_records = []

        mongodbFind = ""#"{ _idNcbi: { $eq:"

        if (searchCriteria.searchByNCBIidCriteria != None):
            mongodbFind += searchCriteria.searchByNCBIidCriteria

        #mongodbFind += "} }"

        collection_from_connectionInstance = self._connectionReference[
            self._collectionName].Constants.MONGODB_COLLECTION_UNFILTERED_NCBI

        collection_cursor = collection_from_connectionInstance.find({"_idNcbi": {"$eq": mongodbFind}})

        for document in collection_cursor:
            single_ncbi_record = NucleotidesFromNCBI()

            single_ncbi_record.idNcbi = document["_idNcbi"]
            single_ncbi_record.sequence = document["_sequence"]

            list_ncbi_records.extend(single_ncbi_record)

        if (len(list_ncbi_records) == 0):
            list_ncbi_records = None

        return list_ncbi_records

client = MongoClient('localhost', 27017)
db = client[str("SparkTest")]

testDao = NCBIDao(db, "test01")

criteria = NCBIsearcher()
criteria.searchByNCBIidCriteria = "001"

listResult = testDao.searchNCBISequence(criteria)


from ntsparksearch.SubsequenceMatcher.ISubSequenceSparkMatcher import ISubSequenceSparkMatcher
from ntsparksearch.Common.NCBItoMongoDAO import NCBItoMongoDAO
from Bio import SeqUtils
from ntsparksearch.Common.Constants import Constants
from pymongo import MongoClient
from copy import deepcopy
import findspark
findspark.init(Constants.SPARK_HOME)
from pyspark.sql import SparkSession


class SubSequenceSparkMatcherBS(ISubSequenceSparkMatcher):
    def get_dict_from_unfiltered_with_sequences(self) -> dict:

        try:

            dict_with_genes = {}

            mongo_dao_retriever = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            dict_with_genes = mongo_dao_retriever.get_all_ncbi_objects_as_dict()

            deep_copy_dict_with_genes = deepcopy(dict_with_genes)

            for id_ncbi, sequence in deep_copy_dict_with_genes.items():
                if sequence is None:
                    del dict_with_genes[id_ncbi]

            return dict_with_genes

        except Exception as error:
            print('Caught exception getting unfiltered sequences (to dict):' + repr(error))

    def filter_sequences_by_sequence_string_to_dict(self, sequence_to_filter: str, remove_previous_result: str) -> dict:

        def map_locator_Spark(x, subsequence):
            return len(SeqUtils.nt_search(x[1], subsequence)) > 1

        try:
            spark_session = SparkSession \
            .builder \
            .appName("ntsparksearch") \
            .config("spark.driver.memory", "4g") \
            .config("spark.driver.maxResultSize", "3g") \
            .config("spark.executor.memory", "3g").getOrCreate()

            dict_to_filter = self.get_dict_from_unfiltered_with_sequences()

            if remove_previous_result == "y":
                self.delete_filtered_collection()

            sc = spark_session.sparkContext

            list_of_list_of_genes = [[k, v] for k, v in dict_to_filter.items()]
            list_of_list_of_genes_rdd = sc.parallelize(list_of_list_of_genes)

            list_of_list_of_genes_filtered = list_of_list_of_genes_rdd. \
                filter(lambda element: map_locator_Spark(element, sequence_to_filter)). \
                collect()

            dict_now_filtered = {gene_id: sequence for (gene_id, sequence) in list_of_list_of_genes_filtered}

            return dict_now_filtered

        except Exception as error:
            print('Caught exception trying to filter the collection:' + repr(error))

    def insert_filtered_dict_in_filtered_collection(self, filtered_dict: dict) -> None:

        try:

            mongo_dao_manager = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            mongo_dao_manager.insert_ncbi_document_from_non_object_dict(filtered_dict)

        except Exception as error:
            print('Caught exception while inserting filtered sequences in mongo collection' + repr(error))

    def get_list_of_ids_from_mongo_filtered(self) -> list:
        list_of_just_ids = None

        try:
            list_of_just_ids = []

            mongo_dao_retriever_filtered = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            list_of_ncbi_objets = mongo_dao_retriever_filtered.get_all_ncbi_objects_as_list()

            for ncbi_object in list_of_ncbi_objets:
                ncbi_id = ncbi_object.idNcbi
                list_of_just_ids.append(ncbi_id)

            return list_of_just_ids

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list (filtered): ' + repr(error))

    def delete_filtered_collection(self) -> None:

        try:

            mongo_dao_retriever_filtered = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            mongo_dao_retriever_filtered.delete_collection()

        except Exception as error:
            print('Caught exception when removing filtered collection' + repr(error))

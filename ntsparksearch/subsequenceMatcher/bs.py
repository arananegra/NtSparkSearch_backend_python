from ntsparksearch.subsequenceMatcher.api import ISubSequenceSparkMatcher
from ntsparksearch.common.dao import NCBItoMongoDAO
from Bio import SeqUtils
from ntsparksearch.common.util import Constants
from pymongo import MongoClient
import findspark
findspark.init("/Users/alvarogomez/spark-2.1.0")
from pyspark.sql import SparkSession
# import os,sys
# os.environ['SPARK_HOME'] = "/Users/alvarogomez/spark-2.1.0"
# sys.path.append("/Users/alvarogomez/spark-2.1.0/python")
# sys.path.append("/Users/alvarogomez/spark-2.1.0/python/python/lib/py4j-0.10.4-src.zip")
# try:
#     from pyspark.sql import SparkSession
#     print("imported spark lib")
# except ImportError as e:
#     print ("Error importing Spark Modules", e)


class SubSequenceSparkMatcherBS(ISubSequenceSparkMatcher):
    def get_dict_from_unfiltered_with_sequences(self) -> dict:

        try:

            dict_with_genes = {}

            mongo_dao_retriever = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            dict_with_genes = mongo_dao_retriever.get_all_ncbi_objects_as_dict()

            for id_ncbi, sequence in dict_with_genes.items():
                if sequence is (None or ""):
                    del dict_with_genes[id_ncbi]

            return dict_with_genes

        except Exception as error:
            print('Caught exception getting unfiltered sequences (to dict):' + repr(error))

    def filter_sequences_by_sequence_string_to_dict(self, sequence_to_filter: str) -> dict:

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

    def obtain_list_of_ids_from_mongo_filtered(self) -> list:
        list_of_just_ids = None

        try:
            list_of_just_ids = []

            mongo_dao_retriever = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            list_of_ncbi_objets = mongo_dao_retriever.get_all_ncbi_objects_as_list()

            for ncbi_object in list_of_ncbi_objets:
                ncbi_id = ncbi_object.idNcbi
                list_of_just_ids.append(ncbi_id)

            return list_of_just_ids

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list (filtered): ' + repr(error))

import findspark
from Bio import SeqUtils
from flask_security import current_user

from ntsparksearch.Common.Constants import Constants
from ntsparksearch.GeneHandler.GeneHandlerBS import GeneHandlerBS
from ntsparksearch.SubsequenceMatcher.ISubSequenceSparkMatcher import ISubSequenceSparkMatcher

findspark.init(Constants.SPARK_HOME)
# from pyspark.sql import SparkSession
from pyspark import SparkContext, SparkConf


class SubSequenceSparkMatcherBS(ISubSequenceSparkMatcher):
    def filter_sequences_by_sequence_string_to_dict_from_custom_list_of_genes(self, sequence_to_filter: str,
                                                                              list_of_ids_to_filter: list) -> dict:

        def map_locator_Spark(x, subsequence):
            return len(SeqUtils.nt_search(x[1], subsequence)) > 1

        try:
            gene_handler_BS = GeneHandlerBS(Constants.MONGODB_DB_INITIAL + str(current_user.id))
            conf = SparkConf().setAppName("ntsparksearch").setMaster("local[8]")

            sc = SparkContext(conf=conf)

            dict_to_filter = gene_handler_BS.get_dict_of_genes_object_from_list_of_ids_with_sequences(
                list_of_ids_to_filter)
            gene_handler_BS.delete_filtered_collection()

            list_of_list_of_genes = [[k, v] for k, v in dict_to_filter.items()]
            list_of_list_of_genes_rdd = sc.parallelize(list_of_list_of_genes)

            list_of_list_of_genes_filtered = list_of_list_of_genes_rdd. \
                filter(lambda element: map_locator_Spark(element, sequence_to_filter)). \
                collect()

            dict_now_filtered = {gene_id: sequence for (gene_id, sequence) in list_of_list_of_genes_filtered}
            return dict_now_filtered

        except Exception as error:
            print('Caught exception trying to filter the collection:' + repr(error))

        finally:
            sc.stop()

    def filter_sequences_by_sequence_string_to_dict(self, sequence_to_filter: str) -> dict:

        def map_locator_Spark(x, subsequence):
            return len(SeqUtils.nt_search(x[1], subsequence)) > 1

        try:
            conf = SparkConf().setAppName("ntsparksearch").setMaster("local[8]")
            # spark_session = SparkSession \
            #     .builder \
            #     .appName("ntsparksearch") \
            #     .master("spark://172.16.239.10:7077") \
            #     .config("spark.driver.memory", "2024m") \
            #     .config("spark.driver.maxResultSize", "2024m") \
            #     .config("spark.executor.memory", "2024m").getOrCreate()

            sc = SparkContext(conf=conf)
            # sc.setLogLevel("ERROR")
            # gene_handler_bs = GeneHandlerBS()

            gene_handler_BS = GeneHandlerBS(Constants.MONGODB_DB_INITIAL + str(current_user.id))

            dict_to_filter = gene_handler_BS.get_dict_from_unfiltered_with_sequences()
            gene_handler_BS.delete_filtered_collection()

            list_of_list_of_genes = [[k, v] for k, v in dict_to_filter.items()]
            list_of_list_of_genes_rdd = sc.parallelize(list_of_list_of_genes)

            list_of_list_of_genes_filtered = list_of_list_of_genes_rdd. \
                filter(lambda element: map_locator_Spark(element, sequence_to_filter)). \
                collect()

            dict_now_filtered = {gene_id: sequence for (gene_id, sequence) in list_of_list_of_genes_filtered}
            return dict_now_filtered

        except Exception as error:
            print('Caught exception trying to filter the collection:' + repr(error))

        finally:
            sc.stop()

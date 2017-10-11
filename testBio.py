from Bio import SeqUtils
from pyspark import SparkContext, SparkConf


# print(SeqUtils.nt_search("GATCAGCAG", "CAG"))


def filter_sequences_by_sequence_string_to_dict(sequence_to_filter: str) -> dict:
    def map_locator_Spark(dict_element, subsequence):
        return [dict_element[0], dict_element[1], SeqUtils.nt_search(dict_element[1], subsequence)]

    def list_filter(pre_filter_dict_element):
        return len(pre_filter_dict_element[2]) > 1

    try:
        conf = SparkConf().setAppName("testSpark").setMaster("local[8]")
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

        # gene_handler_BS = GeneHandlerBS(Constants.MONGODB_DB_INITIAL + str(current_user.id))

        dict_to_filter = {"01": "CAGCA", "02": "TTTTTT"}
        # gene_handler_BS.delete_filtered_collection()

        list_of_list_of_genes = [[k, v] for k, v in dict_to_filter.items()]
        list_of_list_of_genes_rdd = sc.parallelize(list_of_list_of_genes)

        list_of_list_of_genes_filtered = list_of_list_of_genes_rdd. \
            map(lambda element: map_locator_Spark(element, sequence_to_filter)). \
            filter(lambda element: list_filter(element)). \
            collect()

        # dict_now_filtered = {gene_id: sequence for (gene_id, sequence) in list_of_list_of_genes_filtered}
        return list_of_list_of_genes_filtered

    except Exception as error:
        print('Caught exception trying to filter the collection:' + repr(error))

    finally:
        sc.stop()


print(filter_sequences_by_sequence_string_to_dict("NNN"))

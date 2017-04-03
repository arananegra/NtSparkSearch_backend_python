# script for testing
import re
import pickle
from pyspark.sql import SparkSession
from Bio import SeqUtils

if __name__ == "__main__":
    """
        Usage: Spark sequence search testing
    """
    spark = SparkSession\
        .builder\
        .appName("PythonTestingSpark")\
        .config("spark.driver.memory","4g").config("spark.driver.maxResultSize", "3g").config("spark.executor.memory", "3g").getOrCreate()

    sc = spark.sparkContext

    def load_obj(name):
        with open(name, 'rb') as f:
            return pickle.load(f)

    def mapLocatorSpark(x, subsequence):
        return len(SeqUtils.nt_search(x[1],subsequence))>1


    def sequenceRemovalSpark(x):
        return [x[0]]

    genesWithSequenceDictHardcore = load_obj("/Users/alvarogomez/test01_tfg.pickle")

    sub = 'ANAT'

    sub2 = 'AGTTCNCAGCATCANCATCANCATCANCATCA'

    sub3 = 'TTNCNNNAA'

    listOfGenes = [[k,v] for k, v in genesWithSequenceDictHardcore.items()]

    listOfGenesRDD = sc.parallelize(listOfGenes)

    listOfGenesFilter = listOfGenesRDD.filter(lambda element: mapLocatorSpark(element, sub3))

    listOfGenesJustIds = listOfGenesFilter.flatMap(sequenceRemovalSpark).collect()

    print(listOfGenesJustIds)

    print(len(listOfGenesJustIds))

    sc.stop()




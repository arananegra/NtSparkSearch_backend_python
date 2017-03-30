# script for testing
from bioSpark.common.domain import NucleotidesFromNCBI as ncbiNucleotides
import bioSpark.common.pymongoSimpleConnection as mongoConnect

test = ncbiNucleotides(idNcbi = "001", sequence="ACGACT")
test2 = ncbiNucleotides(idNcbi = "002", sequence="GTCAGCATG")


sequencesCollection = mongoConnect.getCollectionFromDb(collectionName="test01", mongoHost='localhost',
                                                       connectionPort=27017, mongoDb="SparkTest")
sequencesCollection.insert(test.__dict__)
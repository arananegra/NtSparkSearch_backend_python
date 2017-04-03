# script for testing
from bioSpark.common.domain import NucleotidesFromNCBI as ncbiNucleotides
import bioSpark.common.pymongoSimpleConnection as mongoConnect

test = ncbiNucleotides()
test2 = ncbiNucleotides()

test.idNcbi = "001"
test.sequence = "GTCAGCATG"
test2.idNcbi = "002"
test2.sequence = "ACGACT"

sequencesCollection = mongoConnect.getCollectionFromDb(collectionName="test01", mongoHost='localhost',
                                                       connectionPort=27017, mongoDb="SparkTest")
sequencesCollection.insert(test.__dict__)
sequencesCollection.insert(test2.__dict__)

cursor = sequencesCollection.find({})

dictionaryWithSeqs = {}

for document in cursor:
    dictionaryWithSeqs[document["_idNcbi"]] = document["_sequence"]

print(dictionaryWithSeqs)


hola = "hola"
hola += " y adios"

print(hola)
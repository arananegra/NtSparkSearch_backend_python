# script for testing
from bioSpark.common.domain import IdAndSequenceFromNCBIforMongoId
from bioSpark.common.pymongoConnection import dc
import bioSpark.common.pymongoSimpleConnection as mongoConnect

test = IdAndSequenceFromNCBIforMongoId(id="001", sequence="ACGACT")
test2 = IdAndSequenceFromNCBIforMongoId(id="002", sequence="ACGACT")


# using the 'singleton' pattern with classes --> pymongoConnection

# If the document does not specify an _id field, then MongoDB will add the _id field
# and assign a unique ObjectId for the document before inserting.
# Most drivers create an ObjectId and insert the _id field, but
# the mongod will create and populate the _id if the driver or application does not.

# create connection
dc.connect('localhost', 27017)
dc.connect()


# using just modules --> pymongoSimpleConnection 

# from that connection, create (or use) collection and perform operations

dc.sequencesTest.insert(test.__dict__)
dc.test01.insert(test.__dict__)
dc.sequencesTest.insert(test2.__dict__)

sequencesCollection = mongoConnect.getCollectionFromDb(collectionName="test01", mongoHost='localhost',
                                                       connectionPort=27017, mongoDb="SparkTest")

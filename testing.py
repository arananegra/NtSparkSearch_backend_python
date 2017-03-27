# script for testing
from bioSpark.common.domain import NucleotidesFromNCBIusingIdOfSequenceInMongo

from bioSpark.common.pymongoConnection import dc

# client = MongoClient('localhost', 27017)
#
# db = client.SparkTestDB
# sequencesCollection = db.sequences

test = NucleotidesFromNCBIusingIdOfSequenceInMongo(id="001", sequence="ACGACT")
test2 = NucleotidesFromNCBIusingIdOfSequenceInMongo(id="002", sequence="ACGACT")


# If the document does not specify an _id field, then MongoDB will add the _id field
# and assign a unique ObjectId for the document before inserting.
# Most drivers create an ObjectId and insert the _id field, but
# the mongod will create and populate the _id if the driver or application does not.

# print(test.__dict__)


# create connection
dc.connect('localhost',27017)

# from that connection, create (or use) collection and perform operations

dc.sequencesTest.insert(test.__dict__)
dc.sequencesTest.insert(test2.__dict__)



#sequencesCollection.insert(test.__dict__)

# class MongoCon(object):
#     __db = None
#
#     @classmethod
#     def get_connection(client):
#         if client.__db is None:
#             client.__db = Connection()
#         return client.__db

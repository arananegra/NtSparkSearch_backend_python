from pymongo import MongoClient


def getDbConnection(dbName, host, port):
    client = MongoClient(host, port)
    db = client[str(dbName)]
    return db


def getCollectionFromDb(mongoDb, collectionName, mongoHost, connectionPort):
    db = getDbConnection(host= mongoHost, port = connectionPort, dbName = mongoDb)
    return db[str(collectionName)]

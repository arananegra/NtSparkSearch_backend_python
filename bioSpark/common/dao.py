from abc import ABCMeta, abstractmethod

class genericMongoDao:
    __metaclass__ = ABCMeta

    @abstractmethod
    def obtainCollection(self, mongodbName, mongoCollection, host, port):
        """
        Obtain a mongodb collection using pymongoSimpleConnection.py
        :param mongodbName: Database name
        :param mongoCollection: Collection name
        :param host: host (ip)
        :param port
        :return: A collection with the matching parameters or a new one 
        if it does not exists
        """

    @abstractmethod
    def obtainCursor(self, mongodbName, mongoCollection, host, port):
        """
        Obtain a mongodb cursor using the obtainCollection function
        :param mongodbName: Database name
        :param mongoCollection: Collection name
        :param host: host (ip)
        :param port
        :return: A cursor over the collection (if the collection does not
        exists, this function will create a new collection and then, it
        will return a cursor (over that collection)
        """
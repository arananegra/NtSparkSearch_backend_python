from pymongo import MongoClient

class _DataContext(object):
    _client = None
    _db = None

    def connect(self, host='localhost', port=27017, dbName = 'SparkSequences_db'):
        self._client = MongoClient(host, port)
        self._db = self._client[dbName]

    @property
    def db(self):
        if self._db is None:
            raise Exception('Not connected')
        return self._db

    def __getattr__(self, item):
        return self.db[item]

dc = _DataContext()

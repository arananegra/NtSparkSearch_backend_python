import redis
from rq import Worker, Queue, Connection

from ntsparksearch.Common.Constants import Constants

listen = ['high', 'default', 'low']

redis_url = 'redis://' + Constants.REDIS_SERVER + ':6379'

conn = redis.from_url(redis_url)

if __name__ == '__main__':
    with Connection(conn):
        worker = Worker(map(Queue, listen))
        worker.work()

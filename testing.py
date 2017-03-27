# script for testing

# bs 
from mongoengine import *

connect('SparkTest', host='localhost', port=27017)


class Sequences(Document):
    idNcbi = StringField(required=True)
    sequence = StringField()

#Sequences(idNcbi='01', sequence = 'CTGAGCAT').save()
#Sequences(idNcbi='02', sequence = 'CGTAGCTAG').save()

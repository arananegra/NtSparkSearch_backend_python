from flask import Flask
from ntsparksearch.RestApi.GeneRetrieverService import *
from ntsparksearch.RestApi.SubSequenceMatcherService import *
from celery import current_app
from celery.bin import worker

app = Flask(__name__)

app.register_blueprint(GeneRetrieverService_endpoints, url_prefix='/generetriever')
app.register_blueprint(SubSequenceMatcherService_endpoints, url_prefix='/genefilter')

@app.route("/")
def hello():
    return "Hello World!"

if __name__ == "__main__":
    app.run()

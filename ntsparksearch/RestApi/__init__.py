from flask import Flask
from flask_cors import CORS

from ntsparksearch.Common.Constants import Constants
from ntsparksearch.RestApi.GeneHandlerService import GeneHandlerService_endpoints
from ntsparksearch.RestApi.SubSequenceMatcherService import SubSequenceMatcherService_endpoints
from ntsparksearch.RestApi.AsyncDownloader import rq


def create_app():
    app = Flask(__name__)
    app.register_blueprint(GeneHandlerService_endpoints, url_prefix='/genehandler')
    app.register_blueprint(SubSequenceMatcherService_endpoints, url_prefix='/genefilter')
    app.config['RQ_REDIS_URL'] = 'redis://' + Constants.REDIS_SERVER + ':6379'
    CORS(app)
    rq.init_app(app)
    return app

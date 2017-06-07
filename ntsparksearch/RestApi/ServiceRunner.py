from flask import Flask

from ntsparksearch.Common.Constants import Constants
from ntsparksearch.RestApi.GeneHandlerService import GeneHandlerService_endpoints
from ntsparksearch.RestApi.SubSequenceMatcherService import SubSequenceMatcherService_endpoints


def create_app():
    app = Flask(__name__)
    app.register_blueprint(GeneHandlerService_endpoints, url_prefix='/genehandler')
    app.register_blueprint(SubSequenceMatcherService_endpoints, url_prefix='/genefilter')
    app.config['RQ_REDIS_URL'] = 'redis://' + Constants.REDIS_SERVER + ':6379'
    from ntsparksearch.RestApi.AsyncDownloader import rq
    rq.init_app(app)
    return app


app = create_app()


@app.route("/")
def hello():
    return "Hello World!"


if __name__ == "__main__":
    app.run(use_reloader=False,
            threaded=True, host='0.0.0.0')

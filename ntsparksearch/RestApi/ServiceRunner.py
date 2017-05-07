from flask import Flask
from ntsparksearch.RestApi.GeneRetrieverService import *
from ntsparksearch.RestApi.SubSequenceMatcherService import *


def create_app():
    app = Flask(__name__)
    app.register_blueprint(GeneRetrieverService_endpoints, url_prefix='/generetriever')
    app.register_blueprint(SubSequenceMatcherService_endpoints, url_prefix='/genefilter')
    from ntsparksearch.RestApi.SubSequenceMatcherService import rq
    rq.init_app(app)
    return app


app = create_app()


@app.route("/")
def hello():
    return "Hello World!"


if __name__ == "__main__":
    app.run(use_reloader=False,
            threaded=True, host='0.0.0.0')

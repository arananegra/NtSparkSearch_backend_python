import json

from flask import Flask, Response
from flask_cors import CORS
from gevent import monkey
from gevent.wsgi import WSGIServer

from ntsparksearch.Common.Constants import Constants
from ntsparksearch.RestApi.GeneHandlerService import GeneHandlerService_endpoints
from ntsparksearch.RestApi.SubSequenceMatcherService import SubSequenceMatcherService_endpoints

monkey.patch_all()


def create_app():
    app = Flask(__name__)
    app.register_blueprint(GeneHandlerService_endpoints, url_prefix='/genehandler')
    app.register_blueprint(SubSequenceMatcherService_endpoints, url_prefix='/genefilter')
    app.config['RQ_REDIS_URL'] = 'redis://' + Constants.REDIS_SERVER + ':6379'
    CORS(app)
    from ntsparksearch.RestApi.AsyncDownloader import rq
    rq.init_app(app)
    return app


app = create_app()
data = {"user2_proximity": 3, "Wifi_1": -80, "Wifi_2": -40, "Wifi_3": -40,
        "thermostat": 18, "light": 0, "hour_of_day": 0, "user3_proximity": 3,
        "user1_proximity": 1, "day_of_week": 1, "security": 0, "minute_of_hour": 9,
        "Act_1": 1, "Act_2": 0, "Act_3": 0}
json_data = json.dumps(data)


@app.route("/")
def hello():
    response = Response(json.dumps(data), mimetype='application/json')
    return response, Constants.OK


if __name__ == "__main__":
    # app.run(use_reloader=False,
    #         threaded=True, host='0.0.0.0')
    WSGIServer((
        "0.0.0.0",  # str(HOST)
        5000,  # int(PORT)
    ), app.wsgi_app).serve_forever()

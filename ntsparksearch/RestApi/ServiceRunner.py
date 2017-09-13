import json

from flask import Response, render_template, redirect
from flask_security import auth_token_required, current_user
from gevent import monkey
from gevent.wsgi import WSGIServer

from ntsparksearch.Common.Constants import Constants
from ntsparksearch.RestApi import app

monkey.patch_all()

data = {"user2_proximity": 3, "Wifi_1": -80, "Wifi_2": -40, "Wifi_3": -40,
        "thermostat": 18, "light": 0, "hour_of_day": 0, "user3_proximity": 3,
        "user1_proximity": 1, "day_of_week": 1, "security": 0, "minute_of_hour": 9,
        "Act_1": 1, "Act_2": 0, "Act_3": 0}
json_data = json.dumps(data)


@app.route("/")
@auth_token_required
def hello():
    print("CORREO sacado a partir del token " + current_user.email)
    print("TOKEN sacado a partir del token " + current_user.get_auth_token())
    response = Response(json.dumps(data), mimetype='application/json')
    return response, Constants.OK


@app.route('/login', methods=['GET', 'POST'])
def login():
    return render_template('/security/login_user.html')


@app.route('/register', methods=['GET', 'POST'])
def register():
    return render_template('/security/register_user.html')


@app.route('/redirect-to-frontend', methods=['GET'])
def redirect_to_frontend():
    return redirect("http://0.0.0.0:3002/")


if __name__ == "__main__":
    # app.run(use_reloader=False,
    #         threaded=True, host='0.0.0.0')
    WSGIServer((
        "0.0.0.0",  # str(HOST)
        5000,  # int(PORT)
    ), app.wsgi_app).serve_forever()

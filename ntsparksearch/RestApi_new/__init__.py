from ntsparksearch.RestApi_new.blueprint import example
from ntsparksearch.RestApi_new.extensions import mail
from flask import Flask
from ntsparksearch.RestApi_new import settings


def create_app(settings=settings):
    ret_val = Flask(__name__)
    ret_val.config.from_object(settings)
    # initialize extensions...
    mail.init_app(ret_val)
    # register blueprints...
    ret_val.register_blueprint(example)

    return ret_val
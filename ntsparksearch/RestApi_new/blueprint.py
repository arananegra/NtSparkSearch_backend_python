from flask import Blueprint
from ntsparksearch.RestApi_new.tasks import handle_notification


example = Blueprint("example", __name__)


@example.route('/')
def notify():
    ret_val = "Hello World!"
    handle_notification.apply_async([ret_val])

    return ret_val
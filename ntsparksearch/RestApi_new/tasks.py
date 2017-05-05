from ntsparksearch.RestApi_new import celery
from ntsparksearch.RestApi_new.extensions import mail
from flask_mail import Message


@celery.task
def handle_notification(body):
    message = Message("Test")
    message.body = body

    with mail.connect() as connection:
        connection.send(message)
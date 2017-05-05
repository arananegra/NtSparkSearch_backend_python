from ntsparksearch.RestApi_new import *
from ntsparksearch.RestApi_new import settings
from celery import Celery


celery = Celery()
celery.config_from_object(settings)
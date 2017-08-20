from flask import Flask
from flask_cors import CORS
from flask_mail import Mail
from flask_mongoengine import MongoEngine
from flask_security import Security, MongoEngineUserDatastore, \
    UserMixin, RoleMixin

from ntsparksearch.Common.Constants import Constants
from ntsparksearch.RestApi.GeneHandlerService import GeneHandlerService_endpoints
from ntsparksearch.RestApi.SubSequenceMatcherService import SubSequenceMatcherService_endpoints

app = Flask(__name__)
app.register_blueprint(GeneHandlerService_endpoints, url_prefix='/genehandler')
app.register_blueprint(SubSequenceMatcherService_endpoints, url_prefix='/genefilter')
app.config['RQ_REDIS_URL'] = 'redis://' + Constants.REDIS_SERVER + ':6379'
app.config['DEBUG'] = True
app.config['SECRET_KEY'] = 'super-secret'

app.config.update(
    MAIL_SERVER='smtp.gmail.com',
    MAIL_PORT=587,
    MAIL_USE_SSL=False,
    MAIL_USERNAME=Constants.MAIL_USER,
    MAIL_PASSWORD=Constants.MAIL_PASS,
    MAIL_USE_TLS=True,
    DEFAULT_MAIL_SENDER='ntsparksearch@gmail.com')

app.config['SECURITY_PASSWORD_HASH'] = 'sha256_crypt'
app.config['SECURITY_PASSWORD_SALT'] = "wjfiwjfiwamcw3214awcq932"
app.config['SECURITY_REGISTERABLE'] = True
app.config['SECURITY_EMAIL_SENDER'] = "no-reply@localhost."
app.config['SECURITY_REGISTER_URL'] = '/register'
app.config['SECURITY_LOGIN_USER_TEMPLATE'] = 'security/login_user.html'
app.config['SECURITY_REGISTER_USER_TEMPLATE'] = 'security/register_user.html'
app.config['SECURITY_RESET_PASSWORD_TEMPLATE'] = 'security/reset_password.html'
app.config['SECURITY_CHANGE_PASSWORD_TEMPLATE'] = 'security/change_password.html'

app.config['SECURITY_TRACKABLE'] = True
app.config['WTF_CSRF_ENABLED'] = False

app.config['MONGODB_DB'] = Constants.MONGODB_DB_NAME
app.config['MONGODB_HOST'] = Constants.MONGODB_HOST
app.config['MONGODB_PORT'] = Constants.MONGODB_PORT
db = MongoEngine(app)
CORS(app)
mail = Mail(app)
from ntsparksearch.RestApi.AsyncDownloader import rq

rq.init_app(app)


class Role(db.Document, RoleMixin):
    name = db.StringField(max_length=80, unique=True)
    description = db.StringField(max_length=255)


class User(db.Document, UserMixin):
    email = db.StringField(max_length=255)
    password = db.StringField(max_length=255)
    active = db.BooleanField(default=True)
    confirmed_at = db.DateTimeField()
    current_login_at = db.DateTimeField()
    current_login_ip = db.StringField(max_length=45)
    login_count = db.IntField()
    roles = db.ListField(db.ReferenceField(Role), default=[])


user_datastore = MongoEngineUserDatastore(db, User, Role)
security = Security(app, user_datastore)

from flask import Blueprint

SubSequenceMatcherService_endpoints = Blueprint('SubSequenceMatcherService', __name__)

@SubSequenceMatcherService_endpoints.route("/account")
def accountList():
    return "list of accounts"
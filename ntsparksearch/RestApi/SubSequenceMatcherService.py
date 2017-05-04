from flask import Blueprint, abort

SubSequenceMatcherService_endpoints = Blueprint('SubSequenceMatcherService', __name__)

@SubSequenceMatcherService_endpoints.route("/account")
def accountList():
    return "list of accounts"


@SubSequenceMatcherService_endpoints.route("/failed")
def test():
    return abort(401)

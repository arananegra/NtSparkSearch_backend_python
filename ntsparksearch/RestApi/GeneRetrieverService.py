import os
import json
from flask import Blueprint, request, Response, redirect, url_for, send_from_directory, jsonify
from werkzeug.utils import secure_filename
from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS
from ntsparksearch.Common.Constants import Constants

GeneRetrieverService_endpoints = Blueprint('GeneRetrieverService', __name__)


# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in Constants.ALLOWED_EXTENSIONS


# Route that will process the file upload
@GeneRetrieverService_endpoints.route('/upload_file', methods=['POST', 'GET'])
def upload_file():
    if request.method == 'POST':
        # Get the name of the uploaded file
        file = request.files['file']
        # Check if the file is one of the allowed types/extensions
        if file and allowed_file(file.filename):
            # Make the filename safe, remove unsupported chars
            filename = secure_filename(file.filename)
            # Move the file form the temporal folder to
            # the upload folder we setup
            file.save(os.path.join(Constants.UPLOAD_FOLDER, filename))
            # Redirect the user to the uploaded_file route, which
            # will basicaly show on the browser the uploaded file
            # return redirect(url_for('uploaded_file',
            #                   filename=filename))
            return "file uploaded succesfully"
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    '''


@GeneRetrieverService_endpoints.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(Constants.UPLOAD_FOLDER,
                               filename)


@GeneRetrieverService_endpoints.route('/ounfiltered', methods=['GET'])
def obtainUnfiltered():
    retriever_BS = GeneRetrieverBS()
    list_of_genes = retriever_BS.get_list_of_ids_from_mongo()
    return Response(json.dumps(list_of_genes), mimetype='application/json')

    # if __name__ == '__main__':
    #     app.run(
    #         #host="0.0.0.0",
    #         port=int("5000"),
    #         debug=True
    #     )

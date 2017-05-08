import os
import json
from flask import Blueprint, request, Response, redirect, url_for, send_from_directory, jsonify
from ntsparksearch.RestApi.AsyncDownloader import gene_downloader_async
from werkzeug.utils import secure_filename
from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS
from ntsparksearch.Common.Constants import Constants
from ntsparksearch.EmailProcess.EmailSender import EmailSender

GeneRetrieverService_endpoints = Blueprint('GeneRetrieverService', __name__)
retriever_BS = GeneRetrieverBS()

email_manager = EmailSender()

def allowed_file(filename, extension):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in extension


@GeneRetrieverService_endpoints.route('/unfiltered', methods=['GET'])
def obtainUnfiltered():
    list_of_genes = retriever_BS.get_list_of_ids_from_mongo()
    return Response(json.dumps(list_of_genes), mimetype='application/json')


@GeneRetrieverService_endpoints.route('/upload-excel', methods=['POST', 'GET'])
def upload_excel_file_and_download_genes():

    try:
        json_with_gene_ids_sequence_to_filter_and_mail = request.get_json()
        email_receiver = json_with_gene_ids_sequence_to_filter_and_mail["email"]

    except Exception:
        print("Warning: Email account not set")
        pass

    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename, Constants.EXCEL_EXTENSION):
            filename = secure_filename(file.filename)
            file.save(os.path.join(Constants.UPLOAD_FOLDER, filename))

            retriever_BS.insert_in_collection_from_excel(Constants.UPLOAD_FOLDER + filename)
            list_of_genes_without_sequence = retriever_BS.get_list_of_ids_from_mongo_without_sequence()

            if list_of_genes_without_sequence is not None:
                if email_receiver is not None:
                    email_manager.receivers = email_receiver
                    email_manager.send_email_download_initialize(list_of_genes_without_sequence)

                gene_downloader_async.queue(list_of_genes_without_sequence, email_receiver)

        return Response(), Constants.POST_WAIT

    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    '''


# @GeneRetrieverService_endpoints.route('/download-list', methods=['POST'])
# def upload_excel_file():
#     if request.method == 'POST':
#         # insertar lista en unfiltered y obtener de las que no tenga secuencias
#         list_of_gene_ids = request.get_json()
#
#         dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
#             list_of_genes_without_sequence)
#
#



@GeneRetrieverService_endpoints.route('/upload-fasta', methods=['POST', 'GET'])
def upload_fasta_file():
    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename, Constants.FASTA_EXTENSION):
            filename = secure_filename(file.filename)
            file.save(os.path.join(Constants.UPLOAD_FOLDER, filename))

            retriever_BS.insert_in_collection_from_fasta(Constants.UPLOAD_FOLDER + filename)
            return Response(), Constants.POST_CREATED
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    '''

# para recuperar un archivo
@GeneRetrieverService_endpoints.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(Constants.UPLOAD_FOLDER,
                               filename)

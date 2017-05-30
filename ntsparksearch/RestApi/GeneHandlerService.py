import os
import json
import string
import random
from flask import Blueprint, request, Response, send_from_directory, send_file, make_response
from ntsparksearch.RestApi.AsyncDownloader import gene_downloader_async_from_list
from werkzeug.utils import secure_filename
from ntsparksearch.GeneHandler.GeneHandlerBS import GeneHandlerBS
from ntsparksearch.Common.Constants import Constants
from ntsparksearch.EmailProcess.EmailSender import EmailSender

GeneHandlerService_endpoints = Blueprint('GeneHandlerService', __name__)
gene_handler_BS = GeneHandlerBS()
email_manager = EmailSender()


def allowed_file(filename, extension):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in extension


@GeneHandlerService_endpoints.route('/unfiltered', methods=['GET'])
def obtain_unfiltered():
    list_of_genes = gene_handler_BS.get_list_of_ids_from_mongo_unfiltered()
    return Response(json.dumps(list_of_genes), mimetype='application/json')


@GeneHandlerService_endpoints.route('/upload-excel', methods=['POST', 'GET'])
def upload_excel_file_and_download_genes():
    try:
        email_receiver = request.args.getlist("email")

    except Exception:
        print("Warning: Email account not set")
        pass

    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename, Constants.EXCEL_EXTENSION):
            filename = secure_filename(file.filename)
            file.save(os.path.join(Constants.UPLOAD_FOLDER, filename))

            gene_handler_BS.insert_in_unfiltered_collection_from_excel(Constants.UPLOAD_FOLDER + filename)
            list_of_genes_without_sequence = gene_handler_BS.get_list_of_ids_from_mongo_unfiltered_without_sequence()

            if list_of_genes_without_sequence is not None:
                if email_receiver is not None:
                    email_manager.receivers = email_receiver
                    email_manager.send_email_download_initialize(list_of_genes_without_sequence)

                gene_downloader_async_from_list.queue(list_of_genes_without_sequence, email_receiver)

        return Response(), Constants.OK_WAIT

    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    '''


@GeneHandlerService_endpoints.route('/upload-fasta', methods=['POST', 'GET'])
def upload_fasta_file():
    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename, Constants.FASTA_EXTENSION):
            filename = secure_filename(file.filename)
            file.save(os.path.join(Constants.UPLOAD_FOLDER, filename))

            gene_handler_BS.insert_in_unfiltered_collection_from_fasta(Constants.UPLOAD_FOLDER + filename)
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


@GeneHandlerService_endpoints.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(Constants.UPLOAD_FOLDER,
                               filename)


def id_generator(size=6, chars=string.ascii_lowercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


@GeneHandlerService_endpoints.route('/download-fasta-unfiltered')
def download_fasta_file_unfiltered():
    random_string_name = id_generator()
    gene_handler_BS.export_unfiltered_genes_collection_to_fasta(Constants.OUTPUT_FOLDER
                                                                + random_string_name)

    response = make_response(send_file(Constants.OUTPUT_FOLDER + random_string_name + "." + Constants.FASTA_EXTENSION,
                                       mimetype='text/plain'))

    response.headers["Content-Disposition"] = "attachment; filename=" + \
                                              random_string_name \
                                              + "." + Constants.FASTA_EXTENSION
    return response


@GeneHandlerService_endpoints.route('/download-fasta-filtered')
def download_fasta_file_filtered():
    random_string_name = id_generator()
    gene_handler_BS.export_filtered_genes_collection_to_fasta(Constants.OUTPUT_FOLDER
                                                              + random_string_name)

    response = make_response(send_file(Constants.OUTPUT_FOLDER + random_string_name + "." + Constants.FASTA_EXTENSION,
                                       mimetype='text/plain'))

    response.headers["Content-Disposition"] = "attachment; filename=" + \
                                              random_string_name \
                                              + "." + Constants.FASTA_EXTENSION
    return response


@GeneHandlerService_endpoints.route('/download-id-unfiltered')
def download_id_file_unfiltered():
    random_string_name = id_generator()
    gene_handler_BS.export_unfiltered_genes_collection_to_file_with_just_ids(Constants.OUTPUT_FOLDER
                                                                             + random_string_name)

    response = make_response(send_file(Constants.OUTPUT_FOLDER + random_string_name + "." + Constants.ID_EXTENSION,
                                       mimetype='text/plain'))

    response.headers["Content-Disposition"] = "attachment; filename=" + \
                                              random_string_name \
                                              + "." + Constants.ID_EXTENSION
    return response


@GeneHandlerService_endpoints.route('/download-id-filtered')
def download_id_file_filtered():
    random_string_name = id_generator()
    gene_handler_BS.export_filtered_genes_collection_to_file_with_just_ids(Constants.OUTPUT_FOLDER
                                                                           + random_string_name)

    response = make_response(send_file(Constants.OUTPUT_FOLDER + random_string_name + "." + Constants.ID_EXTENSION,
                                       mimetype='text/plain'))

    response.headers["Content-Disposition"] = "attachment; filename=" + \
                                              random_string_name \
                                              + "." + Constants.ID_EXTENSION
    return response


@GeneHandlerService_endpoints.route('/delete-unfiltered', methods=['DELETE'])
def delete_unfiltered_collection():
    gene_handler_BS.delete_unfiltered_collection()
    return Response(), Constants.OK


@GeneHandlerService_endpoints.route('/delete-filtered', methods=['DELETE'])
def delete_filtered_collection():
    gene_handler_BS.delete_filtered_collection()
    return Response(), Constants.OK

from flask import Blueprint, Response, request
import json
from ntsparksearch.SubsequenceMatcher.SubSequenceSparkMatcherBS import SubSequenceSparkMatcherBS
from ntsparksearch.GeneHandler.GeneHandlerBS import GeneHandlerBS
from ntsparksearch.Common.Constants import Constants
from ntsparksearch.EmailProcess.EmailSender import EmailSender
from ntsparksearch.RestApi.AsyncDownloader import gene_downloader_async_from_list

SubSequenceMatcherService_endpoints = Blueprint('SubSequenceMatcherService', __name__)

subsequence_matcher_BS = SubSequenceSparkMatcherBS()
gene_handler_BS = GeneHandlerBS()
email_manager = EmailSender()


@SubSequenceMatcherService_endpoints.route('/sparkmatchall', methods=["GET"])
def spark_matcher():
    sequence_to_filter = request.args.get("sequence")
    dict_filtered_with_spark = subsequence_matcher_BS. \
        filter_sequences_by_sequence_string_to_dict(sequence_to_filter)
    gene_handler_BS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)
    dict_filtered_with_ones = {x: 1 for x in dict_filtered_with_spark}

    list_of_genes_pass_filter = list(dict_filtered_with_ones.keys())

    list_of_unfiltered_ids_from_mongo = gene_handler_BS.get_list_of_ids_from_mongo_unfiltered()

    for gene_unfiltered_id in list_of_unfiltered_ids_from_mongo:
        if gene_unfiltered_id not in list_of_genes_pass_filter:
            dict_filtered_with_ones[gene_unfiltered_id] = 0

    return Response(json.dumps(dict_filtered_with_ones), mimetype='application/json'), Constants.OK


@SubSequenceMatcherService_endpoints.route('/genes-checker', methods=["GET"])
def genes_in_unfiltered_checker():
    try:
        list_of_ids_from_request = request.args.getlist("geneIds")
        sequence_to_filter = request.args.get("sequence")

        list_of_genes_not_in_unfiltered = gene_handler_BS.check_gene_id_list_existance_on_unfiltered_from_list(
            list_of_ids_from_request)

        if len(list_of_genes_not_in_unfiltered) == 0:
            dict_filtered_with_spark = subsequence_matcher_BS. \
                filter_sequences_by_sequence_string_to_dict(sequence_to_filter)
            gene_handler_BS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)
            dict_filtered_with_ones = {x: 1 for x in dict_filtered_with_spark}

            list_of_genes_pass_filter = list(dict_filtered_with_ones.keys())

            list_of_unfiltered_ids_from_mongo = gene_handler_BS.get_list_of_ids_from_mongo_unfiltered()

            for gene_unfiltered_id in list_of_unfiltered_ids_from_mongo:
                if gene_unfiltered_id not in list_of_genes_pass_filter:
                    dict_filtered_with_ones[gene_unfiltered_id] = 0

            return Response(json.dumps(dict_filtered_with_ones), mimetype='application/json'), Constants.OK

        else:

            return Response(), Constants.OK_WAIT

    except Exception as ex:
        print(ex)


@SubSequenceMatcherService_endpoints.route('/genes-downloader', methods=["GET"])
def genes_downloader():
    try:
        list_of_ids_from_request = request.args.getlist("geneIds")
        email_receiver = request.args.getlist("email")

        list_of_genes_not_in_unfiltered = gene_handler_BS.check_gene_id_list_existance_on_unfiltered_from_list(
            list_of_ids_from_request)

        if list_of_genes_not_in_unfiltered != 0:

            if email_receiver is not None:
                email_manager.receivers = email_receiver
                email_manager.send_email_download_initialize(list_of_genes_not_in_unfiltered)

            gene_downloader_async_from_list.queue(list_of_genes_not_in_unfiltered, email_receiver)

            return Response(), Constants.OK_WAIT

        else:

            return Response(), Constants.SERVER_ERROR

    except Exception as ex:
        print(ex)
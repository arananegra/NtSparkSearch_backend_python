from flask import Blueprint, Response, request
import json
from ntsparksearch.SubsequenceMatcher.SubSequenceSparkMatcherBS import SubSequenceSparkMatcherBS
from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS
from ntsparksearch.Common.Constants import Constants
from ntsparksearch.EmailProcess.EmailSender import EmailSender
from ntsparksearch.RestApi.AsyncDownloader import gene_downloader_async

SubSequenceMatcherService_endpoints = Blueprint('SubSequenceMatcherService', __name__)

subsequence_matcher_BS = SubSequenceSparkMatcherBS()
retriever_BS = GeneRetrieverBS()
email_manager = EmailSender()


@SubSequenceMatcherService_endpoints.route('/sparkmatchall', methods=["GET"])
def spark_matcher():
    json_with_gene_ids_sequence_to_filter_and_mail = request.get_json()
    sequence_to_filter = json_with_gene_ids_sequence_to_filter_and_mail["sequence"]

    dict_filtered_with_spark = subsequence_matcher_BS. \
        filter_sequences_by_sequence_string_to_dict(sequence_to_filter)

    subsequence_matcher_BS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)
    dict_filtered_with_ones = {x: 1 for x in dict_filtered_with_spark}
    return Response(json.dumps(dict_filtered_with_ones), mimetype='application/json'), Constants.GET_OK


@SubSequenceMatcherService_endpoints.route('/sparkmatchlist', methods=["POST"])
def spark_matcher_from_list():
    json_with_gene_ids_sequence_to_filter_and_mail = request.get_json()

    try:
        email_receiver = json_with_gene_ids_sequence_to_filter_and_mail["email"]

    except Exception:
        print("Warning: Email account not set")
        pass

    try:
        list_of_ids = json_with_gene_ids_sequence_to_filter_and_mail["geneIds"]
        sequence_to_filter = json_with_gene_ids_sequence_to_filter_and_mail["sequence"]

    except Exception:
        print("ERROR: A list of gene ids and sequences must be provided")
        return Response(), Constants.POST_BAD_REQUEST

    list_of_genes_already_filtered = subsequence_matcher_BS.get_list_of_ids_from_mongo_filtered()

    if list_of_genes_already_filtered is not None:

        genes_have_been_already_filtered_before = True

        for possible_gene_id_to_download_or_filter in list_of_ids:
            if possible_gene_id_to_download_or_filter not in list_of_genes_already_filtered:
                genes_have_been_already_filtered_before = False

        if genes_have_been_already_filtered_before == True:
            dict_with_genes_already_filtered = subsequence_matcher_BS.get_dict_from_filtered_with_sequences()
            dict_filtered_with_ones = {x: 1 for x in dict_with_genes_already_filtered}
            return Response(json.dumps(dict_filtered_with_ones), mimetype='application/json'), Constants.POST_CREATED

    list_of_ids_which_exists_on_NCBI = retriever_BS.check_gene_id_existance_on_ncbi_from_list(list_of_ids)
    retriever_BS.insert_in_collection_from_list_of_ids(list_of_ids_which_exists_on_NCBI)
    list_of_genes_without_sequence = retriever_BS.get_list_of_ids_from_mongo_without_sequence()

    if len(list_of_genes_without_sequence) != 0:

        if email_receiver is not None:
            email_manager.receivers = email_receiver
            email_manager.send_email_download_initialize(list_of_genes_without_sequence)

        gene_downloader_async.queue(list_of_genes_without_sequence, email_receiver)

        return Response(), Constants.POST_WAIT

    else:

        dict_filtered_with_spark = subsequence_matcher_BS. \
            filter_sequences_by_sequence_string_to_dict(sequence_to_filter)
        subsequence_matcher_BS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)
        dict_filtered_with_ones = {x: 1 for x in dict_filtered_with_spark}

        return Response(json.dumps(dict_filtered_with_ones), mimetype='application/json'), Constants.POST_CREATED

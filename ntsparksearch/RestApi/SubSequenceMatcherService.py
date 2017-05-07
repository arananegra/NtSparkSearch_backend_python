from flask import Blueprint, Response, request
import json
from flask_rq2 import RQ
from ntsparksearch.SubsequenceMatcher.SubSequenceSparkMatcherBS import SubSequenceSparkMatcherBS
from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS
from ntsparksearch.Common.Constants import Constants
from ntsparksearch.EmailProcess.EmailSender import EmailSender

SubSequenceMatcherService_endpoints = Blueprint('SubSequenceMatcherService', __name__)

subsequence_matcher_BS = SubSequenceSparkMatcherBS()
retriever_BS = GeneRetrieverBS()
email_manager = EmailSender()
rq = RQ()


@rq.job
def gene_downloader_async(list_of_genes_without_sequence, email_receiver):
    dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
        list_of_genes_without_sequence)
    retriever_BS.update_genes_from_dict(dict_of_genes_complete)
    email_manager.receivers = email_receiver
    email_manager.send_email_download_finished(list_of_genes_without_sequence)


@SubSequenceMatcherService_endpoints.route('/sparkmatch', methods=["POST"])
def spark_matcher():
    json_with_gene_ids_sequence_to_filter_and_mail = request.get_json()
    list_of_ids = json_with_gene_ids_sequence_to_filter_and_mail["geneIds"]
    sequence_to_filter = json_with_gene_ids_sequence_to_filter_and_mail["sequence"]
    email_receiver = json_with_gene_ids_sequence_to_filter_and_mail["email"]

    list_of_genes_without_sequence = retriever_BS.get_list_of_ids_from_mongo_without_sequence()

    if list_of_genes_without_sequence is None:
        retriever_BS.insert_in_collection_from_list_of_ids(list_of_ids)

    if len(list_of_genes_without_sequence) != 0 or list_of_genes_without_sequence is None:
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

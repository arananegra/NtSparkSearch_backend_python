from flask import Blueprint, abort, Response, current_app, request
import redis
import json
from ntsparksearch.SubsequenceMatcher.SubSequenceSparkMatcherBS import SubSequenceSparkMatcherBS
from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS
from ntsparksearch.Common.Constants import Constants

SubSequenceMatcherService_endpoints = Blueprint('SubSequenceMatcherService', __name__)

subsequence_matcher_BS = SubSequenceSparkMatcherBS()
retriever_BS = GeneRetrieverBS()


@SubSequenceMatcherService_endpoints.route('/sparkmatch', methods=["POST"])
def spark_matcher():
    json_with_gene_ids_and_sequence_to_filter = request.get_json()
    list_of_ids = json_with_gene_ids_and_sequence_to_filter["geneIds"]
    sequence_to_filter = json_with_gene_ids_and_sequence_to_filter["sequence"]

    retriever_BS.insert_in_collection_from_list_of_ids(list_of_ids)

    list_of_genes_without_sequence = retriever_BS.get_list_of_ids_from_mongo_without_sequence()

    if list_of_genes_without_sequence is not None:
        # inicio de proceso asincrono
        # TODO: lanzar mensaje de COMIENZO de descarga --> se han encontrado genes que hay que descargar

        dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
            list_of_genes_without_sequence)

        # TODO: lanzar mensaje de FIN de descarga --> se han encontrado genes que hay que descargar
        retriever_BS.update_genes_from_dict(dict_of_genes_complete)
        return Response(), Constants.POST_WAIT

    else:

        dict_filtered_with_spark = subsequence_matcher_BS. \
            filter_sequences_by_sequence_string_to_dict(sequence_to_filter)
        subsequence_matcher_BS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)
        dict_filtered_with_ones = {x: 1 for x in dict_filtered_with_spark}
        return Response(json.dumps(dict_filtered_with_ones), mimetype='application/json'), Constants.POST_CREATED

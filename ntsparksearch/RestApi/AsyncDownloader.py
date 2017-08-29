from flask_rq2 import RQ
from flask_security import current_user

from ntsparksearch.Common.Constants import Constants
from ntsparksearch.EmailProcess.EmailSender import EmailSender
from ntsparksearch.GeneHandler.GeneHandlerBS import GeneHandlerBS

email_manager = EmailSender()
rq = RQ()


@rq.job
def gene_downloader_async_from_list(list_of_genes_without_sequence, email_receiver, user_id_str):
    try:
        gene_handler_BS = GeneHandlerBS(Constants.MONGODB_DB_INITIAL + user_id_str)
        gene_handler_BS.insert_in_unfiltered_collection_from_list_of_ids(list_of_genes_without_sequence)
        dict_of_genes_complete = gene_handler_BS.download_sequences_from_list_as_dict_from_NCBI(
            list_of_genes_without_sequence)
        gene_handler_BS.update_genes_from_dict(dict_of_genes_complete)

        if email_receiver is not None:
            email_manager.receivers = email_receiver
            email_manager.send_email_download_finished(list_of_genes_without_sequence)
    except Exception as ex:
        print(ex)

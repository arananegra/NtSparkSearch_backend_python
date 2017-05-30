from flask_rq2 import RQ
from ntsparksearch.GeneHandler.GeneHandlerBS import GeneHandlerBS
from ntsparksearch.EmailProcess.EmailSender import EmailSender

gene_handler_bs = GeneHandlerBS()
email_manager = EmailSender()
rq = RQ()


@rq.job
def gene_downloader_async_from_list(list_of_genes_without_sequence, email_receiver):

    gene_handler_bs.insert_in_unfiltered_collection_from_list_of_ids(list_of_genes_without_sequence)
    dict_of_genes_complete = gene_handler_bs.download_sequences_from_list_as_dict_from_NCBI(
        list_of_genes_without_sequence)
    gene_handler_bs.update_genes_from_dict(dict_of_genes_complete)
    if email_receiver is not None:
        email_manager.receivers = email_receiver
        email_manager.send_email_download_finished(list_of_genes_without_sequence)

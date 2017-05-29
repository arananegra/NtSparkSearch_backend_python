from flask_rq2 import RQ
from ntsparksearch.GeneHandler.GeneHandlerBS import GeneHandlerBS
from ntsparksearch.EmailProcess.EmailSender import EmailSender

retriever_BS = GeneHandlerBS()
email_manager = EmailSender()
rq = RQ()


@rq.job
def gene_downloader_async_from_list(list_of_genes_without_sequence, email_receiver):

    retriever_BS.insert_in_collection_from_list_of_ids(list_of_genes_without_sequence)
    dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
        list_of_genes_without_sequence)
    retriever_BS.update_genes_from_dict(dict_of_genes_complete)
    if email_receiver is not None:
        email_manager.receivers = email_receiver
        email_manager.send_email_download_finished(list_of_genes_without_sequence)

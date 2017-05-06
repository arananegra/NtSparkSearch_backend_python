from ntsparksearch.DownloadScheduler.QueueBS import QueueBS
from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS


class DownloadServiceScheduler(object):
    def downloader_schedulirazed(self, list_of_genes_without_sequence):
        retriever_BS = GeneRetrieverBS()

        dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
            list_of_genes_without_sequence)

from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS
import schedule
import time


def downloader_schedulirazed():
    # retriever_BS = GeneRetrieverBS()
    #
    # dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
    #     list_of_genes_without_sequence)
    #
    # retriever_BS.update_genes_from_dict(dict_of_genes_complete)

    schedule.every(2).seconds.do(tonteria)

    while True:
        schedule.run_continuously()
        time.sleep(1)

def tonteria():
    print("pintalo")
from wsgiref import validate

from ntsparksearch.NCBIretriever.api import INCBIretriever
from ntsparksearch.NCBIretriever.dao import FileRetriverDAO
from ntsparksearch.common.dao import NCBItoMongoDAO
from ntsparksearch.common.domain import NCBIsearcher
from ntsparksearch.common.domain import NucleotidesFromNCBI
from ntsparksearch.common.util import Constants
from progressbar import ProgressBar
from Bio import Entrez
from pymongo import MongoClient


class NCBIretrieverBS(INCBIretriever):
    def insert_in_collection_from_excel(self, file_path: str, sheet: str, column_name: str) -> None:

        try:

            if file_path is not None:
                file_retriever_and_mongo_manager = FileRetriverDAO(
                    MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                    Constants.MONGODB_DB_NAME,
                    Constants.MONGODB_COLLECTION_UNFILTERED,
                    file_path)

                list_of_genes_just_ids = file_retriever_and_mongo_manager.get_list_of_genes_from_xlrd(sheet,
                                                                                                      column_name)

                for ncbi_object_just_id in list_of_genes_just_ids:

                    criteria = NCBIsearcher()
                    criteria.search_by_NCBI_id_criteria = ncbi_object_just_id.idNcbi

                    if file_retriever_and_mongo_manager.search_ncbi_objects_and_return_as_list(criteria) is None:
                        file_retriever_and_mongo_manager.insert_ncbi_document_from_object(ncbi_object_just_id)
            else:
                print("No file attached to this operation")

        except Exception as error:
            print('Caught exception when inserting all data from excel: ' + repr(error))

    def obtain_list_of_ids_from_mongo(self) -> list:
        list_of_just_ids = None

        try:
            list_of_just_ids = []

            mongo_dao_retriever = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            list_of_ncbi_objets = mongo_dao_retriever.get_all_ncbi_objects_as_list()

            for ncbi_object in list_of_ncbi_objets:
                ncbi_id = ncbi_object.idNcbi
                list_of_just_ids.append(ncbi_id)

            return list_of_just_ids

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list: ' + repr(error))

    def update_genes_from_dict(self, dict_of_genes: dict) -> None:

        try:

            if dict_of_genes is not None:
                mongo_dao_retriever = NCBItoMongoDAO(
                    MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                    Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

                for id_ncbi, sequence in dict_of_genes.items():
                    ncbi_object_to_update = NucleotidesFromNCBI()
                    ncbi_object_to_update.idNcbi = id_ncbi
                    ncbi_object_to_update.sequence = sequence
                    mongo_dao_retriever.update_ncbi_element_from_object(ncbi_object_to_update, False)

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list: ' + repr(error))

    def download_sequences_from_list_as_dict(self, list_of_genes: list) -> dict:

        pbar = ProgressBar()
        Entrez.email = "notfunny@notanemail.org"
        dict_id_and_sequences = {}
        for gene_id in pbar(list_of_genes):
            try:
                handle1 = Entrez.efetch(db="gene", id=gene_id, retmode="xml", validate=False)
            except Exception:
                print("gene " + gene_id + " not found")
                continue
            try:
                record = Entrez.read(handle1)
                geneLoci = record[0]["Entrezgene_locus"]
                geneRegion = geneLoci[0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
            except (KeyError, Exception):
                continue
            startPos = int(geneRegion["Seq-interval_from"]) + 1
            endPos = int(geneRegion["Seq-interval_to"]) + 1
            intervalId = geneRegion["Seq-interval_id"]["Seq-id"]["Seq-id_gi"]
            strandSense = geneRegion["Seq-interval_strand"]["Na-strand"].attributes["value"]
            strandSense = 2 if strandSense.lower() == "minus" else 1
            handle1.close()

            handle2 = Entrez.efetch(db="nucleotide", id=intervalId, rettype="fasta", retmode="text", seq_start=startPos,
                                    seq_stop=endPos, strand=strandSense)
            fasta = handle2.read()
            fastaString = str(fasta)
            fastaWithoutId = '\n'.join(fastaString.split('\n')[1:])
            fastaOneLine = fastaWithoutId.replace('\n', '')

            dict_id_and_sequences[gene_id] = fastaOneLine
            handle2.close()

        return dict_id_and_sequences

# test_bs = NCBIretrieverBS()

# test_bs.insert_in_collection_from_excel(0, "gene_id")

# list_of_genes = test_bs.obtain_list_of_ids_from_mongo()
#
# test_bs.download_sequences_from_list_as_dict(list_of_genes)

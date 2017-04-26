from ntsparksearch.NCBIretriever.INCBIretriever import INCBIretriever
from ntsparksearch.NCBIretriever.FileRetriverDAO import FileRetriverDAO
from ntsparksearch.common.NCBItoMongoDAO import NCBItoMongoDAO
from ntsparksearch.common.NucleotidesFromNCBIDTO import NCBIsearcher
from ntsparksearch.common.NucleotidesFromNCBIDTO import NucleotidesFromNCBIDTO
from ntsparksearch.common.Constants import Constants
from progressbar import ProgressBar
import os.path
from Bio import Entrez
from pymongo import MongoClient


class NCBIretrieverBS(INCBIretriever):
    def insert_in_collection_from_excel(self, file_path: str, sheet: str, column_name: str) -> None:

        try:

            if os.path.isfile(file_path) is True:
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
                raise OSError('The file does not exists')

        except Exception as error:
            print('Caught exception when inserting all data from excel: ' + repr(error))

    def insert_in_collection_from_fasta(self, file_path: str) -> None:

        try:

            if os.path.isfile(file_path) is True:
                file_retriever_and_mongo_manager = FileRetriverDAO(
                    MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                    Constants.MONGODB_DB_NAME,
                    Constants.MONGODB_COLLECTION_UNFILTERED,
                    file_path)

                list_of_ncbi_objects_fasta = file_retriever_and_mongo_manager.get_list_of_ncbi_objects_from_multi_fasta()

                for ncbi_object in list_of_ncbi_objects_fasta:

                    criteria = NCBIsearcher()
                    criteria.search_by_NCBI_id_criteria = ncbi_object.idNcbi

                    if file_retriever_and_mongo_manager.search_ncbi_objects_and_return_as_list(criteria) is None:
                        file_retriever_and_mongo_manager.insert_ncbi_document_from_object(ncbi_object)
            else:
                raise OSError('The file does not exists')

        except Exception as error:
            print('Caught exception when inserting all data from fasta: ' + repr(error))

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
            print('Caught exception when getting all ids from mongo as list (unfiltered): ' + repr(error))

    def obtain_list_of_ids_from_mongo_without_sequence(self) -> list:
        list_of_just_ids = None

        try:
            list_of_just_ids = []

            mongo_dao_retriever = NCBItoMongoDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            list_of_ncbi_objets = mongo_dao_retriever.get_all_ncbi_objects_as_list()

            for ncbi_object in list_of_ncbi_objets:
                if ncbi_object.sequence is None:
                    ncbi_id = ncbi_object.idNcbi
                    list_of_just_ids.append(ncbi_id)

            return list_of_just_ids

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list without sequence (unfiltered): ' + repr(
                error))

    def update_genes_from_dict(self, dict_of_genes: dict) -> None:

        try:

            if dict_of_genes is not None:
                mongo_dao_retriever = NCBItoMongoDAO(
                    MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                    Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

                for id_ncbi, sequence in dict_of_genes.items():
                    ncbi_object_to_update = NucleotidesFromNCBIDTO()
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
                record = Entrez.read(handle1)
                gene_loci = record[0]["Entrezgene_locus"]
                gene_region = gene_loci[0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
                start_pos = int(gene_region["Seq-interval_from"]) + 1
                end_pos = int(gene_region["Seq-interval_to"]) + 1
                interval_id = gene_region["Seq-interval_id"]["Seq-id"]["Seq-id_gi"]
                strand_sense = gene_region["Seq-interval_strand"]["Na-strand"].attributes["value"]
                strand_sense = 2 if strand_sense.lower() == "minus" else 1
                handle1.close()

                handle2 = Entrez.efetch(db="nucleotide", id=interval_id, rettype="fasta", retmode="text",
                                        seq_start=start_pos,
                                        seq_stop=end_pos, strand=strand_sense)
                fasta = handle2.read()
                fasta_string = str(fasta)
                fasta_without_id = '\n'.join(fasta_string.split('\n')[1:])
                fasta_one_line = fasta_without_id.replace('\n', '')

                dict_id_and_sequences[gene_id] = fasta_one_line
                handle2.close()

            except Exception as error:
                print('\nCaught exception when trying to download sequence from NCBI at gene ' +
                      gene_id + " : " + repr(error) + " This sequence will be skipped")
                continue

        return dict_id_and_sequences

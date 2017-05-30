from ntsparksearch.GeneHandler.IGeneHandler import IGeneHandler
from ntsparksearch.GeneHandler.GeneHandlerDAO import GeneHandlerDAO
from ntsparksearch.Common.GeneDAO import GeneDAO
from ntsparksearch.Common.GeneDTO import GeneSearcher
from ntsparksearch.Common.GeneDTO import GeneDTO
from ntsparksearch.Common.Constants import Constants
from progressbar import ProgressBar
from copy import deepcopy
import os.path
from Bio import Entrez
from Bio import SeqIO
from pymongo import MongoClient


class GeneHandlerBS(IGeneHandler):
    def get_dict_from_filtered_with_sequences(self) -> dict:
        try:
            dict_with_genes = {}

            mongo_dao_retriever = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            dict_with_genes = mongo_dao_retriever.get_all_gene_objects_as_dict()

            deep_copy_dict_with_genes = deepcopy(dict_with_genes)

            for id_gene, sequence in deep_copy_dict_with_genes.items():
                if sequence is None:
                    del dict_with_genes[id_gene]

            return dict_with_genes

        except Exception as error:
            print('Caught exception getting filtered sequences (to dict):' + repr(error))

    def get_dict_from_unfiltered_with_sequences(self) -> dict:
        dict_with_genes = None

        try:
            dict_with_genes = {}

            mongo_dao_retriever = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            dict_with_genes = mongo_dao_retriever.get_all_gene_objects_as_dict()

            deep_copy_dict_with_genes = deepcopy(dict_with_genes)

            if dict_with_genes is not None:

                for id_gene, sequence in deep_copy_dict_with_genes.items():
                    if sequence is None:
                        del dict_with_genes[id_gene]

            return dict_with_genes

        except Exception as error:
            print('Caught exception getting unfiltered sequences (to dict):' + repr(error))

    def get_list_of_ids_from_mongo_unfiltered(self) -> list:
        list_of_just_ids = None

        try:
            list_of_just_ids = []

            mongo_dao_retriever = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            list_of_gene_objets = mongo_dao_retriever.get_all_gene_objects_as_list()

            if list_of_gene_objets is not None:
                for gene_object in list_of_gene_objets:
                    gene_id = gene_object.gene_id
                    list_of_just_ids.append(gene_id)

            return list_of_just_ids

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list (unfiltered): ' + repr(error))

    def get_list_of_ids_from_mongo_filtered(self) -> list:
        list_of_just_ids = None

        try:
            list_of_just_ids = []

            mongo_dao_retriever_filtered = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            list_of_gene_objets = mongo_dao_retriever_filtered.get_all_gene_objects_as_list()

            if list_of_gene_objets is not None:
                for gene_object in list_of_gene_objets:
                    gene_id = gene_object.gene_id
                    list_of_just_ids.append(gene_id)

            return list_of_just_ids

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list (filtered): ' + repr(error))

    def get_list_of_ids_from_mongo_unfiltered_without_sequence(self) -> list:
        list_of_just_ids = None

        try:
            list_of_just_ids = []

            mongo_dao_retriever = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            list_of_gene_objects = mongo_dao_retriever.get_all_gene_objects_as_list()

            if list_of_gene_objects is not None:

                for gene_object in list_of_gene_objects:
                    if gene_object.sequence is None:
                        gene_id = gene_object.gene_id
                        list_of_just_ids.append(gene_id)

            return list_of_just_ids

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list without sequence (unfiltered): ' + repr(
                error))

    def get_list_of_ids_from_mongo_unfiltered_with_sequence(self) -> list:
        list_of_just_ids_of_genes_with_sequence = None

        try:
            list_of_just_ids_of_genes_with_sequence = []

            mongo_dao_retriever = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            list_of_gene_objets = mongo_dao_retriever.get_all_gene_objects_as_list()

            if list_of_gene_objets is not None:

                for gene_object in list_of_gene_objets:
                    if gene_object.sequence is not None:
                        gene_id = gene_object.gene_id
                        list_of_just_ids_of_genes_with_sequence.append(gene_id)

            return list_of_just_ids_of_genes_with_sequence

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list with sequence (unfiltered): ' + repr(
                error))

    def insert_in_unfiltered_collection_from_excel(self, file_path: str, sheet="0", column_name="gene_id") -> None:

        try:
            if os.path.isfile(file_path) is True:
                file_retriever_and_mongo_manager = GeneHandlerDAO(
                    MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                    Constants.MONGODB_DB_NAME,
                    Constants.MONGODB_COLLECTION_UNFILTERED,
                    file_path)

                list_of_genes_just_ids = file_retriever_and_mongo_manager.get_list_of_genes_from_xlrd(sheet,
                                                                                                      column_name)

                for ncbi_object_just_id in list_of_genes_just_ids:

                    criteria = GeneSearcher()
                    criteria.search_by_gene_id_criteria = ncbi_object_just_id.gene_id

                    if file_retriever_and_mongo_manager.search_gene_objects_and_return_as_list(criteria) is None:
                        file_retriever_and_mongo_manager.insert_gene_document_from_object(ncbi_object_just_id)
            else:
                raise OSError('The file does not exists')

        except Exception as error:
            print('Caught exception when inserting all data from excel: ' + repr(error))

    def insert_in_unfiltered_collection_from_fasta(self, file_path: str) -> None:

        try:
            if os.path.isfile(file_path) is True:
                file_retriever_and_mongo_manager = GeneHandlerDAO(
                    MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                    Constants.MONGODB_DB_NAME,
                    Constants.MONGODB_COLLECTION_UNFILTERED,
                    file_path)

                list_of_gene_objects_fasta = file_retriever_and_mongo_manager.get_list_of_gene_objects_from_multi_fasta()

                for gene_object in list_of_gene_objects_fasta:

                    criteria = GeneSearcher()
                    criteria.search_by_gene_id_criteria = gene_object.gene_id

                    if file_retriever_and_mongo_manager.search_gene_objects_and_return_as_list(criteria) is None:
                        file_retriever_and_mongo_manager.insert_gene_document_from_object(gene_object)
            else:
                raise OSError('The file does not exists')

        except Exception as error:
            print('Caught exception when inserting all data from fasta: ' + repr(error))

    def insert_in_unfiltered_collection_from_list_of_ids(self, list_of_gene_ids: list) -> None:

        try:
            mongo_dao_retriever = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME,
                Constants.MONGODB_COLLECTION_UNFILTERED)

            for single_id in list_of_gene_ids:

                criteria = GeneSearcher()
                criteria.search_by_gene_id_criteria = single_id

                gene = GeneDTO()
                gene.gene_id = single_id

                if mongo_dao_retriever.search_gene_objects_and_return_as_list(criteria) is None:
                    mongo_dao_retriever.insert_gene_document_from_object(gene)

        except Exception as error:
            print('Caught exception when inserting all data from list of ids: ' + repr(error))

    def insert_filtered_dict_in_filtered_collection(self, filtered_dict: dict) -> None:

        try:
            mongo_dao_manager = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            mongo_dao_manager.insert_gene_document_from_non_object_dict(filtered_dict)

        except Exception as error:
            print('Caught exception while inserting filtered sequences in mongo collection' + repr(error))

    def delete_gene_document_by_gene_id(self, gene_id: str) -> None:
        try:
            mongo_dao_retriever = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME,
                Constants.MONGODB_COLLECTION_UNFILTERED)

            criteria = GeneSearcher()
            criteria.search_by_gene_id_criteria = gene_id

            if mongo_dao_retriever.search_gene_objects_and_return_as_list(criteria) is not None:
                mongo_dao_retriever.delete_gene_document_by_gene_id(gene_id)

        except Exception as error:
            print('Caught exception at delete operation: ' + repr(error))

    def export_unfiltered_genes_collection_to_fasta(self, fasta_name: str) -> None:

        try:
            mongo_dao_retriever = GeneDAO(MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                                          Constants.MONGODB_DB_NAME,
                                          Constants.MONGODB_COLLECTION_UNFILTERED)

            list_of_seqrecords = mongo_dao_retriever.get_list_of_seqrecords_from_collection()

            if list_of_seqrecords is None:
                print(Constants.MSG_WARNING_UNFILTERED_COLLECTION_EMPTY)
                list_of_seqrecords = []

            SeqIO.write(list_of_seqrecords, fasta_name + "." + Constants.FASTA_EXTENSION, "fasta")

        except Exception as error:
            print('Caught exception when exporting all data to fasta from unfiltered collection: ' + repr(error))

    def export_filtered_genes_collection_to_fasta(self, fasta_name: str) -> None:

        try:
            mongo_dao_retriever = GeneDAO(MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                                          Constants.MONGODB_DB_NAME,
                                          Constants.MONGODB_COLLECTION_FILTERED)

            list_of_seqrecords = mongo_dao_retriever.get_list_of_seqrecords_from_collection()

            if list_of_seqrecords is None:
                print(Constants.MSG_WARNING_FILTERED_COLLECTION_EMPTY)
                list_of_seqrecords = []

            SeqIO.write(list_of_seqrecords, fasta_name + "." + Constants.FASTA_EXTENSION, "fasta")

        except Exception as error:
            print('Caught exception when exporting all data to fasta from filtered collection: ' + repr(error))

    def export_unfiltered_genes_collection_to_file_with_just_ids(self, file_name: str) -> None:

        try:
            list_of_unfiltered_genes_ids = self.get_list_of_ids_from_mongo_unfiltered()

            if list_of_unfiltered_genes_ids is None:
                print(Constants.MSG_WARNING_UNFILTERED_COLLECTION_EMPTY)
                list_of_unfiltered_genes_ids = []

            file_with_ids = open(file_name + '.' + Constants.ID_EXTENSION, 'w')

            file_with_ids.write(
                Constants.MSG_FILE_NUMBER_OF_ELEMENTS + " " + str(len(list_of_unfiltered_genes_ids)) + "\n")

            for id in list_of_unfiltered_genes_ids:
                file_with_ids.write("%s\n" % id)

        except Exception as error:
            print('Caught exception when exporting all ids to file from unfiltered collection: ' + repr(error))

    def export_filtered_genes_collection_to_file_with_just_ids(self, file_name: str) -> None:

        try:
            list_of_filtered_genes_ids = self.get_list_of_ids_from_mongo_filtered()

            if list_of_filtered_genes_ids is None:
                print(Constants.MSG_WARNING_FILTERED_COLLECTION_EMPTY)
                list_of_filtered_genes_ids = []

            file_with_ids = open(file_name + '.' + Constants.ID_EXTENSION, 'w')

            file_with_ids.write(
                Constants.MSG_FILE_NUMBER_OF_ELEMENTS + " " + str(len(list_of_filtered_genes_ids)) + "\n")

            for id_gene in list_of_filtered_genes_ids:
                file_with_ids.write("%s\n" % id_gene)

        except Exception as error:
            print('Caught exception when exporting all ids to file from filtered collection: ' + repr(error))

    def update_genes_from_dict(self, dict_of_genes: dict) -> None:

        try:
            if dict_of_genes is not None:
                mongo_dao_retriever = GeneDAO(
                    MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                    Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

                for id_ncbi, sequence in dict_of_genes.items():
                    ncbi_object_to_update = GeneDTO()
                    ncbi_object_to_update.gene_id = id_ncbi
                    ncbi_object_to_update.sequence = sequence
                    mongo_dao_retriever.update_gene_element_from_object(ncbi_object_to_update, True)

        except Exception as error:
            print('Caught exception when getting all ids from mongo as list: ' + repr(error))

    def delete_unfiltered_collection(self) -> None:

        try:
            mongo_dao_retriever_unfiltered = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_UNFILTERED)

            mongo_dao_retriever_unfiltered.delete_collection()

        except Exception as error:
            print('Caught exception when removing unfiltered collection' + repr(error))

    def delete_filtered_collection(self) -> None:

        try:
            mongo_dao_retriever_filtered = GeneDAO(
                MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                Constants.MONGODB_DB_NAME, Constants.MONGODB_COLLECTION_FILTERED)

            mongo_dao_retriever_filtered.delete_collection()

        except Exception as error:
            print('Caught exception when removing filtered collection' + repr(error))

    def check_gene_id_list_existance_on_unfiltered_from_list(self, list_of_genes: list) -> list:

        try:
            list_of_unfiltered_genes_from_mongo = self.get_list_of_ids_from_mongo_unfiltered()

            if list_of_unfiltered_genes_from_mongo is None:
                list_of_unfiltered_genes_from_mongo = []

            list_of_genes_not_in_unfiltered = []

            for gene in list_of_genes:
                if gene not in list_of_unfiltered_genes_from_mongo:
                    list_of_genes_not_in_unfiltered.append(gene)

            return list_of_genes_not_in_unfiltered

        except Exception as error:
            print('Caught exception at checking existence of genes from unfiltered collection ' + repr(error))

    def check_gene_id_list_existance_on_ncbi_from_list(self, list_of_genes: list) -> list:
        pbar = ProgressBar()
        Entrez.email = Constants.MAIL_SENDER
        list_of_available_genes = []
        print("Checking if all genes are available")
        for gene_id in pbar(list_of_genes):

            try:
                handle1 = Entrez.efetch(db="gene", id=gene_id, retmode="xml", validate=False)
                handle1.close()

                list_of_available_genes.append(gene_id)

            except Exception as error:
                print('\nGene ' +
                      gene_id + ' does not exists on NCBI : ' +
                      repr(error))
                continue

        return list_of_available_genes

    def download_sequences_from_list_as_dict_from_NCBI(self, list_of_genes: list) -> dict:
        pbar = ProgressBar()
        Entrez.email = Constants.MAIL_SENDER
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
                      gene_id + " : " + repr(error) + " This sequence will be removed from results")

                continue

        return dict_id_and_sequences

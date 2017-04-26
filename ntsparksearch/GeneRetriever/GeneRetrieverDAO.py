import xlrd
from ntsparksearch.Common.NucleotidesFromNCBIDTO import NucleotidesFromNCBIDTO
from ntsparksearch.Common.NCBItoMongoDAO import NCBItoMongoDAO
from pymongo import MongoClient
from Bio.SeqIO import parse


class GeneRetrieverDAO(NCBItoMongoDAO):
    def __init__(self, client_reference: MongoClient, database_name: str, collection_name: str, file_path: str):
        super(GeneRetrieverDAO, self).__init__(client_reference, database_name, collection_name)
        self._file_path = file_path

    def get_list_of_genes_from_xlrd(self, sheet: str, column_name: str) -> list:

        try:
            excel = xlrd.open_workbook(self._file_path)
            sh = excel.sheet_by_index(int(sheet))
            rows = []
            list_ncbi_records = []

            col = int()
            for cx in range(sh.ncols):
                if sh.cell_value(rowx=0, colx=cx) == str(column_name):
                    col = cx

            for rx in range(sh.nrows):
                if sh.cell_value(rx, colx=col) != str(column_name):
                    valor = sh.cell_value(rx, colx=col)
                    rows.append(str(valor))

            rows_no_repeated = tuple(set(rows))

            for id_ncbi in rows_no_repeated:
                gene_id = str(id_ncbi)
                gene_id = gene_id[:-2]
                ncbi_record_only_id = NucleotidesFromNCBIDTO()
                ncbi_record_only_id.idNcbi = gene_id
                list_ncbi_records.append(ncbi_record_only_id)

            return list_ncbi_records

        except Exception as error:
            print('Caught exception while reading from excel file: ' + repr(error))

    def get_list_of_ncbi_objects_from_multi_fasta(self) -> list:

        try:
            list_ncbi_records = []

            for record in parse(self._file_path, "fasta"):
                gene_id = str(record.id)
                gene_sequence = str(record.seq)

                ncbi_single_record = NucleotidesFromNCBIDTO()
                ncbi_single_record.idNcbi = gene_id
                ncbi_single_record.sequence = gene_sequence

                list_ncbi_records.append(ncbi_single_record)

            return list_ncbi_records

        except Exception as error:
            print('Caught exception while reading from fasta file: ' + repr(error))

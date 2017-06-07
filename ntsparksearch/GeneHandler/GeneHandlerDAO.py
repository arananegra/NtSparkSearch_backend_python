import xlrd
from Bio.SeqIO import parse
from pymongo import MongoClient

from ntsparksearch.Common.GeneDAO import GeneDAO
from ntsparksearch.Common.GeneDTO import GeneDTO


class GeneHandlerDAO(GeneDAO):
    def __init__(self, client_reference: MongoClient, database_name: str, collection_name: str, file_path: str):
        super(GeneHandlerDAO, self).__init__(client_reference, database_name, collection_name)
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
                ncbi_record_only_id = GeneDTO()
                ncbi_record_only_id.gene_id = gene_id
                list_ncbi_records.append(ncbi_record_only_id)

            return list_ncbi_records

        except Exception as error:
            print('Caught exception while reading from excel file: ' + repr(error))

    def get_list_of_gene_objects_from_multi_fasta(self) -> list:
        try:
            list_gene_records = []

            for record in parse(self._file_path, "fasta"):
                gene_id = str(record.id)
                gene_sequence = str(record.seq)

                gene_single_record = GeneDTO()
                gene_single_record.gene_id = gene_id
                gene_single_record.sequence = gene_sequence

                list_gene_records.append(gene_single_record)

            return list_gene_records

        except Exception as error:
            print('Caught exception while reading from fasta file: ' + repr(error))

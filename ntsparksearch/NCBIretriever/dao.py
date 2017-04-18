import xlrd
from ntsparksearch.common.domain import NucleotidesFromNCBI
from ntsparksearch.common.dao import NCBItoMongoDAO
from pymongo import MongoClient


class FileRetriverDAO(NCBItoMongoDAO):
    def __init__(self, client_reference: MongoClient, database_name: str, collection_name: str, file_path: str):
        super(FileRetriverDAO, self).__init__(client_reference, database_name, collection_name)
        self._file_path = file_path

    def get_list_of_genes_from_xlrd(self, sheet: int, column_name: str) -> list:

        try:
            excel = xlrd.open_workbook(self._file_path)
            sh = excel.sheet_by_index(sheet)
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
                ncbi_record_only_id = NucleotidesFromNCBI()
                ncbi_record_only_id.idNcbi = gene_id
                list_ncbi_records.append(ncbi_record_only_id)

            return list_ncbi_records

        except Exception as error:
            print('Caught exception while reading from excel file: ' + repr(error))


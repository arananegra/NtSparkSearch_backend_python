import xlrd
from bioSpark.common.domain import NucleotidesFromNCBI
from bioSpark.common.dao import NCBItoMongoDAO
from pymongo import MongoClient


class FileRetriverDAO(NCBItoMongoDAO):
    def __init__(self, client_reference: MongoClient, database_name: str, collection_name: str, file_path: str):
        super(FileRetriverDAO, self).__init__(client_reference, database_name, collection_name)
        self._file_path = file_path

    def obtain_list_of_genes_from_xlrd(self, sheet: str, column_name: str):

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

            for id_ncbi in rows:
                ncbi_record_only_id = NucleotidesFromNCBI()
                ncbi_record_only_id.idNcbi = id_ncbi
                list_ncbi_records.append(ncbi_record_only_id)

            return list_ncbi_records

        except Exception as error:
            print('Caught exception while reading from excel file: ' + repr(error))

    # crear supers o no??
    def insert_ncbi_documents_without_sequence_from_list(self, list_of_just_ids: list):
        super().insert_ncbi_document_from_list_of_objects(list_of_just_ids)

#    def update_record_to_add_sequence(self, id_to_updte):


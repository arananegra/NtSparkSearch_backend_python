import xlrd
from bioSpark.common.domain import NucleotidesFromNCBI


class FileRetriverDAO(object):
    def __init__(self):
        self._file_path = None

    def _set_file_path(self, path) -> None:
        self._file_path = path

    file_path = property(fset=_set_file_path)

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
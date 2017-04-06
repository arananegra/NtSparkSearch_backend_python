import xlrd


class FileRetriverDAO(object):
    def __init__(self, file_path: str):
        self._file_path = file_path

    def obtain_list_of_genes_from_xlrd(self, sheet: str, column_name: str):

        try:
            excel = xlrd.open_workbook(self._file_path)
            sh = excel.sheet_by_index(sheet)
            rows = []
            col = int()
            for cx in range(sh.ncols):
                if sh.cell_value(rowx=0, colx=cx) == str(column_name):
                    col = cx

            for rx in range(sh.nrows):
                if sh.cell_value(rx, colx=col) != str(column_name):
                    valor = sh.cell_value(rx, colx=col)
                    rows.append(str(valor))
            return rows

        except Exception as error:
            print('Caught exception at reading from excel file: ' + repr(error))

import argparse

from ntsparksearch.common.util import Constants
from ntsparksearch.common.dao import NCBItoMongoDAO
from ntsparksearch.common.domain import NucleotidesFromNCBI
from ntsparksearch.NCBIretriever.bs import NCBIretrieverBS


class App(object):
    """App class containing the main method"""

    @staticmethod
    def main():

        try:

            parser = argparse.ArgumentParser()

            parser.add_argument(Constants.COMMAND_DOWNLOAD_FROM_EXCEL,
                                metavar=(Constants.ARG_EXCEL_FILE_PATH,
                                         Constants.ARG_EXCEL_SHEET_NUMBER,
                                         Constants.ARG_EXCEL_COLUMN_NAME),
                                nargs=3, type=str,
                                help=Constants.HELP_COMMAND_DOWNLOAD_FROM_EXCEL)

            parser.add_argument(Constants.COMMAND_OBTAIN_ALL_IDS_FROM_UNFILTERED,
                                action='store_true',
                                help=Constants.HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED)

            args = parser.parse_args()

            if args.downloadGenesFromExcel:
                retrieverBS = NCBIretrieverBS()
                retrieverBS.insert_in_collection_from_excel(args.downloadGenesFromExcel[0],
                                                            args.downloadGenesFromExcel[1],
                                                            args.downloadGenesFromExcel[2])

                list_of_genes_empty = retrieverBS.obtain_list_of_ids_from_mongo()
                dict_of_genes_complete = retrieverBS.download_sequences_from_list_as_dict(list_of_genes_empty)
                retrieverBS.update_genes_from_dict(dict_of_genes_complete)
                print("Operation finished")


            if args.obtainUnfiltered:
                retrieverBS = NCBIretrieverBS()
                list_of_genes = retrieverBS.obtain_list_of_ids_from_mongo()
                print(list_of_genes)

        except (ValueError, OSError) as err:
            print(Constants.MSG_ERROR_INPUT, err)
        except Exception as err:
            print(Constants.MSG_ERROR_UNEXPECTED, err)


if __name__ == '__main__':
    App.main()

import argparse

from bioSpark.common.util import Constants
from bioSpark.common.dao import NCBItoMongoDAO
from bioSpark.common.domain import NucleotidesFromNCBI
from bioSpark.NCBIretriever.bs import NCBIretrieverBS


class App(object):
    """App class containing the main method"""

    @staticmethod
    def main():

        try:

            parser = argparse.ArgumentParser()

            parser.add_argument(Constants.COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED,
                                nargs=0, type=list,
                                help=Constants.HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED)

        # metavar=(Constants.ARG_EXCEL_FILE_PATH,
        #                                  Constants.ARG_EXCEL_SHEET_NUMBER,
        #                                  Constants.ARG_EXCEL_COLUMN_NAME)

            args = parser.parse_args()
            if args.obtainUnfiltered:
                test_bs = NCBIretrieverBS()
                list_of_genes = test_bs.obtain_list_of_ids_from_mongo()
                print(list_of_genes)

        except (ValueError, OSError) as err:
            print(Constants.MSG_ERROR_INPUT, err)
        except Exception as err:
            print(Constants.MSG_ERROR_UNEXPECTED, err)






if __name__ == '__main__':
    App.main()



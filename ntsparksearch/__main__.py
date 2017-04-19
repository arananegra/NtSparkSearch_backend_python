import argparse
from ntsparksearch.common.util import Constants
from ntsparksearch.subsequenceMatcher.bs import SubSequenceSparkMatcherBS
from ntsparksearch.NCBIretriever.bs import NCBIretrieverBS



class App(object):
    """App class containing the main method"""

    @staticmethod
    def main():

        try:

            parser = argparse.ArgumentParser()

            parser.add_argument(Constants.COMMAND_OBTAIN_ALL_IDS_FROM_UNFILTERED,
                                action='store_true',
                                help=Constants.HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED)

            parser.add_argument(Constants.COMMAND_DOWNLOAD_FROM_EXCEL,
                                metavar=(Constants.ARG_EXCEL_FILE_PATH,
                                         Constants.ARG_EXCEL_SHEET_NUMBER,
                                         Constants.ARG_EXCEL_COLUMN_NAME),
                                nargs=3, type=str,
                                help=Constants.HELP_COMMAND_DOWNLOAD_FROM_EXCEL)

            parser.add_argument(Constants.COMMAND_EXACT_SUB_MATCH_SPARK,
                                metavar=Constants.ARG_EXCEL_SEQUENCE_TO_FETCH,
                                nargs=1, type=str,
                                help=Constants.HELP_COMMAND_EXACT_SUB_MATCH_SPARK)

            args = parser.parse_args()

            if args.obtainUnfiltered:
                retrieverBS = NCBIretrieverBS()
                list_of_genes = retrieverBS.obtain_list_of_ids_from_mongo()
                print(list_of_genes)

            if args.downloadGenesFromExcel:
                retrieverBS = NCBIretrieverBS()
                retrieverBS.insert_in_collection_from_excel(args.downloadGenesFromExcel[0],
                                                            args.downloadGenesFromExcel[1],
                                                            args.downloadGenesFromExcel[2])

                list_of_genes_empty = retrieverBS.obtain_list_of_ids_from_mongo()
                dict_of_genes_complete = retrieverBS.download_sequences_from_list_as_dict(list_of_genes_empty)
                retrieverBS.update_genes_from_dict(dict_of_genes_complete)
                print("Operation finished")

                # TODO: cuando se realiza la descarga de un nuevo excel, se vuelve a hacer la descarga completa
                # a partir de la base de datos, en la proxima iteracion, descargar solo los elementos nuevos
                # a no ser que se indique lo contrario

            if args.sparkseqmatch:
                subsequence_matcherBS = SubSequenceSparkMatcherBS()

                dict_filtered_with_spark = subsequence_matcherBS. \
                    filter_sequences_by_sequence_string_to_dict(args.sparkseqmatch[0])

                subsequence_matcherBS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)

        except (ValueError, OSError) as err:
            print(Constants.MSG_ERROR_INPUT, err)
        except Exception as err:
            print(Constants.MSG_ERROR_UNEXPECTED, err)


if __name__ == '__main__':
    App.main()

import argparse
from ntsparksearch.Common.Constants import Constants
from ntsparksearch.SubsequenceMatcher.SubSequenceSparkMatcherBS import SubSequenceSparkMatcherBS
from ntsparksearch.GeneRetriever.GeneRetrieverBS import GeneRetrieverBS


class App(object):
    """App class containing the main method"""

    @staticmethod
    def main():

        try:

            parser = argparse.ArgumentParser()

            parser.add_argument(Constants.COMMAND_OBTAIN_ALL_IDS_FROM_UNFILTERED,
                                action='store_true',
                                help=Constants.HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED)

            parser.add_argument(Constants.COMMAND_OBTAIN_ALL_IDS_FROM_FILTERED,
                                action='store_true',
                                help=Constants.HELP_COMMAND_OBTAIN_ALL_SEQUENCES_FILTERED)

            parser.add_argument(Constants.COMMAND_REMOVE_UNFILTERED_COLLECTION,
                                action='store_true',
                                help=Constants.HELP_COMMAND_REMOVE_UNFILTERED)

            parser.add_argument(Constants.COMMAND_REMOVE_FILTERED_COLLECTION,
                                action='store_true',
                                help=Constants.HELP_COMMAND_REMOVE_FILTERED)

            parser.add_argument(Constants.COMMAND_DOWNLOAD_FROM_EXCEL,
                                metavar=(Constants.ARG_EXCEL_FILE_PATH,
                                         Constants.ARG_EXCEL_SHEET_NUMBER,
                                         Constants.ARG_EXCEL_COLUMN_NAME),
                                nargs=3, type=str,
                                help=Constants.HELP_COMMAND_DOWNLOAD_FROM_EXCEL)

            parser.add_argument(Constants.COMMAND_IMPORT_FROM_FASTA,
                                metavar=Constants.ARG_FASTA_FILE_PATH,
                                nargs=1, type=str,
                                help=Constants.HELP_COMMAND_IMPORT_FROM_FASTA)

            parser.add_argument(Constants.COMMAND_EXACT_SUB_MATCH_SPARK,
                                metavar=(Constants.ARG_SEQUENCE_TO_FETCH,
                                         Constants.ARG_REMOVE_PREVIOUS_RESULT),
                                nargs=2, type=str,
                                help=Constants.HELP_COMMAND_EXACT_SUB_MATCH_SPARK)

            args = parser.parse_args()

            if args.obtainUnfiltered:
                retriever_BS = GeneRetrieverBS()
                list_of_genes = retriever_BS.get_list_of_ids_from_mongo()

                if list_of_genes is None:
                    print("The unfiltered collection of genes is empty")

                else:
                    print(list_of_genes)

            if args.obtainFiltered:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()
                list_of_genes_filtered = subsequence_matcher_BS.get_list_of_ids_from_mongo_filtered()

                if list_of_genes_filtered is None:
                    print("The filtered collection of genes is empty")

                else:
                    print(list_of_genes_filtered)

            if args.removeUnfiltered:
                retriever_BS = GeneRetrieverBS()
                retriever_BS.delete_unfiltered_collection()
                print("Operation finished")

            if args.removeFiltered:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()
                subsequence_matcher_BS.delete_filtered_collection()
                print("Operation finished")

            if args.downloadGenesFromExcel:
                retriever_BS = GeneRetrieverBS()
                retriever_BS.insert_in_collection_from_excel(args.downloadGenesFromExcel[0],
                                                             args.downloadGenesFromExcel[1],
                                                             args.downloadGenesFromExcel[2])

                list_of_genes_empty = retriever_BS.get_list_of_ids_from_mongo_without_sequence()

                if list_of_genes_empty is not None:
                    print("Downloading the content of the unfiltered collection of genes")
                    dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
                        list_of_genes_empty)
                    retriever_BS.update_genes_from_dict(dict_of_genes_complete)
                else:
                    print("WARNING: The unfiltered collection is empty")
                print("Operation finished")

            if args.retrieveFromFasta:
                retriever_BS = GeneRetrieverBS()
                retriever_BS.insert_in_collection_from_fasta(args.retrieveFromFasta[0])
                print("Operation finished")

            if args.sparkSeqMatch:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()

                dict_filtered_with_spark = subsequence_matcher_BS. \
                    filter_sequences_by_sequence_string_to_dict(args.sparkSeqMatch[0], args.sparkSeqMatch[1])

                subsequence_matcher_BS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)
                print("Operation finished")

        except (ValueError, OSError) as err:
            print(Constants.MSG_ERROR_INPUT, err)
        except Exception as err:
            print(Constants.MSG_ERROR_UNEXPECTED, err)


if __name__ == '__main__':
    App.main()

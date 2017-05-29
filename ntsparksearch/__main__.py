import argparse
from ntsparksearch.Common.Constants import Constants
from ntsparksearch.SubsequenceMatcher.SubSequenceSparkMatcherBS import SubSequenceSparkMatcherBS
from ntsparksearch.GeneHandler.GeneHandlerBS import GeneHandlerBS


class App(object):
    """App class containing the main method"""

    @staticmethod
    def main():

        try:

            parser = argparse.ArgumentParser()

            parser.add_argument(Constants.COMMAND_OBTAIN_ALL_IDS_FROM_UNFILTERED,
                                action='store_true',
                                help=Constants.HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED)

            parser.add_argument(Constants.COMMAND_EXPORT_FASTA_FROM_UNFILTERED,
                                metavar=Constants.ARG_FASTA_FILE_EXPORT,
                                nargs=1,
                                help=Constants.HELP_COMMAND_EXPORT_FASTA_FROM_UNFILTERED)

            parser.add_argument(Constants.COMMAND_EXPORT_TEXTFILE_ID_FROM_UNFILTERED,
                                metavar=Constants.ARG_TEXT_FILE_WITH_IDS_EXPORT,
                                nargs=1,
                                help=Constants.HELP_COMMAND_EXPORT_TEXTFILE_ID_FROM_UNFILTERED)

            parser.add_argument(Constants.COMMAND_OBTAIN_ALL_IDS_FROM_FILTERED,
                                action='store_true',
                                help=Constants.HELP_COMMAND_OBTAIN_ALL_SEQUENCES_FILTERED)

            parser.add_argument(Constants.COMMAND_EXPORT_FASTA_FROM_FILTERED,
                                metavar=Constants.ARG_FASTA_FILE_EXPORT,
                                nargs=1,
                                help=Constants.HELP_COMMAND_EXPORT_FASTA_FROM_FILTERED)

            parser.add_argument(Constants.COMMAND_EXPORT_TEXTFILE_ID_FROM_FILTERED,
                                metavar=Constants.ARG_TEXT_FILE_WITH_IDS_EXPORT,
                                nargs=1,
                                help=Constants.HELP_COMMAND_EXPORT_TEXTFILE_ID_FROM_FILTERED)

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
                                metavar=Constants.ARG_FASTA_FILE_IMPORT,
                                nargs=1, type=str,
                                help=Constants.HELP_COMMAND_IMPORT_FROM_FASTA)

            parser.add_argument(Constants.COMMAND_EXACT_SUB_MATCH_SPARK,
                                metavar=(Constants.ARG_SEQUENCE_TO_FETCH,
                                         Constants.ARG_REMOVE_PREVIOUS_RESULT),
                                nargs=2, type=str,
                                help=Constants.HELP_COMMAND_EXACT_SUB_MATCH_SPARK)

            args = parser.parse_args()

            if args.obtainUnfiltered:
                retriever_BS = GeneHandlerBS()
                list_of_genes = retriever_BS.get_list_of_ids_from_mongo()

                if list_of_genes is None:
                    print(Constants.MSG_WARNING_UNFILTERED_COLLECTION_EMPTY)

                else:
                    print(list_of_genes)

            if args.exportUnfilteredFasta:
                retriever_BS = GeneHandlerBS()
                retriever_BS.export_unfiltered_genes_collection_to_fasta(
                    Constants.OUTPUT_FOLDER + args.exportUnfilteredFasta[0])
                print(Constants.MSG_PROCESS_FINISHED)

            if args.exportUnfilteredId:
                retriever_BS = GeneHandlerBS()
                retriever_BS.export_unfiltered_genes_collection_to_file_with_just_ids(
                    Constants.OUTPUT_FOLDER + args.exportUnfilteredId[0])
                print(Constants.MSG_PROCESS_FINISHED)

            if args.obtainFiltered:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()
                list_of_genes_filtered = subsequence_matcher_BS.get_list_of_ids_from_mongo_filtered()

                if list_of_genes_filtered is None:
                    print(Constants.MSG_WARNING_FILTERED_COLLECTION_EMPTY)

                else:
                    print(list_of_genes_filtered)

            if args.exportFilteredFasta:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()
                subsequence_matcher_BS.export_filtered_genes_collection_to_fasta(
                    Constants.OUTPUT_FOLDER + args.exportFilteredFasta[0])
                print(Constants.MSG_PROCESS_FINISHED)

            if args.exportFilteredId:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()
                subsequence_matcher_BS.export_filtered_genes_collection_to_file_with_just_ids(
                    Constants.OUTPUT_FOLDER + args.exportFilteredId[0])
                print(Constants.MSG_PROCESS_FINISHED)

            if args.removeUnfiltered:
                retriever_BS = GeneHandlerBS()
                retriever_BS.delete_unfiltered_collection()
                print(Constants.MSG_PROCESS_FINISHED)

            if args.removeFiltered:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()
                subsequence_matcher_BS.delete_filtered_collection()
                print(Constants.MSG_PROCESS_FINISHED)

            if args.downloadGenesFromExcel:
                retriever_BS = GeneHandlerBS()
                retriever_BS.insert_in_collection_from_excel(Constants.INPUT_FOLDER + args.downloadGenesFromExcel[0],
                                                             args.downloadGenesFromExcel[1],
                                                             args.downloadGenesFromExcel[2])

                list_of_genes_without_sequence = retriever_BS.get_list_of_ids_from_mongo_without_sequence()

                if list_of_genes_without_sequence is not None:
                    print(Constants.MSG_PROCESS_DOWNLOADING_GENES)
                    dict_of_genes_complete = retriever_BS.download_sequences_from_list_as_dict_from_NCBI(
                        list_of_genes_without_sequence)
                    retriever_BS.update_genes_from_dict(dict_of_genes_complete)
                else:
                    print(Constants.MSG_WARNING_UNFILTERED_COLLECTION_EMPTY)
                print(Constants.MSG_PROCESS_FINISHED)

            if args.retrieveFromFasta:
                retriever_BS = GeneHandlerBS()
                retriever_BS.insert_in_collection_from_fasta(Constants.INPUT_FOLDER + args.retrieveFromFasta[0])
                print(Constants.MSG_PROCESS_FINISHED)

            if args.sparkSeqMatch:
                subsequence_matcher_BS = SubSequenceSparkMatcherBS()

                dict_filtered_with_spark = subsequence_matcher_BS. \
                    filter_sequences_by_sequence_string_to_dict(args.sparkSeqMatch[0], args.sparkSeqMatch[1])

                subsequence_matcher_BS.insert_filtered_dict_in_filtered_collection(dict_filtered_with_spark)
                print(Constants.MSG_PROCESS_FINISHED)

        except (ValueError, OSError) as err:
            print(Constants.MSG_ERROR_INPUT, err)
        except Exception as err:
            print(Constants.MSG_ERROR_UNEXPECTED, err)


if __name__ == '__main__':
    App.main()

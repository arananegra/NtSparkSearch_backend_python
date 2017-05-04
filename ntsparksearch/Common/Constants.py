import os


class Constants(object):
    # Database parameters
    MONGODB_DB_NAME = 'bioSpark'
    # url in docker-compose
    #MONGODB_HOST = 'database'
    MONGODB_HOST = 'localhost'
    MONGODB_PORT = 27017
    MONGODB_COLLECTION_FILTERED = 'filteredGenes'
    MONGODB_COLLECTION_UNFILTERED = 'unfilteredGenes'

    # Collection field names
    GENE_ID = "_gene_id"
    SEQUENCE = "_sequence"

    # indexes
    GENE_ID_INDEX = "_gene_id_text"

    # file paths
    OUTPUT_FOLDER = "output/"
    INPUT_FOLDER = "input/"

    # Env variables
    SPARK_HOME = os.environ['SPARK_HOME']

    # Main args
    ARG_EXCEL_FILE_PATH = "<EXCEL file path>"
    ARG_EXCEL_SHEET_NUMBER = "<EXCEL sheet number>"
    ARG_EXCEL_COLUMN_NAME = "<EXCEL column name>"
    ARG_FASTA_FILE_EXPORT = "<FASTA file name to export>"
    ARG_TEXT_FILE_WITH_IDS_EXPORT = "<Text file name with ids to export>"
    ARG_FASTA_FILE_IMPORT = "<FASTA file name to import>"
    ARG_SEQUENCE_TO_FETCH = "<Nucleotide Sequence to fecth>"
    ARG_REMOVE_PREVIOUS_RESULT = "<Remove all data from previous executions (y/n)>"

    # Commands
    COMMAND_OBTAIN_ALL_IDS_FROM_UNFILTERED = "--obtainUnfiltered"
    COMMAND_EXPORT_FASTA_FROM_UNFILTERED = "--exportUnfilteredFasta"
    COMMAND_EXPORT_TEXTFILE_ID_FROM_UNFILTERED = "--exportUnfilteredId"
    COMMAND_OBTAIN_ALL_IDS_FROM_FILTERED = "--obtainFiltered"
    COMMAND_EXPORT_FASTA_FROM_FILTERED = "--exportFilteredFasta"
    COMMAND_EXPORT_TEXTFILE_ID_FROM_FILTERED = "--exportFilteredId"
    COMMAND_REMOVE_UNFILTERED_COLLECTION = "--removeUnfiltered"
    COMMAND_REMOVE_FILTERED_COLLECTION = "--removeFiltered"
    COMMAND_DOWNLOAD_FROM_EXCEL = "--downloadGenesFromExcel"
    COMMAND_IMPORT_FROM_FASTA = "--retrieveFromFasta"
    COMMAND_EXACT_SUB_MATCH_SPARK = "--sparkSeqMatch"

    # Process mssg
    MSG_PROCESS_FINISHED = "Operation finished"
    MSG_PROCESS_DOWNLOADING_GENES = "Downloading the content of the unfiltered collection of genes"

    # Error mssg
    MSG_ERROR_INPUT = "INPUT NOT VALID:"
    MSG_ERROR_UNEXPECTED = "UNEXPECTED ERROR:"

    # Warning mssg
    MSG_WARNING_UNFILTERED_COLLECTION_EMPTY = "WARNING: The unfiltered collection is empty"
    MSG_WARNING_FILTERED_COLLECTION_EMPTY = "WARNING: The filtered collection is empty"

    # Help
    HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED = "Obtain ids from unfiltered collection" \
                                                   " of sequences"

    HELP_COMMAND_EXPORT_FASTA_FROM_UNFILTERED = "Export a fasta file from the unfiltered collection"

    HELP_COMMAND_EXPORT_TEXTFILE_ID_FROM_UNFILTERED = "Export a text file from the unfiltered collection" \
                                                      " with just the ids"

    HELP_COMMAND_OBTAIN_ALL_SEQUENCES_FILTERED = "Obtaing ids from filtered collection" \
                                                 " of sequences"

    HELP_COMMAND_EXPORT_FASTA_FROM_FILTERED = "Export a fasta file from the filtered collection"

    HELP_COMMAND_EXPORT_TEXTFILE_ID_FROM_FILTERED = "Export a text file from the filtered collection" \
                                                      " with just the ids"

    HELP_COMMAND_REMOVE_UNFILTERED = "Removes the whole unfiltered collection of genes"

    HELP_COMMAND_REMOVE_FILTERED = "Removes the whole filtered collection of genes"

    HELP_COMMAND_DOWNLOAD_FROM_EXCEL = "Obtain ids from an excel file and downloads" \
                                       " all the sequences from the NCBI"

    HELP_COMMAND_IMPORT_FROM_FASTA = "Obtain ids and sequences from a FASTA file and import them" \
                                     " to the database"

    HELP_COMMAND_EXACT_SUB_MATCH_SPARK = "From the sequences loaded in the database," \
                                         " performs an exact subsequence search using" \
                                         " the provided sequence (ambiguity codes are" \
                                         " allowed)"

    UPLOAD_FOLDER = '/Users/alvarogomez/Desktop/uploads'
    ALLOWED_EXTENSIONS = {'txt','jpg', 'fasta'}
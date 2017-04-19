class Constants(object):
    AMBIGUOUS_DNA_VALUES = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "M": "AC",
        "R": "AG",
        "W": "AT",
        "S": "CG",
        "Y": "CT",
        "K": "GT",
        "V": "ACG",
        "H": "ACT",
        "D": "AGT",
        "B": "CGT",
        "X": "GATC",
        "N": "GATC",
    }

    # Database parameters
    MONGODB_DB_NAME = 'bioSpark'
    # url in docker-compose
    # MONGODB_HOST = 'database'
    MONGODB_HOST = 'localhost'
    MONGODB_PORT = 27017
    MONGODB_COLLECTION_FILTERED = 'ncbifiltered'
    MONGODB_COLLECTION_UNFILTERED = 'ncbiunfiltered'

    # Main args
    ARG_EXCEL_FILE_PATH = "<EXCEL file path>"
    ARG_EXCEL_SHEET_NUMBER = "<EXCEL sheet number>"
    ARG_EXCEL_COLUMN_NAME = "<EXCEL column name>"
    ARG_EXCEL_SEQUENCE_TO_FETCH = "Nucleotide Sequence to fecth"

    # Commands
    COMMAND_OBTAIN_ALL_IDS_FROM_UNFILTERED = "--obtainUnfiltered"
    COMMAND_OBTAIN_ALL_IDS_FROM_FILTERED = "--obtainFiltered"
    COMMAND_DOWNLOAD_FROM_EXCEL = "--downloadGenesFromExcel"
    COMMAND_EXACT_SUB_MATCH_SPARK = "--sparkseqmatch"

    # Error mssg
    MSG_ERROR_INPUT = "INPUT NOT VALID:"
    MSG_ERROR_UNEXPECTED = "UNEXPECTED ERROR:"

    # Help
    HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED = "Obtaing ids from unfiltered collection" \
                                                   "of sequences"

    HELP_COMMAND_OBTAIN_ALL_SEQUENCES_FILTERED = "Obtaing ids from filtered collection" \
                                                   "of sequences"
    HELP_COMMAND_DOWNLOAD_FROM_EXCEL = "Obtain ids from an excel file and downloads" \
                                       "all the sequences from the NCBI"

    HELP_COMMAND_EXACT_SUB_MATCH_SPARK = "From the sequences loaded in the database," \
                                         "performs a exact subsequence search using" \
                                         "the provided sequence (ambiguity codes are" \
                                         "allowed)"

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

    MONGODB_DB_NAME = 'bioSpark'
    # url in docker-compose
    #MONGODB_HOST = 'database'
    MONGODB_HOST = 'localhost'
    MONGODB_PORT = 27017
    MONGODB_COLLECTION_FILTERED = 'ncbifiltered'
    MONGODB_COLLECTION_UNFILTERED = 'ncbiunfiltered'

    # Main args
    ARG_EXCEL_FILE_PATH = "<EXCEL file path>"
    ARG_EXCEL_SHEET_NUMBER = "<EXCEL sheet number>"
    ARG_EXCEL_COLUMN_NAME = "<EXCEL column name>"
    # Commands
    COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED = "--obtainUnfiltered"
    COMMAND_DOWNLOAD_FROM_EXCEL = "--download-ncbi"

    # Error mssg

    MSG_ERROR_INPUT = "INPUT NOT VALID:"
    MSG_ERROR_UNEXPECTED = "UNEXPECTED ERROR:"
    MSG_ERROR_DNA_SEQUENCE = "The DNA sequence should only contain the symbols 'A', 'C', 'G', and 'T'"

    # Help

    HELP_COMMAND_OBTAIN_ALL_SEQUENCES_UNFILTERED = "Obtaing all info from unfiltered collection" \
                                                   "of sequences"

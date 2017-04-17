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
    MONGODB_HOST = 'database'
    #MONGODB_HOST = 'localhost'
    MONGODB_PORT = 27017
    # MONGODB_COLLECTION_UNFILTERED_NCBI = 'test01'
    MONGODB_COLLECTION_FILTERED = 'ncbifiltered'
    MONGODB_COLLECTION_UNFILTERED = 'ncbiunfiltered'

    REG_EXP_DNA_SEQUENCE = "[ACGT]+"

    # Error mssg

    MSG_ERROR_DNA_SEQUENCE = "The DNA sequence should only contain the symbols 'A', 'C', 'G', and 'T'"

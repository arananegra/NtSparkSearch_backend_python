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

    MONGODB_LOCALHOST = 'localhost'
    MONGODB_PORT = 27017

    REG_EXP_DNA_SEQUENCE = "[ACGT]+"

    # File open mode

    MODE_OPEN_FILE_W = 'w'
    MODE_OPEN_FILE_R = 'r'
    MODE_OPEN_FILE_R_W = 'r+'

    # Error mssg

    MSG_ERROR_DNA_SEQUENCE = "The DNA sequence should only contain the symbols 'A', 'C', 'G', and 'T'"

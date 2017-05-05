import os


class Constants(object):
    # Database parameters
    MONGODB_DB_NAME = 'bioSpark'
    # url in docker-compose
    # MONGODB_HOST = 'database'
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
    OUTPUT_FOLDER = "/Users/alvarogomez/Desktop/uploads/"
    INPUT_FOLDER = "input/"

    # Env variables
    SPARK_HOME = os.environ['SPARK_HOME']

    # Warning mssg
    MSG_WARNING_UNFILTERED_COLLECTION_EMPTY = "WARNING: The unfiltered collection is empty"
    MSG_WARNING_FILTERED_COLLECTION_EMPTY = "WARNING: The filtered collection is empty"

    # Services constants
    UPLOAD_FOLDER = OUTPUT_FOLDER
    EXCEL_EXTENSION = "xlsx"
    FASTA_EXTENSION = "fasta"

    # Service responses
    POST_CREATED = 201

    # Mail messages
    MAIL_SENDER = "arananegrayeye@gmail.com"

    MESSAGE_DOWNLOAD_INITIALIZE = """From: <from@fromdomain.com>
        To: To Person <to@todomain.com>
        Subject: SMTP e-mail test

        This is a test e-mail message.
        """

import os


class Constants(object):
    # Database parameters
    MONGODB_USERS_DB = 'ntSparkUsersMeta'
    MONGODB_DB_INITIAL = 'ntSparkGenes_'
    # url in docker-compose
    MONGODB_HOST = 'database'
    #MONGODB_HOST = 'localhost'
    MONGODB_PORT = 27017
    MONGODB_COLLECTION_FILTERED = 'filteredGenes'
    MONGODB_COLLECTION_UNFILTERED = 'unfilteredGenes'

    # Gene collection field names
    GENE_ID = "_gene_id"
    SEQUENCE = "_sequence"

    # Queue collection field names
    REDIS_SERVER = "redis-server"
    #REDIS_SERVER = "localhost"

    # indexes
    GENE_ID_INDEX = "_gene_id_text"

    # file paths
    OUTPUT_FOLDER = "/usr/spark-2.2.0/output/"
    #OUTPUT_FOLDER = "/Users/alvarogomez/Desktop/uploads/"
    INPUT_FOLDER = "usr/spark-2.2.0/input/"

    # Env variables
    SPARK_HOME = os.environ['SPARK_HOME']

    # Warning mssg
    MSG_WARNING_UNFILTERED_COLLECTION_EMPTY = "WARNING: The unfiltered collection is empty"
    MSG_WARNING_FILTERED_COLLECTION_EMPTY = "WARNING: The filtered collection is empty"

    # File mssg
    MSG_FILE_NUMBER_OF_ELEMENTS = "The number of item(s) in this file is (are)"

    # Services extensions
    UPLOAD_FOLDER = OUTPUT_FOLDER
    EXCEL_EXTENSION = "xlsx"
    FASTA_EXTENSION = "fasta"
    ID_EXTENSION = "txt"

    # Service responses
    OK = 200
    SERVER_ERROR = 500
    BAD_REQUEST = 400
    POST_CREATED = 201
    OK_WAIT = 202
    POST_BAD_REQUEST = 202

    # Service constants
    EMAIL_SERVICE_PARAMETER_NAME_CONSTANT = "email"
    GENE_ID_SERVICE_PARAMETER_NAME_CONSTANT = "gene_id"
    SEQUENCE_SERVICE_PARAMETER_NAME_CONSTANT = "sequence"

    # Mail messages
    MAIL_SENDER = "ntsparksearch@gmail.com"
    MAIL_HOST = 'smtp.gmail.com:587'
    MAIL_USER = "ntsparksearcher"
    MAIL_PASS = 'ntsparksearcher1234'

    MSG_DOWNLOAD_INITIALIZE = "\r\n".join([
        "Subject: Your download request has been accepted",
        "",
        "When the download process is finished you will receive another email. This is the list of genes required:\n\n"])

    MSG_DOWNLOAD_FINISHED = "\r\n".join([
        "Subject: Your download is already finished",
        "",
        "You can now check subsequences about these genes: \n\n"])

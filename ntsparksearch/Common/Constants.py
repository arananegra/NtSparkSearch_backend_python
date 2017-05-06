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
    MONGODB_COLLECTION_QUEUE = 'queueGenes'
    MONGODB_COLLECTION_ASCENDING_ORDER = 1


    # Gene collection field names
    GENE_ID = "_gene_id"
    SEQUENCE = "_sequence"

    # Queue collection field names
    QUEUE_ENTRY_DATE = "_entry_date_in_collection"
    QUEUE_LIST_OF_GENES = "_gene_id_list_to_download"

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
    POST_WAIT = 202

    # Mail messages
    MAIL_SENDER = "ntsparksearch@gmail.com"
    MAIL_HOST = 'smtp.gmail.com:587'
    MAIL_USER = "ntsparksearcher"
    MAIL_PASS = 'ntsparksearcher1234'

    MSG_DOWNLOAD_INITIALIZE = "\r\n".join([
        "Subject: Your download request has been accepted",
        "",
        "This is the list of genes required: When the download process is finished you will receive another email\n"])

    MSG_DOWNLOAD_FINISHED = "\r\n".join([
        "Subject: Your download is already finished",
        "",
        "This is the list of required genes that are already downloaded: \n"])

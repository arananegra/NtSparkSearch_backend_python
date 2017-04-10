from bioSpark.NCBIretriever.api import INCBIretriever
from bioSpark.NCBIretriever.dao import FileRetriverDAO
from bioSpark.common.domain import NCBIsearcher
from bioSpark.common.util import Constants
from pymongo import MongoClient


class NCBIretrieverBS(INCBIretriever):
    def __init__(self, file_path: str):
        self._file_path = file_path

    def insert_in_collection_from_excel(self, sheet: int, column_name: str):
        file_retriever_and_mongo_manager = FileRetriverDAO(MongoClient(Constants.MONGODB_HOST, Constants.MONGODB_PORT),
                                                           Constants.MONGODB_DB_NAME,
                                                           Constants.MONGODB_COLLECTION_UNFILTERED,
                                                           self._file_path)

        list_of_genes_just_ids = file_retriever_and_mongo_manager.obtain_list_of_genes_from_xlrd(sheet, column_name)

        for ncbi_object_just_id in list_of_genes_just_ids:

            criteria = NCBIsearcher()
            criteria.search_by_NCBI_id_criteria = ncbi_object_just_id.idNcbi

            if file_retriever_and_mongo_manager.search_ncbi_objects_and_return_as_list(criteria) == None:
                file_retriever_and_mongo_manager.insert_ncbi_document_from_object(ncbi_object_just_id)

    def download_sequences_from_list_to_dict(self, list_of_genes: list):
        return None

#
# test_from_excel = NCBIretrieverBS("/Users/alvarogomez/DEG.xlsx")
#
# test_from_excel.insert_in_collection_from_excel(0, "gene_id")


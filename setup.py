from setuptools import setup

setup(name='Nucleotide Spark Exact Searcher',
      version='1.0.0',
      description='This application is the final degree project of Álvaro Gómez Jáuregui. It uses Spark to'
                  'perform exact searching through a set of genes imported from excel or fasta files. This program also'
                  'performs the download of the sequences from the NCBI and management of a mongodb database'
                  'to store the data',
      author='Álvaro Gómez Jáuregui',
      author_email='alvarogj@lcc.uma.es',
      maintainer='Álvaro Gómez Jáuregui',
      maintainer_email='alvarogj@lcc.uma.es',
      license='GNU GENERAL PUBLIC LICENSE',
      packages=['ntsparksearch', 'ntsparksearch.GeneHandler', 'ntsparksearch.SubsequenceMatcher',
                'ntsparksearch.Common', 'ntsparksearch.RestApi', 'ntsparksearch.EmailProcess'],
      zip_safe=False
      )

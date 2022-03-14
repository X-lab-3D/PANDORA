Database Generation
===================


.. automodule:: PANDORA.Database

This module contains all the tools to compute the template structure and reference sequence databases. The main class used for the data generation is ``PANDORA.Database.Database``.
Through this class you can specify the type of structures you want to download (MHCI and/or II).

Basic building Example:

>>> from PANDORA.Database import Database
>>>
>>> db = Database.Database()
>>> db.construct_database(save='./database/full_db.pkl', download=True)
>>>

A database can be easily loaded and edited.
Custom structures (i.e. not coming from IMGT) can be also added to the database by using ``Database.add_structure()`` function.
To do so, you will need to read your structure into a Biopython structure object:

>>> from PANDORA.Database import Database
>>>
>>> #Load the database
>>> db = Database.load('./database/full_db.pkl')
>>>
>>> db.add_structure(id='0000', allele_type=['HLA-A*02:01'],
                      peptide = 'AAALLLAAA', MHC_class = 'I',
                      chain_seq = [], anchors = [],
                      pdb_path = './PDBs/pMHCI/0000.pdb', pdb = False)

Database
----------------------------------------

.. automodule:: PANDORA.Database.Database
    :members:
    :undoc-members:
    :show-inheritance:

Database_functions
----------------------------------------

.. automodule:: PANDORA.Database.Database_functions
    :noindex:
    :members:
    :undoc-members:

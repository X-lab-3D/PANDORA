Wrapper Module
===============


.. automodule:: PANDORA.Wrapper

This module is used to run PANDORA for large number of cases.
It allows the user to insert all the cases in a csv/tsv file, add various information to them (template, anchor positions, M chain seq).

Wrapper Example:

>>> from PANDORA.Pandora import Pandora
>>> from PANDORA.Database import Database
>>> from PANDORA.Wrapper import Wrapper
>>>
>>> ## A. Load pregenerated database of all pMHC PDBs as templates
>>> db = Database.load('pandora_Database')
>>>
>>> ## B. Create the wrapper object
>>> wrap =  Wrapper()
>>>
>>> ## C. Create all Target Objects based on peptides in the .tsv file
>>> wrap.create_targets('datafile.tsv', IDs_col=0,
>>>                      peptides_col=1, allele_col=2, db)
>>>
>>> ## C. Perform modelling
>>> wrap.run_pandora(num_cores=24)


Wrapper
----------------------------------------

.. automodule:: PANDORA.Wrapper.Wrapper
    :members:
    :undoc-members:
    :show-inheritance:

run_model
----------------------------------------

.. automodule:: PANDORA.Wrapper.run_model
    :noindex:
    :members:
    :undoc-members:

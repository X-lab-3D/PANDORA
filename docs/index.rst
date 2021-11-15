.. Pandora documentation master file, created by
   sphinx-quickstart on Mon Feb 26 14:44:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pandora
========

`Pandora`_ is a fast, modularized and highly flexible anchor-restrained modelling pipeline 
for generating peptide-MHC structures.

PANDORA contains useful APIs for pre-processing pMHC template data, identify the most adequate
template and adapt the modelling to the specifc case.

.. _`Pandora`: https://github.com/X-lab-3D/PANDORA

**DeepRank highlights**:

- Predefined atom-level and residue-level PPI feature types
   - *e.g. atomic density, vdw energy, residue contacts, PSSM, etc.*
- Predefined target types
   - *e.g. binary class, CAPRI categories, DockQ, RMSD, FNAT, etc.*
- Flexible definition of both new features and targets
- 3D grid feature mapping
- Efficient data storage in HDF5 format
- Support both classification and regression (based on PyTorch)

Tutorial
--------

.. toctree::
   :maxdepth: 3

   tutorial

API Reference
-------------

.. toctree::
   :maxdepth: 3

   Documentation

Indices
-------

* :ref:`genindex`
* :ref:`modindex`

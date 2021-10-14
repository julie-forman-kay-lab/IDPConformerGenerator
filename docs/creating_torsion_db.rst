===================================
Creating the torsion angle database
===================================

Why do you need this DB?
------------------------

IDPConfGen builds conformers by pulling torsion angles from a pool. This pool of
angles is connected to the residue identity where these angles have its origin.
So we need to provide idpconfgen a specific file for it to operate. This
"database", as we call it, is in fact a text-based JSON file. Therefore, it can
be visually inspected and manually corrected, if needed (though this is never
necessary). We decided for an open-text format for traceability, openness, and
reproducibility purposes.

Building the IDPConfGen torsion angle DB is the first step you need to do
before generating conformers. However, you **don't** need to generate this DB
every time you need to run idpconfgen. Once create, you can keep using it until
your project has different requirements and you need to extend it with some
new structures, or specifically remove some.

You can also distribute the DB you created among your collaborators or save it
for the record of your protocols or data. However, process of creating the DB is
reproducible. So, having the initial input and the series of commands would
recreate the same DB entirely.

How to build the DB
-------------------

IDPConformerGenerator has a series of command-line interfaces (CLIs) that must be used
in tandem to create the torsion angle database that idpconfgen requires to build
conformers. These CLIs provide the user with the flexible means to create this
DB at will and according to the project specificity's.

First step) the list of PDB chains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step for creating the DB is to define which PDB chains you want to use
for the database. Normally, these are "non-redundant PDBs under certain angstrom
resolution", which can be a list of about 20,000 PDBs. Defining such a list
manually is impossible. But we don't have to do it. `Dr. Dunbrack hosts the PICES
server`_ where different lists of non-redundant PDBs of different parameters are
periodically created. Our preference is to use the::

    cullpdb_pc90.0_res0.0-2.2_len40-10000_R0.25_Xray_d*

At the time of writing, the latest was `d2021_10_01`_.

.. _Dr. Dunbrack hosts the PICES server: https://dunbrack.fccc.edu/pisces/
.. _d2021_10_01: https://dunbrack.fccc.edu/pisces/download/cullpdb_pc90.0_res0.0-2.2_len40-10000_R0.25_Xray_d2021_10_01_chains29602

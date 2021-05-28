IDPConfGen Example
==================

To build the torsion angle database you need to provide IDPConfGen with a list
of PDB/mmCIF files. Our advice is that you use a culled list of your choice from
the `Dr. Dunbrack PISCES database <http://dunbrack.fccc.edu/PISCES.php>`_.

The PDB chain id list should have the format of the :code:`cull100` file in this
folder, which emulates the format provided by PISCES. Actually, the only column
that is read by IDPConfGen is the first column. No header lines are allowed.

::

    12E8H       221  XRAY        1.900    0.22    0.27  
    16PKA       415  XRAY        1.600    0.19    0.23  
    1A05A       358  XRAY        2.000    0.20    0.28  
    1A0JA       223  XRAY        1.700    0.17    0.21  
    1A12A       413  XRAY        1.700    0.19    0.22  
    1A1XA       108  XRAY        2.000    0.21    0.25  

The first three alphanumberic characters are the PDBID codes. The forth (or
more) are the PDB chain identifier. mmCIF files have chain ids of several
characters.

Run the following command to down the PDB files::

    idpconfgen pdbdl cull100 -u -n 2 -d pdbs.tar

You can inspect all options of this subclient with::

    idpconfgen pdbdl -h

This execution will create a :code:`pdbs.tar` file with the parser PDBs. Those
PDBs contain only the information needed for IDPConfGen. Unnecessary chains or
residues were removed.

We now proceed to the identification of secondary structure elements. For
this you need to have DSSP program installed. IDPConfGen is built to be modular.
You can also use any other program to calculate secondary structure but you
would need to implement the respective parser. The DSSP parser is implemented.
To install DSSP follow these instructions: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/issues/48.

The following command will operate on the :code:`pdbs.tar` file and will create
temporary files and a result file with the DSSP information::

    idpconfgen sscalc mkdssp pdbs.tar -m 3 -rd

You will see that the files :code:`sscalc.json` and :code:`sscalc_splittled.tar`
were created. :code:`sscalc.json` matches the sequence information with that of
the secondary structure identity. :code:`sscalc_splitted.tar` contains the PDB
chains split into continuous chains.

Now we need to calculate the torsion angles. There are several options available
in the command line but these are good defaults::

    idpconfgen torsions sscalc_splittled.tar -sc sscalc.json -o idpconfgen_database.json

This will create the final file :code:`idpconfgen_database.json`. This is the
file IDPConfGen needs to generate conformers. This is the torsion angle database
file. If you open it, you will see it is a regular human-readable :code:`json` file.

Finally, to generate conformers you will use the :code:`build` interface. The
build interface has several parameters that can be use to fine tune the
conformer construction protocol. You can read deeper instructions in the
documentations and in the client help. The following is a good default::

    idpconfgen build -db idpconfgen_database.json -seq EGAAGAASS -nc 10 -dr L+ -et 0 -xp 1 1 1 1 -rs 0

After some time you will see 10 conformers in the folder.

All IDPConfGen operations can be distributed over multiple cores. Use the flag
:code:`-n` to indicate the number of cores you wish to use.

This is all you need to know for a basic usage of IDPConfGen. Now you can use a
larger PDB database to feed the conformational sampling.

We look forward to your feedback,
The IDPConfGen Team

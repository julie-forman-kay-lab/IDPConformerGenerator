IDPConfGen Small Peptide Test
=============================

.. start-description

.. note::
    The purpose of this module is to test the installation of IDPConformerGenerator
    on a small example peptide. The database you will procure here will not be used
    further cases, much less in practice due to the database housing only 100 PDB IDs.

To build the torsion angle database you need to provide IDPConfGen with a list
of PDB/mmCIF files. Our advice is that you use a culled list of your choice from
the `Dunbrack PISCES database <http://dunbrack.fccc.edu/PISCES.php>`_.

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

The first three alphanumeric characters are the PDBID codes. The forth (or
more) are the PDB chain identifier. mmCIF files have chain IDs of several
characters.

Run the following command to download the PDB files::

    idpconfgen pdbdl cull100 -u -n 2 -d pdbs.tar

You can inspect all options of this (and any other) subclient with :code:`-h`::

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

    idpconfgen sscalc pdbs.tar -m 3 -rd -cmd <DSSP EXEC>

You will see that the files :code:`sscalc.json` and :code:`sscalc_splittled.tar`
were created. :code:`sscalc.json` matches the sequence information with that of
the secondary structure identity. :code:`sscalc_splitted.tar` contains the PDB
chains split into continuous chains.

Now we need to calculate the torsion angles. There are several options available
in the command line but these are good defaults::

    idpconfgen torsions sscalc_splitted.tar -sc sscalc.json -o idpconfgen_database.json

This will create the final file :code:`idpconfgen_database.json`. This is the
file IDPConfGen needs to generate conformers. This is the torsion angle database
file. If you open it, you will see it is a regular human-readable :code:`json` file.

Finally, to generate conformers you will use the :code:`build` interface. The
build interface has several parameters that can be use to fine-tune the
conformer construction protocol. You can read deeper instructions in the
documentation and client help. The following is a good default that uses 
the FASPR method for adding side chains::

    idpconfgen build -db idpconfgen_database.json -seq EGAAGAASS -nc 10 --dhelix --dstrand -et 'pairs' -rs 0

After some time you will see 10 conformers in the folder.
Please note that searching for loops is enabled by default for :code:`--dloop`.
Appending :code:`--dhelix --dstrand` will extend sampling to alpha-helicies and 
beta-strands in addition to loops. For more information on usage, please view :code:`idpconfgen build -h`.

If you would like to use a program to assign secondary structure bias based on NMR chemical shifts (e.g. CheSPI or δ2D) to
employ probabilistic custom secondary structure sampling (CSSS), the :code:`probs8_[ID].txt` output
from CheSPI or .TXT output from δ2D would have to be standardized into a user-editable text file indicating the
probability of secondary structures (based on DSSP codes) on a per residue basis.
The following example will process CheSPI output and assign probabilities to L/H/E based on H/G/I/E/ /T/S/B
structures on a per residue basis::

    idpconfgen csssconv -p8 probs8_ex.txt -o csss_ex.json

For simplicity, secondary structures from CheSPI and δ2D are grouped into L/H/E as defined by idpconfgen.
If you do not want this grouping feature, please build the database above without :code:`-rd` 
and run :code:`csssconv` with :code:`--full` to avoid grouping. 
To build with the CSSS file, :code:`-csss` would have to point to the CSSS.JSON file::

    idpconfgen build -db idpconfgen_database.json -seq EGAAGAASS -nc 10 -csss csss_ex.json -et 'pairs' -rs 0 --dloop-off

After some time you will see 10 conformers in the folder with the probabilistic CSSS.

All IDPConfGen operations can be distributed over multiple cores. Use the flag
:code:`-n` to indicate the number of cores you wish to use. Appending only :code:`-n`
will use all available CPU threads except for one.

This is all you need to know for a basic usage of IDPConfGen. Now you can use a
larger PDB database to feed the conformational sampling.

.. end-description

We look forward to your feedback,
The IDPConfGen Team

Building Folded Regions - N-term. and C-term. IDRs (FLDR/S)
===========================================================

.. start-description

The following example will walk you through building a IDR tails on folded regions
using the Folded Local Domain Region Sampling (FLDR/S ``fldrs`` subclient).

For this exercise, we will be constructing the N-IDR and C-IDR tails on the
human CNOT7 deadenylase. Please enter the example :code:`example/cnot7_example` folder
where you will find the FASTA sequence: ``cnot7.fasta``, and a PDB of the folded region
from PDB ID 4GMJ: ``cnot7_4gmj_fld_11-263.pdb``.

Steps from now will assume you're in the working director of ``example/cnot7_example``
and have already compiled your preferred reusable IDPConformerGenerator database. For
instructions on the database, please visit the previous exercise "A Real Case Scenario".

Aside from the ``build`` subclient, we will be using the ``fldrs`` subclient to build
disordered tails on folded domains like so::

    idpconfgen fldrs \
        -db <PATH TO DATABASE.JSON> \
        -seq cnot7.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        -fld cnot7_4gmj_fld_11-263.pdb \
        -of ./cnot7_fldrs_L+_faspr \
        -n

The ``fldrs`` subclient with automatically detect the N-IDR and C-IDR tail based on mismatches
in the primary sequence of the ``.fasta`` file (or input sequence from ``-seq``) and the PDB
file of the folded domain.

Any other parameters will only impact the disordered regions generated. Additional settings
include ``-tol`` and ``-kt``, where the former sets a tolerance for clashes between the
disordered tail(s) and folded domain while the latter acts as a switch to retain the
disordered tail(s) generated in the building process. By default, disordered tails are
deleted after full length conformers are generated.

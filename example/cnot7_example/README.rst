Modeling Disordered Region Tails on a Folded Domain
===================================================

.. start-description

The following example will walk you through building N-terminal and C-terminal
IDR tails on folded regions using the Local Disordered Region Sampling (LDRS
``ldrs`` subclient).

For this exercise, we will be constructing the tails on the human CNOT7
deadenylase protein.  Please enter the example :code:`example/cnot7_example`
folder where you will find the complete CNOT7 sequence: ``cnot7.fasta``, and a
PDB of the folded region from PDB ID 4GMJ: ``4GMJ_CNOT7.pdb``.

.. note::
    If your input PDB has phosphorylated residues such as phosphorylated
    threonine and serine, please change the three-letter code in the PDB file
    indicating the residue label to the non-modified version. For example:
    ``TPO`` phosphorylated threonine will become ``THR`` and ``SEP``
    phosphorylated serine will become ``SER``.

    Using the ``resre`` subclient can help you with this.

Steps from now will assume you're in the ``example/cnot7_example`` directory and
have already created your preferred reusable IDPConformerGenerator *database*.
For instructions on the database, please visit the previous exercise "A Real
Case Scenario".

To generate the disordered terminal tails on CNOT7 run the following command::

    idpconfgen ldrs \
        -db <PATH TO DATABASE.JSON> \
        -seq cnot7.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        -fld 4GMJ_CNOT7.pdb \
        -of ./cnot7_ldrs_L+_faspr \
        -n

The ``ldrs`` subclient with automatically detect the N-IDR and C-IDR tail based on mismatches
in the primary sequence of the ``.fasta`` file (or input sequence from ``-seq``) and the PDB
file of the folded domain. This command took approximately 3 minutes on a single workstation with
64 GB DDR4 RAM and 50 CPU threads (``-n 50``) clocked at 3.0 GHz.

To check your outputs against what is to be expected for this tutorial section. Please click
`here <https://www.dropbox.com/sh/6j9ahb4r2od45kh/AAAqPWyMoS9cZQiiaWQrpv7Ua?dl=0>`_
and download the archive named ``cnot7_ldrs_example.zip``.

.. note::
    Sidechain clashes may appear if you use the FASPR method for packing on sidechains
    above. To guarantee no sidechain clashes, we recommended either lowering the steric-clash
    tolerance using the ``-tol`` flag above or generating backbone-only conformers first
    then packing sidechains later with MCSCE as described below.

To generate backbone-only IDR tails on CNOT7 then pack sidechains on the IDRs with MCSCE.
We will be using agnostic secondary structure sampling here with ``--dany``.::

    idpconfgen ldrs \
        -db <PATH TO DATABASE.JSON> \
        -seq cnot7.fasta \
        --dloop-off \
        --dany \
        -etbb 100 \
        -dsd \
        -nc 1000 \
        -fld 4GMJ_CNOT7.pdb \
        -of ./cnot7_ldrs_ANY_bb \
        -n
    
    idpconfgen resre \
        ./cnot7_ldrs_ANY_bb/ \
        -of ./cnot7_ldrs_ANY_bb_resre \
        -pt 126:HIP,157:HIP,225:HIP,249:HIP,258:HIP, \
        -n

    mkdir cnot7_ldrs_ANY_mcsce

    mcsce \
        ./cnot7_ldrs_ANY_bb_resre \
        64 \
        -w \
        -o ./cnot7_ldrs_ANY_mcsce \
        -l ./mcsce_log \
        -s \
        -m simple \
        -f 12-262

.. note::
    You can access the MCSCE software `here <https://github.com/THGLab/MCSCE>`_
    to ignore folded regions and add post-translational modifications during the
    sidechain packing process.
    
    If you run into an error with ``mcsce`` and your input PDB has histines labled as ``HIS``,
    please change the three-letter code in the PDB file to ``HIP`` to account for all
    protonation states.
    
    Using the ``resre`` subclient like so above can help you with this.

Any other parameters will only impact the disordered regions generated. Additional settings
include ``-tol`` and ``-kt``, where the former sets a tolerance for clashes between the
disordered tail(s) and folded domain while the latter acts as a switch to retain the
disordered tail(s) generated in the building process. By default, disordered tail only
conformers are deleted after full length conformers are generated.

.. end-description

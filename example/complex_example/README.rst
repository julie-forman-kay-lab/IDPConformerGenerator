Modeling Disordered Regions in a Multi-Chain Protein Complex
============================================================

.. start-description

The following example will walk you through building intrinsically disordered
regions on multiple chains of a protein complex using the Local Disordered
Region Sampling (LDRS ``ldrs`` subclient).

For this exercise, we will be constructing a combination of the three cases of
N-IDR, L-IDR, and C-IDR on both chain A and chain B of the crystal structure of the
D1D2 sub-complex from the human SNRNP core domain. Please enter the example 
:code:`example/complex_example` folder where you will find the complete set of
sequences: ``D1D2.fasta``, and a PDB of the complex from the
RCSB PDB ID 1B34: ``1B34.pdb``.

Steps from now will assume you're in the ``example/complex_example`` directory and
have already created your preferred reusable IDPConformerGenerator *database*.
For instructions on the database, please visit the previous exercise "A Real
Case Scenario".

Due to the automated process of multi-chain detection and building we can generate
a set of 10 structures with a single command::

    idpconfgen ldrs \
        -db <PATH TO DATABASE.JSON> \
        -seq D1D2.fasta \
        -nc 10 \
        -fld 1B34.pdb \
        --dloop-off \
        --dany \
        -of ./D1D2_ldrs_ANY_faspr \
        -n

The ``ldrs`` subclient with automatically detect the all IDRs and their corresponding
chains based on sequence similarity and mismatches in the primary sequence of the ``.fasta``
file and the PDB file of the folded domain. This command took approximately an hour on a
single workstation with 64 GB DDR4 RAM and 10 CPU threads (``-n 10``) clocked at 3.0 GHz.

To check your outputs against what is to be expected for this tutorial section. Please click
`here <https://www.dropbox.com/sh/6j9ahb4r2od45kh/AAAqPWyMoS9cZQiiaWQrpv7Ua?dl=0>`_
and download the archive named ``d1d2_complex_ldrs_example.zip``.

.. end-description

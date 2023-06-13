Modeling Disordered Regions Within Folded Domains
=================================================

.. start-description

The following example will walk you through building an 86-residue-long IDR
connecting two folded domains. By now, we expect you to be familiar with the
Local Disordered Region Sampling (LDRS ``ldrs`` subclient). If not, please visit
the ``cnot7_example``.

For this exercise, we will construct the intrinsically disordered region from
residues 49-135 on the STAS domain of SLC26A9 (PDB ID 7CH1, UniProt A0A7K9KBT8).
In the :code:`example/slc26a9_example` folder you will find the complete FASTA
sequence: ``SLC26A9_STAS.fasta``, and a PDB of the folded region from PDB ID
4GMJ: ``7CH1_SLC26A9.pdb``.

.. note::
    Modeling an IDR between folded regions may take a while depending on various
    factors such as the length of the IDR to model, the distance between the
    chain breaks, the location of the chain break, and the presence of folded
    parts that restrict the growth of the chain.


To continue the tutorial, navigate to the ``example/slc26a9_example`` directory.
Ensure you have already compiled your preferred reusable IDPConformerGenerator
database. For instructions on the database, please visit the previous exercise,
"A Real Case Scenario".

We will be using the ``ldrs`` subclient to model 50 conformations of the
intrinsically disordered region::

    idpconfgen ldrs \
        -db <PATH TO DATABASE.JSON> \
        -seq slc26a9_example.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 50 \
        --dloop-off \
        --dany \
        -fld 7CH1_SLC26A9.pdb \
        -of ./slc26a9_ldrs_ANY_faspr \
        -n

From the `.fasta` file, the ``ldrs`` subclient will automatically identify the
N-IDR, the C-IDR, and any IDRs missing between folded domains; and construct
those.

Advanced LDRS Usage
-------------------

For those more proficient in Python, the modularity of LDRS allows users access
to **advanced features** that we describe here. Such features include profiting
from better parallelization that allows modeling longer IDRs in a shorter time
(e.g., 221 residues in one or two days). Here, we explain how to use the
modularity of IDPConformerGenerator to exploit its total capacity for modeling
IDRs by writing two new Python scripts that import IDPConfGen machinery.

The logic behind the LDRS subclient for modeling and IDR connecting two folded
domains assumes that we have an N-IDR-like case at the C-terminal region of the
chain break and a C-IDR-like case at the N-terminal region of the other chain
break. Thus, when defining the IDR sequence in the fast file given to the `-seq`
parameter, we need to provide two overlapping residues at each. Those will be
"QK" and "LA" in this example.

We have already prepared the IDR sequence for this example; see the
``SLC26A9_IDR.fasta`` file. You can perform sequence alignment between the
``SLC26A9_IDR.fasta`` file and the ``slc26a9_example.fasta`` file to visualize
where the IDR fits within the whole protein sequence.

Here is a brief overview of what we will do to speed up the process of closing
the chain break with an all-atom IDR model. We show two scripts you can use as
templates for using the LDRS features of IDPConformgerGenerator.

1. Generate 10,000 backbone-only structures (more = better sampling) of the IDR:

    .. code-block:: bash
        idpconfgen build \
            -db <PATH TO DATABASE.JSON> \
            -seq SLC26A9_IDR.fasta \
            -etbb 100 \
            -dsd \
            -nc 10000 \
            --dloop-off \
            --dany \
            -of ./slc26a9_idr_ANY_bb \
            -n

The idea is to have IDPConformerGenerator generate a large library of IDPs that
may represent the IDR to model. We are exploiting IDPConformerGenerator's speed
and diversity for generating conformers. The time rate-limiting step here is the
``next-seeker`` protocol (see scripts), where we have to compare all of the
structures in ``slc26a9_cterm`` to all structures in ``slc26a9_nterm`` to find
our candidates for sidechain addition.

2. Create the necessary folders for the script to run:

    .. code-block:: bash
        mkdir slc26a9_cterm
        mkdir slc26a9_nterm
        mkdir slc26a9_matches
        mkdir slc26a9_sidechains
        mkdir slc26a9_results

3. Change the paths in the script ``slc26a9_shortcut.py`` and run it; always use
the ``idpconfgen`` Python environment.

4. Model the sidechains onto the conformers generated previously in the
B
``matches`` folder using `MC-SCE software<https://github.com/THGLab/MCSCE>`_:

    .. code-block:: bash
        mcsce ./slc26a9_matches 64 -w -s -o ./slc26a9_sidechains -l ./mcsce_log

5. Use ``psurgeon()`` in ``slc26a9_stitching.py`` script to attach the all-atom
IDR models to the folded domain.

To further save time, especially on a computing cluster, we can split the
conformers in the ``cterm`` and ``nterm`` folders and run jobs in parallel or
request more workers. Please note that this shortcut is not a memory-intensive
task, so 8 GB of RAM is sufficient to run the ``next-seeker`` protocol.

.. end-description

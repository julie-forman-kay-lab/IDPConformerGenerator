A Real Case Scenario
====================

.. start-description

.. note::
    The purpose of this module is to use IDPConformerGenerator as you would with
    a real protein of interest (unfolded drkN SH3 is presented here). The database
    procured from this module can be re-usable for future cases.

The example with a small peptide in the :code:`example` folder is a good way to
get introduced to IDPConfGen. Although building other IDP conformer ensembles
use the same workflow as the previous example, we will go over more detailed
usage examples with a well studied IDP, the unfolded state of the drkN SH3 domain.

Chemical shift data for the unfolded state of the drkN SH3 domain (BMRB ID: 25501) has been already processed with
δ2D and CheSPI and secondary structure propensity calculations can be found in 
:code:`example/drksh3_example` as :code:`drk_d2D.txt` and :code:`probs8_25501_unfolded.txt`
respectively for δ2D and CheSPI output.

An extensive culled list is in ``cull.tar``. Unpack it with::

    tar -xf cull.tar

This is the same culled list used in the IDPConfGen `main publication <link-to-DOI>`_.
However feel free to choose your own from the `Dunbrack PISCES database
<http://dunbrack.fccc.edu/PISCES.php>`_.

Steps from now on will assume you're in the working directory of :code:`example/drksh3_examples`.

To initialize the database if you do not already have one, we must download the PDB files from our culled list
(can be found in the supplemental package in the IDPConformerGenerator paper)::

    idpconfgen pdbdl cullpdb_pc90_res2.0_R0.25_d201015_chains24003 -u -n -d pdbs.tar

Next we will create temporary files storing the secondary structure information for each
PDB file downloaded. Later to be processed for their torsion angles::

    idpconfgen sscalc pdbs.tar -rd -n -cmd <DSSP EXEC>

Please note that since IDPConfGen is a toolkit, many of these modules can be used with
custom folders or ``.tar`` files.

Finally, torsion angles are extracted and the database we will use for future calculations
can be created with the :code:`torsions` subclient::

    idpconfgen torsions sscalc_splitted.tar -sc sscalc.json -n -o idpconfgen_database.json

Now we're ready to construct multiple conformer ensembles for the unfolded states of the drkN SH3 domain. To build 100 conformers,
sampling only the loop region, the default limits to the backbone and side chain L-J energy potentials are 
be 100 kJ and 250 kJ respectively, using default fragment sizes, no substitutions, and to have side chains added with FASPR::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -nc 100 \
        -of ./drk_L+_nosub_faspr \
        -n

:code:`idpconfgen` is deterministic. Therefore, the random seed defines the sampling progression - 
read :ref:`here <F.A.Q.s>` for more information.

To switch the side chain building algorithm to MC-SCE (recommended), you would first have to install MC-SCE.
Please re-visit the :ref:`installation <Installation>` page to get MC-SCE set up. Here's the following example::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -nc 100 \
        -scm mcsce \
        -of ./drk_L+_nosub_mcsce \
        -n

.. note::
    Running MC-SCE within IDPConformerGenerator can be memory (RAM) intensive.
    Consider running with a lower number of CPU threads using the `-n` flag if
    necessary.

The defaults for :code:`--mcsce-n_trials` is 16 while using the :code:`--mcsce-mode exhaustive`, however
we recommend trials larger or equal to 100 for smaller conformer pools. In this exercise, we will be using the
default MC-SCE side chain building mode :code:`simple`.

If you're encountering an error with MC-SCE running interally through IDPConformerGenerator,
we recommend you to generate backbones first, then pack sidechains after. For example, these would be the commands,
required to generate backbones first and then sidechains.::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -nc 100 \
        -dsd \
        -of ./drk_L+_nosub_bb \
        -n
    
    # Make the output folder for MC-SCE
    mkdir ./drk_L+_nosub_mcsce

    # Run MC-SCE in the idpconfgen environment because it's already installed
    mcsce ./drk_L+_nosub_bb 100 \
        -w \
        -o ./drk_L+_nosub_mcsce \
        -s \
        -m simple \
        -l drk_L+_nosub_mcsce_log

As stated in the :code:`idpconfgen build -h`, sampling using other secondary structure
parameters required :code:`--dloop` to be turned off :code:`--dloop-off`. For example, if we'd like to 
sample only helices and extended strands::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -nc 100 \
        -et 'pairs' \
        --dstrand \
        --dhelix \
        --dloop-off \
        -of ./drk_H+E+_nosub \
        -n

For sampling loops, helices, and strands, we would specify :code:`--dhelix --dstrand`
where :code:`--dloop` is turned on by default. However, sampling without biasing for secondary structure
can be done with :code:`--dany --dloop-off`.

To sample using custom secondary structure sampling (CSSS) a CSSS database (.JSON) file needs
to be created specifying the secondary structure probabilities for each residue. This can be
done using the :code:`makecsss` module if chemical shift data is not readily available, if you'd
like to edit a pre-existing CSSS.JSON, or create a new file. Here's an example for making a 
custom CSSS.JSON file that samples only helices for residues 15-25 of the unfolded state of the drkN SH3 domain
and loops for everything else::

    idpconfgen makecsss -cp 1-14 L 1.0|15-25 H 1.0|26-59 L 1.0 -o cust_csss_drk.json

If chemical shift files are readily available, consider using CheSPI or δ2D to generate the CSSS.JSON.
δ2D predictions have been included in the :code:`example/drksh3_ex_resources` folder as :code:`drk_d2D.txt`.
CheSPI :code:`probs8_*` predictions have been included in the :code:`example/drksh3_ex_resources` folder
as :code:`probs8_25501_unfolded.txt`.

To convert output from δ2D to CSSS, use the :code:`csssconv` subclient with flag :code:`-d2D`::

    idpconfgen csssconv -d2D drk_d2D.txt -o csss_drk_d2D.json

To convert output from CheSPI to CSSS, use the :code:`csssconv` subclient with flag :code:`-p8`::

    idpconfgen csssconv -p8 probs8_25501_unfolded.txt -o csss_drk_chespi.json

The outputted :code:`csss_*.json` files will be used for the :code:`-csss` flag in the :code:`build` subclient.
For example, constructing 100 conformers for the unfolded state of the drkN SH3 domain using the δ2D predictions and the same settings for
energy and MC-SCE as above::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -nc 100 \
        -csss csss_drk_d2D.json \
        --dloop-off \
        -et 'pairs' \
        -of ./drk_CSSSd2D_nosub \
        -n

The default fragment size probabilities for building are (1, 1, 3, 3, 2) for fragment sizes of (1, 2, 3, 4, 5) respectively.
To change this, we would have to create a :code:`.TXT` file with two columns, the first specifying what fragment sizes
from lowest to highest, the second specifying their relative probabilities. We have provided an example in
:code:`example/drksh3_ex_resources` as :code:`customFragments.txt`. To use these custom fragment size probabilities with CSSS::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -nc 100 \
        -xp customFragments.txt \
        -csss csss_drk_d2D.txt \
        --dloop-off \
        -et 'pairs' \
        -of ./drk_fragN_CSSSd2D_nosub \
        -n

Finally, to expand torsion angle sampling beyond the residue identity, we can provide a residue tolerance map using the :code:`-urestol` flag in the
:code:`build` subclient. For this example, we will be using columns 5, 3, and 2 from the `EDSSMat50 <https://www.nature.com/articles/s41598-019-52532-8>`_
substitution matrix::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -nc 100 \
        --dany \
        --dloop-off \
        -urestol '{"R":"RK","D":"DE","C":"CY","C":"CW","Q":"QH","E":"ED","H":"HYQ","I":"IVM","I":"IL","K":"KR","M":"MI","M":"MVL","F":"FY","F":"FWL","W":"WYFC","Y":"YF","Y":"YC","Y":"YWH"}' \
        -et 'pairs' \
        -of ./drk_ANY_sub532 \
        -n

Please note for the above run, we are sampling the torsion angle database disregarding secondary structure
with the :code:`--dany` flag.

Hopefully this more in-depth realistic example with the unfolded state of the drkN SH3 domain has provided you with the utilities and usage examples
to explore IPDConfGen more with your custom protein systems.


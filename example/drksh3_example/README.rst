A real case scenario with drkSH3
================================

.. start-description

The example with a small peptide in the :code:`example` folder is a good way to get introduced
to IDPConfGen. Although building other IDP conformer ensembles use the same workflow as the 
one provided in :code:`example/README.rst`, we will go over more detailed usage examples with
a well studied IDP, the drkN SH3 domain.

Chemical shift data for the unfolded state of the drkN SH3 domain (BMRB ID: 25501) has been already processed with
δ2D and CheSPI and secondary structure propensity calculations can be found in 
:code:`example/drksh3_ex_resources` as :code:`drk_d2D.txt` and :code:`probs8_25501_unfolded.txt`
respectively for δ2D and CheSPI output.

We have also provided a culled list of PDB IDs in the :code:`example` folder. 
This is the same culled list used in the IDPConfGen `paper <link-to-DOI>`_. However feel free to choose your
own from the the `Dunbrack PISCES database <http://dunbrack.fccc.edu/PISCES.php>`_.

Steps from now on will assume you're in the working directory of :code:`example/drksh3_ex_resources`.

To initialize the database if you do not already have one, we must download the PDB files from our culled list
(can be found in the supplemental package in the IDPConformerGenerator paper)::

    idpconfgen pdbdl cullpdb_pc90_res2.0_R0.25_d201015_chains24003 -u -n -d pdbs.tar

Next we will create temporary files storing the secondary structure information for each
PDB file downloaded. Later to be processed for their torsion angles::

    idpconfgen sscalc pdbs.tar -rd -n -cmd <DSSP EXEC>

Please note that since IDPConfGen is a toolkit, many of these modules can be used with
custom folders or .tar files.

Finally, torsion angles are extracted and the database we will use for future calculations
can be created with the :code:`torsions` subclient::

    idpconfgen torsions sscalc_splitted.tar -sc sscalc.json -n -o idpconfgen_database.json

Now we're ready to construct multiple conformer ensembles for the unfolded states of the drkN SH3 domain. To build 100 conformers,
sampling only the loop region, limiting the backbone and side chain L-J energy potentials to 
be 100 kJ and 250 kJ respectively, using default chunk sizes, no substitutions, and to have side chains added with FASPR::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        -rs 0 \
        -et 'pairs' \
        -of ./drk_L+_nosub_faspr \
        -n

:code:`idpconfgen` is deterministic. Therefore, the random seed defines the sampling progression - 
read :ref:`here <F.A.Q.s>` for more information.

To switch the side chain building algorithm to MC-SCE (recommended), you would first have to install MC-SCE.
Please re-visit the :ref:`installation <Installation>` page to get MC-SCE set up. Here's the following example::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        -rs 0 \
        -et 'pairs' \
        -scm mcsce \
        -of ./drk_L+_nosub_mcsce \
        -n

The defaults for :code:`--mcsce-n_trials` is 16 while using the :code:`--mcsce-mode exhaustive`, however
we recommend trials larger than 100 for smaller conformer pools. In this exercise, we will be using the
default MC-SCE side chain building mode :code:`simple`.

As stated in the :code:`idpconfgen build -h`, sampling using other secondary structure
parameters required :code:`--dloop` to be turned off :code:`--dloop-off`. For example, if we'd like to 
sample only helices and extended strands::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        -rs 0 \
        -et 'pairs' \
        --dstrand \
        --dhelix \
        --dloop-off \
        -scm mcsce \
        -of ./drk_H+E+_nosub_mcsce \
        -n

For sampling loops, helices, and strands, we would specify :code:`--dhelix --dstrand`
where :code:`--dloop` is turned on by default. However, sampling without biasing for secondary structure
can be done with :code:`--dany --dloop-off`.

To sample using custom secondary structure sampling (CSSS) a CSSS database (.JSON) file needs
to be created specifying the secondary structure probabilities for each residue. This can be
done using the :code:`makecsss` module if chemical shift data is not readily available, if you'd
like to edit a pre-existing CSSS.JSON, or create a new file. Here's an example for making a 
custom CSSS.JSON file that samples only helices for residues 15-25 of drkN SH3 and loops for everything else::

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
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        -csss csss_drk_d2D.json \
        --dloop-off \
        -et 'pairs' \
        -scm mcsce \
        -of ./drk_CSSSd2D_nosub_mcsce \
        -n

The default chunk size probabilities for building are (1, 1, 3, 3, 2) for chunk sizes of (1, 2, 3, 4, 5) respectively.
To change this, we would have to create a :code:`.TXT` file with two columns, the first specifying what chunk sizes
from lowest to highest, the second specifying their relative probabilities. We have provided an example in
:code:`example/drksh3_ex_resources` as :code:`customChunk.txt`. To use these custom chunk size probabilities with CSSS::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        -xp customChunk.txt \
        -csss csss_drk_d2D.txt \
        --dloop-off \
        -et 'pairs' \
        -scm mcsce \
        -of ./drk_chunkN_CSSSd2D_nosub_mcsce \
        -n

Finally, to expand torsion angle sampling beyond the residue identity, we can provide a residue tolerance map using the :code:`-subs` flag in the
:code:`build` subclient. For this example, we will be using columns 5, 3, and 2 from the `EDSSMat50 <https://www.nature.com/articles/s41598-019-52532-8>`_
substitution matrix::

    idpconfgen build \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 100 \
        --dany \
        --dloop-off \
        -subs '{"R":"RK","D":"DE","C":"CY","C":"CW","Q":"QH","E":"ED","H":"HYQ","I":"IVM","I":"IL","K":"KR","M":"MI","M":"MVL","F":"FY","F":"FWL","W":"WYFC","Y":"YF","Y":"YC","Y":"YWH"}' \
        -et 'pairs' \
        -scm mcsce \
        -of ./drk_ANY_sub532_mcsce \
        -n

Please note for the above run, we are sampling the torsion angle database disregarding secondary structure
with the :code:`--dany` flag.

Hopefully this more in-depth realistic example with the unfolded state of the drkN SH3 domain has provided you with the utilities and usage examples
to explore IPDConfGen more with your custom protein systems.


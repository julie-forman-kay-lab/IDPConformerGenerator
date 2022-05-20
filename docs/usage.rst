=====
Usage
=====

Command-line
------------

To execute :code:`idpconfgen` command-line, run :code:`idpconfgen` in your
terminal window, after :ref:`installation <Installation>`::

    idpconfgen

or::

    idpconfgen -h

both will output the help menu.

Please note that all subclients have the :code:`-h` option to show help information.

:code:`idpconfgen` has several interfaces that perform different functions.
However, there is a sequence of interfaces that need to be used to prepare the
local torsion angle database and the files needed to build conformers. After
these operations executed, you will end up with a single :code:`json` file that
you can use to feed the build calculations. The other files are safe to be
removed.

The :code:`example` folder contains instructions to setup
IDPConformerGenerator. We review them here as well:

.. include:: ../example/README.rst


IDPConfGen Library
--------------------

To use IDP Conformer Generator in your project, import it as a library::

    import idpconfgen


Getting Started with IDPConfGen
-------------------------------

The example with a small peptide in the :code:`example` folder is a good way to get introduced
to IDPConfGen. Although building other IDP conformer ensembles use the same workflow as the 
one provided in :code:`example/README.rst`, we will go over more detailed usage examples with
a well studied IDP, the drkN SH3 domain.

Chemical shift data for drkN SH3 (BMRB ID: 25501) has been already processed with
delta2D and CheSPI and secondary structure propensity calculations can be found in 
:code:`example/drksh3_ex_resources` as :code:`drk_d2D.txt` and :code:`probs8_25501_unfolded.txt`
respectively for delta2D and CheSPI output.

We have also provided a culled list of PDB IDs in the :code:`example` folder. 
This is the same culled list used in the IDPConfGen `paper <link-to-DOI>`_. However feel free to choose your
own from the the `Dr. Dunbrack PISCES database <http://dunbrack.fccc.edu/PISCES.php>`_.

Steps from now on will assume you're in the working directory of :code:`example/drksh3_ex_resources`.

To initialize the database if you do not already have one, we must download the PDB files from our culled list::

    idpconfgen pdbdl cullpdb_pc90_res2.0_R0.25_d201015_chains24003 -u -n -d pdbs.tar

Next we will create temporary files storing the secondary structure information for each
PDB file downloaded. Later to be processed for their torsion angles.::

    idpconfgen sscalc pdbs.tar -rd -n -cmd <DSSP EXEC>

Please note that since IDPConfGen is a toolkit, many of these modules can be used with
custom folders or .tar files.

Finally, torsion angles are extracted and the database we will use for future calculations
can be created with the :code:`torsions` subclient::

    idpconfgen torsions sscalc_splitted.tar -sc sscalc.json -n -o idpconfgen_database.json

Now we're ready to construct multiple conformer ensembles for drkN SH3. To build 100 conformers,
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
read :ref:`here <Frequently Asked Questions (F.A.Q.s)>` for more information.

To switch the side chain building algorithm to MCSCE (recommended), you would first have to install MCSCE.
Please re-visit the :code:`docs/installation.rst` to get MCSCE set up. Here's the following example::

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
default MCSCE side chain building mode `simple`.

As stated in the :code:`idpconfgen build -h`, sampling using other secondary structure
parameters required :code:`--dloop` to be turned off `--dloop-off`. For example, if we'd like to 
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
custom CSSS.JSON file that samples only helices for residues 15-25 of drkN SH3 and loops for everything else.::

    idpconfgen makecsss -cp 1-14 L 1.0|15-25 H 1.0|26-59 L 1.0 -o cust_csss_drk.json

If chemical shift files are readily available, consider using CheSPI or delta2D to generate the CSSS.JSON.
delta2D predictions have been included in the :code:`example/drksh3_ex_resources` folder as `drk_d2D.txt`.
CheSPI `probs8_*` predictions have been included in the :code:`example/drksh3_ex_resources` folder
as `probs8_25501_unfolded.txt`.

To convert output from delta2D to CSSS, use the :code:`csssconv` subclient with flag :code:`-d2D`::

    idpconfgen csssconv -d2D drk_d2D.txt -o csss_drk_d2D.json

To convert output from CheSPI to CSSS, use the :code:`csssconv` subclient with flag :code:`-p8`::

    idpconfgen csssconv -p8 probs8_25501_unfolded.txt -o csss_drk_chespi.json

The outputted `csss_*.json` files will be used for the :code:`-csss` flag in the :code:`build` subclient.
For example, constructing 100 conformers for drkN SH3 using the delta2D predictions and the same settings for
energy and MCSCE as above::

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
:code:`example/drksh3_ex_resources` as `customChunk.txt`. To use these custom chunk size probabilities with CSSS.::

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
substitution matrix.::

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

Hopefully this more in-depth realistic example with drkN SH3 has provided you with the utilities and usage examples
to explore IPDConfGen more with your custom protein systems.

Exploring IDPConfGen Analyses Functions
---------------------------------------

Our vision for IDPConformerGenerator as a platform includes the analysis of your database and the PDBs generated by IDPConfGen.
To get started, the :code:`stats` subclient is a quick way to check how many hits for different sequence fragment matches you will
find in the database for your protein system of choice. It is also possible to include different secondary structure filters as well
as amino-acid substitutions to get a more accurate representation the number of hits in the database for your system.::

    idpconfgen stats \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        --dloop-off \
        --dany \
        -op drk_any \
        -of ./drk_any_dbStats

Another tool to investiagte the database is the :code:`search` subclient. To use this, you will need a tarball or folder of raw PDBs required
from the :code:`fetch` subclient. The :code:`search` function goes through the PDB headers to find keywords of your choice and returns the
number of hits and their associated PDBIDs in .JSON format.::

    idpconfgen fetch \
        ../cull100 \
        -d ./cull100pdbs/ \
        -u \
        -n

    idpconfgen search \
        -fpdb ./cull100pdbs/ \
        -kw 'drk,beta' \
        -n

After generating conformer ensembles with IDPConfGen, it is possible to do some basic plotting with the integrated plotting flags
in the :code:`torsions` and :code:`sscalc` subclients. For :code:`torsions`, you can choose to plot either omega, phi, or psi dihedral
angle distributions in a scatter plot format. For :code:`sscalc`, fractional secondary structure will be plotted in terms of DSSP codes
as well as fractions from the alpha, beta, or other regions of the Ramachandran space for your conformers of choice. The following example
plots the psi angle distributions and the fractional secondary structure of the `drk_CSSSd2D_nosub_mcsce` ensemble generated in the previous
module.::

    idpconfgen torsions \
        ./drk_CSSSd2D_nosub_mcsce \
        -deg \
        -n \
        --plot type=psi xlabel=drk_residues
    
To plot the fractional Ramachandran space information.::

    idpconfgen torsions \
        ./drk_CSSSd2D_nosub_mcsce \
        -deg \
        -n \
        --ramaplot filename=fracDrkRama.png colors=['o', 'b', 'k']

To plot the fractional secondary structure information.::

    idpconfgen sscalc \
        ./drk_CSSSd2D_nosub_mcsce \
        -u \
        -rd \
        -n \
        --plot filename=dssp_reduced_drk_.png

To see which plotting parameters can be modified, please refer to `src/idpconfgen/plotfuncs.py`. We have given a short list of modifyable parameters here.::

    --plot title=<TITLE> title_fs=<TITLE FONT SIZE> xlabel=<X-AXIS LABEL> xlabel_fs=<X-AXIS LABEL FONT SIZE> colors=<LIST_OF_COLORS>

Exploring MC-SCE and Int2Cart Integrations
------------------------------------------

Integrating the functions from our collaborators at the `Head-Gordon Lab <https://thglab.berkeley.edu/>_`,
IDPConformerGenerator has the ability to build with bond geometries derived from a recurrent neural network
machine learning model `Int2Cart <https://github.com/THGLab/int2cart>_`. Furthermore, as we introduced
the `MC-SCE <https://github.com/THGLab/MCSCE>_` method for building sidechains in the previous modules,
we would like to provide some examples on changing the default sidechain settings.

To use the Int2Cart method for bond geometries, the :code:`-bgeo_int2cart` flag needs to be defined during
the building stage.::

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
        -bgeo_int2cart \
        -of ./drk_CSSSd2D_nosub_int2cart_mcsce \
        -n

To change the number of trials for MC-SCE to increase success rate.::

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
        --mcsce-n_trials 32 \
        -of ./drk_CSSSd2D_nosub_32_trials_mcsce \
        -n


How to Efficiently Set Jobs up for HPC Clusters
-----------------------------------------------

Using the :code:`sethpc` subclient, users can generate bash scripts for SLURM managed
systems. Due to architecture of Python's multiprocessing module, IDPConformerGenerator is
unable to utilize the resources of multiple nodes on HPC clusters. However, with :code:`sethpc`,
users are able to request multiple nodes per job and :code:`sethpc` will automatically generate
the SBATCH scripts needed, along with an :code:`all*.sh` and :code:`cancel*.sh` script to run/cancel
all of the jobs generated with ease.

Please note that your queuing priority will not change requesting 5 nodes per job or 1 node per 5 jobs.

If multiple nodes are requested, at the end of all jobs, the :code:`merge` subclient can be run to
merge all of the conformers generated into one folder with the option of modify the naming-pattern
for each structure. Please see below for an example of running :code:`sethpc` and :code:`merge`.

To request 3 nodes to generate 512,000 structures of drkN SH3 with 10 hours per node.::
    
    idpconfgen sethpc \
        -des ./drk_hpc_jobs/ \
        --account def-username \
        --job-name drk_hpc \
        --nodes 3 \
        --ntasks-per-node 32 \
        --mem 16g \
        --time-per-node 0-10:00:00 \
        --mail-user your@email.com \
        -db idpconfgen_database.json \
        -seq drksh3.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 512000 \
        -csss csss_drk_d2D.json \
        --dloop-off \
        -et 'pairs' \
        -scm mcsce \
        -bgeo_int2cart \
        -of /scratch/user/drk/ \
        -n 32 \
        -rs 12 

To merge all of the folders created by the multi-node jobs.::

    idpconfgen merge \
        -tgt /scratch/user/drk/ \
        -des /scratch/user/drk/drk_CSSSd2D_nosub_multiple_mcsce \
        -pre drk_confs \
        -del
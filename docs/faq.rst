====================================
Frequently Asked Questions (F.A.Q.s)
====================================

How long does the database generation take?
-------------------------------------------
Although it is possible to generate the torsion angle and secondary structure database
for IDPConformerGenerator on the login node on ComputeCanada, we strongly advise against this.
It is recommended that the generation of this database be done locally or as a submitted job on HPC clusters.  

As per ComputeCanada recommendations, processes taking up more than 10 CPU-minutes and 
4 GB of RAM should be submitted as jobs via the slurm job scheduler. Furthermore, due to
the quantity of PDB files downloaded, many read-write requests will slow down the login node
for other users as well as create errors in file I/O processing in IDPConfGen. For reference,
generating a database using 24,002 PDB files takes 1 hour and 13 minutes (using 31 workers)
as a job on Graham@ComputeCanada. 

Generating this same database locally (using 63 workers of an AMD Threadripper 2990wx), takes
only 37 minutes (this is approximately half the time compared to 31 workers on Graham@ComputeCanada
due to majority of the time is taken up by the PDB downloading process). It is important to note that
this database can then be transferred onto HPC clusters and can be reused for all future calculations.
Users do not need to generate specific input files for each case of IDP conformer generation
(unless CSSS is used, in that case standard NMR STAR files can be read by CheSPI and the output of
CheSPI/delta2D can be converted to CSSS.JSON files by IDPConfGen’s :code:`csssconv` subclient.

What does it mean that IDPConfGen is deterministic?
---------------------------------------------------
Reproducibility is one of the most important pillars in scientific research. Thus, we've ensured
reproducibility by implementing a :code:`--random-seed` flag while building. Therefore, generating
ensembles with the same database file and building parameters will generate the same set of conformers.

The random seed parameter is also helpful for appending to incomplete conformer pools. For example,
if your target was 3000 conformers and the requested job-time was not enough and only 2000 were generated,
1000 more unique conformers can be generated simply by using a different integer for the :code:`-rs` flag.

What forcefield does IDPConfGen use by default to generate PDB files?
---------------------------------------------------------------------
IDPConfGen uses the Amberff14SB forcefield by default from `OpenMM <https://github.com/openmm/openmmforcefields>`_.
The forcefield can be changed via the :code:`--forcefield` flag during building.


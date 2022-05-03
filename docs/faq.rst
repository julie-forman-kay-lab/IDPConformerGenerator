====================================
Frequently Asked Questions (F.A.Q.s)
====================================

What are the best options for number of conformers and cores to use?
--------------------------------------------------------------------

It is best to have the number of conformers/trials equal to an integer multiple of the number of cores/workers.
For example, if you choose to use 32 cores via :code:`-n 32`, it's ideal to have :code:`-nc` as 32, 64, 128, etc.

For hyperthreadded systems, :code:`-n` will use the maximum number of THREADS - 1. For the best experience however,
we recommend running on the maximum number of physical cores for hyperthreadded CPUs.
For example, a Threadripper 2990wx has 32 cores and 64 threads, it's best to use :code:`-n 32` to avoid overloading
the system.

How long does the database generation take?
-------------------------------------------

Although it is possible to generate the torsion angle and secondary structure database for
IDPConformerGenerator on any system, including HPC clusters, we recommend local systems with the
largest number of workers and the fastest read-write times. Due to the quantity of PDB files downloaded,
many read-write requests will slow down the login node for other users on HPC resources as well as create
errors in file I/O processing in IDPConfGen. For reference, generating a database using 24,002 PDB files
takes 1 hour 13 minutes (using 31 workers) as a job on Graham@ComputeCanada, an HPC resource. Generating
this same database locally (using 63 workers of an AMD Threadripper 2990wx) takes only 37 minutes
(approximately half the time compared to 31 workers on Graham@ComputeCanada) since the bulk of the time
is the PDB downloading process. Importantly, this database can then be transferred onto HPC clusters and
can be reused for all future calculations.

What does it mean that IDPConfGen is deterministic?
---------------------------------------------------

Reproducibility is one of the most important pillars in scientific research. Thus, we've ensured
reproducibility by implementing a :code:`--random-seed` flag while building. Therefore, generating
ensembles with the same database file and building parameters on the same processing system
will generate the same set of conformers.

The random seed parameter is also helpful for appending to incomplete conformer pools. For example,
if your target was 3000 conformers and the requested job-time was not enough and only 2000 were generated,
1000 more unique conformers can be generated simply by using a different integer for the :code:`-rs` flag.

What forcefield does IDPConfGen use by default to generate PDB files?
---------------------------------------------------------------------

IDPConfGen uses the Amberff14SB forcefield by default from `OpenMM <https://github.com/openmm/openmmforcefields>`_.
The forcefield can be changed via the :code:`--forcefield` flag during building.


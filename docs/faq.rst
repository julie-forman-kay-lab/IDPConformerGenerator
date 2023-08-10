=======
F.A.Q.s
=======

What are the recommended hardware specifications?
-------------------------------------------------

IDPConformerGenerator was built with computational efficiency in mind. However,
having more CPU cores and more RAM would benefit larger protein systems (e.g.
for systems greater than 350 residues, having more than 16 cores and 64GB of RAM
is suggested). However, additional resources might be required to use certain
third-party integrations. For example, `Int2Cart
<https://github.com/THGLab/int2cart>`_, required high performance CUDA
compatible GPUs.

IDPConformerGenerator can also run on HPC clusters easily. Please refer to the :ref:`installation <Installation>`
instructions to set-up IDPConformerGenerator using :code:`conda` or :code:`virtualenv`.

How long does the database generation take?
-------------------------------------------

Although it is possible to generate the torsion angle and secondary structure database for
IDPConformerGenerator on any system, including HPC clusters, we recommend local systems with the
largest number of workers and the fastest read-write times. Due to the quantity of PDB files downloaded,
many read-write requests will slow down the login node for other users on HPC resources as well as create
errors in file I/O processing in IDPConfGen. For reference, generating a database using 24,002 PDB files
takes 1 hour 13 minutes (using 31 workers) as a job on the Graham cluster, an HPC resource of the Digital 
Research Alliance of Canada (DRAC), an HPC resource. Generating this same database locally (using 63 workers
of an AMD Threadripper 2990wx) takes only 37 minutes (approximately half the time compared to 31 workers on 
the Graham cluster (DRAC)) since the bulk of the time is the PDB downloading process. Importantly,
this database can then be transferred onto HPC clusters and can be reused for
all future calculations, as well as shared between collaborators.

What are the best options for combination of conformers and cores to use?
-------------------------------------------------------------------------

It is best to have the number of conformers/trials equal to an integer multiple of the number of cores/workers.
For example, if you choose to use 32 cores via :code:`-n 32`, it's ideal to have :code:`-nc` as 32, 64, 128, etc.
However, this is not mandatory and IDPConfGen will accommodate your `-nc` and
`-n` request accordingly.

For hyperthreadded systems, :code:`-n` will use the maximum number of THREADS - 1. For the best experience however,
we recommend running on the maximum number of physical cores for hyperthreadded CPUs.
For example, a Threadripper 2990wx has 32 cores and 64 threads, it's best to use :code:`-n 32` to avoid overloading
the system.

What's the difference between MC-SCE and FASPR?
-----------------------------------------------

IDPConfGen incorporates `FASPR <https://github.com/tommyhuangthu/FASPR>`_
to natively build sidechains on the conformer's backbone.

A part from the algorithmic differences between FASPR and MC-SCE (please refer
to their publications for details), FASPR does not add hydrogens to the
sidechains, so other tools are required to add hydrogens if they are desired.
MC-SCE does add hydrogens to the sidechains it generates.

Although MC-SCE is slower, it produces structures with no steric clashes, as it
uses the same clash quality controls as IDPConfGen. Furthermore, it also has
many settings that can be tweaked for flexible in silico experiments as opposed
to the default settings of FASPR.

How can I optimize the MC-SCE settings for large proteins?
----------------------------------------------------------

For larger proteins such as the Tau fragment (441 residues), we recommend
generating backbone-only conformers and use MC-SCE as a stand-alone program
afterwards.

During the backbone generation stage (i.e. :code:`idpconfgen build ... -dsd`), speed could be
improved by setting a higher backbone energy threshold (:code:`-etbb`). We recommend 250 as a minimum.

For the MC-SCE sidechain step, we recommend doing a small benchmark with 128 trials to see what the median
number of trials MC-SCE requires for a successful sidechain addition.


What does it mean that IDPConformerGenerator is deterministic?
--------------------------------------------------------------

Reproducibility is one of the most important pillars in scientific research. Thus, we've ensured
reproducibility by implementing a :code:`--random-seed` flag while building. Therefore, generating
ensembles with the same database file and building parameters on the same processing system
will generate the same set of conformers.

The random seed parameter is also helpful for appending to incomplete conformer pools. For example,
if your target was 3000 conformers and the requested job-time was not enough and only 2000 were generated,
1000 more unique conformers can be generated simply by using a different integer for the :code:`-rs` flag.

What forcefield does IDPConformerGenerator use by default to generate PDB files?
--------------------------------------------------------------------------------

IDPConfGen uses the Amberff14SB forcefield by default from `OpenMM <https://github.com/openmm/openmmforcefields>`_.
The forcefield can be changed via the :code:`--forcefield` flag, but no other is
currently implemented. For more advanced information see
https://github.com/julie-forman-kay-lab/IDPConformerGenerator/blob/master/src/idpconfgen/core/data/README.md.


Out-of-memory (OOM) error when running the ``count_clashes`` function 
---------------------------------------------------------------------

If you are scripting the ``count_clashes`` function or using an exceptionally large template structure
(e.g. greater than 450,000 atoms), you could run into a ``numpy.core._exceptions._ArrayMemoryError``.

To avoid this, try reducing the number of cores/multiple processes with the ``--ncores`` flag. Or
Decrease the size of your template structure, or request more RAM when submitting a job to a cluster.

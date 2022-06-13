IDPConformerGenerator
=======================

.. start-description

IDPConformerGenerator is a flexible, modular platform for generating ensembles
of disordered protein states that builds conformers by sampling backbone torsion
angles of relevant sequence fragments extracted from protein structures in the
RCSB Protein Data Bank.

IDPConformerGenerator can efficiently build large and diverse conformer pools of
disordered proteins, with user defined options enabling variable fractional
population of secondary structures, including matching those assigned based on
NMR chemical shift data. These conformer pools are intended to be utilized as
input for further approaches to match experimental data, such as re-weighting or
sub-setting algorithms.

.. end-description

How to Install
--------------

Installation instructions are described in ``docs/installation.rst``.

How to Use
----------

Usage instructions are described in the file ``docs/usage.rst``. See also
the example in the ``example/`` folder.

How to cite
-----------

.. start-citing

If you use IDPConformerGenerator, please cite::

    João M.C. Teixeira, Zi Hao Liu, Ashley Namini, Jie Li, Robert M. Vernon, Mickaël Krzeminski, Alaa A. Shamandy, Oufan Zhang, Mojtaba Haghighatlari, Lei Yu, Teresa Head-Gordon, Julie D. Forman-Kay
    IDPConformerGenerator: A Flexible Software Suite for Sampling Conformational Space of Disordered Protein States.
    bioRxiv 2022.05.28.493726 (2022). doi:10.1101/2022.05.28.493726.

.. end-citing

Version
-------
v0.6.7

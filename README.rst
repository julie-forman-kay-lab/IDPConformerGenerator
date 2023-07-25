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

    IDPConformerGenerator: A Flexible Software Suite for Sampling the Conformational Space of Disordered Protein States
    João M. C. Teixeira, Zi Hao Liu, Ashley Namini, Jie Li, Robert M. Vernon, Mickaël Krzeminski, Alaa A. Shamandy, Oufan Zhang, Mojtaba Haghighatlari, Lei Yu, Teresa Head-Gordon, and Julie D. Forman-Kay
    The Journal of Physical Chemistry A 2022 126 (35), 5985-6003
    DOI: 10.1021/acs.jpca.2c03726

.. end-citing

Version
-------
v0.7.0

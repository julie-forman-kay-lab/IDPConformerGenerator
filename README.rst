IDPConformerGenerator
=======================

.. start-description

.. image:: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/actions/workflows/tests.yml/badge.svg?branch=main
    :target: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/actions/workflows/tests.yml
    :alt: Test Status

.. image:: https://readthedocs.org/projects/idpconformergenerator/badge/?version=latest
    :target: https://idpconformergenerator.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/badge/idpconfgen-10.1021%2Facs.jpca.2c03726-blue
    :target: https://doi.org/10.1021/acs.jpca.2c03726
    :alt: IDPConformerGenerator Publication

.. image:: https://img.shields.io/badge/LDRS-10.1093%2Fbioinformatics%2Fbtad739-blue
    :target: https://doi.org/10.1093/bioinformatics/btad739
    :alt: LDRS Publication

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

**Note to users:** IDR conformers generated with ``ldrs``, specifically processed
using the ``align_coords()`` function in ``ldrs_helper.py`` prior to v0.7.17
may have the wrong stereochemistry. This bug has since been fixed. Thank you for
your understanding and we apologize for any inconvenience this has caused.

.. end-description

First Steps
-----------

To get a first glance on IDPConformerGenerator, read through our `first steps <https://idpconformergenerator.readthedocs.io/en/latest/first_steps.html>`_
documentation page.

How to Install
--------------

Installation instructions are described on the `installation page <https://idpconformergenerator.readthedocs.io/en/latest/installation.html>`_.

How to Use
----------

Usage instructions are described in the `usage page <https://idpconformergenerator.readthedocs.io/en/latest/usage.html>`_. See also
tutorial examples in the ``example/`` folder or by following different sections on the usage page.

How to cite
-----------

.. start-citing

If you use IDPConformerGenerator, please always cite its original publication::

    IDPConformerGenerator: A Flexible Software Suite for Sampling the Conformational Space of Disordered Protein States
    João M. C. Teixeira, Zi Hao Liu, Ashley Namini, Jie Li, Robert M. Vernon, Mickaël Krzeminski, Alaa A. Shamandy, Oufan Zhang, Mojtaba Haghighatlari, Lei Yu, Teresa Head-Gordon, and Julie D. Forman-Kay
    The Journal of Physical Chemistry A 2022 126 (35), 5985-6003
    DOI: 10.1021/acs.jpca.2c03726

If you use the Local Disordered Region Sampling (LDRS) module, please also cite::

    Zi Hao Liu, João M C Teixeira, Oufan Zhang, Thomas E Tsangaris, Jie Li, Claudiu C Gradinaru, Teresa Head-Gordon,
    Julie D Forman-Kay, Local Disordered Region Sampling (LDRS) for ensemble modeling of proteins with experimentally undetermined
    or low confidence prediction segments, Bioinformatics, Volume 39, Issue 12, December 2023, btad739,
    https://doi.org/10.1093/bioinformatics/btad739

.. end-citing

Version
-------
v0.7.30

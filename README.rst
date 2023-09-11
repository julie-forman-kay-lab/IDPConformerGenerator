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

.. image:: https://img.shields.io/badge/LDRS-10.1101%2F2023.07.25.550520-blue
    :target: https://doi.org/10.1101/2023.07.25.550520
    :alt: LDRS Preprint

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

If you use the Local Disordered Region Sampling (LDRS) module, please cite::

    Local Disordered Region Sampling (LDRS) for Ensemble Modeling of Proteins with Experimentally Undetermined or Low Confidence Prediction Segments
    Zi Hao Liu, João M.C. Teixeira, Oufan Zhang, Thomas E. Tsangaris, Jie Li, Claudiu C. Gradinaru, Teresa Head-Gordon, Julie D. Forman-Kay
    bioRxiv 2023.07.25.550520; doi: https://doi.org/10.1101/2023.07.25.550520

.. end-citing

Version
-------
v0.7.5

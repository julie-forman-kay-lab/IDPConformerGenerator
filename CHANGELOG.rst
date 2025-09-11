Versioning
==========

This project follows strictly `Semantic Versioning 2.0 <https://semver.org/#semantic-versioning-200>`_ for version control. 

All additions to the ``master`` branch are done by PR followed by its respective version increment.
While in version ``0``, minor and patch upgrades converge in the ``patch`` number while major upgrades are reflected in the ``minor`` number.

Changelog
=========

v0.7.31 (2025-09-11)
------------------------------------------------------------

* Add new sublcient ``contacts`` to build and evaluate extended idpconfgen database
* Add new subclient ``complex`` to analyze sequences against the database and build dynamic complexes
* Add ``BioPython`` as a dependency for calulation purposes (Shrake-Rupley SASA method)
* Optimize ``next_seeker`` protocol to identify Linker-IDR solutions
* Add a new property ``residues_splitted`` in ``libstructure.py`` to get the residue numbers in all chains
* Add new functions ``split_consecutive_groups()`` and ``update_chars_lower()`` in ``libparse.py``
* New library ``libmultichain.py`` added to process multiple chains for building complexes
* Added pKa, pKb, and pKx values to known residues in ``definitions.py``
* Added plotting capability for generated contact heatmaps and automated font-size calculation in ``plotfuncs.py``
* Catch unique KeyError when calculating torsion angles in ``cli_torsions.py``

v0.7.30 (2025-07-11)
------------------------------------------------------------

* Emergency fix for required packaged since setuptools 80.0.0 is not supported

v0.7.29 (2025-07-10)
------------------------------------------------------------

* Minor update the environment requirements specifying Python version to be between 3.8 and 3.12.

v0.7.28 (2025-03-11)
------------------------------------------------------------

* Corrected a bug in issue #287 regarding building with ``ldrs`` for 3 or more chains

v0.7.27 (2025-02-14)
------------------------------------------------------------

* Adds ``--force-long`` flag into ``build`` client to force long-method for sequences longer than 100 AA

v0.7.26 (2025-01-14)
------------------------------------------------------------

* Update documentation to clarify DSSPv3 installation
* Update usage documentation for recognition of PDB files for ``ldrs``
* Addresses issue #283

v0.7.25 (2024-05-09)
------------------------------------------------------------

* Update user install instructions to use only ``miniconda3``.

v0.7.24 (2024-05-07)
------------------------------------------------------------

* Added first steps page to documentation.

v0.7.23 (2024-03-27)
------------------------------------------------------------

* Fixes bug in ``inter_chain_cc()`` for ``ldrs`` so conformers generated with ``-scm mcsce`` also work

v0.7.22 (2024-03-22)
------------------------------------------------------------

* Fixes issue where clash-checking IDR combinations was not being performed between different chains

v0.7.21 (2024-03-13)
------------------------------------------------------------

* When a given sequence overlaps with sequence from folded template in ``ldrs``, ignore building on that chain
* Fixes issue #269

v0.7.20 (2024-02-01)
------------------------------------------------------------

* Add a note for users in the ``README.rst`` regarding conformers generated with ``ldrs`` prior to version 0.7.17

v0.7.19 (2024-01-29)
------------------------------------------------------------

* Update Sphinx requirement to v5 to fix ``docs`` workflow and ``tox`` error

v0.7.18 (2024-01-25)
------------------------------------------------------------

* Correct case where multiple chains are given in a template and only some need ``ldrs``
* Update documentation to clarify sequence formatting example for complexes

v0.7.17 (2024-01-20)
------------------------------------------------------------

* Fix bug where torsion angles from ``ldrs`` would be swapped sometimes due to mirroring in rotation matrix

v0.7.16 (2024-01-11)
------------------------------------------------------------

* Update citation for ``ldrs`` module to Bioinformatics paper

v0.7.15 (2024-01-05)
------------------------------------------------------------

* Bug-fix for ``ldrs`` file processing during ``psurgeon`` stage for single IDR cases

v0.7.14 (2024-01-04)
------------------------------------------------------------

* Correct reference hyperlink in usage documentation

v0.7.13 (2024-01-03)
------------------------------------------------------------

* Update documentation for clarity
* Update example folder with processing AlphaFold structures

v0.7.12 (2023-11-17)
------------------------------------------------------------

* Update documentation for DSSP due to clarity with version requirements
* Set default of ``-cmd`` in ``sscalc`` to ``mkdssp``

v0.7.11 (2023-11-02)
------------------------------------------------------------

* Update citing information.

v0.7.10 (2023-10-27)
------------------------------------------------------------

* Correct reading CIF files without `#` end block.
* Fix minor bugs for reading CIF templates in ``ldrs``

v0.7.9 (2023-10-27)
------------------------------------------------------------

* Update `libstructure.Structure` documentation.

v0.7.8 (2023-10-02)
------------------------------------------------------------

* Bug-fix in ``ldrs`` for input sequence recognition

v0.7.7 (2023-09-26)
------------------------------------------------------------

* Bug-fix MC-SCE integration with force-field topology generation

v0.7.6 (2023-09-25)
------------------------------------------------------------

* Update documentation for installation and usage

v0.7.5 (2023-09-11)
------------------------------------------------------------

* Fix tests for ``cli_build.py`` that arose from #245

v0.7.4 (2023-09-11)
------------------------------------------------------------

* Add ``--long`` and ``--long-ranges`` to the ``build`` subclient to build long (300+ AA) IDPs much quicker
* Change default backbone energy ``-etbb`` in ``build`` to 100.0
* Change default sidechain energy ``-etss`` in ``build`` to 250.0

v0.7.3 (2023-09-11)
------------------------------------------------------------

* Add select ionic radii from CRC Handbook of Chemistry and Physics, 82nd Ed
* Update ``merge`` module to use multiprocessing
* Increase efficiency of clash-checking algorithm in ``ldrs_helper.py``
* Acceptance and automated processing of proteins with a membrane/bilayer
* Automated clash-checking between all combinations of IDRs built
* Added ability for multi-chain protein detection and processing
* Added ability to build IDRs on multi-chain complexes
* Change default backbone energy ``-etbb`` in ``ldrs`` to 100.0
* Change default sidechain energy ``-etss`` in ``ldrs`` to 250.0

v0.7.2 (2023-07-27)
------------------------------------------------------------

* Correct documentation formatting reference to MC-SCE

v0.7.1 (2023-07-27)
------------------------------------------------------------

* Update documentation formatting for LDRS
* Update README citations
* Temporarily removed MacOS GitHub actions tests

v0.7.0 (2023-07-25)
------------------------------------------------------------

* Added new module ``ldrs`` to build IDRs in the context of a folded region
* Added new module ``resre`` to rename certain residue names of multiple PDB files
* Minor bug-fixes to numpy references (#231)

v0.6.17 (2023-06-12)
------------------------------------------------------------

* Bug-fix for installation and tox build/pr

v0.6.16 (2022-11-08)
------------------------------------------------------------

* Update citation to J Phys Chem A reference

v0.6.15 (2022-09-26)
------------------------------------------------------------

* Added new subclient ``bgeodb`` to append to the database exact bgeos
* Implement ``exact`` bgeo-strategy
* Add documentation regarding new features with ``exact`` and ``bgeodb``
* PR #220

v0.6.14 (2022-07-13)
------------------------------------------------------------

* Add tests to ``cli_build.main`` (#225)

v0.6.13 (2022-07-13)
------------------------------------------------------------

* Correct how some loops are written due to redundancies (#226)

v0.6.12 (2022-07-07)
------------------------------------------------------------

* lint files

v0.6.11 (2022-06-14)
------------------------------------------------------------

* Add option to build backbone conformers with fixed bond angle and
  distances. (#217)

v0.6.10 (2022-06-14)
------------------------------------------------------------

* Implement boxplot plotting feature for bond angle distributions
* Append documentation for ``bgeo`` subclient
* Ability to use ``degrees`` for ``bgeo`` subclient
* Ability to change the name of the ouput file for ``bgeo`` subclient
* PR #219

v0.6.9 (2022-06-13)
------------------------------------------------------------

* Correct CLI help for the ``--mem`` flag in the ``sethpc`` subclient (#218)

v0.6.8 (2022-06-13)
------------------------------------------------------------

* Add bioRxiv citation

v0.6.7 (2022-05-31)
------------------------------------------------------------

* update actions from ``master`` to ``main``
* update install instructions for MC-SCE in GRAHAM (#215)

v0.6.6 (2022-05-27)
------------------------------------------------------------

* Add ``pr`` env to ``tox``
* Update CI workflows
* Update ReadTheDocs python version to 3.8
* Dropped python 3.7 after Numpy

v0.6.5 (2022-05-25)
------------------------------------------------------------

* Correct typo bugs in ``fastext`` and ``bgeo``
* General lints

v0.6.4 (2022-05-25)
------------------------------------------------------------

* Re-licensed to Apache-2.0

v0.6.3 (2022-05-25)
------------------------------------------------------------

v0.6.2 (2022-05-25)
------------------------------------------------------------

* Update usage instructions for `bgeo`

v0.6.1 (2022-05-25)
------------------------------------------------------------

* updated GRAHAM install instruction (#207)

v0.6.0 (2022-05-24)
------------------------------------------------------------

* Add bond geometry option to build with `Int2Cart` software
* PR #203

v0.5.1 (2022-05-24)
------------------------------------------------------------

* add plot functions to ``sscalc`` and ``torsions``
* PR #198

v0.5.0 (2022-05-24)
------------------------------------------------------------

* Add residue tolerance matrices: EDSS50
* Update/improve parameters to residue tolerance options
* PR #183

v0.4.10 (2022-05-23)
------------------------------------------------------------

* Add documentation RTD format
* Add documentation for several features and examples
* PR #171

v0.4.9 (2022-05-23)
------------------------------------------------------------

* Add ``sethpc`` client.
* Add ``merge`` client.
* PR #202

v0.4.8 (2022-05-23)
------------------------------------------------------------

* Add ``stats`` client
* Add ``search`` client
* PR #200

v0.4.7 (2022-05-23)
------------------------------------------------------------

* update CI methods
* PR #205

v0.4.6 (2022-04-22)
------------------------------------------------------------

v0.4.5 (2022-04-21)
------------------------------------------------------------

v0.4.4 (2022-03-29)
------------------------------------------------------------

* Fixes MC-SCE integration when sidechain packing fails
* Corrects MC-SCE installation
* #190

v0.4.3 (2022-03-26)
------------------------------------------------------------

v0.4.2 (2022-03-20)
------------------------------------------------------------

v0.4.1 (2022-03-17)
------------------------------------------------------------

* Adds support for single residues when not specified. Addresses #184

v0.4.0 (2022-03-15)
------------------------------------------------------------

* Integrates the MC-SCE protocol in the building process as part of the
  sidechain packing method options.

v0.3.3 (2022-03-14)
------------------------------------------------------------

* removes assert in 0.3.2

v0.3.2 (2022-03-14)
------------------------------------------------------------

* improves regex creation to avoid silent bugs in possible parallel
  futures

v0.3.1 (2022-03-13)
------------------------------------------------------------

* incorporates `G` in `H` when treating DSSP with reduced labels

v0.3.0 (2022-03-13)
------------------------------------------------------------

* see #168
* Revisited the whole regex sampling machinery during conformer building
* A initial major part for preparing the regex database was dropped
* applied multiprocessing to the regex database preparation steps
* updated the `cli_build` API with 4 new command options
* dropped using regex in the `cli_build` command line

v0.2.6 (2022-03-13)
------------------------------------------------------------

* corrected `sscalc` from * input in command-line #175

v0.2.5 (2022-03-11)
------------------------------------------------------------

* Implemented capacity to read PDBs with names different from cull #167

v0.2.4 (2022-03-11)
------------------------------------------------------------

* implemented support for N-terminal Proline residues #166

v0.2.3 (2022-03-08)
------------------------------------------------------------

* corrected energy.log #162

v0.2.2 (2022-03-07)
------------------------------------------------------------

* incorporated `libfuncpy` internally

v0.2.1 (2022-03-03)
------------------------------------------------------------

v0.2.0 (2022-02-10)
------------------------------------------------------------

v0.1.0 (2021-07-24)
------------------------------------------------------------

* Implements energy calculation to individual pairs. Energy threshold
  can now be compared to `pairs` or `whole`.

v0.0.24 (2021-07-01)
------------------------------------------------------------

* Corrects `make_folder` function in `cli_build`.

v0.0.23 (2021-07-01)
------------------------------------------------------------

* Added libfuncpy to requirements.yml

v0.0.22 (2021-06-30)
------------------------------------------------------------

* Users are now able to fully configure the size of fragments and probabilities,
    via the flag `-xp` that expects a two column file.

v0.0.21 (2021-06-28)
------------------------------------------------------------

* Now build prints log to terminal.
* improved other minor logging issues

v0.0.20 (2021-06-21)
------------------------------------------------------------

* Decoupled ``energy-threshold`` parameters. Now Backbone and sidechains,
    can be configured separately.

v0.0.19 (2021-06-14)
------------------------------------------------------------

* Saves a table with energy values per conformer.
* Crash reports now saved in execution folder (CLI build).

v0.0.18 (2021-06-10)
------------------------------------------------------------

* Improves sampling of multiple secondary structure regexes.
    Now, when given multiple regex, angle sampling will be biased towards
    the number of occurrences in each regex.

v0.0.17 (2021-06-10)
------------------------------------------------------------

* Corrects bug in Coulomb formula

v0.0.16 (2021-06-09)
------------------------------------------------------------

* Add output-folder option for the ``build`` interface

v0.0.15 (2021-06-09)
------------------------------------------------------------

* corrected typo in example/ commands

v0.0.14 (2021-06-05)
------------------------------------------------------------

* Users can now select single residue fragment size
* ``-xp`` parameter was updated with checks and completion

v0.0.13 (2021-05-28)
------------------------------------------------------------

* Added usage example and documentation.

v0.0.12 (2021-05-28)
------------------------------------------------------------

* Corrects path suffix evaluation in ``cli_torsions.py``

v0.0.11 (2021-05-28)
------------------------------------------------------------

* corrects var name bug in ProgressBar

v0.0.10 (2021-05-27)
------------------------------------------------------------

* Implements residue substitution/tolerance during conformer build

v0.0.9 (2021-05-27)
------------------------------------------------------------

* user can now define the fragment size selection probabilities

v0.0.8 (2021-05-09)
------------------------------------------------------------

* Expands try:catch to avoid index error when restarting conformer

v0.0.7 (2021-05-09)
------------------------------------------------------------

* saves version number to file before running a client

v0.0.6 (2021-04-20)
------------------------------------------------------------

* additional functions for logging
* add logging to build and other parts

v0.0.5 (2021-04-19)
------------------------------------------------------------

* added ``--energy-threshold`` flag to control energy threshold after sidechain addition

v0.0.4 (2021-04-19)
------------------------------------------------------------

* ``builder`` CLI now accepts ``.fasta`` files.

v0.0.3 (2021-04-19)
------------------------------------------------------------

* added matplotlib in requirements.yml as dependency

v0.0.2 (2021-04-03)
------------------------------------------------------------

* corrects variable name in ``libbuild`` that was breaking sidechain
    construction.

v0.0.1 (2021-04-02)
------------------------------------------------------------

* added CI integration files

v0.0.0
------

* Any development previous to version 0.0.1 is registered in PRs up to #102.

Versioning
==========

This project follows strictly `Semantic Versioning 2.0 <https://semver.org/#semantic-versioning-200>`_ for version control. 

All additions to the ``master`` branch are done by PR followed by its respective version increment.
While in version ``0``, minor and patch upgrades converge in the ``patch`` number while major upgrades are reflected in the ``minor`` number.

Changelog
=========

* Improve variable names related to ``int2cart`` in ``build`` that were
  prone to bugs.

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

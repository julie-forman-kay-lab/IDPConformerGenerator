Versioning
==========

This project follows strictly `Semantic Versioning 2.0 <https://semver.org/#semantic-versioning-200>`_ for version control. 

All additions to the ``master`` branch are done by PR followed by its respective version increment.
While in version ``0``, minor and patch upgrades converge in the ``patch`` number while major upgrades are reflected in the ``minor`` number.

Changelog
=========

* Users are now able to fully configure the size of chunks and probabilities,
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

* Users can now select single residue chunk size
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

* user can now define the chunk size selection probabilities

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

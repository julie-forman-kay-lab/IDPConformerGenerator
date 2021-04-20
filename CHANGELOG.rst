Versioning
==========

This project follows strictly `Semantic Versioning 2.0 <https://semver.org/#semantic-versioning-200>`_ for version control. 

All additions to the ``master`` branch are done by PR followed by its respective version increment.
While in version ``0``, minor and patch upgrades converge in the ``patch`` number while major upgrades are reflected in the ``minor`` number.

Changelog
=========

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

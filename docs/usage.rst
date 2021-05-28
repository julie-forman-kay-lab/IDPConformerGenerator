=====
Usage
=====

Command-line
------------

To execute :code:`idpconfgen` command-line just run :code:`idpconfgen` in your
terminal window, after :ref:`installation <Installation>`::

    idpconfgen

or::

    idpconfgen -h

both will output the help menu.

:code:`idpconfgen` has several interfaces that perform different functions.
However, there is a sequence of interfaces that need to be used to prepare the
local database of torsion angles and the files needed to build conformers. After
these operations executed, you will end up with a single :code:`json` file that
you can use to feed the build calculations. The other files are safe to be
removed.

The :code:`example` folder contains instructions to setup
IDPConformerGenerator. We review them here as well:

.. include:: ../example/README.rst


IDP Conf Gen Library
--------------------

To use IDP Conformer Generator in your project, import it as a library::

    import idpconfgen

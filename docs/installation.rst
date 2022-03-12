============
Installation
============

To install IDPConformerGenerator (idpconfgen for short):

From source
-----------

Clone from the official repository::

    git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator

And navigate to the new ``IDPConformerGenerator`` folder::

    cd IDPConformerGenerator

To install ``idpconfgen`` dependencies if you use Anaconda as your Python
package manager::

    conda env create -f requirements.yml
    conda activate idpconfgen

To install ``idpconfgen`` dependencies with ``pip``::

    virtualenv idpconfgenenv --python=3.9
    source idpconfgenvenv/bin/activate
    pip install -r requirements.txt

Install ``idpconfgen`` in development mode in order for your installation to be
always up-to-date with the repository::

    python setup.py develop --no-deps

**Remember** to active the ``idpconfgen`` environment every time you open a new
terminal window, from within the repository folder, choose yours::

    conda activate idpconfgen
    source idpconfgenvenv/bin/activate


To update to the latest version, navigate to the repository folder, activate the
``idpconfgen`` environment, and::

    git pull
    python setup.py develop --no-deps

Your installation will become up to date with the latest developments.

From PyPI
---------

To install ``idpconfgen`` from PyPI, first create a new environment::

    virtualenv idpconfgenenv --python=3.9
    source idpconfgenvenv/bin/activate

Install ``idpconfgen``::

    pip install IDPConformerGenerator

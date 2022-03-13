============
Installation
============
IDPConformerGenerator v0.3.X has been tested to work with Ubuntu 18.04 LTS and 20.04 LTS as well as on WSL2.0 and Graham@ComputeCanada.
Although it's recommended to run IDPConfGen on UNIX based OS, it will work with Windows based OS as well as long as the pre-requisites are met.

**Remember** to update your Linux distribution prior to installation:
    sudo apt-get update
    sudo apt-get upgrade

Pre-installation Reqirements
----------------------------
An up-to-date version of anaconda3, pip, gcc, java, csh, perl, automake, boostlib, and dssp are required

    sudo apt install python3-pip

    sudo apt-get install build-essential

    sudo apt-get install default-jdk

    sudo apt-get install csh

    sudo perl -MCPAN -e'install "LWP::Simple"'
    sudo apt-get install libswitch-perl

    sudo apt-get install automake

    sudo apt-get install libboost-all-dev

    sudo apt install dssp

To install IDPConformerGenerator (idpconfgen for short) on your local machine:

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


To install IDPConformerGenerator (idpconfgen for short) on your Graham@ComputeCanada:

From source
-----------
Log-in and make sure you're in the /home directory

    cd

Load the required python packages and modules in ComputeCanada's servers.
    
    module load scipy-stack dssp boost

Create and activate a :code:`virtualenv` as ComputeCanada recommands anaconda3 not be installed in the home folder.

    virtualenv --no-download idpconfgen
    source idpconfgen/bin/activate

For the first time installation, install dependencies manually using :code:`pip`.
Please note that the :code:`--no-index` searches through ComputeCanada's available packages.
If they're not available, it will install from the web.

    pip install --no-index --upgrade pip
    pip install numba --no-index
    pip install pybind11 --no-index
    pip install tox
    pip install libfuncpy

We are ready to clone from source and installation from here will be similar to local.

    git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator
    cd IDPConformerGenerator

Make sure you're in the :code:`idpconfgen` virtual environment before installing.

    python setup.py develop --no-deps
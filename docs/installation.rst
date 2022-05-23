============
Installation
============

IDPConformerGenerator v0.4.X has been tested to work with Ubuntu 18.04 LTS and 20.04 LTS as well as on WSL2.0 and the Graham cluster, 
an HPC resource of the Digital Research Alliance of Canada (DRAC). Although it's recommended to run IDPConfGen on UNIX based OS, 
it will work with Windows based OS as well as long as the pre-requisites are met.

**Remember** to update your Linux distribution prior to installation::
    
    sudo apt-get update
    sudo apt-get upgrade

Pre-installation Reqirements
----------------------------

(Required) An up-to-date version of anaconda3, and pip.
(Required) For DSSP installation, please refer to: https://github.com/julie-forman-kay-lab/IDPConformerGenerator/issues/48

(Recommended) For MC-SCE installation, please refer to: https://github.com/THGLab/MCSCE

(Optional) To use CSSS with CheSPI, please refer to: https://github.com/protein-nmr/CheSPI
(Optional) To use CSSS with δ2D, please refer to: https://github.com/carlocamilloni/d2D

To install IDPConformerGenerator (idpconfgen or IDPConfGen for short) on your local machine:

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

    # Installation with Anaconda
    conda activate idpconfgen

    # Installation with pip
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


From source on the Graham Cluster (DRAC)
----------------------------------------

Log-in and make sure you're in the /home directory::

    cd

Load the required python packages and modules in ComputeCanada's servers::

    module load scipy-stack dssp boost

Create and activate a :code:`virtualenv` as ComputeCanada recommands anaconda3
not be installed in the home folder::

    virtualenv --no-download idpconfgen
    source idpconfgen/bin/activate

For the first time installation, install dependencies manually using :code:`pip`.
Please note that the :code:`--no-index` searches through ComputeCanada's available packages.
If they're not available, it will install from the web::

    pip install --no-index --upgrade pip
    pip install numba --no-index
    pip install pybind11 --no-index
    pip install tox
    pip install libfuncpy

We are ready to clone from source and installation from here will be similar to
local::

    git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator
    cd IDPConformerGenerator

Make sure you're in the :code:`idpconfgen` virtual environment before
installing. Install with::

    python setup.py develop --no-deps

Install MC-SCE
--------------

IDPConformerGenerator can integrate MC-SCE to generate sidechains on top of the
backbone conformers it generates, on the fly. For that you need to install MC-SCE
on top of the `idpconfgen` Python environment. First, install IDPConfGen as
described above. Next, follow these steps::

    # clone MC-SCE, navigate to a folder of your preference
    git clone https://github.com/THGLab/MCSCE

    # Install MC-SCE on top of idpconfgen
    cd MCSCE
    conda env update --file requirements.yml --name idpconfgen

    # deactivate the environment and come back
    conda deactivate
    conda activate idpconfgen

    # install MC-SCE within the `idpconfgen` environment
    python setup.py develop --no-deps

    # navigate back to the idpconfgen github folder and re-run
    python setup.py develop --no-deps

Now, if you choose the flag :code:`-scm mcsce`, IDPConfGen will use MC-SCE to build
sidechains as backbone conformers are generated. You will see :code:`idpconfgen build -h`
has a specific group of parameters dedicated to MC-SCE, you can explore those as well.

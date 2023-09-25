Installation
============

IDPConformerGenerator uses only Python-based APIs for which we expect it to run
native on any system Python can run, as long as the third-party installation
requirements are met.

We tested IDPConfGen on Ubuntu 18.04 LTS and 20.04 LTS as well as on WSL2.0 and
the Graham cluster, an HPC resource of the Digital Research Alliance of Canada
(DRAC).


Follow the steps below to install IDPConformerGenerator (``idpconfgen`` or
``IDPConfGen`` for short) on your local machine:

From source
-----------

Clone from the official repository::

    git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator

And navigate to the new ``IDPConformerGenerator`` folder::

    cd IDPConformerGenerator

Run the following commands to install ``idpconfgen`` dependencies if you use
Anaconda as your Python package manager::

    conda env create -f requirements.yml
    conda activate idpconfgen

.. note::
    If you don't use Anaconda to manage your Python installations, you can use
    ``virtualenv`` and the ``requirements.txt`` file following the commands:

    | ``virtualenv idpcgenv --python=3.9``
    | ``source venv/bin/activate``
    | ``pip install -r requirements.txt``

    If you have difficulties installing ``idpconfgen``, raise an Issue in the
    main GitHub repository, and we will help you.

Install ``idpconfgen`` in development mode in order for your installation to be
always up-to-date with the repository::

    python setup.py develop --no-deps

.. note::
    The above applies also if you used ``virtualenv`` instead of ``conda``.

**Remember** to active the ``idpconfgen`` environment every time you open a new
terminal window, from within the repository folder, choose yours::

    # Installation with Anaconda
    conda activate idpconfgen

    # Installation with virtualenv
    source idpcgenv/bin/activate


To update to the latest version, navigate to the repository folder, activate the
``idpconfgen`` python environment as described above, and run the commands::

    git pull

    # if you used anaconda to create the python environment, run:
    conda env update -f requirements.yml

    # if you used venv to create the python environment, run:
    pip install -r requirements.txt  --upgrade

    python setup.py develop --no-deps

Your installation will become up to date with the latest developments.

From source on the Graham Cluster (DRAC)
----------------------------------------

Log-in and make sure you are in your user home directory::

    cd

Load the required python packages and modules on Graham's servers::

    module load scipy-stack dssp boost

Create and activate a ``virtualenv`` as DRAC recommends ``anaconda3``
not be installed in the home folder::

    virtualenv --no-download idpcgenv
    source idpcgenv/bin/activate

For the first time installation, install dependencies manually using :code:`pip`.
Please note that the :code:`--no-index` searches through DRAC's available packages.
If they're not available, it will install from the web::

    pip install --no-index --upgrade pip
    pip install numba --no-index
    pip install pybind11 --no-index
    pip install tox

We are ready to clone from source and installation from here will be similar to
local::

    git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator
    cd IDPConformerGenerator

Make sure you're in the :code:`idpcgenv` virtual environment before
installing. Install with::

    python setup.py develop --no-deps

When you login again to your cluster account remember to reactive the
``idpconfgen`` environment before using ``idpconfgen``::

    cd
    source idpcgenv/bin/activate

Installing third-party software
-------------------------------

Some functionalities of ``idpconfgen`` require third-party software. These
are not mandatory to install unless you want to use such operations.

DSSP
````

IDPConfGen uses `DSSP <https://github.com/cmbi/dssp>`_ to calculate secondary
structures. However, you only need DSSP if you are generated the database from
scratch. If you use a prepared database JSON file you don't need to install
DSSP.

To install DSSP, please refer to:
https://github.com/julie-forman-kay-lab/IDPConformerGenerator/issues/48

Install MC-SCE
``````````````

IDPConformerGenerator can integrate MC-SCE to generate sidechains on top of the
backbone conformers it generates, on the fly. For that you need to install
MC-SCE on top of the `idpconfgen` Python environment. First, install IDPConfGen
as described above. Next, follow these steps::

    # ensure you are in the parent IDPConformerGenerator GitHub folder
    # clone and enter the MC-SCE GitHub repository
    git clone https://github.com/THGLab/MCSCE
    cd MCSCE

    # Make sure you're in the idpconfgen environment then
    # install the additional dependencies using pip
    pip install tensorflow
    pip install tqdm
    pip install pathos

    # Install MC-SCE on top of IDPConformerGenerator
    python setup.py develop --no-deps

    # Navigate back to the IDPConformerGenerator GitHub folder and install
    # `idpconfgen` again if needed
    cd ../IDPConformerGenerator
    python setup.py develop --no-deps

Now, if you choose the flag :code:`-scm mcsce` in ``idpconfgen build`` command,
IDPConfGen will use MC-SCE to build sidechains as backbone conformers are
generated. You will see :code:`idpconfgen build -h` has a specific group of
parameters dedicated to MC-SCE, you can explore those as well.

For installation on a cluster via virtualenv, dependencies need to be manually installed
as the following for MC-SCE::

    # ensure you're in the idpcgenv and the IDPConformerGenerator GitHub folder
    git clone https://github.com/THGLab/MCSCE

    # MC-SCE also requires numba and tox but that's already handled in previous steps
    pip install tensorflow --no-index
    pip install keras --no-index
    pip install tqdm --no-index
    pip install pathos --no-index

    # cd into the MCSCE GitHub folder and install MC-SCE
    cd MCSCE
    python setup.py develop --no-deps

    # cd back into the IDPConformerGenerator GitHub folder and install idpconfgen on top of MC-SCE
    cd ..
    python setup.py develop --no-deps

Install Int2Cart
````````````````

IDPConformerGenerator can use Int2Cart on the fly to optimize bond geometries
of the backbones calculated. For this feature, you must have a CUDA compatible
GPU as well as install Int2Cart on top of the ``idpconfgen`` Python environment.
First, install IDPConfGen as described above. Next, follow these steps. Please
note that these steps are the same if you have installed idpconfgen through ``virtualenv``::

    # ensure you are in the IDPConformerGenerator GitHub folder

    # Install a pre-requisite of Int2Cart: sidechainnet
    git clone https://github.com/THGLab/sidechainnet
    cd sidechainnet
    pip install -e .
    cd ..

    # Install Int2Cart
    git clone https://github.com/THGLab/int2cart
    cd int2cart
    pip install -e .
    pip install pyyaml
    cd ..

    # you should be back in the IDPConformerGenerator GitHub folder


Running Int2Cart on the Graham cluster requires GPU allocations and ``module load cuda``.
Otherwise, installation is the same within the ``idpconfgen`` virtualenv.

Troubleshooting Int2Cart installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If IDPConfGen is still giving you an error that Int2Cart is not installed, please test this import
in the ``idpconfgen`` environment::

    python
    >>> from modelling.models.builder import BackboneBuilder

If you receieve this error: ``ImportError: TensorBoard logging requires TensorBoard version 1.15 or above``,
do the following::

    pip install tensorboard==1.15.0

CheSPI
``````

To use CSSS via the ``idpconfgen csssconv`` command you need CheSPI. Please
refer to https://github.com/protein-nmr/CheSPI to install CheSPI.

δ2D
```

The use δ2D via the ``idpconfgen csssconv`` command you need δ2D.
Please refer to https://github.com/carlocamilloni/d2D.

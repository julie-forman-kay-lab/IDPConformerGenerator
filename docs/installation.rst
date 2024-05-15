Installation
============

IDPConformerGenerator uses only Python-based APIs for which we expect it to run
native on any system Python can run, as long as the third-party installation
requirements are met.

Please note that `SPyCi-PDB <https://github.com/julie-forman-kay-lab/SPyCi-PDB>`_ and
`X-EISDv2 <https://github.com/THGLab/X-EISDv2>`_ can be installed on top of the ``idpconfgen``
Python environment. It is actually recommended since they both have IDPConfGen as a dependency.

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

.. note::

    The ``requirements.yml`` describe the Python dependencies of
    IDPConformerGenerator. If you are skilled managing python environments you
    can go on your own. Otherwise, you can calmly follow our install
    instructions.

At the end of the installation process, you will have a ``miniconda3``
directory inside the ``IDPConformerGenerator`` main directory where the whole
installation is placed. If your ever want to delete ``IDPCG`` from your computer,
simply delete the ``IDPConformerGenerator`` directory.

**To install IDPConfGen**, run the following three commands. Wait until one
finishes before running the second one::

    ./install_miniconda3.sh
    source activate.sh
    ./install_deps.sh

Once this finishes, ``idpconfgen`` is ready to be used. Go to the :ref:`usage
<Usage>` and continue from there.

**Remember** to active the ``idpconfgen`` environment every time you open a
new terminal window. For that, navigate to the ``IDPConformerGenerator``
repository folder and ``source`` the ``activate.sh`` file::

    cd path/to/my/IDPConformerGenerator
    source activate.sh

Update
------

To update to the latest version, open a new terminal window, and navigate to the
``IDPConformerGenerator`` source folder. Remove the ``miniconda3`` environment::

    rm -rf miniconda3

Update the source to the latest version::

    git pull

Reinstall the project and it's dependencies. Run the following commands, one
after the other, wait for them to finish before running the next one::

    ./install_miniconda3.sh
    source activate.sh
    ./install_deps.sh

Your installation will become up to date with the latest developments.  If you
had installed MC-SCE, Int2Cart, SPyCi-PDB, or X-EISDv2,  you need to reinstall
them again in the ``idpconfgen`` environment.

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

Please note we are only compatible with DSSP versions 2 and 3. If you have
installed DSSP version 4 (check by using the command ``mkdssp --version``) please
refer `to this issue <https://github.com/julie-forman-kay-lab/IDPConformerGenerator/issues/48>`_
for a proper re-installation after removing DSSP version 4.

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

Installing back-calculators and reweighting protocols
-----------------------------------------------------

Both SPyCi-PDB and X-EISDv2 have been developed in-house with considerations
for protein structural ensembles in mind. We recommend to install both of
these packages on-top of the ``idpconfgen`` environment for streamlined usage.

Install SPyCi-PDB
`````````````````

Clone the SPyCi-PDB repository to the parent directory of where IDPConformerGenerator was cloned::
    
    git clone https://github.com/julie-forman-kay-lab/SPyCi-PDB

Activate the ``idpconfgen`` environment and install the missing dependencies::

    pip install pandas
    pip install natsort

Move into the SPyCi-PDB directory and install on top of IDPConfGen::

    cd SPyCi-PDB
    python setup.py develop --no-deps

.. note::

    For the usage of all the back-calculators, please refer to the installation
    directions documented for SPyCi-PDB that can be found `here <https://spyci-pdb.readthedocs.io/en/stable/installation.html>`_.

    The publication for SPyCi-PDB can be found `here <https://joss.theoj.org/papers/10.21105/joss.04861>`_.

Install X-EISDv2
````````````````

Clone the X-EISDv2 repository to the parent directory of where IDPConformerGenerator was cloned::
    
    git clone https://github.com/THGLab/X-EISDv2

Activate the ``idpconfgen`` environment and install the missing dependencies.
You can skip this step if you've already installed SPyCi-PDB::

    pip install pandas
    pip install natsort

Move into the X-EISDv2 directory and install on top of IDPConfGen::

    cd X-EISDv2
    python setup.py develop --no-deps

.. note::

    Usage directions for X-EISDv2 can be found within the command-line interface
    by using the ``-h`` command. For example: ``xeisd -h``, ``xeisd score -h``.

    The publication for X-EISD can be found `here <https://pubs.acs.org/doi/10.1021/jacs.6b00351>`_.
    The original X-EISD repository can be found `here <https://github.com/THGLab/X-EISD>`_.

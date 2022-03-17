============
Installation
============

How to install the current beta-version:

1. Clone this repository::

    git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator

And navigate to the new ``IDPConformerGenerator`` folder::

    cd IDPConformerGenerator

2. Install the required dependencies manually, by creating a dedicated environment in your Anaconda::

    conda env create -f requirements.yml

If you don't Anaconda to manage your Python installations, and have difficulties
installing ``IDPConfGen``, raise an Issue in the main GitHub repository, and we
will help you.

3. Activate the new conda environment::

    conda activate idpconfgen

4. Install IDPConfGen in development mode in order for your installation to be
always up-to-date with the repository::

    python setup.py develop --no-deps

5. To update to the latest version::

    git pull


Install MCSCE
-------------

IDPConformerGenerator can integrate MCSCE to generate sidechains on top of the
backbone conformers it generates, on the fly. For that you need to install MCSCE
on top of the `idpconfgen` Python environment. First, install IDPConfGen as
described above. Next, follow these steps::

    # clone MCSCE, navigate to a folder of your preference
    git clone https://github.com/THGLab/MCSCE

    # Install MCSCE on top of idpconfgen
    cd MCSCE
    conda env update --file requirements.yml --name idpconfgen

    # deactivate the environment and come back
    conda deactivate
    conda activate idpconfgen

    # navigate back to the idpconfgen github folder and re-run
    python setup-py develop --no-deps

Now, if you choose the flag `-scm mcsce`, IDPConfGen will use MCSCE to build
sidechains as backbone conformers are generated. You will see `idpconfgen build
-h` has a specific group of parameters dedicated to MCSCE, you can explore those
as well.

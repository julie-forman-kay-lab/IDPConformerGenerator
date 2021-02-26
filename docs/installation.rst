============
Installation
============

How to install the current beta-version:

1. Clone this repository
2. Install the require dependencies manually::

    conda env create -f requirements.yml

3. Activate the new conda environment::

    conda activate idpconfgen

4. Install IDPConfGen in development mode in order for your installation to be
always up-to-date with the repository::

    python setup.py develop --no-deps



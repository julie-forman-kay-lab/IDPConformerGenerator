AlphaFlex (AFX-IDPCG) Protocols
===============================
Scripts for building and analysis of AFX-IDPCG ensembles can be found in the ``.zip`` file within this directory. Scripts are archived to prevent
interference with the main IDPConformerGenerator package.

AFX-IDPCG Python scripts must be run within the ``idpconfgen`` Python environment created during the installation of IDPConformerGenerator.

The ``idpconfgen`` conda environment will need to be updated with the following packages::

    conda install -c conda-forge pdbfixer
    pip install pdb-tools
    pip install chardet==5.2.0
    pip install tqdm==4.67.1
    pip install biopython==1.85
    pip install mdtraj==1.11.0

The AlphaFlex database (``AlphaFlex_database_Jul2024.json``) is a standard JSON dictionary database where each key is a UniProt ID and the values correspond
to IDR boundaries (``idrs``), mean PAEs between folded (F) and disordered (D) regions (``mean_pae``), and any interactions
between folded domains where the mean PAE is less than 15 Angstroms is documented in ``interactions``.

The ``AF2_9606_HUMAN_v4_num_residues.json`` is a database that contains the total length of the amino-acid residues from the
AlphaFold2 9606 Human v4 database.

Resources
---------
- AlphaFlex Manuscript (pre-print): https://doi.org/10.1101/2025.11.24.690279
- AlphaFlex Zenodo Repository: https://zenodo.org/records/17684898
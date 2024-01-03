Processing Low-Confidence Predicted Residues
============================================

.. start-description

Although the ``ldrs`` subclient accepts any PDB or mmCIF file to be used
as a template. A script has been prepared in this folder ``/remove_lowconfidence_residues.py``
to automate the removal of low-confidence predicted residues.

Sample thresholds have been given within the Python script that uses the ``idpconfgen``
environment based on the source of each structure prediction algorithm. For example,
a threshold of 70 needs to be applied to AlphaFold structures and 0.7 for ESMFold structures.

Please edit the ``input_file``, ``output_file`` and ``threshold`` variables in the script
before running with ``python remove_lowconfidence_residues.py`` in the ``idpconfgen`` Python
environment.

.. note::
    It is sometimes preferable to remove residues manually using a molecular viewer
    such as PyMOL to avoid ending up with short segments of confident, yet unwanted
    residues.

    For example, 2 "confident" residues in between 10 unconfident residues should be
    removed as well for optimal performance.

.. end-description
